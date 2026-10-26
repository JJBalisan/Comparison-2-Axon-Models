"""Coincidence-detection windows comparable to Myoga et al 2014.

Myoga et al 2014: Nat Commun 5:3790.

Their protocol: two identical EPSGs at relative times stepped in 20 us, peak
conductance set just above G_t (50% spike probability for coincident inputs),
spike probability measured per delay, and the width of that curve at half
maximum reported (221 us without inhibition, 35 C, adult gerbil).

With a fixed input size G, a model spikes at delay dt exactly when its threshold
there is at or below G. With symmetric trial-to-trial noise in the effective
threshold, probability is 50% where threshold(dt) == G, so the 50%-probability
delays follow from a finely resolved threshold curve. That equals the width at
half maximum when probability peaks near 100%, which is how Myoga et al set up
their protocol. If noise keeps the peak lower (e.g. 85%), the half-maximum width
comes out wider than this estimate (~10% in the two-compartment model).
`probability_trials` checks this with noisy runs.

A run spikes when the axon rises `criterion` above the soma: 10 mV (FACTOR) by
default, as in the MATLAB EPSG-pair sweeps, or a fraction of the model's own
spike height, criterion=Relative(0.5) (issue #21; see spike_height).

This is new analysis code, not part of the MATLAB port: it uses the predictive
search of _bisect (bisection helped by interpolation) rather than
BinarySearch.m's halving search.

Event lists (issue #78). The searches above drive an EPSG pair; passing
`events=` drives them with the "Events" stimulus instead, and the protocols
at the end of the module (fibre counts, pre-pulses, Casarez-style ITD sweeps,
are built on it. `events` is an event
*source*, a callable (delay, amp) -> tuple of synaptic.Event, or a plain
event tuple, which scales with the searched amplitude; the search runs it
with I = 1 and varies amp, so a source can hold some events fixed (an IPSG, a
pre-pulse) while it scales the rest. Calls without `events` are untouched.
"""

import dataclasses
from collections.abc import Sequence
from concurrent import futures
from typing import Any, NamedTuple

import numpy as np
from scipy import optimize

from msoaxon import (
    _bisect,
    _common,
    _dispatch,
    _fast,
    _parallel,
    measure,
    synaptic,
)

START = _common.START
T_END = _common.T_END
# spike = axon rises this far above soma (as in the EPSGpair sweeps)
FACTOR = _common.FACTOR
HEIGHT_AT = 2.0  # spike_height's run, as a multiple of the 10 mV threshold
HEIGHT_REL_TOL = 1e-3  # that threshold's tolerance: the height is flat there


@dataclasses.dataclass(frozen=True)
class Relative:
    """A spike criterion as a fraction of the model's own spike height.

    criterion=Relative(0.5) asks for axon - soma to reach half of
    spike_height() for the same configuration, instead of FACTOR mV.

    Attributes:
        fraction: The fraction of the spike height, strictly between 0
            and 1.
    """

    fraction: float

    def __post_init__(self) -> None:
        """Checks the fraction.

        Raises:
            ValueError: If the fraction is not strictly between 0 and 1.
        """
        if not 0 < self.fraction < 1:
            raise ValueError(
                f"Relative needs a fraction strictly between 0 and 1, got "
                f"{self.fraction!r}"
            )


# a float is absolute [mV]; a Relative is a fraction of spike_height()
Criterion = float | Relative


def _respond(
    model: str,
    I: float,
    delay: float,
    node: int,
    v0: float,
    epsg_tau: tuple[float, float],
    site: tuple[str, int, int | None],
    model_kw: dict[str, Any],
    start: float = START,
    *,
    factor: float = FACTOR,
    source: synaptic.Source | None = None,
) -> tuple[bool, float]:
    """Runs one EPSG pair, `delay` apart, or one run of an event source.

    Returns (spiked, peak axon - soma); site = (stim, input_node,
    input_node2); factor is the resolved spike criterion [mV]. A source
    replaces the pair (site and start are then unused) and is run with
    I = 1 and the amplitude passed to it.
    """
    if source is not None:
        return _dispatch.response(
            model,
            "Events",
            0.0,
            0.0,
            1.0,
            node,
            T_END,
            v0,
            1,
            synaptic.event_syn(source, delay, I, T_END),
            factor=factor,
            **model_kw,
        )
    stim, input_node, input_node2 = site
    syn = synaptic.SynParams(t_end=T_END, epsg_tau=tuple(epsg_tau))
    if input_node2 is not None:
        model_kw = {**model_kw, "input_node2": input_node2}
    return _dispatch.response(
        model,
        stim,
        start,
        start + delay,
        I,
        node,
        T_END,
        v0,
        input_node,
        syn,
        factor=factor,
        **model_kw,
    )


def _spikes(
    model: str,
    I: float,
    delay: float,
    node: int,
    v0: float,
    epsg_tau: tuple[float, float],
    site: tuple[str, int, int | None],
    model_kw: dict[str, Any],
    start: float = START,
    *,
    factor: float = FACTOR,
    source: synaptic.Source | None = None,
) -> bool:
    return _respond(
        model,
        I,
        delay,
        node,
        v0,
        epsg_tau,
        site,
        model_kw,
        start,
        factor=factor,
        source=source,
    )[0]


def _source(events: synaptic.EventsLike | None) -> synaptic.Source | None:
    """Returns the event source of a call's events, or None for a pair."""
    return None if events is None else synaptic.as_source(events)


def _check_delay(source: synaptic.Source | None, delay: float) -> None:
    """Rejects a delay for an event list that has none."""
    if isinstance(source, synaptic.Fixed) and delay != 0.0:
        raise ValueError(
            f"a fixed event tuple has no delay (got {delay}); pass an event "
            "source such as synaptic.PairSource to vary it"
        )


def threshold(
    model: str,
    delay: float,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    guess: float = 50.0,
    ceiling: float = 2000.0,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
    criterion: Criterion = FACTOR,
    events: synaptic.EventsLike | None = None,
) -> float:
    """Finds the smallest EPSG-pair amplitude that spikes at this delay.

    Bisection helped by interpolating the peak axon - soma of the runs
    that don't spike (_bisect.smallest_crossing): ~10 runs instead of
    ~15, same guarantee.

    Args:
        model: "multi" or "two".
        delay: Time from the first EPSG onset to the second [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords, e.g. morph/mem, or r1/tau_est
            for "two".
        rel_tol: The relative tolerance of the threshold.
        guess: The first amplitude tried.
        ceiling: The largest amplitude worth trying.
        stim: "EPSGpair" (both EPSGs at input_node) or, for the
            multi-compartment model, "EPSGbilateral" (the second at
            input_node2).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".
        criterion: The spike criterion: axon - soma [mV], or a
            Relative fraction of spike_height().
        events: Drive with the "Events" stimulus instead of a pair: an
            event source (callable (delay, amp) -> tuple of
            synaptic.Event) or an event tuple whose peaks are scaled by
            the amplitude. The search then finds the smallest amp that
            spikes; `delay` goes to the source (a plain tuple has none:
            pass 0.0), and stim, input_node, input_node2 and epsg_tau are
            ignored. A Relative criterion's spike height is still
            measured with the EPSG pair.

    Returns:
        The threshold amplitude, to rel_tol, or inf if nothing up to
        `ceiling` spikes.

    Raises:
        ValueError: If rel_tol or guess is not > 0 (checked before any
            run, a Relative criterion's included), or if events is a
            plain tuple and delay is not 0.
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
    _bisect.check_positive("rel_tol", rel_tol)
    _bisect.check_positive("guess", guess)
    source = _source(events)
    _check_delay(source, delay)
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
    factor = resolve_criterion(
        criterion,
        model,
        node,
        v0,
        epsg_tau,
        model_kw,
        stim=stim,
        input_node=input_node,
        input_node2=input_node2,
    )
    try:
        return _bisect.smallest_crossing(
            lambda I: _respond(
                model,
                I,
                delay,
                node,
                v0,
                epsg_tau,
                site,
                model_kw,
                factor=factor,
                source=source,
            ),
            guess,
            ceiling,
            rel_tol,
            factor,
        )
    except _bisect.FiresWithoutInputError as e:
        config = _dispatch.describe(
            model,
            node,
            v0,
            model_kw.get("mem"),
            delay=delay,
            stim=stim,
            input_node=input_node,
            input_node2=input_node2,
            model_kw=sorted(model_kw),
            **({} if source is None else {"events": type(source).__name__}),
        )
        raise _bisect.FiresWithoutInputError(
            f"coincidence.threshold({config}): {e}"
        ) from e


def _threshold_task(task: tuple[tuple[Any, ...], dict[str, Any]]) -> float:
    args, kw = task
    return threshold(*args, **kw)


def threshold_curve(
    model: str,
    delays: Sequence[float] | np.ndarray,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
    criterion: Criterion = FACTOR,
    events: synaptic.EventsLike | None = None,
) -> np.ndarray:
    """Computes threshold() at each delay, in parallel.

    Scripts that call this in parallel need a __main__ guard.

    Args:
        model: "multi" or "two".
        delays: The delays [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        rel_tol: The relative tolerance of each threshold.
        workers: See _parallel.map_tasks.
        executor: See _parallel.map_tasks.
        stim: "EPSGpair" (both EPSGs at input_node) or, for the
            multi-compartment model, "EPSGbilateral" (the second at
            input_node2).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".
        criterion: The spike criterion (see threshold()), resolved here
            once rather than in every worker.
        events: An event source or tuple (see threshold()); the delays
            go to the source.

    Returns:
        An array of thresholds, one per delay.

    Raises:
        ValueError: If rel_tol is not > 0, or events is a plain tuple
            and a delay is not 0.
    """
    _bisect.check_positive("rel_tol", rel_tol)
    site = {"stim": stim, "input_node": input_node, "input_node2": input_node2}
    site["criterion"] = resolve_criterion(
        criterion, model, node, v0, epsg_tau, model_kw or {}, **site
    )
    if events is not None:
        site["events"] = events = synaptic.as_source(events)
        for d in delays:
            _check_delay(events, float(d))
    tasks = [
        (
            (
                model,
                float(d),
                node,
                v0,
                tuple(epsg_tau),
                model_kw or {},
                rel_tol,
            ),
            site,
        )
        for d in delays
    ]
    return np.array(
        _parallel.map_tasks(_threshold_task, tasks, workers, executor)
    )


def half_width(
    delays: Sequence[float] | np.ndarray,
    thresholds: Sequence[float] | np.ndarray,
    margin: float,
) -> float:
    """Returns the full width at half maximum of the spike-probability curve.

    Inputs are set `margin` (fraction) above the coincident threshold;
    probability is 50% where the threshold curve crosses that level. The
    EPSGs are identical and either share one input site or sit at
    mirror-image sites (EPSGbilateral on the two identical dendrites),
    so the curve is symmetric in the delay and the width is twice the
    crossing delay.

    Args:
        delays: The delays of the threshold curve, from 0 [ms].
        thresholds: The threshold at each delay.
        margin: How far above the coincident threshold the input is, as
            a fraction.

    Returns:
        The width [ms], or nan if the crossing lies beyond the delays
        given.
    """
    delays, thresholds = (
        np.asarray(delays, float),
        np.asarray(thresholds, float),
    )
    level = (1 + margin) * thresholds[0]
    above = np.nonzero(thresholds > level)[0]
    if len(above) == 0:
        return np.nan
    j = above[0]
    d = np.interp(level, thresholds[j - 1 : j + 1], delays[j - 1 : j + 1])
    return 2 * d


def window(
    model: str,
    margin: float,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    delay_tol: float = 1e-4,
    max_delay: float = 1.0,
    threshold0: float | None = None,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
    criterion: Criterion = FACTOR,
    events: synaptic.EventsLike | None = None,
) -> tuple[float, float]:
    """Finds the coincidence-window width at `margin` directly.

    Myoga et al held the input size fixed and varied the delay, and this
    does the same: find the coincident threshold th0 (one bisection,
    skipped if threshold0 is given), fix the amplitude at
    (1 + margin) * th0, and bisect over delay, to delay_tol, for where
    that input stops spiking. The width is twice that delay.

    It measures the same thing as
    half_width(delays, threshold_curve(...), margin), but needs ~15 +
    ~14 runs instead of ~15 per delay of the curve, and has no linear
    interpolation between grid points. It assumes a single boundary:
    once the threshold curve rises above (1 + margin) * th0 it stays
    above. Every saved curve does (their only dips are <0.1% of th0, on
    the plateau at about twice th0).

    Args:
        model: "multi" or "two".
        margin: How far above th0 the input is, as a fraction.
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        rel_tol: The relative tolerance of th0.
        delay_tol: The tolerance of the boundary delay [ms].
        max_delay: The largest delay tried [ms].
        threshold0: th0 if already known, to skip its search.
        stim: "EPSGpair" (both EPSGs at input_node) or, for the
            multi-compartment model, "EPSGbilateral" (the second at
            input_node2).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".
        criterion: The spike criterion (see threshold()). A given
            threshold0 must have been found with the same one.
        events: An event source (see threshold()) to measure instead of
            the pair, or None. With events the window is measured by
            window_edges (the delay is signed, and its reference is the
            lowest threshold over delay, th0 below), and threshold0 may
            not be given.

    Returns:
        (width, th0): the width [ms], nan if the input still spikes at
        max_delay or if th0 is not finite (nothing fires up to the
        ceiling, or a given threshold0 is inf or nan), and the
        coincident threshold (with events: the lowest threshold over
        delay).

    Raises:
        ValueError: If delay_tol is not > 0, or rel_tol is not > 0 when
            th0 is searched for (checked before any run), or if
            threshold0 is given with events.
    """
    _bisect.check_positive("delay_tol", delay_tol)
    source = _source(events)
    if source is not None:
        if threshold0 is not None:
            raise ValueError(
                "threshold0 does not apply to an event source: the "
                "reference threshold is the lowest over delay"
            )
        edges = window_edges(
            model,
            margin,
            node,
            v0,
            epsg_tau,
            model_kw,
            rel_tol,
            delay_tol,
            max_delay,
            events=source,
            criterion=criterion,
        )
        return edges.width, edges.th0
    if threshold0 is None:
        _bisect.check_positive("rel_tol", rel_tol)
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
    factor = resolve_criterion(
        criterion,
        model,
        node,
        v0,
        epsg_tau,
        model_kw,
        stim=stim,
        input_node=input_node,
        input_node2=input_node2,
    )
    th0 = (
        threshold0
        if threshold0 is not None
        else threshold(
            model,
            0.0,
            node,
            v0,
            epsg_tau,
            model_kw,
            rel_tol,
            stim=stim,
            input_node=input_node,
            input_node2=input_node2,
            criterion=factor,
        )
    )
    if not np.isfinite(th0):
        return np.nan, th0  # no threshold, so no window
    amp = (1 + margin) * th0

    def respond(delay: float) -> tuple[bool, float]:
        return _respond(
            model, amp, delay, node, v0, epsg_tau, site, model_kw, factor=factor
        )

    spiked, peak = respond(max_delay)
    if spiked:
        return np.nan, th0
    # amp > th0, so it spikes at delay 0; the peak falls smoothly through the
    # criterion as the delay grows, so the same interpolation as in threshold()
    # applies
    far, near = _bisect.crossing(
        respond,
        max_delay,
        0.0,
        lambda _no, _yes: delay_tol,
        factor,
        [(max_delay, peak)],
    )
    return far + near, th0  # twice the midpoint of the final bracket


class WindowEdges(NamedTuple):
    """The two edges of a coincidence window on a signed delay.

    Attributes:
        lo: The delay below which the input stops spiking [ms].
        hi: The delay above which it stops spiking [ms].
        best: The delay of the lowest threshold [ms].
        th0: The lowest threshold, the reference of the margin.
    """

    lo: float
    hi: float
    best: float
    th0: float

    @property
    def width(self) -> float:
        """The full width, hi - lo [ms]."""
        return self.hi - self.lo

    @property
    def centre(self) -> float:
        """The window's midpoint, (hi + lo) / 2 [ms]."""
        return (self.hi + self.lo) / 2


def best_delay(
    model: str,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    span: float = 0.5,
    delay_tol: float = 2e-3,
    *,
    events: synaptic.EventsLike | None = None,
    criterion: Criterion = FACTOR,
) -> tuple[float, float]:
    """Finds the delay where an event source's threshold is lowest.

    A bounded scalar minimisation of threshold() over delay (Brent's
    method, scipy), to delay_tol. The curve is smooth and has one
    minimum over the +-0.5 ms where the inputs overlap.

    Args:
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        rel_tol: The relative tolerance of each threshold.
        span: The delays searched are -span to +span [ms].
        delay_tol: The tolerance of the delay [ms].
        events: The event source (see threshold()).
        criterion: The spike criterion (see threshold()).

    Returns:
        (delay, threshold) at the minimum; (nan, inf) if the threshold
        at delay 0 is not finite.
    """
    factor = resolve_criterion(criterion, model, node, v0, epsg_tau, model_kw)

    def th(d: float) -> float:
        return threshold(
            model,
            float(d),
            node,
            v0,
            epsg_tau,
            model_kw,
            rel_tol,
            criterion=factor,
            events=events,
        )

    if not np.isfinite(th(0.0)):
        return np.nan, np.inf
    res = optimize.minimize_scalar(
        th,
        bounds=(-span, span),
        method="bounded",
        options={"xatol": delay_tol},
    )
    return float(res.x), float(res.fun)


def window_edges(
    model: str,
    margin: float,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    delay_tol: float = 1e-4,
    max_delay: float = 1.0,
    reference: tuple[float, float] | None = None,
    *,
    events: synaptic.EventsLike,
    criterion: Criterion = FACTOR,
) -> WindowEdges:
    """Finds both edges of the coincidence window of an event source.

    window() fixes the amplitude at (1 + margin) times the coincident
    threshold and bisects over delay for where the input stops spiking,
    counting on the identical EPSGs making the curve symmetric about
    delay 0. An inhibition timed against one EPSG breaks that, so here the
    delay is signed (the second EPSG at start + delay, see PairSource), the
    reference is the lowest threshold over delay (best_delay), the
    amplitude is (1 + margin) times it, and each edge is bisected outwards
    from that delay. The window's width is hi - lo and its centre, the
    optimum Myoga et al 2014 report shifting, their midpoint. It assumes a
    single window: the input spikes between the edges and not beyond.

    Args:
        model: "multi" or "two".
        margin: How far above the reference threshold the input is, as
            a fraction.
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        rel_tol: The relative tolerance of the thresholds.
        delay_tol: The tolerance of each edge [ms].
        max_delay: The farthest each edge is looked for from the best
            delay [ms].
        reference: (best delay, threshold) if already known, to skip
            best_delay.
        events: The event source (see threshold()), with a signed delay.
        criterion: The spike criterion (see threshold()).

    Returns:
        The WindowEdges; its edges are nan where the input still spikes
        at max_delay from the best delay, both are nan if it does not
        spike at the best delay (a reference threshold too low), and
        all of lo, hi and best are nan if nothing fires.

    Raises:
        ValueError: If delay_tol or rel_tol is not > 0.
    """
    _bisect.check_positive("delay_tol", delay_tol)
    _bisect.check_positive("rel_tol", rel_tol)
    source = synaptic.as_source(events)
    model_kw = model_kw or {}
    factor = resolve_criterion(criterion, model, node, v0, epsg_tau, model_kw)
    best, th0 = reference or best_delay(
        model,
        node,
        v0,
        epsg_tau,
        model_kw,
        rel_tol,
        events=source,
        criterion=factor,
    )
    nan = float("nan")
    if not np.isfinite(th0):
        return WindowEdges(nan, nan, nan, th0)
    amp = (1 + margin) * th0
    site = ("EPSGpair", 1, None)

    def run(delay: float) -> tuple[bool, float]:
        return _respond(
            model,
            amp,
            delay,
            node,
            v0,
            epsg_tau,
            site,
            model_kw,
            factor=factor,
            source=source,
        )

    if not run(best)[0]:  # a given reference may be too low to spike at all
        return WindowEdges(nan, nan, best, th0)

    def edge(sign: float) -> float:
        def respond(away: float) -> tuple[bool, float]:
            return run(best + sign * away)

        spiked, peak = respond(max_delay)
        if spiked:
            return nan
        far, near = _bisect.crossing(
            respond,
            max_delay,
            0.0,
            lambda _no, _yes: delay_tol,
            factor,
            [(max_delay, peak)],
        )
        return best + sign * (far + near) / 2

    return WindowEdges(edge(-1.0), edge(1.0), best, th0)


def trial_spikes(
    model: str,
    event_lists: Sequence[Sequence[synaptic.Event]],
    node: int = 3,
    v0: float = -68.0,
    model_kw: dict[str, Any] | None = None,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    *,
    criterion: Criterion = FACTOR,
) -> np.ndarray:
    """Runs one trial per event list, in parallel, and says which spiked.

    The events are run as given (I = 1): the caller draws each trial's
    onsets, peaks and compartments (see synaptic.fibre_events and
    bilateral_events). Scripts that call this in parallel need a __main__
    guard.

    Args:
        model: "multi" or "two".
        event_lists: The trials' event tuples.
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        model_kw: Further model keywords (see threshold()).
        workers: See _parallel.map_tasks.
        executor: See _parallel.map_tasks.
        criterion: The spike criterion (see threshold()), resolved here
            once rather than in every worker.

    Returns:
        A boolean array, True where the trial spiked.

    Raises:
        ValueError: If an event list is empty or an event is out of
            range (see synaptic.check_events), before it runs.
    """
    model_kw = model_kw or {}
    factor = resolve_criterion(criterion, model, node, v0, model_kw=model_kw)
    tasks = [
        (model, tuple(ev), node, v0, model_kw, factor) for ev in event_lists
    ]
    hits = _parallel.map_tasks(
        _event_trial, tasks, workers, executor, chunksize=16
    )
    return np.array(hits, bool)


def _event_trial(args: tuple[Any, ...]) -> bool:
    model, events, node, v0, model_kw, factor = args
    return _dispatch.spikes(
        model,
        "Events",
        0.0,
        0.0,
        1.0,
        node,
        T_END,
        v0,
        1,
        synaptic.SynParams(t_end=T_END, events=events),
        factor=factor,
        **model_kw,
    )


def _trial(args: tuple[Any, ...]) -> bool:
    model, amp, t1, t2, node, v0, epsg_tau, site, model_kw, factor = args
    first, second = min(t1, t2), max(t1, t2)
    return _spikes(
        model,
        amp,
        second - first,
        node,
        v0,
        epsg_tau,
        site,
        model_kw,
        start=first,
        factor=factor,
    )


def probability_trials(
    model: str,
    delays: Sequence[float] | np.ndarray,
    amplitude: float,
    n_trials: int = 200,
    amp_cv: float = 0.03,
    jitter: float = 0.015,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    seed: int = 0,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
    criterion: Criterion = FACTOR,
    events: synaptic.EventsLike | None = None,
) -> np.ndarray:
    """Measures spike probability per delay from noisy trials.

    This checks half_width's shortcut. Noise per trial: the pair's
    amplitude scaled by N(1, amp_cv), and each EPSG's onset jittered by
    N(0, jitter) ms independently. For EPSGbilateral the earlier EPSG is
    always at input_node.

    With events the trial is the event source's events at that delay and
    amplitude, with the same noise: every onset moves by N(0, jitter) ms.
    The delays are signed, and onsets are kept at or after 0.

    Args:
        model: "multi" or "two".
        delays: The delays [ms].
        amplitude: The mean EPSG-pair amplitude.
        n_trials: Trials per delay.
        amp_cv: Coefficient of variation of the amplitude.
        jitter: Standard deviation of each onset [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        seed: Seed of the noise generator.
        workers: See _parallel.map_tasks.
        executor: See _parallel.map_tasks.
        stim: "EPSGpair" (both EPSGs at input_node) or, for the
            multi-compartment model, "EPSGbilateral" (the second at
            input_node2).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".
        criterion: The spike criterion (see threshold()), resolved here
            once rather than in every worker.
        events: An event source or tuple (see threshold()), or None.

    Returns:
        An array with the fraction of trials that spiked, per delay.
    """
    source = _source(events)
    if source is not None:
        return _probability_events(
            model,
            delays,
            source,
            amplitude,
            n_trials,
            amp_cv,
            jitter,
            node,
            v0,
            model_kw,
            seed,
            workers,
            executor,
            criterion,
        )
    site = (stim, input_node, input_node2)
    factor = resolve_criterion(
        criterion,
        model,
        node,
        v0,
        epsg_tau,
        model_kw or {},
        stim=stim,
        input_node=input_node,
        input_node2=input_node2,
    )
    rng = np.random.default_rng(seed)
    tasks = []
    for d in delays:
        for _ in range(n_trials):
            amp = amplitude * rng.normal(1, amp_cv)
            t1 = START + rng.normal(0, jitter)
            t2 = START + d + rng.normal(0, jitter)
            tasks.append(
                (
                    model,
                    amp,
                    t1,
                    t2,
                    node,
                    v0,
                    tuple(epsg_tau),
                    site,
                    model_kw or {},
                    factor,
                )
            )
    hits = _parallel.map_tasks(_trial, tasks, workers, executor, chunksize=16)
    return np.array(hits, float).reshape(len(delays), n_trials).mean(axis=1)


def _probability_events(
    model: str,
    delays: Sequence[float] | np.ndarray,
    source: synaptic.Source,
    amplitude: float,
    n_trials: int,
    amp_cv: float,
    jitter: float,
    node: int,
    v0: float,
    model_kw: dict[str, Any] | None,
    seed: int,
    workers: int | None,
    executor: futures.Executor | None,
    criterion: Criterion,
) -> np.ndarray:
    """Does probability_trials for an event source."""
    rng = np.random.default_rng(seed)
    trials = []
    for d in delays:
        _check_delay(source, float(d))
        for _ in range(n_trials):
            amp = amplitude * rng.normal(1, amp_cv)
            trials.append(
                tuple(
                    dataclasses.replace(
                        e,
                        t=max(e.t + rng.normal(0, jitter), 0.0),
                    )
                    for e in source(float(d), amp)
                )
            )
    hits = trial_spikes(
        model,
        trials,
        node,
        v0,
        model_kw,
        workers,
        executor,
        criterion=criterion,
    )
    return hits.reshape(len(delays), n_trials).mean(axis=1)


# --- the spike height, for a relative criterion -------------------------------

# spike_height() results, keyed on _height_key(); failures are not stored
_HEIGHTS: dict[tuple[Any, ...], float] = {}


def _frozen(name: str, value: object) -> object:
    """Returns a hashable stand-in for one model keyword's value.

    A morphology or membrane (a dataclass whose ndarray fields are
    unhashable) stands in as its `key`, which names every setting it was built
    from, plus a membrane's morph_key.

    Raises:
        TypeError: If the value is neither keyed nor hashable.
    """
    key = getattr(value, "key", None)
    if key is not None:
        return ("key", key, getattr(value, "morph_key", None))
    try:
        hash(value)
    except TypeError:
        raise TypeError(
            f"spike_height can't cache on model keyword {name}={value!r}: "
            "it is neither hashable nor a keyed morphology or membrane"
        ) from None
    return value


def _height_key(
    model: str,
    node: int,
    v0: float,
    epsg_tau: tuple[float, float],
    model_kw: dict[str, Any],
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
) -> tuple[Any, ...]:
    """Names everything that shapes spike_height()'s runs, hashably."""
    if model == "multi":  # the mode the runs will use, not the one asked for
        model_kw = {**model_kw, "mode": _fast.effective(model_kw.get("mode"))}
    return (
        model,
        int(node),
        float(v0),
        tuple(float(x) for x in epsg_tau),
        stim,
        int(input_node),
        input_node2,
        tuple(sorted((k, _frozen(k, v)) for k, v in model_kw.items())),
    )


def spike_height(
    model: str,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
) -> float:
    """Measures the model's spike height, for a Relative criterion.

    The height is the largest axon - soma of one coincident EPSG pair at
    HEIGHT_AT (2x) the 10 mV threshold, run to T_END with no spike stop.
    It is flat there (within 2.5% from 1.5x to 3x threshold), so the
    10 mV threshold that seeds it doesn't bias it: 38.4 mV for the
    45-compartment model, 38.0 for TwoCpt.m, 20.2 with Goldwyn's
    calibration (model EPSG). It costs one threshold search plus one run,
    and is cached per configuration: the model, node, v0, EPSG kinetics,
    stimulus sites and every model keyword (a morphology or membrane by
    its key, not by identity), since variants such as dendrites or a KLT
    block change the spike.

    Args:
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        stim: "EPSGpair" or "EPSGbilateral" (see threshold()).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".

    Returns:
        The peak axon - soma [mV].

    Raises:
        _bisect.FiresWithoutInputError: If the cell fires with no input.
        ValueError: If no EPSG pair up to threshold()'s ceiling spikes,
            so there is no spike to measure.
        TypeError: If a model keyword can't be part of the cache key.
    """
    model_kw = model_kw or {}
    site = {"stim": stim, "input_node": input_node, "input_node2": input_node2}
    key = _height_key(model, node, v0, epsg_tau, model_kw, **site)
    if key in _HEIGHTS:
        return _HEIGHTS[key]
    try:
        th = threshold(
            model,
            0.0,
            node,
            v0,
            epsg_tau,
            model_kw,
            HEIGHT_REL_TOL,
            criterion=FACTOR,
            **site,
        )
    except _bisect.FiresWithoutInputError as e:
        raise _bisect.FiresWithoutInputError(
            f"coincidence.spike_height: {e}"
        ) from e
    if not np.isfinite(th):
        config = _dispatch.describe(
            model, node, v0, model_kw.get("mem"), **site
        )
        raise ValueError(
            f"coincidence.spike_height({config}): no EPSG pair up to the "
            "ceiling spikes, so there is no spike height"
        )
    syn = synaptic.SynParams(t_end=T_END, epsg_tau=tuple(epsg_tau))
    kw = model_kw
    if input_node2 is not None:
        kw = {**model_kw, "input_node2": input_node2}
    t, x = _dispatch.run_model(
        model,
        stim,
        START,
        START,
        HEIGHT_AT * th,
        node,
        T_END,
        v0,
        input_node,
        syn,
        **kw,
    )
    axon = _dispatch.axon_column(model, node)
    _HEIGHTS[key] = float((x[:, axon] - x[:, 0]).max())
    return _HEIGHTS[key]


def resolve_criterion(
    criterion: Criterion,
    model: str,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    *,
    stim: str = "EPSGpair",
    input_node: int = 1,
    input_node2: int | None = None,
) -> float:
    """Turns a spike criterion into mV for this configuration.

    Args:
        criterion: axon - soma [mV], or a Relative fraction of
            spike_height().
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSG (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        stim: "EPSGpair" or "EPSGbilateral" (see threshold()).
        input_node: Where the (first) EPSG enters, 1-indexed.
        input_node2: Where the second EPSG enters, for "EPSGbilateral".

    Returns:
        The criterion [mV]: a float unchanged, a Relative times the
        spike height.

    Raises:
        _bisect.FiresWithoutInputError: See spike_height().
        ValueError: See spike_height().
    """
    if not isinstance(criterion, Relative):
        return criterion
    return criterion.fraction * spike_height(
        model,
        node,
        v0,
        epsg_tau=epsg_tau,
        model_kw=model_kw,
        stim=stim,
        input_node=input_node,
        input_node2=input_node2,
    )


# --- protocols on event lists (issue #78) -------------------------------------

# Casarez et al 2026: about 125 synapses per side, release probability 0.45.
# Their unitary conductance is not in the notes this was built from;
# CASAREZ_G is this repo's choice: the smallest round peak at which the
# default sweep is tuned in two_cpt (1.1 nS: probability 1.0 at ITD 0, 0 at
# +-0.5 ms; 0.8 nS never fires, 1.4 nS fires over +-0.3 ms). It is a
# calibration to this model, not a measurement.
CASAREZ_N = 125
CASAREZ_P = 0.45
CASAREZ_G = 1.1  # [nS]


def bilateral_probability(
    model: str,
    delays: Sequence[float] | np.ndarray,
    n_trials: int = 50,
    repeats: int = 1,
    node: int = 3,
    v0: float = -68.0,
    model_kw: dict[str, Any] | None = None,
    seed: int = 0,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    *,
    criterion: Criterion = FACTOR,
    **drive: Any,
) -> np.ndarray:
    """Measures spike probability per delay for two sides of random fibres.

    Each trial draws synaptic.bilateral_events afresh (jittered onsets,
    random release and peaks) with the second side `delay` after the
    first, and runs it. The defaults are the Couchman et al 2010 input:
    4 fibres of 37 nS per side with a 0.32 ms Gaussian jitter (the
    excitatory jitter Myoga et al 2014 used for their probability
    curves). One seeded generator draws everything, so a seed reproduces
    the result whatever the worker count.

    Args:
        model: "multi" or "two".
        delays: The delays of the second side [ms], signed.
        n_trials: Trials per delay and repeat.
        repeats: Independent repeats of the whole sweep, each a full set
            of n_trials per delay: the spread over repeats is the
            distribution of a width measured from n_trials trials.
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        model_kw: Further model keywords (see threshold()).
        seed: Seed of the generator.
        workers: See _parallel.map_tasks.
        executor: See _parallel.map_tasks.
        criterion: The spike criterion (see threshold()).
        **drive: Further synaptic.bilateral_events arguments: n, g,
            g_sd, jitter, release_p, nodes, tau, rate, duration. The
            defaults here are n=4, g=37, jitter=0.32.

    Returns:
        An array of spike probabilities, shape (repeats, len(delays)).
    """
    drive = {"n": 4, "g": 37.0, "jitter": 0.32, **drive}
    rng = np.random.default_rng(seed)
    trials = [
        synaptic.bilateral_events(float(d), rng, **drive)
        for _ in range(repeats)
        for d in delays
        for _ in range(n_trials)
    ]
    hits = trial_spikes(
        model,
        trials,
        node,
        v0,
        model_kw,
        workers,
        executor,
        criterion=criterion,
    )
    return hits.reshape(repeats, len(delays), n_trials).mean(axis=2)


def fibre_window_trials(
    model: str,
    delays: Sequence[float] | np.ndarray,
    n_fibres: int = 2,
    jitter: float = 0.32,
    n_trials: int = 50,
    repeats: int = 5,
    **kw: Any,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    """Measures the half-width of the fibre-drive window over repeats.

    The Myoga et al 2014 probability curve, built from random fibres
    (bilateral_probability) instead of a threshold search: each repeat
    gives a spike-probability curve over delay and measure.curve_summary
    its half-width, so the repeats give the distribution of the width one
    n_trials-trial experiment would report.

    Args:
        model: "multi" or "two".
        delays: The delays of the second side [ms], signed and increasing.
        n_fibres: Fibres per side, 37 nS each unless g is given. The
            default 2 makes 4 fibres in all, Couchman et al 2010's
            lower bound; 3 per side already fire the cell alone.
        jitter: The standard deviation of each fibre's onset [ms].
        n_trials: Trials per delay and repeat.
        repeats: Independent repeats.
        **kw: Further bilateral_probability arguments (g, g_sd, seed,
            workers, executor, node, v0, model_kw, criterion, ...).

    Returns:
        (probability, summaries): probability of shape (repeats,
        len(delays)), and measure.curve_summary of each repeat's curve
        against the delays.
    """
    prob = bilateral_probability(
        model,
        delays,
        n_trials,
        repeats,
        n=n_fibres,
        jitter=jitter,
        **kw,
    )
    return prob, [measure.curve_summary(delays, row) for row in prob]


def itd_sweep(
    model: str,
    itds: Sequence[float] | np.ndarray | None = None,
    n_trials: int = 20,
    n_syn: int = CASAREZ_N,
    release_p: float = CASAREZ_P,
    g: float = CASAREZ_G,
    jitter: float = 0.1,
    **kw: Any,
) -> tuple[np.ndarray, dict[str, float]]:
    """Runs a Casarez-style trial-based ITD sweep.

    Casarez et al 2026: about 125 synapses per side, release probability
    0.45, 400 trials per ITD, ITDs of +-0.5 ms in 10 us steps. The
    defaults are a small version (11 ITDs, 20 trials) so that it runs in
    a minute; pass itds and n_trials for the real one. The sides enter
    nodes[0] and nodes[1] (both the soma by default: a model without
    dendrites has one input site, so the ITD tuning is then coincidence
    alone). Each trial is the second side's synapses `itd` after the
    first's, each synapse jittered, released with probability
    release_p, and peaking at g.

    Args:
        model: "multi" or "two".
        itds: The ITDs [ms], increasing; -0.5 to 0.5 in 0.1 steps by
            default.
        n_trials: Trials per ITD.
        n_syn: Synapses per side.
        release_p: Release probability of each.
        g: The unitary peak [nS].
        jitter: The standard deviation of each onset [ms].
        **kw: Further bilateral_probability arguments (seed, workers,
            executor, node, v0, model_kw, criterion, nodes, g_sd, ...).

    Returns:
        (probability, summary): the spike probability at each ITD and
        measure.curve_summary of the curve (its peak, half-width,
        centre and centroid).
    """
    if itds is None:
        itds = np.round(np.linspace(-0.5, 0.5, 11), 4)
    prob = bilateral_probability(
        model,
        itds,
        n_trials,
        n=n_syn,
        release_p=release_p,
        g=g,
        jitter=jitter,
        **kw,
    )[0]
    return prob, measure.curve_summary(itds, prob)


def prepulse_thresholds(
    model: str,
    pre_dts: Sequence[float] = (-1.0, -0.75, -0.5),
    g_pre: float = 37.0,
    delay: float = 0.0,
    node: int = 3,
    v0: float = -68.0,
    epsg_tau: tuple[float, float] = synaptic.EPSG_TAU,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = 1e-4,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    *,
    criterion: Criterion = FACTOR,
) -> tuple[float, np.ndarray]:
    """Finds the EPSG-pair threshold after a preceding EPSG.

    Franken et al 2015: an EPSP at 1.00 / 0.75 ms before the test raises
    the EPSG that fires the cell half the time by 7.2 / 6.6 nS, because
    the pre-pulse activates Kv1. The pre-pulse is one EPSG of peak g_pre
    at each of pre_dts relative to the pair's first EPSG, with the pair's
    kinetics, at the pair's node; the search scales the pair only.

    Args:
        model: "multi" or "two".
        pre_dts: The pre-pulse onsets relative to the pair [ms]
            (negative: before).
        g_pre: The pre-pulse peak [nS]; the default is one Couchman et
            al 2010 fibre.
        delay: The pair's delay [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        epsg_tau: The EPSGs' (rise, decay) time constants [ms].
        model_kw: Further model keywords (see threshold()).
        rel_tol: The relative tolerance of each threshold.
        workers: See _parallel.map_tasks.
        executor: See _parallel.map_tasks.
        criterion: The spike criterion (see threshold()), resolved here
            once rather than in every worker.

    Returns:
        (th0, thresholds): the pair's threshold amplitude with no
        pre-pulse, and with each pre-pulse (so the rise, in nS, is
        thresholds - th0).
    """
    factor = resolve_criterion(criterion, model, node, v0, epsg_tau, model_kw)
    sources = [
        synaptic.PairSource(tau=tuple(epsg_tau), pre=((dt, g_pre),))
        for dt in pre_dts
    ]
    sources.insert(0, synaptic.PairSource(tau=tuple(epsg_tau)))
    tasks = [
        (
            (model, delay, node, v0, tuple(epsg_tau), model_kw or {}, rel_tol),
            {"events": src, "criterion": factor},
        )
        for src in sources
    ]
    th = _parallel.map_tasks(_threshold_task, tasks, workers, executor)
    return th[0], np.array(th[1:])
