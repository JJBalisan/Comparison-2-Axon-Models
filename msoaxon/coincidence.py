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
"""

import dataclasses
from collections.abc import Sequence
from concurrent import futures
from typing import Any

import numpy as np

from msoaxon import _bisect, _common, _dispatch, _fast, _parallel, synaptic

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
) -> tuple[bool, float]:
    """Runs one EPSG pair, `delay` apart.

    Returns (spiked, peak axon - soma); site = (stim, input_node,
    input_node2); factor is the resolved spike criterion [mV].
    """
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
    )[0]


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

    Returns:
        The threshold amplitude, to rel_tol, or inf if nothing up to
        `ceiling` spikes.

    Raises:
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
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

    Returns:
        An array of thresholds, one per delay.
    """
    site = {"stim": stim, "input_node": input_node, "input_node2": input_node2}
    site["criterion"] = resolve_criterion(
        criterion, model, node, v0, epsg_tau, model_kw or {}, **site
    )
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

    Returns:
        (width, th0): the width [ms], nan if the input still spikes at
        max_delay or if th0 is not finite (nothing fires up to the
        ceiling, or a given threshold0 is inf or nan), and the
        coincident threshold.
    """
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
) -> np.ndarray:
    """Measures spike probability per delay from noisy trials.

    This checks half_width's shortcut. Noise per trial: the pair's
    amplitude scaled by N(1, amp_cv), and each EPSG's onset jittered by
    N(0, jitter) ms independently. For EPSGbilateral the earlier EPSG is
    always at input_node.

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

    Returns:
        An array with the fraction of trials that spiked, per delay.
    """
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
    axon = node - 1 if model == "multi" else 1
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
