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

This is new analysis code, not part of the MATLAB port: it uses a standard
bracketing bisection rather than BinarySearch.m's halving search.
"""

import numpy as np

from msoaxon import _bisect, _dispatch, _parallel, synaptic

START = 5.0
T_END = 20.0
FACTOR = (
    10.0  # spike = axon rises this far above soma (as in the EPSGpair sweeps)
)


def _respond(model, I, delay, node, v0, epsg_tau, site, model_kw, start=START):
    """Runs one EPSG pair, `delay` apart.

    Returns (spiked, peak axon - soma); site = (stim, input_node,
    input_node2).
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
        factor=FACTOR,
        **model_kw,
    )


def _spikes(model, I, delay, node, v0, epsg_tau, site, model_kw, start=START):
    return _respond(model, I, delay, node, v0, epsg_tau, site, model_kw, start)[
        0
    ]


def threshold(
    model,
    delay,
    node=3,
    v0=-68.0,
    epsg_tau=synaptic.EPSG_TAU,
    model_kw=None,
    rel_tol=1e-4,
    guess=50.0,
    ceiling=2000.0,
    *,
    stim="EPSGpair",
    input_node=1,
    input_node2=None,
):
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

    Returns:
        The threshold amplitude, to rel_tol, or inf if nothing up to
        `ceiling` spikes.

    Raises:
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
    try:
        return _bisect.smallest_crossing(
            lambda I: _respond(
                model, I, delay, node, v0, epsg_tau, site, model_kw
            ),
            guess,
            ceiling,
            rel_tol,
            FACTOR,
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


def _threshold_task(task):
    args, kw = task
    return threshold(*args, **kw)


def threshold_curve(
    model,
    delays,
    node=3,
    v0=-68.0,
    epsg_tau=synaptic.EPSG_TAU,
    model_kw=None,
    rel_tol=1e-4,
    workers=None,
    executor=None,
    *,
    stim="EPSGpair",
    input_node=1,
    input_node2=None,
):
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

    Returns:
        An array of thresholds, one per delay.
    """
    site = dict(stim=stim, input_node=input_node, input_node2=input_node2)
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


def half_width(delays, thresholds, margin):
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
    model,
    margin,
    node=3,
    v0=-68.0,
    epsg_tau=synaptic.EPSG_TAU,
    model_kw=None,
    rel_tol=1e-4,
    delay_tol=1e-4,
    max_delay=1.0,
    threshold0=None,
    *,
    stim="EPSGpair",
    input_node=1,
    input_node2=None,
):
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

    Returns:
        (width, th0): the width [ms], nan if the input still spikes at
        max_delay, and the coincident threshold.
    """
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
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
        )
    )
    amp = (1 + margin) * th0

    def respond(delay):
        return _respond(model, amp, delay, node, v0, epsg_tau, site, model_kw)

    spiked, peak = respond(max_delay)
    if spiked:
        return np.nan, th0
    # amp > th0, so it spikes at delay 0; the peak falls smoothly through FACTOR
    # as the delay grows, so the same interpolation as in threshold() applies
    far, near = _bisect.crossing(
        respond,
        max_delay,
        0.0,
        lambda n, y: delay_tol,
        FACTOR,
        [(max_delay, peak)],
    )
    return far + near, th0  # twice the midpoint of the final bracket


def _trial(args):
    model, amp, t1, t2, node, v0, epsg_tau, site, model_kw = args
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
    )


def probability_trials(
    model,
    delays,
    amplitude,
    n_trials=200,
    amp_cv=0.03,
    jitter=0.015,
    node=3,
    v0=-68.0,
    epsg_tau=synaptic.EPSG_TAU,
    model_kw=None,
    seed=0,
    workers=None,
    executor=None,
    *,
    stim="EPSGpair",
    input_node=1,
    input_node2=None,
):
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

    Returns:
        An array with the fraction of trials that spiked, per delay.
    """
    site = (stim, input_node, input_node2)
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
                )
            )
    hits = _parallel.map_tasks(_trial, tasks, workers, executor, chunksize=16)
    return np.array(hits, float).reshape(len(delays), n_trials).mean(axis=1)
