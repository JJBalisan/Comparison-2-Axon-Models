"""Coincidence-detection windows comparable to Myoga et al 2014 (Nat Commun 5:3790).

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

from ._bisect import smallest_firing
from ._parallel import map_tasks
from ._dispatch import spikes
from .synaptic import EPSG_TAU, SynParams

START = 5.0
T_END = 20.0
FACTOR = 10.0  # spike = axon rises this far above soma (as in the EPSGpair sweeps)


def _spikes(model, I, delay, node, v0, epsg_tau, site, model_kw, start=START):
    """One EPSG pair, `delay` apart; site = (stim, input_node, input_node2)."""
    stim, input_node, input_node2 = site
    syn = SynParams(t_end=T_END, epsg_tau=tuple(epsg_tau))
    if input_node2 is not None:
        model_kw = {**model_kw, "input_node2": input_node2}
    return spikes(model, stim, start, start + delay, I, node, T_END, v0, input_node, syn,
                  factor=FACTOR, **model_kw)


def threshold(model, delay, node=3, v0=-68.0, epsg_tau=EPSG_TAU, model_kw=None,
              rel_tol=1e-4, guess=50.0, ceiling=2000.0, *, stim="EPSGpair", input_node=1,
              input_node2=None):
    """Smallest EPSG-pair amplitude that spikes at this delay, to rel_tol (bisection).

    stim: "EPSGpair" (both EPSGs at input_node) or, for the multi-compartment
    model, "EPSGbilateral" (the second at input_node2).
    model_kw: further model keywords, e.g. morph/mem, or r1/tau_est for "two".
    Returns inf if nothing up to `ceiling` spikes (see _bisect.smallest_firing).
    """
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
    return smallest_firing(lambda I: _spikes(model, I, delay, node, v0, epsg_tau, site, model_kw),
                           guess, ceiling, rel_tol)


def _threshold_task(task):
    args, kw = task
    return threshold(*args, **kw)


def threshold_curve(model, delays, node=3, v0=-68.0, epsg_tau=EPSG_TAU, model_kw=None,
                    rel_tol=1e-4, workers=None, executor=None, *, stim="EPSGpair",
                    input_node=1, input_node2=None):
    """threshold() at each delay, in parallel (needs a __main__ guard in scripts).

    workers / executor: see _parallel.map_tasks.
    """
    site = dict(stim=stim, input_node=input_node, input_node2=input_node2)
    tasks = [((model, float(d), node, v0, tuple(epsg_tau), model_kw or {}, rel_tol), site)
             for d in delays]
    return np.array(map_tasks(_threshold_task, tasks, workers, executor))


def half_width(delays, thresholds, margin):
    """Full width at half maximum of the spike-probability curve.

    Inputs are set `margin` (fraction) above the coincident threshold; probability
    is 50% where the threshold curve crosses that level. The EPSGs are identical
    and either share one input site or sit at mirror-image sites (EPSGbilateral
    on the two identical dendrites), so the curve is symmetric in the delay and
    the width is twice the crossing delay. Returns nan if the crossing lies beyond
    the delays given.
    """
    delays, thresholds = np.asarray(delays, float), np.asarray(thresholds, float)
    level = (1 + margin) * thresholds[0]
    above = np.nonzero(thresholds > level)[0]
    if len(above) == 0:
        return np.nan
    j = above[0]
    d = np.interp(level, thresholds[j - 1:j + 1], delays[j - 1:j + 1])
    return 2 * d


def window(model, margin, node=3, v0=-68.0, epsg_tau=EPSG_TAU, model_kw=None, rel_tol=1e-4,
           delay_tol=1e-4, max_delay=1.0, threshold0=None, *, stim="EPSGpair", input_node=1,
           input_node2=None):
    """Coincidence-window width [ms] at `margin`, found directly. Returns (width, th0).

    Myoga et al held the input size fixed and varied the delay, and this does the
    same: find the coincident threshold th0 (one bisection, skipped if threshold0 is
    given), fix the amplitude at (1 + margin) * th0, and bisect over delay, to
    delay_tol, for where that input stops spiking. The width is twice that delay.

    It measures the same thing as half_width(delays, threshold_curve(...), margin),
    but needs ~15 + ~14 runs instead of ~15 per delay of the curve, and has no
    linear interpolation between grid points. It assumes a single boundary: once
    the threshold curve rises above (1 + margin) * th0 it stays above. Every saved
    curve does (their only dips are <0.1% of th0, on the plateau at about twice
    th0). width is nan if the input still spikes at max_delay.
    """
    site, model_kw = (stim, input_node, input_node2), model_kw or {}
    th0 = threshold0 if threshold0 is not None else threshold(
        model, 0.0, node, v0, epsg_tau, model_kw, rel_tol, stim=stim, input_node=input_node,
        input_node2=input_node2)
    amp = (1 + margin) * th0

    def fires(delay):
        return _spikes(model, amp, delay, node, v0, epsg_tau, site, model_kw)

    if fires(max_delay):
        return np.nan, th0
    lo, hi = 0.0, max_delay  # amp > th0, so it fires at delay 0
    while hi - lo > delay_tol:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if fires(mid) else (lo, mid)
    return lo + hi, th0  # twice the midpoint of the final bracket


def _trial(args):
    model, amp, t1, t2, node, v0, epsg_tau, site, model_kw = args
    first, second = min(t1, t2), max(t1, t2)
    return _spikes(model, amp, second - first, node, v0, epsg_tau, site, model_kw, start=first)


def probability_trials(model, delays, amplitude, n_trials=200, amp_cv=0.03, jitter=0.015,
                       node=3, v0=-68.0, epsg_tau=EPSG_TAU, model_kw=None, seed=0,
                       workers=None, executor=None, *, stim="EPSGpair", input_node=1,
                       input_node2=None):
    """Spike probability per delay from noisy trials, to check half_width's shortcut.

    Noise per trial: the pair's amplitude scaled by N(1, amp_cv), and each EPSG's
    onset jittered by N(0, jitter) ms independently. stim/input_node/input_node2
    as in threshold(); for EPSGbilateral the earlier EPSG is always at input_node.
    """
    site = (stim, input_node, input_node2)
    rng = np.random.default_rng(seed)
    tasks = []
    for d in delays:
        for _ in range(n_trials):
            amp = amplitude * rng.normal(1, amp_cv)
            t1 = START + rng.normal(0, jitter)
            t2 = START + d + rng.normal(0, jitter)
            tasks.append((model, amp, t1, t2, node, v0, tuple(epsg_tau), site, model_kw or {}))
    hits = map_tasks(_trial, tasks, workers, executor, chunksize=16)
    return np.array(hits, float).reshape(len(delays), n_trials).mean(axis=1)
