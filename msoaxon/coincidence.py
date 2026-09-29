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

from ._parallel import map_tasks
from .multi import mso_axon
from .spiking import count_spikes
from .synaptic import SynParams
from .two import two_cpt

START = 5.0
T_END = 20.0
FACTOR = 10.0  # spike = axon rises this far above soma (as in the EPSGpair sweeps)


def _spikes(model, I, delay, node, v0, epsg_tau, model_kw, start=START):
    """model_kw may carry stim ("EPSGpair" or, for the 45-compartment model,
    "EPSGbilateral"), input_node, and model keywords such as morph/input_node2/mem."""
    f = mso_axon if model == "multi" else two_cpt
    syn = SynParams(t_end=T_END, epsg_tau=tuple(epsg_tau))
    kw = dict(model_kw)
    stim, input_node = kw.pop("stim", "EPSGpair"), kw.pop("input_node", 1)
    t, x = f(stim, start, start + delay, I, node, "active-full", T_END, v0, input_node, syn,
             stop_on_spike=FACTOR, **kw)
    axon = node - 1 if model == "multi" else 1
    return t[-1] < T_END or count_spikes(x[:, 0], x[:, axon], FACTOR) > 0


def threshold(model, delay, node=3, v0=-68.0, epsg_tau=(0.1, 0.18), model_kw=None,
              rel_tol=1e-4, guess=50.0, ceiling=2000.0):
    """Smallest EPSG-pair amplitude that spikes at this delay, to rel_tol (bisection)."""
    model_kw = model_kw or {}
    lo, hi = 0.0, guess
    while not _spikes(model, hi, delay, node, v0, epsg_tau, model_kw):
        lo, hi = hi, hi * 2
        if hi > ceiling:
            return np.inf
    while (hi - lo) > rel_tol * hi:
        mid = (lo + hi) / 2
        if _spikes(model, mid, delay, node, v0, epsg_tau, model_kw):
            hi = mid
        else:
            lo = mid
    return hi


def _threshold_task(args):
    return threshold(*args)


def threshold_curve(model, delays, node=3, v0=-68.0, epsg_tau=(0.1, 0.18), model_kw=None,
                    rel_tol=1e-4, workers=None, executor=None):
    """threshold() at each delay, in parallel (needs a __main__ guard in scripts).

    workers / executor: see _parallel.map_tasks.
    """
    tasks = [(model, float(d), node, v0, tuple(epsg_tau), model_kw or {}, rel_tol)
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


def _trial(args):
    model, amp, t1, t2, node, v0, epsg_tau, model_kw = args
    first, second = min(t1, t2), max(t1, t2)
    return _spikes(model, amp, second - first, node, v0, epsg_tau, model_kw, start=first)


def probability_trials(model, delays, amplitude, n_trials=200, amp_cv=0.03, jitter=0.015,
                       node=3, v0=-68.0, epsg_tau=(0.1, 0.18), model_kw=None, seed=0,
                       workers=None, executor=None):
    """Spike probability per delay from noisy trials, to check half_width's shortcut.

    Noise per trial: the pair's amplitude scaled by N(1, amp_cv), and each EPSG's
    onset jittered by N(0, jitter) ms independently.
    """
    rng = np.random.default_rng(seed)
    tasks = []
    for d in delays:
        for _ in range(n_trials):
            amp = amplitude * rng.normal(1, amp_cv)
            t1 = START + rng.normal(0, jitter)
            t2 = START + d + rng.normal(0, jitter)
            tasks.append((model, amp, t1, t2, node, v0, tuple(epsg_tau), model_kw or {}))
    hits = map_tasks(_trial, tasks, workers, executor, chunksize=16)
    return np.array(hits, float).reshape(len(delays), n_trials).mean(axis=1)
