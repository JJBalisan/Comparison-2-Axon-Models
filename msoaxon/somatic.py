"""Somatic action potential size, measured the way Scott et al 2005 did.

Scott, Mathews & Golding 2005 (J Neurosci 25:7887): spikes evoked by 100 ms
somatic current steps, amplitude measured "relative to the inflection point".
Mature MSO neurons: 17 +/- 2 mV, 5-15 mV near threshold and graded with the
stimulus; dendrotoxin (Kv1/KLT block) raised it from 15 to 37 mV (P20-21).

MSO neurons fire once, at step onset, while the membrane is still charging at
tens of mV/ms. So a fixed dV/dt criterion latches onto the charging, and the
inflection is taken instead as the peak of d2V/dt2 on the soma, searched after
the charging transient and shortly before the axonal spike.
"""

import numpy as np

from .mso_axon import mso_axon
from .two_cpt import two_cpt

START, STOP, T_END = 5.0, 105.0, 110.0
FACTOR = 10.0  # spike = axon rises this far above soma
CHARGING = 0.15  # ms after step onset excluded from the inflection search
LOOKBACK = 0.5  # ms before the spike-detection time searched for the takeoff


def _run(model, I, node, v0, mem, t_end=T_END, **kw):
    if model == "multi":
        morph = None if mem is None else getattr(mem, "morph", None)
        return mso_axon("step", START, STOP, I, node, "active-full", t_end, v0, 1, mem=mem,
                        morph=morph, **kw)
    return two_cpt("step", START, STOP, I, node, "active-full", t_end, v0, 1, **kw)


def fires(model, I, node=3, v0=-68.0, mem=None):
    t, _ = _run(model, I, node, v0, mem, stop_on_spike=FACTOR)
    return t[-1] < T_END


def rheobase(model, node=3, v0=-68.0, mem=None, rel_tol=1e-4, guess=2000.0):
    """Smallest 100 ms somatic step [pA] that fires, by bisection."""
    lo, hi = 0.0, guess
    while not fires(model, hi, node, v0, mem):
        lo, hi = hi, hi * 2
        if hi > 1e6:
            return np.inf
    while hi - lo > rel_tol * hi:
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if fires(model, mid, node, v0, mem) else (mid, hi)
    return hi


def spike_amplitude(model, I, node=3, v0=-68.0, mem=None):
    """Somatic spike amplitude from the inflection point [mV], plus details.

    Returns dict(amplitude, peak_above_rest, t_spike, t_inflection, v_inflection).
    """
    t0, _ = _run(model, I, node, v0, mem, stop_on_spike=FACTOR)
    ts = t0[-1]
    if ts >= T_END:
        return None
    t, x = _run(model, I, node, v0, mem, t_end=ts + 2.0, max_step=0.005)
    grid = np.arange(START + 0.001, ts + 2.0, 0.001)
    V = np.interp(grid, t, x[:, 0])
    d2 = np.gradient(np.gradient(V, grid), grid)
    i_pk = int(np.argmax(np.where(grid > ts - 0.3, V, -np.inf)))
    lo = max(START + CHARGING, ts - LOOKBACK)
    win = np.flatnonzero((grid >= lo) & (np.arange(len(grid)) < i_pk))
    i_inf = win[np.argmax(d2[win])]
    return dict(amplitude=V[i_pk] - V[i_inf], peak_above_rest=V[i_pk] - v0,
                t_spike=ts, t_inflection=grid[i_inf], v_inflection=V[i_inf])
