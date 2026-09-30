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

from ._bisect import smallest_firing
from ._dispatch import run_model, spikes

START, STOP, T_END = 5.0, 105.0, 110.0
FACTOR = 10.0  # spike = axon rises this far above soma
CHARGING = 0.15  # ms after step onset excluded from the inflection search
LOOKBACK = 0.5  # ms before the spike-detection time searched for the takeoff


def fires(model, I, node=3, v0=-68.0, mem=None):
    """Whether a 100 ms somatic step of I [pA] makes compartment `node` spike."""
    return spikes(model, "step", START, STOP, I, node, T_END, v0, 1, mem=mem, factor=FACTOR)


def rheobase(model, node=3, v0=-68.0, mem=None, rel_tol=1e-4, guess=2000.0):
    """Smallest 100 ms somatic step [pA] that fires, by bisection (inf above 1e6)."""
    return smallest_firing(lambda I: fires(model, I, node, v0, mem), guess, 1e6, rel_tol)


def spike_amplitude(model, I, node=3, v0=-68.0, mem=None):
    """Somatic spike amplitude from the inflection point [mV], plus details.

    Returns dict(amplitude, peak_above_rest, t_spike, t_inflection, v_inflection,
    trace), where trace = (t, x) is the finely sampled run up to 2 ms after the
    spike, or None if the step doesn't fire.
    """
    t0, _ = run_model(model, "step", START, STOP, I, node, T_END, v0, 1, mem=mem,
                      stop_on_spike=FACTOR)
    ts = t0[-1]
    if ts >= T_END:
        return None
    t, x = run_model(model, "step", START, STOP, I, node, ts + 2.0, v0, 1, mem=mem,
                     max_step=0.005)
    grid = np.arange(START + 0.001, ts + 2.0, 0.001)
    V = np.interp(grid, t, x[:, 0])
    d2 = np.gradient(np.gradient(V, grid), grid)
    i_pk = int(np.argmax(np.where(grid > ts - 0.3, V, -np.inf)))
    lo = max(START + CHARGING, ts - LOOKBACK)
    win = np.flatnonzero((grid >= lo) & (np.arange(len(grid)) < i_pk))
    i_inf = win[np.argmax(d2[win])]
    return dict(amplitude=V[i_pk] - V[i_inf], peak_above_rest=V[i_pk] - v0,
                t_spike=ts, t_inflection=grid[i_inf], v_inflection=V[i_inf], trace=(t, x))
