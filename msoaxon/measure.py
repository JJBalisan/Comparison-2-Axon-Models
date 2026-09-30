"""Passive and synaptic measurements on the soma, shared by scripts and tests."""

import numpy as np

from ._dispatch import run_model


def soma_on_grid(t, x, v0, t0, t1, dt=0.001):
    """The soma voltage relative to v0, interpolated onto a uniform grid [t0, t1)."""
    g = np.arange(t0, t1, dt)
    return g, np.interp(g, t, x[:, 0]) - v0


def passive_step(model, v0=-68.0, model_type="active-full", I=-2.0, n_grid=200001, **model_kw):
    """Soma input resistance and time constant from a 40 ms current step (5 to 45 ms).

    Returns dict(rin_steady, rin_peak [MOhm], tau [ms]); tau is the time to 63% of
    the steady response. The grid stops at 44.9 ms: the implicit solver's step that
    ends exactly at switch-off is evaluated with the current already off, so that
    one point has decayed a little. model_kw go to run_model (mem, morph, r1, ...).
    """
    t, x = run_model(model, "step", 5, 45, I, 3, 50, v0, 1, model_type=model_type, **model_kw)
    g = np.linspace(5, 44.9, n_grid)
    dv = np.interp(g, t, x[:, 0]) - v0
    tau = g[np.argmax(dv <= dv[-1] * (1 - np.exp(-1)))] - 5
    return dict(rin_steady=dv[-1] / I * 1e3, rin_peak=dv.min() / I * 1e3, tau=tau)
