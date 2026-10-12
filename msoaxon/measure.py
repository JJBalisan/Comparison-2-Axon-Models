"""Passive and synaptic measurements on the soma.

Shared by the scripts and the tests.
"""

import numpy as np

from msoaxon import _dispatch


def soma_on_grid(t, x, v0, t0, t1, dt=0.001):
    """Interpolates the soma voltage onto a uniform grid.

    Args:
        t: Times of the run [ms].
        x: The run's state, soma voltage in column 0.
        v0: The resting potential to subtract [mV].
        t0: Grid start [ms].
        t1: Grid end (excluded) [ms].
        dt: Grid spacing [ms].

    Returns:
        (grid, soma voltage relative to v0 on it).
    """
    g = np.arange(t0, t1, dt)
    return g, np.interp(g, t, x[:, 0]) - v0


def passive_step(
    model, v0=-68.0, model_type="active-full", I=-2.0, n_grid=200001, **model_kw
):
    """Measures soma input resistance and time constant from a current step.

    The step lasts 40 ms (5 to 45 ms). The grid stops at 44.9 ms: the
    implicit solver's step that ends exactly at switch-off is evaluated
    with the current already off, so that one point has decayed a
    little.

    Args:
        model: "multi" or "two".
        v0: Initial potential [mV].
        model_type: Which channels are active.
        I: Step amplitude [pA].
        n_grid: Points of the analysis grid.
        **model_kw: Further run_model keywords (mem, morph, r1, ...).

    Returns:
        dict(rin_steady, rin_peak [MOhm], tau [ms]); tau is the time to
        63% of the steady response.
    """
    t, x = _dispatch.run_model(
        model, "step", 5, 45, I, 3, 50, v0, 1, model_type=model_type, **model_kw
    )
    g = np.linspace(5, 44.9, n_grid)
    dv = np.interp(g, t, x[:, 0]) - v0
    tau = g[np.argmax(dv <= dv[-1] * (1 - np.exp(-1)))] - 5
    return dict(
        rin_steady=dv[-1] / I * 1e3, rin_peak=dv.min() / I * 1e3, tau=tau
    )
