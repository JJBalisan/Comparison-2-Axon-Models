"""Passive and synaptic measurements on the soma.

Shared by the scripts and the tests.
"""

from typing import Any

import numpy as np

from msoaxon import _dispatch


def soma_on_grid(
    t: np.ndarray,
    x: np.ndarray,
    v0: float,
    t0: float,
    t1: float,
    dt: float = 0.001,
) -> tuple[np.ndarray, np.ndarray]:
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
    model: str,
    v0: float = -68.0,
    model_type: str = "active-full",
    I: float = -2.0,
    n_grid: int = 200001,
    **model_kw: Any,
) -> dict[str, float]:
    """Measures soma input resistance and time constant from a current step.

    The step lasts 40 ms (5 to 45 ms). The response is the soma voltage
    minus that of an identical run without the step, so any drift of the
    membrane away from v0 cancels. The responses are microvolts (2 pA
    through a few MOhm), while a membrane that is not rest-balanced
    drifts by about 0.1 mV before the step and 0.01 mV during it; a
    constant baseline (v0, or even the voltage at the step onset) is
    swamped by that. The grid stops at 44.9 ms: the implicit solver's
    step that ends exactly at switch-off is evaluated with the current
    already off, so that one point has decayed a little.

    Assumes the step is small enough to stay in the linear range, and
    that the response is a monotone relaxation to a steady state: the
    last grid point is taken as steady, and tau is the first crossing of
    63% of it. Either sign of I gives positive resistances; the peak is
    the largest deflection in the direction of the step. Resistances
    depend on the operating point, so the default membrane (resting a
    little below v0) reads about 2% above the rest-balanced one.

    Args:
        model: "multi" or "two".
        v0: Initial potential [mV].
        model_type: Which channels are active.
        I: Step amplitude [pA], nonzero.
        n_grid: Points of the analysis grid.
        **model_kw: Further run_model keywords (mem, morph, r1, ...).

    Returns:
        dict(rin_steady, rin_peak [MOhm], tau [ms]); tau is the time to
        63% of the steady response.

    Raises:
        ValueError: If I is zero.
    """
    if I == 0:
        raise ValueError("passive_step needs a nonzero step current I")
    start = 5

    def soma(current: float) -> np.ndarray:
        t, x = _dispatch.run_model(
            model,
            "step",
            start,
            45,
            current,
            3,
            50,
            v0,
            1,
            model_type=model_type,
            **model_kw,
        )
        return np.interp(g, t, x[:, 0])

    g = np.linspace(start, 44.9, n_grid)
    dv = soma(I) - soma(0.0)
    w = dv * np.sign(I)  # the response in the direction of the step
    tau = g[np.argmax(w >= w[-1] * (1 - np.exp(-1)))] - start
    return {
        "rin_steady": dv[-1] / I * 1e3,
        "rin_peak": dv[np.argmax(w)] / I * 1e3,
        "tau": tau,
    }
