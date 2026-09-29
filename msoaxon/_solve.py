"""ode15s stand-in: scipy BDF, split at stimulus breakpoints."""

import numpy as np
from scipy.integrate import solve_ivp


def epsg_unitary(t):
    """Unitary EPSG waveform, peak-normalised by 0.21317."""
    return (1 / 0.21317) * (t > 0) * (np.exp(-t / 0.18) - np.exp(-t / 0.1))


def breakpoints(stim_type, start, stop, t_end):
    """Times where the stimulus switches on or off.

    Restarting the integrator there keeps it from stepping over a narrow EPSG
    onset after a long quiet stretch. ode15s had no explicit equivalent, but it
    only affects accuracy, not the model.
    """
    pts = {start, stop}
    if stim_type in ("ramp", "ramp2"):
        pts |= {5.0, stop + 5.0}
    return sorted(p for p in pts if 0 < p < t_end)


def integrate(rhs, y0, t_end, cuts, rtol, atol, max_step, jac_sparsity=None):
    """Integrate 0..t_end in pieces; return (t, y) with y shaped (n_times, n_states).

    Like ode15s without a tspan vector, the output is every accepted solver step.
    """
    edges = [0.0, *cuts, t_end]
    ts, ys = [], []
    y = np.asarray(y0, dtype=float)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        sol = solve_ivp(rhs, (a, b), y, method="BDF", rtol=rtol, atol=atol,
                        max_step=max_step, jac_sparsity=jac_sparsity)
        if not sol.success:
            raise RuntimeError(f"integration failed on [{a}, {b}]: {sol.message}")
        skip = 0 if i == 0 else 1  # drop the duplicated segment start
        ts.append(sol.t[skip:])
        ys.append(sol.y[:, skip:])
        y = sol.y[:, -1]
    return np.concatenate(ts), np.concatenate(ys, axis=1).T
