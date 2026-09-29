"""ode15s stand-in: scipy BDF, split at stimulus breakpoints."""

import numpy as np
from scipy.integrate import solve_ivp

_QUIET_CACHE = {}
_QUIET_CACHE_SIZE = 64


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


def pre_stimulus_is_quiet(cuts, start, stop):
    """True when no stimulus acts before the first breakpoint, so [0, cuts[0]] can be shared."""
    return bool(cuts) and cuts[0] <= min(start, stop)


def _segment(rhs, a, b, y, rtol, atol, max_step, jac_sparsity):
    sol = solve_ivp(rhs, (a, b), y, method="BDF", rtol=rtol, atol=atol,
                    max_step=max_step, jac_sparsity=jac_sparsity)
    if not sol.success:
        raise RuntimeError(f"integration failed on [{a}, {b}]: {sol.message}")
    return sol.t, sol.y


def _quiet_prefix(key, rhs, b, y0, rtol, atol, max_step, jac_sparsity):
    """Solution on [0, b] with no stimulus, computed once per (model settings, y0, b)."""
    full_key = (key, np.asarray(y0, dtype=float).tobytes(), b, rtol, atol, max_step)
    hit = _QUIET_CACHE.get(full_key)
    if hit is None:
        hit = _segment(rhs, 0.0, b, y0, rtol, atol, max_step, jac_sparsity)
        if len(_QUIET_CACHE) >= _QUIET_CACHE_SIZE:
            _QUIET_CACHE.pop(next(iter(_QUIET_CACHE)))
        _QUIET_CACHE[full_key] = hit
    return hit


def integrate(rhs, y0, t_end, cuts, rtol, atol, max_step, jac_sparsity=None, quiet=None):
    """Integrate 0..t_end in pieces; return (t, y) with y shaped (n_times, n_states).

    Like ode15s without a tspan vector, the output is every accepted solver step.

    quiet: optional (key, rhs_without_stimulus). When given, the first segment
    [0, cuts[0]] is taken from a cache shared by every run with the same key —
    in a threshold sweep that is every run, since only the stimulus changes.
    `key` must capture everything the stimulus-free dynamics depend on.
    """
    edges = [0.0, *cuts, t_end]
    y = np.asarray(y0, dtype=float)
    ts, ys = [], []
    first = 0
    if quiet is not None and cuts:
        t, sol_y = _quiet_prefix(quiet[0], quiet[1], cuts[0], y, rtol, atol, max_step,
                                 jac_sparsity)
        ts.append(t)
        ys.append(sol_y)
        y = sol_y[:, -1]
        first = 1
    for i in range(first, len(edges) - 1):
        t, sol_y = _segment(rhs, edges[i], edges[i + 1], y, rtol, atol, max_step, jac_sparsity)
        skip = 0 if i == 0 else 1  # drop the duplicated segment start
        ts.append(t[skip:])
        ys.append(sol_y[:, skip:])
        y = sol_y[:, -1]
    return np.concatenate(ts), np.concatenate(ys, axis=1).T
