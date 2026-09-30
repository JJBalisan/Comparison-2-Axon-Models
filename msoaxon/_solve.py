"""ode15s stand-in: scipy BDF, split at stimulus breakpoints."""

import numpy as np
from scipy.integrate import solve_ivp

_QUIET_CACHE = {}
_QUIET_CACHE_SIZE = 64


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


def _segment(rhs, a, b, y, rtol, atol, max_step, jac_sparsity, events=None):
    """One solve; returns (t, y, stopped) where stopped means a terminal event fired."""
    sol = solve_ivp(rhs, (a, b), y, method="BDF", rtol=rtol, atol=atol,
                    max_step=max_step, jac_sparsity=jac_sparsity, events=events)
    if not sol.success:
        raise RuntimeError(f"integration failed on [{a}, {b}]: {sol.message}")
    return sol.t, sol.y, sol.status == 1


def spike_event(axon_col, factor):
    """Terminal event: axon compartment rises `factor` mV above the soma.

    The continuous-time version of count_spikes' first detection.
    """
    def event(t, x):
        return x[axon_col] - x[0] - factor
    event.terminal = True
    event.direction = 1
    return event


def _quiet_prefix(key, rhs, b, y0, rtol, atol, max_step, jac_sparsity):
    """Solution on [0, b] with no stimulus, computed once per (model settings, y0, b)."""
    full_key = (key, np.asarray(y0, dtype=float).tobytes(), b, rtol, atol, max_step)
    hit = _QUIET_CACHE.get(full_key)
    if hit is None:
        hit = _segment(rhs, 0.0, b, y0, rtol, atol, max_step, jac_sparsity)[:2]
        if len(_QUIET_CACHE) >= _QUIET_CACHE_SIZE:
            _QUIET_CACHE.pop(next(iter(_QUIET_CACHE)))
        _QUIET_CACHE[full_key] = hit
    return hit


def integrate(rhs, y0, t_end, cuts, rtol, atol, max_step, jac_sparsity=None, quiet=None,
              stop_event=None):
    """Integrate 0..t_end in pieces; return (t, y) with y shaped (n_times, n_states).

    Like ode15s without a tspan vector, the output is every accepted solver step.

    quiet: optional (key, rhs_without_stimulus). When given, the first segment
    [0, cuts[0]] is taken from a cache shared by every run with the same key —
    in a threshold sweep that is every run, since only the stimulus changes.
    `key` must capture everything the stimulus-free dynamics depend on.

    stop_event: optional terminal event; integration ends where it fires, so
    the returned t ends before t_end. Not applied to the cached quiet prefix,
    where there is no input to drive it.
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
        t, sol_y, stopped = _segment(rhs, edges[i], edges[i + 1], y, rtol, atol, max_step,
                                     jac_sparsity, stop_event)
        skip = 0 if i == 0 else 1  # drop the duplicated segment start
        ts.append(t[skip:])
        ys.append(sol_y[:, skip:])
        y = sol_y[:, -1]
        if stopped:
            break
    return np.concatenate(ts), np.concatenate(ys, axis=1).T
