"""ode15s stand-in: scipy BDF, split at stimulus breakpoints."""

from collections.abc import Callable, Hashable
from typing import Any

import numpy as np
from scipy import integrate as scipy_integrate

# f(t, x) -> dx/dt, and a solve_ivp event function g(t, x)
Rhs = Callable[[float, np.ndarray], np.ndarray | list[float]]
Event = Callable[[float, np.ndarray], float]

_QUIET_CACHE = {}
_QUIET_CACHE_SIZE = 64

# msoAxon.m and TwoCptODE.m hardcode 5 ms in the ramp stimuli: "ramp" is gated
# on 5 < t <= stop + 5, and "ramp2" in two_cpt rises from t = 5 whatever
# `start` is. One name for all of those sites (breakpoints, multi, two).
RAMP_RISE = 5.0  # [ms]
# input_end: an EPSG has decayed to ~2e-9 of its peak 20 decay constants after
# its last onset.
EPSG_TAIL = 20
# The searches stop a run early for one of two reasons, and spikes() tells them
# apart by the final state of axon - soma (in units of the criterion `factor`):
# a spike stop ends at factor (spike_event); a settled stop (settled_event)
# needs axon - soma below SETTLED_FRACTION * factor, so it ends at or below
# that. The test sits between them at SPIKE_TEST_FRACTION * factor, away from
# both: the settled event can fire exactly as axon - soma falls through
# SETTLED_FRACTION * factor, so testing there misread three near-threshold ramp
# runs as spikes.
SETTLED_FRACTION = 0.5
SPIKE_TEST_FRACTION = 0.75


def breakpoints(
    stim_type: str, start: float, stop: float, t_end: float
) -> list[float]:
    """Returns the times where the stimulus switches on or off.

    Restarting the integrator there keeps it from stepping over a narrow
    EPSG onset after a long quiet stretch. ode15s had no explicit
    equivalent, but it only affects accuracy, not the model.

    Args:
        stim_type: The stimulus name.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        t_end: End of the run [ms].

    Returns:
        The sorted breakpoints strictly between 0 and t_end.
    """
    pts = {start, stop}
    if stim_type in ("ramp", "ramp2"):
        pts |= {RAMP_RISE, stop + RAMP_RISE}
    return sorted(p for p in pts if 0 < p < t_end)


def input_end(
    stim_type: str,
    start: float,
    stop: float,
    t_end: float,
    epsg_tau: tuple[float, float],
) -> float:
    """Returns the time after which the stimulus no longer drives the cell.

    EPSGs: the last onset plus EPSG_TAIL decay constants, when the conductance
    is below ~2e-9 of its peak (neither model's EPSG ends sooner:
    two_cpt never cuts it off). step, ramp, sine and the synaptic inputs
    are off after `stop`. ramp2 holds its current to t_end, so it never
    ends.

    Args:
        stim_type: The stimulus name.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        t_end: End of the run [ms].
        epsg_tau: The EPSG (rise, decay) time constants [ms].

    Returns:
        The time [ms]; t_end for a stimulus that never ends.
    """
    if stim_type in ("EPSG", "EPSGpair", "EPSGbilateral"):
        last = start if stim_type == "EPSG" else max(start, stop)
        return last + EPSG_TAIL * max(epsg_tau)
    if stim_type in ("step", "ramp", "sine", "Synaptic", "SynapticPair"):
        return max(start, stop)
    return t_end


def pre_stimulus_is_quiet(cuts: list[float], start: float, stop: float) -> bool:
    """Returns whether no stimulus acts before the first breakpoint.

    If so, the solution on [0, cuts[0]] can be shared between runs.

    Args:
        cuts: The breakpoints (see breakpoints()).
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].

    Returns:
        True when the first segment is stimulus-free.
    """
    return bool(cuts) and cuts[0] <= min(start, stop)


def _segment(
    rhs: Rhs,
    a: float,
    b: float,
    y: np.ndarray,
    rtol: float,
    atol: float,
    max_step: float,
    jac_sparsity: Any,
    events: Event | list[Event] | None = None,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Solves one segment.

    Returns (t, y, stopped), where stopped means a terminal event fired.
    """
    sol = scipy_integrate.solve_ivp(
        rhs,
        (a, b),
        y,
        method="BDF",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        jac_sparsity=jac_sparsity,
        events=events,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed on [{a}, {b}]: {sol.message}")
    return sol.t, sol.y, sol.status == 1


def spike_event(axon_col: int, factor: float) -> Event:
    """Returns a terminal event: the axon rises `factor` mV above the soma.

    The continuous-time version of count_spikes' first detection.

    Args:
        axon_col: The state column of the compared compartment.
        factor: The spike criterion: axon - soma [mV].

    Returns:
        An event function for solve_ivp (terminal, rising).
    """

    # solve_ivp calls events as g(t, x); this one needs only x
    def event(t: float, x: np.ndarray) -> float:  # noqa: ARG001
        return x[axon_col] - x[0] - factor

    event.terminal = True
    event.direction = 1
    return event


def settled_event(dv: Rhs, t_off: float, axon_col: int, factor: float) -> Event:
    """Returns a terminal event: after t_off, stop once no spike can follow.

    Fires at the first time after t_off when axon - soma is below
    SETTLED_FRACTION * factor (factor / 2) and every compartment's
    voltage is falling. t_off (input_end) guarantees the input has
    stopped driving the cell; all voltages falling guarantees its
    response has peaked everywhere, a
    stronger condition than axon - soma falling (which also passed the
    benchmark's --check, but stops later). With no input and every
    compartment repolarising, the run's yes/no is decided; the searches,
    which only need that answer, stop here instead of integrating the
    return to rest (scripts/benchmark_search.py measures the saving).

    The event is -1 until then and only switches sign when the condition
    holds, so it never alters the solver's steps: a run that doesn't
    stop is identical to one without the event, and one that does is
    identical up to the stop. spikes() tells the two stops apart by the
    final state: a spike stop ends with axon - soma at factor, a settled
    stop at or below SETTLED_FRACTION * factor (see SPIKE_TEST_FRACTION).

    Only the voltage derivatives are needed, so the models pass a
    function that computes those alone (multi._dv, two._dv), bit for bit
    the leading entries of their right-hand side.

    Args:
        dv: The model's voltage derivatives, f(t, x) -> dV/dt for every
            compartment.
        t_off: When the input stops driving the cell (input_end) [ms].
        axon_col: The state column of the compared compartment.
        factor: The spike criterion: axon - soma [mV].

    Returns:
        An event function for solve_ivp (terminal, rising).
    """

    def event(t: float, x: np.ndarray) -> float:
        if t < t_off or x[axon_col] - x[0] >= factor * SETTLED_FRACTION:
            return -1.0
        return 1.0 if np.all(np.asarray(dv(t, x)) < 0) else -1.0

    event.terminal = True
    event.direction = 1
    return event


def _quiet_prefix(
    key: Hashable,
    rhs: Rhs,
    b: float,
    y0: np.ndarray,
    rtol: float,
    atol: float,
    max_step: float,
    jac_sparsity: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Returns the solution on [0, b] with no stimulus.

    Computed once per (model settings, y0, b) and cached.
    """
    full_key = (
        key,
        np.asarray(y0, dtype=float).tobytes(),
        b,
        rtol,
        atol,
        max_step,
    )
    hit = _QUIET_CACHE.get(full_key)
    if hit is None:
        hit = _segment(rhs, 0.0, b, y0, rtol, atol, max_step, jac_sparsity)[:2]
        if len(_QUIET_CACHE) >= _QUIET_CACHE_SIZE:
            _QUIET_CACHE.pop(next(iter(_QUIET_CACHE)))
        _QUIET_CACHE[full_key] = hit
    return hit


def integrate(
    rhs: Rhs,
    y0: np.ndarray,
    t_end: float,
    cuts: list[float],
    rtol: float,
    atol: float,
    max_step: float,
    jac_sparsity: Any = None,
    quiet: tuple[Hashable, Rhs] | None = None,
    stop_event: Event | list[Event] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrates from 0 to t_end in pieces split at the breakpoints.

    Like ode15s without a tspan vector, the output is every accepted
    solver step.

    Args:
        rhs: The right-hand side, f(t, y).
        y0: The initial state.
        t_end: End of the run [ms].
        cuts: The breakpoints (see breakpoints()).
        rtol: Relative tolerance of the BDF solver.
        atol: Absolute tolerance of the BDF solver.
        max_step: Largest solver step [ms].
        jac_sparsity: Optional sparsity pattern of the Jacobian.
        quiet: Optional (key, rhs_without_stimulus). When given, the
            first segment [0, cuts[0]] is taken from a cache shared by
            every run with the same key; in a threshold sweep that is
            every run, since only the stimulus changes. `key` must
            capture everything the stimulus-free dynamics depend on.
        stop_event: Optional terminal event, or a list of them
            (spike_event and settled_event); integration ends where one
            fires, so the returned t ends before t_end. Not applied to
            the cached quiet prefix, where there is no input to drive it.

    Returns:
        (t, y) with y shaped (n_times, n_states).

    Raises:
        RuntimeError: If the solver fails on a segment.
    """
    edges = [0.0, *cuts, t_end]
    y = np.asarray(y0, dtype=float)
    ts, ys = [], []
    first = 0
    if quiet is not None and cuts:
        t, sol_y = _quiet_prefix(
            quiet[0], quiet[1], cuts[0], y, rtol, atol, max_step, jac_sparsity
        )
        ts.append(t)
        ys.append(sol_y)
        y = sol_y[:, -1]
        first = 1
    for i in range(first, len(edges) - 1):
        t, sol_y, stopped = _segment(
            rhs,
            edges[i],
            edges[i + 1],
            y,
            rtol,
            atol,
            max_step,
            jac_sparsity,
            stop_event,
        )
        skip = 0 if i == 0 else 1  # drop the duplicated segment start
        ts.append(t[skip:])
        ys.append(sol_y[:, skip:])
        y = sol_y[:, -1]
        if stopped:
            break
    return np.concatenate(ts), np.concatenate(ys, axis=1).T
