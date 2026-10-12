"""Somatic action potential size, measured the way Scott et al 2005 did.

Scott, Mathews & Golding 2005 (J Neurosci 25:7887): spikes evoked by 100 ms
somatic current steps, amplitude measured "relative to the inflection point".
Mature MSO neurons: 17 +/- 2 mV, 5-15 mV near threshold and graded with the
stimulus; dendrotoxin (Kv1/KLT block) raised it from 15 to 37 mV (P20-21).

MSO neurons fire once, at step onset, while the membrane is still charging at
tens of mV/ms. So a fixed dV/dt criterion latches onto the charging, and the
inflection is taken instead as the peak of d2V/dt2 on the soma, searched after
the charging transient and shortly before the axonal spike.

A step fires when the axon rises `criterion` above the soma: 10 mV (FACTOR) by
default, or coincidence.Relative(f), a fraction f of the model's spike height
(coincidence.spike_height, measured with a coincident EPSG pair on the same
membrane).
"""

import types
from typing import Any

import numpy as np

from msoaxon import _bisect, _dispatch, coincidence

START, STOP, T_END = 5.0, 105.0, 110.0
FACTOR = 10.0  # spike = axon rises this far above soma
CHARGING = 0.15  # ms after step onset excluded from the inflection search
LOOKBACK = 0.5  # ms before the spike-detection time searched for the takeoff


def _factor(
    criterion: coincidence.Criterion,
    model: str,
    node: int,
    v0: float,
    mem: types.SimpleNamespace | None,
) -> float:
    """Resolves a criterion to mV, a Relative one on this membrane."""
    model_kw = {} if mem is None else {"mem": mem}
    return coincidence.resolve_criterion(
        criterion, model, node, v0, model_kw=model_kw
    )


def fires(
    model: str,
    I: float,
    node: int = 3,
    v0: float = -68.0,
    mem: types.SimpleNamespace | None = None,
    *,
    criterion: coincidence.Criterion = FACTOR,
) -> bool:
    """Returns whether a 100 ms somatic step makes compartment `node` spike.

    Args:
        model: "multi" or "two".
        I: Step amplitude [pA].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None.
        criterion: The spike criterion: axon - soma [mV], or a
            coincidence.Relative fraction of the spike height.

    Returns:
        True if it spikes.
    """
    return _dispatch.spikes(
        model,
        "step",
        START,
        STOP,
        I,
        node,
        T_END,
        v0,
        1,
        mem=mem,
        factor=_factor(criterion, model, node, v0, mem),
    )


def rheobase(
    model: str,
    node: int = 3,
    v0: float = -68.0,
    mem: types.SimpleNamespace | None = None,
    rel_tol: float = 1e-4,
    guess: float = 2000.0,
    *,
    criterion: coincidence.Criterion = FACTOR,
) -> float:
    """Finds the smallest 100 ms somatic step that fires, by bisection.

    Args:
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None.
        rel_tol: The relative tolerance.
        guess: The first amplitude tried [pA].
        criterion: The spike criterion (see fires()).

    Returns:
        The rheobase [pA], or inf if nothing up to 1e6 pA fires.

    Raises:
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
    factor = _factor(criterion, model, node, v0, mem)
    try:
        return _bisect.smallest_firing(
            lambda I: fires(model, I, node, v0, mem, criterion=factor),
            guess,
            1e6,
            rel_tol,
        )
    except _bisect.FiresWithoutInputError as e:
        config = _dispatch.describe(model, node, v0, mem)
        raise _bisect.FiresWithoutInputError(
            f"somatic.rheobase({config}): {e}"
        ) from e


def spike_amplitude(
    model: str,
    I: float,
    node: int = 3,
    v0: float = -68.0,
    mem: types.SimpleNamespace | None = None,
    *,
    criterion: coincidence.Criterion = FACTOR,
) -> dict[str, Any] | None:
    """Measures the somatic spike amplitude from the inflection point.

    Args:
        model: "multi" or "two".
        I: Step amplitude [pA].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None.
        criterion: The spike criterion (see fires()); it also sets
            t_spike, which anchors the inflection search.

    Returns:
        dict(amplitude, peak_above_rest, t_spike, t_inflection,
        v_inflection, trace) in mV and ms, where trace = (t, x) is the
        finely sampled run up to 2 ms after the spike; or None if the
        step doesn't fire.
    """
    t0, _ = _dispatch.run_model(
        model,
        "step",
        START,
        STOP,
        I,
        node,
        T_END,
        v0,
        1,
        mem=mem,
        stop_on_spike=_factor(criterion, model, node, v0, mem),
    )
    ts = t0[-1]
    if ts >= T_END:
        return None
    t, x = _dispatch.run_model(
        model,
        "step",
        START,
        STOP,
        I,
        node,
        ts + 2.0,
        v0,
        1,
        mem=mem,
        max_step=0.005,
    )
    grid = np.arange(START + 0.001, ts + 2.0, 0.001)
    V = np.interp(grid, t, x[:, 0])
    d2 = np.gradient(np.gradient(V, grid), grid)
    i_pk = int(np.argmax(np.where(grid > ts - 0.3, V, -np.inf)))
    lo = max(START + CHARGING, ts - LOOKBACK)
    win = np.flatnonzero((grid >= lo) & (np.arange(len(grid)) < i_pk))
    i_inf = win[np.argmax(d2[win])]
    return {
        "amplitude": V[i_pk] - V[i_inf],
        "peak_above_rest": V[i_pk] - v0,
        "t_spike": ts,
        "t_inflection": grid[i_inf],
        "v_inflection": V[i_inf],
        "trace": (t, x),
    }
