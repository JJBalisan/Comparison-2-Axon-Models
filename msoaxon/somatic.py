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

from typing import Any

import numpy as np
from scipy import interpolate

from msoaxon import _bisect, _common, _dispatch, coincidence, multi

START, STOP, T_END = _common.START, 105.0, 110.0
FACTOR = _common.FACTOR  # spike = axon rises this far above soma
CHARGING = 0.15  # ms after step onset excluded from the inflection search
LOOKBACK = 0.5  # ms before the spike-detection time searched for the takeoff
# Smallest d2V/dt2 at the inflection [mV/ms^2] that counts as a distinct one.
# Set from the multi and two models' scan of step size (issue #38): the peak is
# 0-11 up to 1.25x rheobase (concave throughout, the smooth hump), 60-66 at
# 1.3x, 105-129 at 1.35x, 240-327 at 1.5x and 1600-1950 at 3x. 100 sits between
# 1.3x and 1.35x, 2.4x below the weakest published point (two, 1.5x).
D2_MIN = 100.0


def _merge_kw(
    mem: multi.Membrane | None, model_kw: dict[str, Any] | None
) -> dict[str, Any]:
    """Folds the `mem` shorthand into the model keywords.

    Args:
        mem: Channel overrides from multi.membrane(), or None. A
            shorthand for model_kw={"mem": mem}.
        model_kw: Further model keywords, or None.

    Returns:
        A new dict of the model keywords, with mem among them if given.

    Raises:
        ValueError: If mem is given both as `mem` and in model_kw.
    """
    kw = dict(model_kw or {})
    if mem is not None:
        if kw.get("mem") is not None:
            raise ValueError(
                "mem is given twice: as mem= and in model_kw; use one"
            )
        kw["mem"] = mem
    return kw


def _factor(
    criterion: coincidence.Criterion,
    model: str,
    node: int,
    v0: float,
    model_kw: dict[str, Any],
) -> float:
    """Resolves a criterion to mV, a Relative one on this configuration."""
    return coincidence.resolve_criterion(
        criterion, model, node, v0, model_kw=model_kw
    )


def fires(
    model: str,
    I: float,
    node: int = 3,
    v0: float = -68.0,
    mem: multi.Membrane | None = None,
    *,
    criterion: coincidence.Criterion = FACTOR,
    model_kw: dict[str, Any] | None = None,
) -> bool:
    """Returns whether a 100 ms somatic step makes compartment `node` spike.

    Args:
        model: "multi" or "two".
        I: Step amplitude [pA].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None. A
            shorthand for model_kw={"mem": mem}.
        criterion: The spike criterion: axon - soma [mV], or a
            coincidence.Relative fraction of the spike height.
        model_kw: Further model keywords, as in coincidence.threshold():
            mode, morph, model_type, max_step, or r1/tau_est for "two"
            (Goldwyn 2019: v0=-58, r1=8.5, tau_est=0.34).

    Returns:
        True if it spikes.

    Raises:
        ValueError: If mem is given both as `mem` and in model_kw.
    """
    kw = _merge_kw(mem, model_kw)
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
        factor=_factor(criterion, model, node, v0, kw),
        **kw,
    )


def response(
    model: str,
    I: float,
    node: int = 3,
    v0: float = -68.0,
    mem: multi.Membrane | None = None,
    *,
    criterion: coincidence.Criterion = FACTOR,
    model_kw: dict[str, Any] | None = None,
) -> tuple[bool, float]:
    """Runs a 100 ms somatic step like fires(), and says how close it came.

    Args:
        model: "multi" or "two".
        I: Step amplitude [pA].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None. A
            shorthand for model_kw={"mem": mem}.
        criterion: The spike criterion (see fires()).
        model_kw: Further model keywords (see fires()).

    Returns:
        (spiked, peak): fires()'s answer, and the largest axon - soma of
        the run [mV] (about the criterion if it spiked).

    Raises:
        ValueError: If mem is given both as `mem` and in model_kw.
    """
    kw = _merge_kw(mem, model_kw)
    return _dispatch.response(
        model,
        "step",
        START,
        STOP,
        I,
        node,
        T_END,
        v0,
        1,
        factor=_factor(criterion, model, node, v0, kw),
        **kw,
    )


def rheobase(
    model: str,
    node: int = 3,
    v0: float = -68.0,
    mem: multi.Membrane | None = None,
    rel_tol: float = 1e-4,
    guess: float = 2000.0,
    *,
    criterion: coincidence.Criterion = FACTOR,
    model_kw: dict[str, Any] | None = None,
) -> float:
    """Finds the smallest 100 ms somatic step that fires.

    By the predictive search (_bisect.smallest_crossing), as for the EPSG
    thresholds: in the multi and two models the step's peak axon - soma
    rises smoothly through the criterion (9.96 mV at 0.9999 rheobase,
    10.02 at 1.00005), so a few runs below it predict where it lies,
    10-11 runs instead of bisection's 15 (issue #19). With KLT blocked the
    response is all or none (0.3 mV just below rheobase, a 15 mV spike
    just above): no prediction lands inside the bracket and the search,
    halving down from `guess` (descent 2), makes exactly bisection's runs.

    Args:
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None. A
            shorthand for model_kw={"mem": mem}.
        rel_tol: The relative tolerance.
        guess: The first amplitude tried [pA].
        criterion: The spike criterion (see fires()).
        model_kw: Further model keywords (see fires()).

    Returns:
        The rheobase [pA], or inf if nothing up to 1e6 pA fires.

    Raises:
        _bisect.FiresWithoutInputError: If the cell fires with no input.
        ValueError: If mem is given both as `mem` and in model_kw.
    """
    kw = _merge_kw(mem, model_kw)
    factor = _factor(criterion, model, node, v0, kw)
    try:
        return _bisect.smallest_crossing(
            lambda I: response(
                model, I, node, v0, criterion=factor, model_kw=kw
            ),
            guess,
            1e6,
            rel_tol,
            factor,
            descent=2.0,
        )
    except _bisect.FiresWithoutInputError as e:
        extra = sorted(k for k in kw if k != "mem")
        config = _dispatch.describe(
            model,
            node,
            v0,
            kw.get("mem"),
            **({"model_kw": extra} if extra else {}),
        )
        raise _bisect.FiresWithoutInputError(
            f"somatic.rheobase({config}): {e}"
        ) from e


def spike_amplitude(
    model: str,
    I: float,
    node: int = 3,
    v0: float = -68.0,
    mem: multi.Membrane | None = None,
    *,
    criterion: coincidence.Criterion = FACTOR,
    model_kw: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Measures the somatic spike amplitude from the inflection point.

    The inflection exists only for part of the range of step sizes. In the
    multi and two models with the default criterion that is about 1.35x to 4.5x
    rheobase. Below it the soma is a smooth hump with no upstroke curvature
    (d2V/dt2 never rises above ~0 before the peak). Above it the axon crosses
    the criterion so early that the takeoff falls in the first CHARGING ms,
    which the search skips, and the picker would land on the upstroke next to
    the peak. In both cases the step still fires, so the result is not None,
    but `distinct` is False and the inflection fields are None instead of a
    number from the wrong feature. The limits move with the membrane and the
    criterion: with KLT removed the inflection stays distinct up to 10x.

    The inflection is the argmax of the second derivative of a cubic spline
    through the solver's own samples, evaluated on a 1 us grid. Differentiating
    a linear interpolant of the trace twice (the rule before issue #56) turns
    every solver step into a tooth of d2V/dt2 and the argmax into a lottery
    over the teeth, so a different platform's step placement moved a steep
    takeoff's amplitude by ~1 mV (issue #56).

    The inflection is distinct when all of these hold:
      - the d2V/dt2 peak in the search window is the largest before the
        spike peak, so the takeoff was not left out of the window;
      - that peak is not at either end of the window; and
      - it is at least D2_MIN, a real convex takeoff.

    Args:
        model: "multi" or "two".
        I: Step amplitude [pA].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None. A
            shorthand for model_kw={"mem": mem}.
        criterion: The spike criterion (see fires()); it also sets
            t_spike, which anchors the inflection search.
        model_kw: Further model keywords (see fires()); the finely
            sampled second run keeps its own max_step of 0.005 ms.

    Returns:
        dict(amplitude, distinct, peak_above_rest, t_spike, t_inflection,
        v_inflection, trace) in mV and ms, where trace = (t, x) is the
        finely sampled run up to 2 ms after the spike. `distinct` says
        whether the soma has an inflection; when it is False, amplitude,
        t_inflection and v_inflection are None, while peak_above_rest,
        t_spike and trace are still valid. None if the step doesn't fire.

    Raises:
        ValueError: If mem is given both as `mem` and in model_kw.
    """
    kw = _merge_kw(mem, model_kw)
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
        **{**kw, "stop_on_spike": _factor(criterion, model, node, v0, kw)},
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
        **{**kw, "max_step": 0.005},
    )
    grid = np.arange(START + 0.001, ts + 2.0, 0.001)
    keep = (t > START) & np.r_[True, np.diff(t) > 0]
    spline = interpolate.CubicSpline(t[keep], x[keep, 0])
    V = spline(grid)
    d2 = spline(grid, 2)
    i_pk = int(np.argmax(np.where(grid > ts - 0.3, V, -np.inf)))
    lo = max(START + CHARGING, ts - LOOKBACK)
    win = np.flatnonzero((grid >= lo) & (np.arange(len(grid)) < i_pk))
    result = {
        "amplitude": None,
        "distinct": False,
        "peak_above_rest": V[i_pk] - v0,
        "t_spike": ts,
        "t_inflection": None,
        "v_inflection": None,
        "trace": (t, x),
    }
    if win.size == 0:
        return result
    i_inf = win[np.argmax(d2[win])]
    before_peak = np.arange(len(grid)) < i_pk
    i_takeoff = np.argmax(np.where(before_peak, d2, -np.inf))
    if (
        i_takeoff != i_inf  # the sharpest curvature precedes the window
        or i_inf in (win[0], win[-1])
        or d2[i_inf] < D2_MIN
    ):
        return result
    result.update(
        amplitude=V[i_pk] - V[i_inf],
        distinct=True,
        t_inflection=grid[i_inf],
        v_inflection=V[i_inf],
    )
    return result
