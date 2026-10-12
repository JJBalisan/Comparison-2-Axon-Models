"""Searches for the smallest input that fires.

New analysis code: smallest_firing (plain bisection, used by somatic.rheobase)
and crossing / smallest_crossing (bisection helped by interpolation, used by
coincidence.threshold and coincidence.window). The MATLAB port's own search
(threshold._search, BinarySearch.m's halving search) is separate and must stay
as it is.

Both searches assume zero input doesn't fire. When their first value fires
they check that with one run at zero, and raise FiresWithoutInputError if the
cell fires anyway (otherwise they would halve towards zero: smallest_firing
for ~1000 runs, smallest_crossing for ever).
"""

from collections.abc import Callable, Iterable

import numpy as np

Response = Callable[[float], tuple[bool, float]]


def smallest_firing(
    fires: Callable[[float], bool],
    guess: float,
    ceiling: float,
    rel_tol: float,
) -> float:
    """Finds the smallest x with fires(x), by doubling and then bisection.

    Assumes firing only gets easier as x grows. Doubles from `guess`
    until fires(hi), giving up with inf once the next doubling would
    exceed `ceiling` (so values between the last tested one and the
    ceiling are never tried), then bisects [lo, hi] until
    hi - lo <= rel_tol * hi. If `guess` itself fires, the bracket's
    lower end is zero, so zero is tested once first.

    Args:
        fires: The predicate, fires(x) -> bool.
        guess: The first value tried.
        ceiling: The largest value worth trying.
        rel_tol: The relative width of the final bracket.

    Returns:
        hi, the smallest value seen to fire, or inf.

    Raises:
        FiresWithoutInputError: If fires(0.0), tested only when `guess`
            fires.
    """
    lo, hi = 0.0, guess
    while not fires(hi):
        lo, hi = hi, hi * 2
        if hi > ceiling:
            return np.inf
    if lo == 0.0:  # guess fired: nothing yet shows zero is quiet
        _check_quiet(fires(0.0))
    while hi - lo > rel_tol * hi:
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if fires(mid) else (mid, hi)
    return hi


def crossing(
    respond: Response,
    x_no: float,
    x_yes: float,
    tol: Callable[[float, float], float],
    level: float,
    points: Iterable[tuple[float, float]] = (),
    log: bool = False,
) -> tuple[float, float]:
    """Narrows a bracket [x_no, x_yes] around where a response first spikes.

    Runs that don't spike give their peak below `level`, and the peak
    rises smoothly to `level` at the boundary (B4 of the search-speed
    plan: at the searches' default 10 mV criterion there is no jump
    there). A relative criterion of half the spike height (~19 mV, on
    the upstroke of the 45-compartment and TwoCpt.m spikes) still
    converges, with ~7% more runs (issue #21).

    Each step predicts the boundary from the non-spiking points
    (_predict) and probes there, kept 0.4 tol inside the bracket so that
    a near-exact prediction closes it in two runs. It bisects instead
    when there is no usable prediction, when the prediction falls
    outside the bracket, or when three probes in a row landed on the
    same side without halving the bracket.

    Args:
        respond: respond(x) -> (spiked, peak).
        x_no: A value that doesn't spike; it may be the larger one.
        x_yes: A value that spikes.
        tol: tol(x_no, x_yes) -> the bracket width to stop at.
        level: The peak at which a run spikes (the criterion).
        points: Known non-spiking (x, peak) pairs to start from.
        log: Predict in log(x) rather than x.

    Returns:
        The final (x_no, x_yes), once |x_yes - x_no| <= tol(x_no, x_yes).
    """
    pts = list(points)
    f, finv = (np.log, np.exp) if log else (float, float)
    same, last, mark = 0, None, abs(x_yes - x_no)
    while abs(x_yes - x_no) > tol(x_no, x_yes):
        lo, hi, t = min(x_no, x_yes), max(x_no, x_yes), 0.4 * tol(x_no, x_yes)
        x = _predict(pts, level, f, finv)
        if x is None or not lo < x < hi or (same >= 3 and hi - lo > mark / 2):
            x, same, mark = (lo + hi) / 2, 0, (hi - lo) / 2
        elif hi - lo <= 2 * t:
            x = (lo + hi) / 2
        else:
            x = min(max(x, lo + t), hi - t)
        spiked, peak = respond(x)
        if spiked:
            x_yes = x
        else:
            x_no = x
            pts.append((x, peak))
        same = same + 1 if spiked == last else 1
        if same == 1:
            mark = abs(x_yes - x_no)
        last = spiked
    return x_no, x_yes


def _predict(
    pts: list[tuple[float, float]],
    level: float,
    f: Callable[[float], float],
    finv: Callable[[float], float],
) -> float | None:
    """Predicts where the peak reaches `level` from the nearest points.

    Inverse quadratic interpolation (secant with two points) of f(x)
    against h = log(peak / level), through the three non-spiking points
    with the smallest |h|; one point alone uses the slope measured for
    the EPSG-pair thresholds, dh/dlog x ~ 6 (after the B4 experiments,
    scripts in the search-speed work). None if there is nothing to go
    on.
    """
    near = sorted(
        (abs(h), f(x), h)
        for x, p in pts
        for h in [np.log(max(p, 1e-3) / level)]
    )[:3]
    xs, hs = [x for _, x, _ in near], [h for _, _, h in near]
    if len(near) == 3 and len(set(hs)) == 3:
        (x0, x1, x2), (h0, h1, h2) = xs, hs
        est = (
            x0 * h1 * h2 / ((h0 - h1) * (h0 - h2))
            + x1 * h0 * h2 / ((h1 - h0) * (h1 - h2))
            + x2 * h0 * h1 / ((h2 - h0) * (h2 - h1))
        )
    elif len(near) >= 2 and hs[0] != hs[1]:
        est = xs[0] - hs[0] * (xs[1] - xs[0]) / (hs[1] - hs[0])
    elif len(near) == 1 and f is np.log:
        est = xs[0] - hs[0] / 6.0
    else:
        return None
    return finv(est) if np.isfinite(est) else None


def smallest_crossing(
    respond: Response,
    guess: float,
    ceiling: float,
    rel_tol: float,
    level: float,
) -> float:
    """Does smallest_firing with the bisection replaced by prediction.

    Starts at `guess` like smallest_firing, but finds the bracket by
    prediction (see crossing) rather than by doubling up or halving
    towards zero: while everything spikes it steps down (1.25x first,
    since thresholds usually sit just below the guess, then 2x per
    step), while nothing does it steps up (to the prediction, at most 2x
    per step). It gives up with inf past `ceiling` and keeps the same
    guarantee as smallest_firing; the value differs from
    smallest_firing's within the tolerance. Like smallest_firing, it tests
    zero once if `guess` spikes; that run is not used for prediction.

    Args:
        respond: respond(x) -> (spiked, peak).
        guess: The first value tried.
        ceiling: The largest value worth trying.
        rel_tol: The relative width of the final bracket.
        level: The peak at which a run spikes (the criterion).

    Returns:
        A value that spikes, with a non-spiking one within rel_tol of
        it, or inf.

    Raises:
        FiresWithoutInputError: If respond(0.0) spikes, tested only when
            `guess` spikes.
    """
    pts, x = [], guess
    lo = hi = None
    while lo is None or hi is None:
        if x > ceiling:
            return np.inf
        spiked, peak = respond(x)
        if spiked and lo is None and hi is None:  # only ever the guess
            _check_quiet(respond(0.0)[0])
        if spiked:
            hi = x if hi is None else min(hi, x)
        else:
            lo = x if lo is None else max(lo, x)
            pts.append((x, peak))
        if hi is None:  # nothing spikes yet: go up
            x = min(
                max(_predict(pts, level, np.log, np.exp) or 2 * lo, 1.01 * lo),
                2 * lo,
            )
        elif (
            lo is None
        ):  # everything spikes: come down, 1.25x first, then halving
            x = hi / 1.25 if hi == guess else hi / 2
    _, hi = crossing(
        respond, lo, hi, lambda _, y: rel_tol * y, level, pts, log=True
    )
    return hi


# --- a cell that fires with no input ------------------------------------------

REST_HINT = (
    "A multi-compartment cell built with rebalance_rest=False (the default "
    "of multi.membrane) settles ~0.9 mV away from v0 and, with KLT reduced, "
    "can fire on its own; build it with multi.membrane(v0, ..., "
    "rebalance_rest=True)."
)


class FiresWithoutInputError(ValueError):
    """The cell fires with no input, so no smallest firing input exists."""


def _check_quiet(spiked_at_zero: bool) -> None:
    """Raises FiresWithoutInputError if the run at zero input spiked."""
    if spiked_at_zero:
        raise FiresWithoutInputError(
            "the cell fires with no input (x = 0), so there is no smallest "
            "input that fires. " + REST_HINT
        )
