"""Searches for the smallest input that fires.

New analysis code: smallest_firing (plain bisection, used by somatic.rheobase)
and crossing / smallest_crossing (bisection helped by interpolation, used by
coincidence.threshold and coincidence.window). The MATLAB port's own search
(threshold._search, BinarySearch.m's halving search) is separate and must stay
as it is.
"""

import numpy as np


def smallest_firing(fires, guess, ceiling, rel_tol):
    """Smallest x with fires(x), assuming firing only gets easier as x grows.

    Doubles from `guess` until fires(hi), giving up with inf once the next doubling
    would exceed `ceiling` (so values between the last tested one and the ceiling
    are never tried), then bisects [lo, hi] until hi - lo <= rel_tol * hi.
    Returns hi, the smallest value seen to fire.
    """
    lo, hi = 0.0, guess
    while not fires(hi):
        lo, hi = hi, hi * 2
        if hi > ceiling:
            return np.inf
    while hi - lo > rel_tol * hi:
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if fires(mid) else (mid, hi)
    return hi


def crossing(respond, x_no, x_yes, tol, level, points=(), log=False):
    """Narrow a bracket [x_no, x_yes] around where a response first spikes.

    respond(x) -> (spiked, peak); runs that don't spike give their peak below
    `level`, and the peak rises smoothly to `level` at the boundary (B4 of the
    search-speed plan: at the searches' 10 mV criterion there is no jump there).
    x_no doesn't spike, x_yes does; either may be the larger. Stops once
    |x_yes - x_no| <= tol(x_no, x_yes) and returns the final (x_no, x_yes).

    Each step predicts the boundary from the non-spiking points (_predict) and
    probes there, kept 0.4 tol inside the bracket so that a near-exact prediction
    closes it in two runs. It bisects instead when there is no usable prediction,
    when the prediction falls outside the bracket, or when three probes in a row
    landed on the same side without halving the bracket. points: known
    non-spiking (x, peak) pairs to start from.
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


def _predict(pts, level, f, finv):
    """Where the peak reaches `level`, from the non-spiking points nearest it.

    Inverse quadratic interpolation (secant with two points) of f(x) against
    h = log(peak / level), through the three points with the smallest |h|; one
    point alone uses the slope measured for the EPSG-pair thresholds, dh/dlog x
    ~ 6 (after the B4 experiments, scripts in the search-speed work). None if
    there is nothing to go on.
    """
    near = sorted((abs(h), f(x), h) for x, p in pts for h in [np.log(max(p, 1e-3) / level)])[:3]
    xs, hs = [x for _, x, _ in near], [h for _, _, h in near]
    if len(near) == 3 and len(set(hs)) == 3:
        (x0, x1, x2), (h0, h1, h2) = xs, hs
        est = (x0 * h1 * h2 / ((h0 - h1) * (h0 - h2)) + x1 * h0 * h2 / ((h1 - h0) * (h1 - h2))
               + x2 * h0 * h1 / ((h2 - h0) * (h2 - h1)))
    elif len(near) >= 2 and hs[0] != hs[1]:
        est = xs[0] - hs[0] * (xs[1] - xs[0]) / (hs[1] - hs[0])
    elif len(near) == 1 and f is np.log:
        est = xs[0] - hs[0] / 6.0
    else:
        return None
    return finv(est) if np.isfinite(est) else None


def smallest_crossing(respond, guess, ceiling, rel_tol, level):
    """smallest_firing with the bisection replaced by prediction (see crossing).

    Starts at `guess` like smallest_firing, but finds the bracket by prediction
    rather than by doubling up or halving towards zero: while everything spikes it
    steps down (1.25x first, since thresholds usually sit just below the guess,
    then 2x per step), while nothing does it steps up (to the prediction, at most
    2x per step). Same giving up with inf
    past `ceiling`, and the same guarantee: returns a value that spikes, with a
    non-spiking one within rel_tol of it; the value differs from smallest_firing's
    within that tolerance.
    """
    pts, x = [], guess
    lo = hi = None
    while lo is None or hi is None:
        if x > ceiling:
            return np.inf
        spiked, peak = respond(x)
        if spiked:
            hi = x if hi is None else min(hi, x)
        else:
            lo = x if lo is None else max(lo, x)
            pts.append((x, peak))
        if hi is None:  # nothing spikes yet: go up
            x = min(max(_predict(pts, level, np.log, np.exp) or 2 * lo, 1.01 * lo), 2 * lo)
        elif lo is None:  # everything spikes: come down, 1.25x first, then halving
            x = hi / 1.25 if hi == guess else hi / 2
    _, hi = crossing(respond, lo, hi, lambda n, y: rel_tol * y, level, pts, log=True)
    return hi
