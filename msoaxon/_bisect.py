"""Bracketing bisection for the smallest input that fires.

New analysis code, used by coincidence.threshold and somatic.rheobase. The
MATLAB port's own search (threshold._search, BinarySearch.m's halving search)
is separate and must stay as it is.
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
