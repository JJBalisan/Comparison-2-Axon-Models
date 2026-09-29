"""Spike detection (port of Spiking.m and the inline counters in BinarySearch.m)."""

import math

import numpy as np


def spiking(x, factor, model):
    """Mark the first sample of each excursion where compartment j exceeds the soma by `factor`.

    Faithful to Spiking.m, including `reset` carrying over from one column to the
    next rather than being cleared per column.
    """
    ncols = 2 if model in ("Two", "two") else 45
    diff = x[:, :ncols] - x[:, [0]]
    spikes = np.zeros_like(diff)
    reset = False
    for j in range(ncols):
        for i in range(diff.shape[0]):
            if diff[i, j] > factor and not reset:
                spikes[i, j] = 1
                reset = True
            if diff[i, j] <= 0:
                reset = False
    return spikes


def count_rising_edges(col):
    """Combine_all.m's counter: 0 -> 1 transitions in a spike-marker column."""
    col = np.asarray(col)
    return int(np.sum((col[1:] > 0) & (col[:-1] <= 0)))


def count_spikes(soma, axon, factor):
    """BinarySearch.m's counter: axon above soma by `factor`, re-armed once axon < soma."""
    count, armed = 0, True
    for u1, u2 in zip(soma, axon):
        if u2 - u1 > factor and armed:
            count += 1
            armed = False
        if u2 < u1:
            armed = True
    return count


def matlab_round(x, digits=0):
    """MATLAB round(x, n): half away from zero (Python's round is half-to-even)."""
    scale = 10.0 ** digits
    return math.copysign(math.floor(abs(x) * scale + 0.5), x) * 10.0 ** -digits
