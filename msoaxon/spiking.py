"""Spike detection.

Ports Spiking.m and the inline counters in BinarySearch.m.
"""

import math

import numpy as np


def spiking(x, factor, model):
    """Marks where each compartment first exceeds the soma by `factor`.

    Marks the first sample of each excursion. Faithful to Spiking.m,
    including `reset` carrying over from one column to the next rather
    than being cleared per column.

    Args:
        x: A run's state, voltages in the first columns.
        factor: The spike criterion: compartment - soma [mV].
        model: "Two" or "two" for the two-compartment model (2
            columns); anything else means 45 columns.

    Returns:
        An array like x[:, :ncols] with 1 at each marked sample.
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
    """Counts 0 -> 1 transitions in a spike-marker column (Combine_all.m).

    Args:
        col: One column of spiking()'s markers.

    Returns:
        The number of rising edges.
    """
    col = np.asarray(col)
    return int(np.sum((col[1:] > 0) & (col[:-1] <= 0)))


def count_spikes(soma, axon, factor):
    """Counts spikes the way BinarySearch.m does.

    A spike is the axon above the soma by `factor`; the counter re-arms
    once axon < soma.

    Args:
        soma: Soma voltage per sample [mV].
        axon: Axon voltage per sample [mV].
        factor: The spike criterion: axon - soma [mV].

    Returns:
        The number of spikes.
    """
    count, armed = 0, True
    for u1, u2 in zip(soma, axon):
        if u2 - u1 > factor and armed:
            count += 1
            armed = False
        if u2 < u1:
            armed = True
    return count


def matlab_round(x, digits=0):
    """Rounds like MATLAB's round(x, n): half away from zero.

    Python's round is half-to-even.

    Args:
        x: The value.
        digits: Decimal digits to keep; negative rounds to tens etc.

    Returns:
        The rounded value.
    """
    scale = 10.0**digits
    return math.copysign(math.floor(abs(x) * scale + 0.5), x) * 10.0**-digits
