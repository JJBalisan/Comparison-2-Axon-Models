"""Spiking-threshold search (port of BinarySearch.m)."""

import os
from concurrent.futures import ProcessPoolExecutor

from .mso_axon import mso_axon
from .spiking import count_spikes, matlab_round
from .synaptic import SynParams
from .two_cpt import two_cpt

START = 5.0
T_END = 20.0
V0 = -68.0
MODEL_TYPE = "active-full"
INPUT_NODE = 1


def sweep_setting(stim_type, i, epsg_pair_dt=1 / 25):
    """(stop, syn, I_override) for sweep point i (1-indexed, as in the MATLAB loop)."""
    syn = SynParams(t_end=T_END)
    if stim_type == "step":
        return 15.0, syn, None
    if stim_type in ("ramp", "ramp2"):
        return 5 + i / 10, syn, None
    if stim_type == "sine":
        syn.f = i * 100
        return START + 1000 / (2 * syn.f), syn, None
    if stim_type == "EPSG":
        return 10.0, syn, None
    if stim_type == "EPSGpair":
        return START + (i - 1) * epsg_pair_dt, syn, None  # stop = second EPSG onset
    if stim_type == "Synaptic":
        syn.random_in = 13986
        return 15.0, syn, 0.0
    raise ValueError(f"stimType {stim_type!r} is not supported by the threshold search")


def _search(run, soma_col, axon_col, stim_type, i, factor, max_I, zoom, epsg_pair_dt):
    """The halving search exactly as BinarySearch.m does it; returns the last tested value."""
    location, previous, distance, first = max_I, 0.0, max_I, 0.0
    while distance > zoom and location <= max_I:
        stop, syn, I_override = sweep_setting(stim_type, i, epsg_pair_dt)
        I = location if I_override is None else I_override
        _, x = run(stim_type, START, stop, I, syn)
        spiked = count_spikes(x[:, soma_col], x[:, axon_col], factor) != 0
        tested = location
        first = location
        distance = abs((location - previous) / 2)
        location = location - distance if spiked else location + distance
        previous = tested
    return first


def _run(model, stim, start, stop, I, syn, node):
    f = mso_axon if model == "multi" else two_cpt
    return f(stim, start, stop, I, node, MODEL_TYPE, T_END, V0, INPUT_NODE, syn)


def _point(task):
    """One sweep point for one model; module-level so worker processes can import it."""
    model, stim_type, i, node, factor, max_I, zoom, epsg_pair_dt = task
    axon_col = node - 1 if model == "multi" else 1
    run = lambda stim, start, stop, I, syn: _run(model, stim, start, stop, I, syn, node)
    return _search(run, 0, axon_col, stim_type, i, factor, max_I, zoom, epsg_pair_dt)


def binary_search(stim_type, n_points, node, factor, max_I, zoom=1.0,
                  epsg_pair_dt=1 / 25, rounded=True, workers=None):
    """Return (thresholds_multi, thresholds_two), one value per sweep point.

    epsg_pair_dt: spacing of the second-EPSG delay. BinarySearch.m uses 1/25;
    the EPSGpair_Thresholds*.jpg plots in the repo were made with 0.1 and 11 points.
    rounded: apply BinarySearch.m's rounding (to 1 for EPSGpair, to 10 otherwise).
    workers: processes for the independent sweep points (default: all CPUs; 1 runs
    in-process). Results are identical either way. Scripts that call this with
    workers != 1 need an `if __name__ == "__main__":` guard, since macOS starts
    worker processes by re-importing the main module.

    BinarySearch.m tests location1 (the multi-compartment variable) in the second
    loop's while condition; here each search tests its own location.
    """
    # multi-compartment points first: they take ~10x longer, so starting them early
    # keeps the pool busy while the cheap two-compartment points fill the gaps
    tasks = [(model, stim_type, i, node, factor, max_I, zoom, epsg_pair_dt)
             for model in ("multi", "two") for i in range(1, n_points + 1)]
    workers = min(workers or os.cpu_count() or 1, len(tasks))
    if workers == 1:
        results = [_point(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(_point, tasks))

    digits = 0 if stim_type == "EPSGpair" else -1
    finish = (lambda v: matlab_round(v, digits)) if rounded else (lambda v: v)
    results = [finish(v) for v in results]
    return results[:n_points], results[n_points:]
