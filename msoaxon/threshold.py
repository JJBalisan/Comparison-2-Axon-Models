"""Spiking-threshold search (port of BinarySearch.m)."""

from msoaxon import _dispatch, _parallel, spiking, synaptic

START = 5.0
T_END = 20.0
V0 = -68.0
MODEL_TYPE = "active-full"
INPUT_NODE = 1


def sweep_setting(stim_type, i, epsg_pair_dt=1 / 25):
    """Returns the settings of sweep point i, as BinarySearch.m sets them.

    Args:
        stim_type: The stimulus being swept.
        i: The sweep point, 1-indexed as in the MATLAB loop.
        epsg_pair_dt: Spacing of the second-EPSG delay [ms].

    Returns:
        (stop, syn, I_override); I_override is None unless the stimulus
        fixes I.

    Raises:
        ValueError: If the threshold search doesn't support stim_type.
    """
    syn = synaptic.SynParams(t_end=T_END)
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
        return (
            START + (i - 1) * epsg_pair_dt,
            syn,
            None,
        )  # stop = second EPSG onset
    if stim_type == "Synaptic":
        syn.random_in = 13986
        return 15.0, syn, 0.0
    raise ValueError(
        f"stimType {stim_type!r} is not supported by the threshold search"
    )


def sweep_x(stim_type, i, epsg_pair_dt=1 / 25):
    """Returns what sweep point i varies, for plotting.

    Args:
        stim_type: "ramp", "ramp2", "sine" or "EPSGpair".
        i: The sweep point, 1-indexed.
        epsg_pair_dt: Spacing of the second-EPSG delay [ms].

    Returns:
        The ramp top, the sine frequency, or the delay between the two
        EPSGs.

    Raises:
        ValueError: For a stimulus with no x-axis.
    """
    stop, syn, _ = sweep_setting(stim_type, i, epsg_pair_dt)
    if stim_type in ("ramp", "ramp2"):
        return stop
    if stim_type == "sine":
        return syn.f
    if stim_type == "EPSGpair":
        return (
            (i - 1) * epsg_pair_dt
        )  # stop - START, without the rounding of the subtraction
    raise ValueError(f"no x-axis for stimType {stim_type!r}")


def _search(spiked_at, stim_type, i, factor, max_I, zoom, epsg_pair_dt):
    """Runs the halving search exactly as BinarySearch.m does it.

    Returns the last tested value.
    """
    location, previous, distance, first = max_I, 0.0, max_I, 0.0
    while distance > zoom and location <= max_I:
        stop, syn, I_override = sweep_setting(stim_type, i, epsg_pair_dt)
        I = location if I_override is None else I_override
        # BinarySearch.m counts spikes in the full trace (count_spikes). The run now
        # stops at the first spike, so an early end is that count being nonzero: any
        # sampled crossing is a sign change the terminal event catches
        spiked = spiked_at(stim_type, START, stop, I, syn, factor)
        tested = location
        first = location
        distance = abs((location - previous) / 2)
        location = location - distance if spiked else location + distance
        previous = tested
    return first


def _point(task):
    """Searches one sweep point for one model.

    Module-level so worker processes can import it.
    """
    model, stim_type, i, node, factor, max_I, zoom, epsg_pair_dt = task

    def spiked_at(
        stim, start, stop, I, syn, factor
    ):  # the search only needs yes/no
        return _dispatch.spikes(
            model,
            stim,
            start,
            stop,
            I,
            node,
            T_END,
            V0,
            INPUT_NODE,
            syn,
            model_type=MODEL_TYPE,
            factor=factor,
        )

    return _search(spiked_at, stim_type, i, factor, max_I, zoom, epsg_pair_dt)


def binary_search(
    stim_type,
    n_points,
    node,
    factor,
    max_I,
    zoom=1.0,
    epsg_pair_dt=1 / 25,
    rounded=True,
    workers=None,
    executor=None,
):
    """Runs BinarySearch.m's threshold search for both models.

    BinarySearch.m tests location1 (the multi-compartment variable) in
    the second loop's while condition; here each search tests its own
    location.

    Args:
        stim_type: The stimulus being swept.
        n_points: The number of sweep points.
        node: The compartment compared with the soma, 1-indexed.
        factor: The spike criterion: axon - soma [mV].
        max_I: The search's starting value and ceiling.
        zoom: The search stops once its step is at most this.
        epsg_pair_dt: Spacing of the second-EPSG delay. BinarySearch.m
            uses 1/25; the EPSGpair_Thresholds*.jpg plots in the repo
            were made with 0.1 and 11 points.
        rounded: Apply BinarySearch.m's rounding (to 1 for EPSGpair, to
            10 otherwise).
        workers: Processes for the independent sweep points (default:
            all CPUs; 1 runs in-process). Results are identical either
            way. Scripts that call this with workers != 1 need an
            `if __name__ == "__main__":` guard, since macOS starts worker
            processes by re-importing the main module.
        executor: An existing process pool to use instead (see
            _parallel.map_tasks).

    Returns:
        (thresholds_multi, thresholds_two), one value per sweep point.
    """
    # multi-compartment points first: they take ~10x longer, so starting them early
    # keeps the pool busy while the cheap two-compartment points fill the gaps
    tasks = [
        (model, stim_type, i, node, factor, max_I, zoom, epsg_pair_dt)
        for model in ("multi", "two")
        for i in range(1, n_points + 1)
    ]
    results = _parallel.map_tasks(_point, tasks, workers, executor)

    digits = 0 if stim_type == "EPSGpair" else -1
    finish = (
        (lambda v: spiking.matlab_round(v, digits))
        if rounded
        else (lambda v: v)
    )
    results = [finish(v) for v in results]
    return results[:n_points], results[n_points:]
