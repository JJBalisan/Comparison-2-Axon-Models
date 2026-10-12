"""Stimulus and argument handling shared by both models (from TwoCpt.m)."""

import types

from msoaxon import _solve, constants as C, synaptic

STIM_TYPES = (
    "step",
    "ramp",
    "ramp2",
    "sine",
    "Synaptic",
    "SynapticPair",
    "EPSG",
    "EPSGpair",
)
MODEL_TYPES = (
    "passive",
    "active-KLT",
    "active-H",
    "active-KLT+H",
    "active-sodium",
    "Active-sodium",
    "active-KHT",
    "active-full",
)


def check_args(
    stim_type,
    model_type,
    node,
    input_node,
    min_node,
    stim_types=STIM_TYPES,
    n_max=C.N_CPT,
):
    """Rejects inputs that numpy's negative indexing would accept silently.

    Args:
        stim_type: The stimulus name; must be one of `stim_types`.
        model_type: The channel set; must be one of MODEL_TYPES.
        node: The compared compartment, 1-indexed, min_node..N_CPT.
        input_node: Where the stimulus enters, 1-indexed, 1..n_max.
        min_node: The smallest allowed `node` (2 in two_cpt, where node
            1 would be the soma itself).
        stim_types: The stimulus names the calling model accepts.
        n_max: The number of compartments an input may enter.

    Raises:
        ValueError: If any argument is out of range or unknown.
    """
    if stim_type not in stim_types:
        raise ValueError(f"unknown stimType {stim_type!r}")
    if model_type not in MODEL_TYPES:
        raise ValueError(
            f"unknown model type {model_type!r}; expected one of {MODEL_TYPES}"
        )
    if not min_node <= node <= C.N_CPT:
        raise ValueError(
            f"node must be {min_node}..{C.N_CPT} (1-indexed), got {node}"
        )
    if not 1 <= input_node <= n_max:
        raise ValueError(
            f"input_node must be 1..{n_max} (1-indexed), got {input_node}"
        )


def stimulus(stim_type, start, stop, I, t_end, syn):
    """Bundles the stimulus settings (what TwoCpt.m stored on P).

    Args:
        stim_type: The stimulus name; "sine" and the synaptic inputs
            take extra settings from `syn`.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms] (its meaning depends on stim_type).
        I: Stimulus amplitude (pA or nS, depending on stim_type).
        t_end: End of the run [ms].
        syn: A SynParams with the waveform settings.

    Returns:
        A SimpleNamespace with start, stop, I, t_end and epsg_tau, plus
        f for "sine", and t_syn, g_syn, VsynE and diff for the
        synaptic inputs.
    """
    s = types.SimpleNamespace(
        start=start, stop=stop, I=I, t_end=t_end, epsg_tau=tuple(syn.epsg_tau)
    )
    if stim_type == "sine":
        s.f = syn.f
    if stim_type in ("Synaptic", "SynapticPair"):
        s.t_syn, s.g_syn = synaptic.synaptic(syn)
        s.VsynE, s.diff = syn.VsynE, syn.diff
    return s


def _stops(rhs, axon_col, stop_on_spike, stop_when_settled, stim_type, s, n_v):
    """Returns the terminal events for a run: none, the spike stop, or both.

    n_v: how many leading state entries are voltages (2 for two_cpt, n
    for mso_axon).
    """
    if stop_on_spike is None:
        if stop_when_settled:
            raise ValueError("stop_when_settled needs stop_on_spike")
        return None
    events = [_solve.spike_event(axon_col, stop_on_spike)]
    if stop_when_settled:
        t_off = _solve.input_end(
            stim_type, s.start, s.stop, s.t_end, s.epsg_tau
        )
        if t_off < s.t_end:
            events.append(
                _solve.settled_event(rhs, t_off, axon_col, stop_on_spike, n_v)
            )
    return events
