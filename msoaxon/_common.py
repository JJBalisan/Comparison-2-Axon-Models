"""Stimulus and argument handling shared by both models (from TwoCpt.m)."""

from types import SimpleNamespace

from . import constants as C
from .synaptic import synaptic

STIM_TYPES = ("step", "ramp", "ramp2", "sine", "Synaptic", "SynapticPair", "EPSG", "EPSGpair")
MODEL_TYPES = ("passive", "active-KLT", "active-H", "active-KLT+H", "active-sodium",
               "Active-sodium", "active-KHT", "active-full")


def check_args(stim_type, model_type, node, input_node, min_node, stim_types=STIM_TYPES,
               n_max=C.N_CPT):
    """Reject inputs that numpy's negative indexing would otherwise accept silently."""
    if stim_type not in stim_types:
        raise ValueError(f"unknown stimType {stim_type!r}")
    if model_type not in MODEL_TYPES:
        raise ValueError(f"unknown model type {model_type!r}; expected one of {MODEL_TYPES}")
    if not min_node <= node <= C.N_CPT:
        raise ValueError(f"node must be {min_node}..{C.N_CPT} (1-indexed), got {node}")
    if not 1 <= input_node <= n_max:
        raise ValueError(f"input_node must be 1..{n_max} (1-indexed), got {input_node}")


def stimulus(stim_type, start, stop, I, t_end, syn):
    """Bundle stimulus settings (what TwoCpt.m stored on P)."""
    s = SimpleNamespace(start=start, stop=stop, I=I, t_end=t_end, epsg_tau=tuple(syn.epsg_tau))
    if stim_type == "sine":
        s.f = syn.f
    if stim_type in ("Synaptic", "SynapticPair"):
        s.t_syn, s.g_syn = synaptic(syn)
        s.VsynE, s.diff = syn.VsynE, syn.diff
    return s
