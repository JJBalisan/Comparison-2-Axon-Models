"""One entry point for both models, used by the analysis code.

The analysis functions pick the model with model="multi" / "two".
"""

import types
from typing import Any

import numpy as np

from msoaxon import multi, synaptic, two


def run_model(
    model: str,
    stim: str,
    start: float,
    stop: float,
    I: float,
    node: int,
    t_end: float,
    v0: float,
    input_node: int,
    syn: synaptic.SynParams | None = None,
    model_type: str = "active-full",
    mem: types.SimpleNamespace | None = None,
    morph: types.SimpleNamespace | None = None,
    input_node2: int | None = None,
    mode: str | None = None,
    **kw: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Runs mso_axon (model="multi") or two_cpt (model="two").

    Args:
        model: "multi" or "two".
        stim: The stimulus name (see mso_axon).
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        I: Stimulus amplitude.
        node: The compared compartment, 1-indexed.
        t_end: End of the run [ms].
        v0: Initial (and default leak reversal) potential [mV].
        input_node: Where the stimulus enters, 1-indexed.
        syn: Synaptic waveform settings (SynParams), or None.
        model_type: Which channels are active.
        mem: Channel overrides from multi.membrane(); multi only.
        morph: A multi-compartment morphology; multi only.
        input_node2: Second input site for "EPSGbilateral"; multi only.
        mode: "fast", "exact" or None, as in mso_axon; multi only (the
            two-compartment model has only its numpy path, so it ignores
            mode, and a search over both models can pass it in model_kw).
        **kw: Further model keywords (max_step, stop_on_spike, r1,
            tau_est, ...), passed on unchanged.

    Returns:
        The model's (t, x).

    Raises:
        ValueError: If model is neither "multi" nor "two", or if mem,
            morph or input_node2 is given with model="two" (they exist
            only in the multi-compartment model, so they are not
            silently dropped).
    """
    if model == "multi":
        return multi.mso_axon(
            stim,
            start,
            stop,
            I,
            node,
            model_type,
            t_end,
            v0,
            input_node,
            syn,
            mem=mem,
            morph=morph,
            input_node2=input_node2,
            mode=mode,
            **kw,
        )
    if model != "two":
        raise ValueError(f'model must be "multi" or "two", got {model!r}')
    given = [
        k
        for k, v in (
            ("mem", mem),
            ("morph", morph),
            ("input_node2", input_node2),
        )
        if v is not None
    ]
    if given:
        raise ValueError(
            f"{', '.join(given)} only apply to the multi-compartment model"
        )
    return two.two_cpt(
        stim, start, stop, I, node, model_type, t_end, v0, input_node, syn, **kw
    )


# spikes() stops runs once the input is over and no spike can follow (A1 in the
# search-speed plan); False restores running every non-spiking run to t_end
STOP_WHEN_SETTLED = True


def response(
    model: str,
    stim: str,
    start: float,
    stop: float,
    I: float,
    node: int,
    t_end: float,
    v0: float,
    input_node: int,
    syn: synaptic.SynParams | None = None,
    *,
    factor: float,
    **kw: Any,
) -> tuple[bool, float]:
    """Runs one stimulus and reports whether it spiked and how close it came.

    For a run that didn't spike the peak says how close it came, which
    lets a search interpolate towards the threshold instead of halving
    (_bisect.crossing); for one that did, the run stopped at the
    crossing, so the peak is ~factor. The settled stop keeps the peak:
    it only ends runs whose response has peaked.

    Args:
        model: "multi" or "two".
        stim: The stimulus name.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        I: Stimulus amplitude.
        node: The compartment compared with the soma, 1-indexed.
        t_end: End of the run [ms].
        v0: Initial potential [mV].
        input_node: Where the stimulus enters, 1-indexed.
        syn: Synaptic waveform settings (SynParams), or None.
        factor: The spike criterion: axon - soma [mV].
        **kw: Further run_model keywords.

    Returns:
        (spiked, peak): whether compartment `node` rose `factor` mV
        above the soma before t_end (as in spikes()), and the largest
        axon - soma of the run [mV].
    """
    t, x = run_model(
        model,
        stim,
        start,
        stop,
        I,
        node,
        t_end,
        v0,
        input_node,
        syn,
        stop_on_spike=factor,
        stop_when_settled=STOP_WHEN_SETTLED,
        **kw,
    )
    axon = node - 1 if model == "multi" else 1
    margin = x[:, axon] - x[:, 0]
    return bool(t[-1] < t_end and margin[-1] > 0.75 * factor), float(
        margin.max()
    )


def spikes(
    model: str,
    stim: str,
    start: float,
    stop: float,
    I: float,
    node: int,
    t_end: float,
    v0: float,
    input_node: int,
    syn: synaptic.SynParams | None = None,
    *,
    factor: float,
    **kw: Any,
) -> bool:
    """Returns whether compartment `node` rises `factor` mV above the soma.

    The run stops at that first spike (stop_on_spike) or, once the input
    is over, as soon as no spike can follow (stop_when_settled). Both
    end before t_end, so the final state tells them apart: a spike stop
    ends with axon - soma at `factor`, a settled stop at or below
    factor / 2 (see _solve.settled_event). The test sits at
    0.75 * factor, away from both: the settled event can fire exactly as
    axon - soma falls through factor / 2, so testing at factor / 2
    itself misread three near-threshold ramp runs as spikes.

    Args:
        model: "multi" or "two".
        stim: The stimulus name.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        I: Stimulus amplitude.
        node: The compartment compared with the soma, 1-indexed.
        t_end: End of the run [ms].
        v0: Initial potential [mV].
        input_node: Where the stimulus enters, 1-indexed.
        syn: Synaptic waveform settings (SynParams), or None.
        factor: The spike criterion: axon - soma [mV].
        **kw: Further run_model keywords.

    Returns:
        True if the run spiked before t_end.
    """
    t, x = run_model(
        model,
        stim,
        start,
        stop,
        I,
        node,
        t_end,
        v0,
        input_node,
        syn,
        stop_on_spike=factor,
        stop_when_settled=STOP_WHEN_SETTLED,
        **kw,
    )
    axon = node - 1 if model == "multi" else 1
    return bool(t[-1] < t_end and x[-1, axon] - x[-1, 0] > 0.75 * factor)


# the fields of multi.membrane()'s cache key, in order
_MEMBRANE_KEY = (
    "v0",
    "soma_klt_scale",
    "ais_klt_scale",
    "soma_na_vhalf",
    "rebalance_rest",
    "dendrite_klt_scale",
)


def describe(
    model: str, node: int, v0: float, mem: object = None, **settings: object
) -> str:
    """Names a model configuration, for error messages.

    Args:
        model: "multi" or "two".
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        mem: Channel overrides from multi.membrane(), or None.
        **settings: Further settings to name, in order.

    Returns:
        E.g. "'multi', node=3, v0=-68.0, mem=membrane(v0=-68.0, ...)".
    """
    parts = [repr(model), f"node={node}", f"v0={v0}"]
    parts += [f"{k}={v!r}" for k, v in settings.items()]
    if mem is not None:
        knobs = ", ".join(
            f"{k}={v}" for k, v in zip(_MEMBRANE_KEY, mem.key, strict=True)
        )
        parts.append(f"mem=membrane({knobs})")
    return ", ".join(parts)
