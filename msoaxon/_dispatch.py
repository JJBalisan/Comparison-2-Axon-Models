"""One entry point for both models, used by the analysis code (model="multi" / "two")."""

from .multi import mso_axon
from .two import two_cpt


def run_model(model, stim, start, stop, I, node, t_end, v0, input_node, syn=None,
              model_type="active-full", mem=None, morph=None, input_node2=None, **kw):
    """mso_axon (model="multi") or two_cpt (model="two") with the shared arguments.

    mem, morph and input_node2 exist only in the multi-compartment model, so passing
    them with model="two" raises instead of being dropped. Other keywords (max_step,
    stop_on_spike, r1, tau_est, ...) go to the model unchanged.
    """
    if model == "multi":
        return mso_axon(stim, start, stop, I, node, model_type, t_end, v0, input_node, syn,
                        mem=mem, morph=morph, input_node2=input_node2, **kw)
    if model != "two":
        raise ValueError(f'model must be "multi" or "two", got {model!r}')
    given = [k for k, v in (("mem", mem), ("morph", morph), ("input_node2", input_node2))
             if v is not None]
    if given:
        raise ValueError(f"{', '.join(given)} only apply to the multi-compartment model")
    return two_cpt(stim, start, stop, I, node, model_type, t_end, v0, input_node, syn, **kw)


# spikes() stops runs once the input is over and no spike can follow (A1 in the
# search-speed plan); False restores running every non-spiking run to t_end
STOP_WHEN_SETTLED = True


def response(model, stim, start, stop, I, node, t_end, v0, input_node, syn=None, *, factor,
             **kw):
    """(spiked, peak): whether compartment `node` rose `factor` mV above the soma
    before t_end (see spikes()), and the largest axon - soma of the run [mV].

    For a run that didn't spike the peak says how close it came, which lets a
    search interpolate towards the threshold instead of halving (_bisect.crossing);
    for one that did, the run stopped at the crossing, so the peak is ~factor.
    The settled stop keeps the peak: it only ends runs whose response has peaked.
    """
    t, x = run_model(model, stim, start, stop, I, node, t_end, v0, input_node, syn,
                     stop_on_spike=factor, stop_when_settled=STOP_WHEN_SETTLED, **kw)
    axon = node - 1 if model == "multi" else 1
    margin = x[:, axon] - x[:, 0]
    return bool(t[-1] < t_end and margin[-1] > 0.75 * factor), float(margin.max())


def spikes(model, stim, start, stop, I, node, t_end, v0, input_node, syn=None, *, factor, **kw):
    """True if compartment `node` rises `factor` mV above the soma before t_end.

    The run stops at that first spike (stop_on_spike) or, once the input is over,
    as soon as no spike can follow (stop_when_settled). Both end before t_end, so
    the final state tells them apart: a spike stop ends with axon - soma at
    `factor`, a settled stop at or below factor / 2 (see _solve.settled_event).
    The test sits at 0.75 * factor, away from both: the settled event can fire
    exactly as axon - soma falls through factor / 2, so testing at factor / 2
    itself misread three near-threshold ramp runs as spikes.
    """
    t, x = run_model(model, stim, start, stop, I, node, t_end, v0, input_node, syn,
                     stop_on_spike=factor, stop_when_settled=STOP_WHEN_SETTLED, **kw)
    axon = node - 1 if model == "multi" else 1
    return bool(t[-1] < t_end and x[-1, axon] - x[-1, 0] > 0.75 * factor)
