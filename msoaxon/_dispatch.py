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


def spikes(model, stim, start, stop, I, node, t_end, v0, input_node, syn=None, *, factor, **kw):
    """True if compartment `node` rises `factor` mV above the soma before t_end.

    The run stops at that first spike (stop_on_spike), so it ends early exactly
    when there is one.
    """
    t, _ = run_model(model, stim, start, stop, I, node, t_end, v0, input_node, syn,
                     stop_on_spike=factor, **kw)
    return bool(t[-1] < t_end)
