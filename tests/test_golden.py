"""Golden traces: soma and axon voltage pinned for every stimulus and model type.

Each case runs one model once and samples the soma and the compared axon
compartment at fixed times (linear interpolation between solver steps). The
values were generated from the code as it stood before the Google-style
refactor and are compared at rtol 1e-9, so any change to the arithmetic, the
solver settings, the breakpoints or the pre-stimulus cache shows up here. The
MATLAB quirks kept on purpose (multi.py's table above external_current) are
pinned through these cases rather than tested one by one.

Regenerate only when a change in results is intended:

    uv run python tests/test_golden.py
"""

import json
import pathlib

import numpy as np
import pytest

from msoaxon import SynParams, mso_axon, two_cpt
from msoaxon import multi
from msoaxon import two

GOLDEN = pathlib.Path(__file__).resolve().parent / "golden_traces.json"
TIMES = (1.0, 5.1, 5.3, 5.6, 6.0, 7.0, 9.0, 12.0)
DTX = dict(soma_klt_scale=0, ais_klt_scale=0, rebalance_rest=True)

# (stimType, start, stop, I, t_end, SynParams keywords); active-full, node 3,
# v0 -68 mV, input at the soma unless a case below says otherwise
STIMULI = {
    "step": ("step", 5, 10, 8500, 20, {}),
    "step-subthreshold": ("step", 5, 10, 500, 20, {}),
    "ramp": ("ramp", 5, 5.5, 8500, 20, {}),
    "ramp2": ("ramp2", 5, 5.5, 8500, 20, {}),
    "sine": ("sine", 5, 7.5, 8500, 20, {}),
    "Synaptic": ("Synaptic", 5, 15, 0, 20, dict(random_in=1396)),
    "SynapticPair": ("SynapticPair", 5, 15, 0, 20, dict(diff=1.0)),
    "EPSG": ("EPSG", 5, 10, 108, 20, {}),
    "EPSG-subthreshold": ("EPSG", 5, 10, 26.7, 20, {}),
    "EPSGpair": ("EPSGpair", 5, 5.2, 100, 20, {}),
    "EPSGpair-subthreshold": ("EPSGpair", 5, 5.3, 30, 20, {}),
}
MULTI_TYPES = ("passive", "active-KLT", "active-H", "active-KLT+H",
               "active-sodium", "active-KHT", "active-full")
TWO_TYPES = MULTI_TYPES + ("Active-sodium",)  # TwoCpt.m's capitalised spelling


def _run(model, stim, start, stop, amp, t_end, syn_kw, model_type="active-full",
         node=3, input_node=1, **kw):
    syn = SynParams(t_end=t_end, **syn_kw)
    fn = mso_axon if model == "multi" else two_cpt
    return fn(stim, start, stop, amp, node, model_type, t_end, -68.0, input_node,
              syn, **kw)


def cases():
    """name -> (model, keyword arguments for _run, axon column)."""
    out = {}
    for model, axon in (("multi", 2), ("two", 1)):
        for name, (stim, start, stop, amp, t_end, syn_kw) in STIMULI.items():
            out[f"{model}/{name}"] = (model, dict(
                stim=stim, start=start, stop=stop, amp=amp, t_end=t_end,
                syn_kw=syn_kw), axon)
        for mt in MULTI_TYPES if model == "multi" else TWO_TYPES:
            out[f"{model}/type/{mt}"] = (model, dict(
                stim="EPSG", start=5, stop=10, amp=108, t_end=15, syn_kw={},
                model_type=mt), axon)
        # input into the axon: two_cpt routes it to compartment 2, mso_axon to
        # input_node, except step, which mso_axon always puts into the soma
        for stim, amp in (("step", 8500), ("EPSGpair", 100)):
            out[f"{model}/{stim}-input-node-3"] = (model, dict(
                stim=stim, start=5, stop=5.2 if stim == "EPSGpair" else 10,
                amp=amp, t_end=20, syn_kw={}, input_node=3), axon)
    d = multi.with_dendrites()
    lat, med = int(d.lateral[2]), int(d.medial[2])
    out["multi/EPSGbilateral"] = ("multi", dict(
        stim="EPSGbilateral", start=5, stop=5.2, amp=60, t_end=20, syn_kw={},
        input_node2=3), 2)
    out["multi/dendrites/EPSGbilateral"] = ("multi", dict(
        stim="EPSGbilateral", start=5, stop=5.1, amp=60, t_end=15, syn_kw={},
        input_node=lat, input_node2=med, morph=d), 2)
    out["multi/dendrites/EPSG-klt-from-soma-ra200"] = ("multi", dict(
        stim="EPSG", start=5, stop=10, amp=60, t_end=15, syn_kw={},
        input_node=lat, morph=multi.with_dendrites(conserve_totals=False,
                                                   dendrite_ra=200.0)), 2)
    out["multi/membrane/dendrotoxin"] = ("multi", dict(
        stim="step", start=5, stop=10, amp=1500, t_end=15, syn_kw={},
        mem=multi.membrane(-68.0, **DTX)), 2)
    out["multi/membrane/soma-na-vhalf"] = ("multi", dict(
        stim="EPSGpair", start=5, stop=5.2, amp=60, t_end=15, syn_kw={},
        mem=multi.membrane(-68.0, soma_klt_scale=3, soma_na_vhalf=-77,
                           rebalance_rest=True)), 2)
    goldwyn = {k: v for k, v in two.GOLDWYN_2019.items() if k != "v0"}
    out["two/goldwyn-2019"] = ("two", dict(
        stim="EPSGpair", start=5, stop=5.2, amp=60, t_end=15, syn_kw={},
        **goldwyn), 1)
    out["two/node-5"] = ("two", dict(
        stim="EPSGpair", start=5, stop=5.2, amp=60, t_end=15, syn_kw={},
        node=5), 1)
    return out


CASES = cases()


def sample(name):
    """[t_end reached, soma samples, axon samples] for one case."""
    model, kw, axon = CASES[name]
    t, x = _run(model, **kw)
    times = [s for s in TIMES if s <= t[-1]]
    return {"t_last": float(t[-1]),
            "soma": [float(v) for v in np.interp(times, t, x[:, 0])],
            "axon": [float(v) for v in np.interp(times, t, x[:, axon])]}


@pytest.fixture(scope="module")
def golden():
    return json.loads(GOLDEN.read_text())


def test_golden_covers_every_case(golden):
    assert sorted(golden) == sorted(CASES)


@pytest.mark.parametrize("name", sorted(CASES))
def test_golden_trace(golden, name):
    got, want = sample(name), golden[name]
    assert got["t_last"] == want["t_last"]
    np.testing.assert_allclose(got["soma"], want["soma"], rtol=1e-9, atol=0)
    np.testing.assert_allclose(got["axon"], want["axon"], rtol=1e-9, atol=0)


def test_step_ignores_input_node_in_mso_axon():
    # msoAxon.m injects the step into the soma whatever input_node says
    kw = dict(stim="step", start=5, stop=10, amp=8500, t_end=20, syn_kw={})
    t1, y1 = _run("multi", **kw, input_node=1)
    t2, y2 = _run("multi", **kw, input_node=3)
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


def test_step_follows_input_node_in_two_cpt():
    # TwoCptODE.m does route it: input_node != 1 means the axon compartment
    kw = dict(stim="step", start=5, stop=10, amp=8500, t_end=20, syn_kw={})
    _, x1 = _run("two", **kw, input_node=1)
    _, x3 = _run("two", **kw, input_node=3)
    assert not np.array_equal(x1, x3)


if __name__ == "__main__":
    GOLDEN.write_text(json.dumps({n: sample(n) for n in sorted(CASES)},
                                 indent=1) + "\n")
    print(f"wrote {len(CASES)} cases to {GOLDEN}")
