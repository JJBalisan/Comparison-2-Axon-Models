"""Stimulus branches, input validation, synaptic input, threshold search, plotting, scripts."""

import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from msoaxon import SynParams, mso_axon, two_cpt
from msoaxon.plotting import graphing, threshold_figure
from msoaxon.spiking import spiking
from msoaxon.synaptic import synaptic
from msoaxon.threshold import binary_search

REPO = Path(__file__).resolve().parents[1]
Q = 150 / 128  # EPSGpair search grid

# (stimType, start, stop, I) that make both models fire at node 3
FIRING = {
    "step": (5, 10, 8500),
    "ramp": (5, 5.5, 8500),
    "ramp2": (5, 5.5, 8500),
    "sine": (5, 7.5, 8500),
    "EPSG": (5, 10, 108),
    "EPSGpair": (5, 5, 100),
}


@pytest.mark.parametrize("model", [two_cpt, mso_axon])
@pytest.mark.parametrize("stim", sorted(FIRING))
def test_every_stimulus_fires(model, stim):
    start, stop, I = FIRING[stim]
    t, x = model(stim, start, stop, I, 3, "active-full", 20, -68, 1)
    axon = 1 if model is two_cpt else 2
    assert np.isfinite(x).all() and t[0] == 0 and t[-1] == 20
    assert (x[:, axon] - x[:, 0]).max() > 30


@pytest.mark.parametrize("model", [two_cpt, mso_axon])
@pytest.mark.parametrize("stim", ["Synaptic", "SynapticPair"])
def test_synaptic_inputs_depolarise(model, stim):
    syn = SynParams(t_end=20, random_in=1396, diff=1.0)
    _, x = model(stim, 5, 15, 0, 3, "active-full", 20, -68, 1, syn)
    assert np.isfinite(x).all()
    assert x[:, 0].max() > -67  # noise-driven EPSPs lift the soma above rest


def test_synaptic_pair_near_t_end_stays_finite():
    # t + diff runs past the conductance table; MATLAB's interp1q gave NaN here
    syn = SynParams(t_end=20, diff=3.0)
    for model in (two_cpt, mso_axon):
        _, x = model("SynapticPair", 5, 20, 0, 3, "active-full", 20, -68, 1, syn)
        assert np.isfinite(x).all()


def test_synaptic_conductance_shape_and_seed():
    t, g = synaptic(SynParams(t_end=20, random_in=7))
    assert len(t) == len(g) == 2001 and np.isclose(t[-1], 20)
    assert (g >= 0).all() and g.max() > 0
    np.testing.assert_array_equal(g, synaptic(SynParams(t_end=20, random_in=7))[1])
    assert not np.array_equal(g, synaptic(SynParams(t_end=20, random_in=8))[1])


def test_input_to_axon_compartment():
    # input_node != 1 routes the 2-CPT input into compartment 2
    _, x = two_cpt("EPSGpair", 5, 5, 100, 3, "active-full", 20, -68, 3)
    assert x[:, 1].max() > x[:, 0].max()


@pytest.mark.parametrize("model,node,input_node", [
    (two_cpt, 1, 1),    # node 1 is the soma; would read coupling[-1]
    (two_cpt, 46, 1),
    (two_cpt, 3, 0),
    (mso_axon, 3, 0),   # would inject into compartment 45
    (mso_axon, 3, 46),
])
def test_out_of_range_nodes_are_rejected(model, node, input_node):
    with pytest.raises(ValueError, match="node"):
        model("EPSGpair", 5, 5, 100, node, "active-full", 20, -68, input_node)


@pytest.mark.parametrize("model", [two_cpt, mso_axon])
def test_unknown_names_are_rejected(model):
    with pytest.raises(ValueError, match="model type"):
        model("EPSGpair", 5, 5, 100, 3, "active-ful", 20, -68, 1)
    with pytest.raises(ValueError, match="stimType"):
        model("epsgpair", 5, 5, 100, 3, "active-full", 20, -68, 1)


def test_every_model_type_runs():
    for mt in ("passive", "active-KLT", "active-H", "active-KLT+H", "active-sodium",
               "active-KHT", "active-full"):
        for model in (two_cpt, mso_axon):
            _, x = model("EPSG", 5, 10, 50, 3, mt, 10, -68, 1)
            assert np.isfinite(x).all(), (model.__name__, mt)


def test_spiking_reset_carries_across_columns_mult():
    # column 2 ends above threshold; with the carried-over reset, column 3's first
    # excursion is not marked until its difference drops to <= 0
    x = np.zeros((4, 315))
    x[:, 1] = [0, 50, 50, 50]
    x[:, 2] = [50, 50, 0, 50]
    s = spiking(x, 30, "Mult")
    np.testing.assert_array_equal(s[:, 1], [0, 1, 0, 0])
    np.testing.assert_array_equal(s[:, 2], [0, 0, 0, 1])


def test_epsgpair_thresholds_match_repo_jpg():
    """Delays 0 and 0.3 ms at node 3: both models read 41q and 55q off EPSGpair_Thresholds_CPT3.jpg."""
    multi, two = binary_search("EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, rounded=False)
    assert [round(v / Q) for v in multi] == [41, 55]
    assert [round(v / Q) for v in two] == [41, 55]


def test_binary_search_rounding_and_ceiling():
    multi, two = binary_search("step", 1, 3, 30, 15000)
    assert all(v % 10 == 0 and 0 < v < 15000 for v in multi + two)
    # a stimulus that never fires reports the ceiling, as BinarySearch.m does
    multi, _ = binary_search("EPSG", 1, 3, 30, 15000)
    assert multi == [15000]


@pytest.fixture(scope="module")
def runs():
    syn = SynParams(t_end=20)
    args = ("EPSGpair", 5, 5, 100)
    t1, y = mso_axon(*args, 3, "active-full", 20, -68, 1, syn)
    t2, x = two_cpt(*args, 3, "active-full", 20, -68, 1, syn)
    t3, z = two_cpt(*args, 5, "active-full", 20, -68, 1, syn)
    graph = dict(node=3, tEnd=20, stimType="EPSGpair", inputNode=1, start=5, stop=5,
                 I=100, Syn=syn)
    return graph, (t1, y, t2, x, t3, z)


@pytest.mark.parametrize("kind", ["Contour", "ModelComparison", "ModelComparison2",
                                  "NodeComparison", "MultiCompartment", "Input", "Multigraph"])
def test_graph_types_render(runs, kind):
    graph, data = runs
    fig = graphing({**graph, "Type": kind}, *data)
    assert fig.axes
    matplotlib.pyplot.close(fig)


def test_two_compartment_graph_and_unknown_type(runs):
    graph, (_, _, t2, x, _, _) = runs
    matplotlib.pyplot.close(graphing({**graph, "Type": "TwoCompartment"}, t2, x))
    with pytest.raises(ValueError):
        graphing({**graph, "Type": "Nope"}, t2, x)
    matplotlib.pyplot.close(threshold_figure([0, 1], [1, 2], [1, 2], "t", "x"))


def test_combine_all_script(tmp_path):
    out = tmp_path / "fig.png"
    r = subprocess.run([sys.executable, str(REPO / "scripts" / "combine_all.py"),
                        "--stim", "EPSGpair", "--I", "100", "--out", str(out)],
                       capture_output=True, text=True, env={"MPLBACKEND": "Agg", "PATH": ""})
    assert r.returncode == 0, r.stderr
    assert "multi-compartment 1, two-compartment 1" in r.stdout and out.stat().st_size > 0


def test_parallel_search_matches_serial():
    serial = binary_search("EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, workers=1)
    assert binary_search("EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, workers=2) == serial


@pytest.mark.parametrize("model", [two_cpt, mso_axon])
def test_pre_stimulus_cache_does_not_change_results(model):
    from msoaxon import _solve
    _solve._QUIET_CACHE.clear()
    t1, x1 = model("EPSGpair", 5, 5.3, 70, 3, "active-full", 20, -68, 1)
    assert _solve._QUIET_CACHE  # first run filled it
    t2, x2 = model("EPSGpair", 5, 5.3, 90, 3, "active-full", 20, -68, 1)  # reuses it
    _solve._QUIET_CACHE.clear()
    t3, x3 = model("EPSGpair", 5, 5.3, 90, 3, "active-full", 20, -68, 1)  # recomputes
    np.testing.assert_array_equal(t2, t3)
    np.testing.assert_array_equal(x2, x3)


@pytest.mark.parametrize("model,axon", [(two_cpt, 1), (mso_axon, 2)])
def test_stop_on_spike_truncates_the_full_run(model, axon):
    args = ("EPSGpair", 5, 5.3, 70, 3, "active-full", 20, -68, 1)
    t, x = model(*args)
    te, xe = model(*args, stop_on_spike=10)
    assert te[-1] < 20 and np.isclose(xe[-1, axon] - xe[-1, 0], 10)
    # same accepted steps as the full run until the step where the spike is found
    np.testing.assert_array_equal(te[:-1], t[:len(te) - 1])
    np.testing.assert_array_equal(xe[:-1], x[:len(te) - 1])
    # and the first sample above the factor in the full run comes right after
    first = np.argmax(x[:, axon] - x[:, 0] > 10)
    assert t[first - 1] <= te[-1] <= t[first]


def test_stop_on_spike_leaves_non_spiking_runs_alone():
    args = ("EPSGpair", 5, 5, 40, 3, "active-full", 20, -68, 1)
    t, x = mso_axon(*args)
    te, xe = mso_axon(*args, stop_on_spike=10)
    np.testing.assert_array_equal(t, te)
    np.testing.assert_array_equal(x, xe)
