"""Stimulus branches, input validation and synaptic input.

Also the threshold search, plotting and the combine_all script.
"""

import pathlib
import runpy
import subprocess
import sys

import matplotlib as mpl

mpl.use("Agg")

import numpy as np
import pytest

import msoaxon
from msoaxon import (
    _bisect,
    _dispatch,
    _parallel,
    _solve,
    coincidence,
    plotting,
    spiking,
    synaptic,
    threshold as threshold_mod,
)

REPO = pathlib.Path(__file__).resolve().parents[1]
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


@pytest.mark.parametrize("model", [msoaxon.two_cpt, msoaxon.mso_axon])
@pytest.mark.parametrize("stim", sorted(FIRING))
def test_every_stimulus_fires(model, stim):
    start, stop, I = FIRING[stim]
    t, x = model(stim, start, stop, I, 3, "active-full", 20, -68, 1)
    axon = 1 if model is msoaxon.two_cpt else 2
    assert np.isfinite(x).all() and t[0] == 0 and t[-1] == 20
    assert (x[:, axon] - x[:, 0]).max() > 30


@pytest.mark.parametrize("model", [msoaxon.two_cpt, msoaxon.mso_axon])
@pytest.mark.parametrize("stim", ["Synaptic", "SynapticPair"])
def test_synaptic_inputs_depolarise(model, stim):
    syn = msoaxon.SynParams(t_end=20, random_in=1396, diff=1.0)
    _, x = model(stim, 5, 15, 0, 3, "active-full", 20, -68, 1, syn)
    assert np.isfinite(x).all()
    assert x[:, 0].max() > -67  # noise-driven EPSPs lift the soma above rest


def test_synaptic_pair_near_t_end_stays_finite():
    # t + diff runs past the conductance table; MATLAB's interp1q gave NaN here
    syn = msoaxon.SynParams(t_end=20, diff=3.0)
    for model in (msoaxon.two_cpt, msoaxon.mso_axon):
        _, x = model(
            "SynapticPair", 5, 20, 0, 3, "active-full", 20, -68, 1, syn
        )
        assert np.isfinite(x).all()


def test_synaptic_conductance_shape_and_seed():
    t, g = synaptic.synaptic(msoaxon.SynParams(t_end=20, random_in=7))
    assert len(t) == len(g) == 2001 and np.isclose(t[-1], 20)
    assert (g >= 0).all() and g.max() > 0
    np.testing.assert_array_equal(
        g, synaptic.synaptic(msoaxon.SynParams(t_end=20, random_in=7))[1]
    )
    assert not np.array_equal(
        g, synaptic.synaptic(msoaxon.SynParams(t_end=20, random_in=8))[1]
    )


def test_input_to_axon_compartment():
    # input_node != 1 routes the 2-CPT input into compartment 2
    _, x = msoaxon.two_cpt("EPSGpair", 5, 5, 100, 3, "active-full", 20, -68, 3)
    assert x[:, 1].max() > x[:, 0].max()


@pytest.mark.parametrize(
    "model,node,input_node",
    [
        (msoaxon.two_cpt, 1, 1),  # node 1 is the soma; would read coupling[-1]
        (msoaxon.two_cpt, 46, 1),
        (msoaxon.two_cpt, 3, 0),
        (msoaxon.mso_axon, 3, 0),  # would inject into compartment 45
        (msoaxon.mso_axon, 3, 46),
    ],
)
def test_out_of_range_nodes_are_rejected(model, node, input_node):
    with pytest.raises(ValueError, match="node"):
        model("EPSGpair", 5, 5, 100, node, "active-full", 20, -68, input_node)


@pytest.mark.parametrize("model", [msoaxon.two_cpt, msoaxon.mso_axon])
def test_unknown_names_are_rejected(model):
    with pytest.raises(ValueError, match="model type"):
        model("EPSGpair", 5, 5, 100, 3, "active-ful", 20, -68, 1)
    with pytest.raises(ValueError, match="stimType"):
        model("epsgpair", 5, 5, 100, 3, "active-full", 20, -68, 1)


def test_every_model_type_runs():
    for mt in (
        "passive",
        "active-KLT",
        "active-H",
        "active-KLT+H",
        "active-sodium",
        "active-KHT",
        "active-full",
    ):
        for model in (msoaxon.two_cpt, msoaxon.mso_axon):
            _, x = model("EPSG", 5, 10, 50, 3, mt, 10, -68, 1)
            assert np.isfinite(x).all(), (model.__name__, mt)


def test_spiking_reset_carries_across_columns_mult():
    # column 2 ends above threshold; with the carried-over reset, column 3's
    # first excursion is not marked until its difference drops to <= 0
    x = np.zeros((4, 315))
    x[:, 1] = [0, 50, 50, 50]
    x[:, 2] = [50, 50, 0, 50]
    s = spiking.spiking(x, 30, "Mult")
    np.testing.assert_array_equal(s[:, 1], [0, 1, 0, 0])
    np.testing.assert_array_equal(s[:, 2], [0, 0, 0, 1])


# EPSGpair thresholds on the grid of the repo's EPSGpair_Thresholds_CPT*.jpg
# (making_threshold_graphs.py --only EPSGpair --jpg-grid), delays 0-0.6 ms.
# The halving search only returns multiples of Q, so its output is pinned
# exactly: the closest true threshold to any probe is 8.5e-4 relative away,
# against ~1e-14 drift across platforms and modes. The jpg column is the
# MATLAB value read off the plot, kept as history; 20 of 28 match and the
# other 8 are one output (2Q) away (PYTHON.md, "Validation").
# node, model, delay [ms], threshold / Q, jpg / Q, match
JPG_GRID = [
    (3, "multi", 0.0, 41, 41, True),
    (3, "multi", 0.1, 43, 43, True),
    (3, "multi", 0.2, 47, 47, True),
    (3, "multi", 0.3, 55, 55, True),
    (3, "multi", 0.4, 71, 71, True),
    (3, "multi", 0.5, 81, 81, True),
    (3, "multi", 0.6, 81, 83, False),  # the plateau, 94.9 vs 97.3
    (3, "two", 0.0, 41, 41, True),
    (3, "two", 0.1, 41, 41, True),
    (3, "two", 0.2, 47, 45, False),
    (3, "two", 0.3, 55, 55, True),
    (3, "two", 0.4, 75, 73, False),
    (3, "two", 0.5, 81, 81, True),
    (3, "two", 0.6, 81, 81, True),
    (5, "multi", 0.0, 39, 39, True),
    (5, "multi", 0.1, 39, 41, False),
    (5, "multi", 0.2, 43, 43, True),
    (5, "multi", 0.3, 51, 51, True),
    (5, "multi", 0.4, 65, 65, True),
    (5, "multi", 0.5, 77, 77, True),
    (5, "multi", 0.6, 77, 77, True),
    (5, "two", 0.0, 37, 37, True),
    (5, "two", 0.1, 39, 39, True),
    (5, "two", 0.2, 41, 43, False),  # one step low at four delays
    (5, "two", 0.3, 49, 51, False),
    (5, "two", 0.4, 67, 69, False),
    (5, "two", 0.5, 73, 75, False),
    (5, "two", 0.6, 75, 75, True),
]


def test_jpg_grid_table_records_the_matlab_comparison():
    assert len(JPG_GRID) == 28
    assert all(match == (q == jpg) for *_, q, jpg, match in JPG_GRID)
    assert sum(match for *_, match in JPG_GRID) == 20
    assert all(abs(q - jpg) in (0, 2) for *_, q, jpg, _ in JPG_GRID)


@pytest.mark.slow
@pytest.mark.parametrize("node", [3, 5])
def test_epsgpair_jpg_grid_thresholds(node):
    script = runpy.run_path(
        str(REPO / "scripts" / "making_threshold_graphs.py")
    )
    _, kw = script["epsgpair_grid"](jpg_grid=True)
    _, factor, max_I, *_ = script["SWEEPS"]["EPSGpair"]
    rows = [r for r in JPG_GRID if r[0] == node]
    n = len(rows) // 2  # the jpg grid's first 7 points
    multi, two = threshold_mod.binary_search(
        "EPSGpair", n, node, factor, max_I, workers=4, **kw
    )
    got = [
        (node, model, round(d, 1), v / Q)
        for model, values in (("multi", multi), ("two", two))
        for d, v in zip(
            (
                threshold_mod.sweep_x("EPSGpair", i, kw["epsg_pair_dt"])
                for i in range(1, n + 1)
            ),
            values,
            strict=True,
        )
    ]
    assert got == [r[:4] for r in rows]


def test_binary_search_rounding_and_ceiling():
    multi, two = threshold_mod.binary_search("step", 1, 3, 30, 15000, workers=1)
    assert all(v % 10 == 0 and 0 < v < 15000 for v in multi + two)
    # a stimulus that never fires reports the ceiling, as BinarySearch.m does
    multi, _ = threshold_mod.binary_search("EPSG", 1, 3, 30, 15000, workers=1)
    assert multi == [15000]


@pytest.fixture(scope="module")
def runs():
    syn = msoaxon.SynParams(t_end=20)
    args = ("EPSGpair", 5, 5, 100)
    t1, y = msoaxon.mso_axon(*args, 3, "active-full", 20, -68, 1, syn)
    t2, x = msoaxon.two_cpt(*args, 3, "active-full", 20, -68, 1, syn)
    t3, z = msoaxon.two_cpt(*args, 5, "active-full", 20, -68, 1, syn)
    graph = {
        "node": 3,
        "tEnd": 20,
        "stimType": "EPSGpair",
        "inputNode": 1,
        "start": 5,
        "stop": 5,
        "I": 100,
        "Syn": syn,
    }
    return graph, (t1, y, t2, x, t3, z)


@pytest.mark.parametrize(
    "kind",
    [
        "Contour",
        "ModelComparison",
        "ModelComparison2",
        "NodeComparison",
        "MultiCompartment",
        "Input",
        "Multigraph",
    ],
)
def test_graph_types_render(runs, kind):
    graph, data = runs
    fig = plotting.graphing({**graph, "Type": kind}, *data)
    assert fig.axes
    mpl.pyplot.close(fig)


def test_two_compartment_graph_and_unknown_type(runs):
    graph, (_, _, t2, x, _, _) = runs
    mpl.pyplot.close(
        plotting.graphing({**graph, "Type": "TwoCompartment"}, t2, x)
    )
    with pytest.raises(ValueError):
        plotting.graphing({**graph, "Type": "Nope"}, t2, x)
    mpl.pyplot.close(
        plotting.threshold_figure([0, 1], [1, 2], [1, 2], "t", "x")
    )


def test_combine_all_script(tmp_path):
    out = tmp_path / "fig.png"
    r = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "combine_all.py"),
            "--stim",
            "EPSGpair",
            "--I",
            "100",
            "--out",
            str(out),
        ],
        capture_output=True,
        text=True,
        env={"MPLBACKEND": "Agg", "PATH": ""},
        check=False,
    )
    assert r.returncode == 0, r.stderr
    assert (
        "multi-compartment 1, two-compartment 1" in r.stdout
        and out.stat().st_size > 0
    )


@pytest.mark.parametrize("model", [msoaxon.two_cpt, msoaxon.mso_axon])
def test_pre_stimulus_cache_does_not_change_results(model):

    _solve._QUIET_CACHE.clear()
    t1, x1 = model("EPSGpair", 5, 5.3, 70, 3, "active-full", 20, -68, 1)
    assert _solve._QUIET_CACHE  # first run filled it
    t2, x2 = model(
        "EPSGpair", 5, 5.3, 90, 3, "active-full", 20, -68, 1
    )  # reuses it
    _solve._QUIET_CACHE.clear()
    t3, x3 = model(
        "EPSGpair", 5, 5.3, 90, 3, "active-full", 20, -68, 1
    )  # recomputes
    np.testing.assert_array_equal(t2, t3)
    np.testing.assert_array_equal(x2, x3)


@pytest.mark.parametrize(
    "model,axon", [(msoaxon.two_cpt, 1), (msoaxon.mso_axon, 2)]
)
def test_stop_on_spike_truncates_the_full_run(model, axon):
    args = ("EPSGpair", 5, 5.3, 70, 3, "active-full", 20, -68, 1)
    t, x = model(*args)
    te, xe = model(*args, stop_on_spike=10)
    assert te[-1] < 20 and np.isclose(xe[-1, axon] - xe[-1, 0], 10)
    # same accepted steps as the full run until the step where the spike is
    # found
    np.testing.assert_array_equal(te[:-1], t[: len(te) - 1])
    np.testing.assert_array_equal(xe[:-1], x[: len(te) - 1])
    # and the first sample above the factor in the full run comes right after
    first = np.argmax(x[:, axon] - x[:, 0] > 10)
    assert t[first - 1] <= te[-1] <= t[first]


def test_stop_on_spike_leaves_non_spiking_runs_alone():
    args = ("EPSGpair", 5, 5, 40, 3, "active-full", 20, -68, 1)
    t, x = msoaxon.mso_axon(*args)
    te, xe = msoaxon.mso_axon(*args, stop_on_spike=10)
    np.testing.assert_array_equal(t, te)
    np.testing.assert_array_equal(x, xe)


def test_model_modules_are_importable():
    # the re-exported functions used to shadow their own submodules
    from msoaxon import multi, two  # noqa: PLC0415 (the import is the test)

    assert multi.mso_axon is msoaxon.mso_axon and two.two_cpt is msoaxon.two_cpt


def test_parallel_runs_match_serial():

    serial = threshold_mod.binary_search(
        "EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, workers=1
    )
    assert (
        threshold_mod.binary_search(
            "EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, workers=2
        )
        == serial
    )
    curve = coincidence.threshold_curve("two", [0.0, 0.2], workers=1)
    with _parallel.process_pool(2) as pool:  # one pool reused across calls
        assert (
            threshold_mod.binary_search(
                "EPSGpair", 2, 3, 10, 150, epsg_pair_dt=0.3, executor=pool
            )
            == serial
        )
        np.testing.assert_array_equal(
            coincidence.threshold_curve("two", [0.0, 0.2], executor=pool), curve
        )


def test_smallest_firing():

    assert _bisect.smallest_firing(
        lambda x: x >= 3.0, 1.0, 100.0, 1e-6
    ) == pytest.approx(3.0, rel=1e-6)
    assert (
        _bisect.smallest_firing(lambda x: x >= 3.0, 1.0, 100.0, 1e-6) >= 3.0
    )  # returns a firing value
    assert _bisect.smallest_firing(lambda x: False, 1.0, 100.0, 1e-3) == np.inf
    # the ceiling bounds the doubling: 64 is tested, 128 exceeds 100, so 90 is
    # never reached
    assert (
        _bisect.smallest_firing(lambda x: x >= 90.0, 1.0, 100.0, 1e-3) == np.inf
    )


def test_run_model_rejects_multi_only_keywords_for_two():

    with pytest.raises(ValueError, match="mem only apply"):
        _dispatch.run_model(
            "two", "step", 5, 10, 100.0, 3, 10, -68.0, 1, mem=object()
        )
    with pytest.raises(ValueError, match='"multi" or "two"'):
        _dispatch.run_model("Two", "step", 5, 10, 100.0, 3, 10, -68.0, 1)


def _smooth_response(threshold, level=10.0):
    """A stand-in for an EPSG-pair search: the peak rises convexly and smoothly
    through `level` at `threshold`, like the models' axon - soma at 10 mV."""
    calls = []

    def respond(x):
        calls.append(x)
        peak = level * (x / threshold) ** 12
        return peak >= level, min(peak, level)

    return respond, calls


@pytest.mark.parametrize("threshold", [3.7, 47.72, 48.0, 190.0, 1500.0])
def test_smallest_crossing_keeps_bisections_guarantee(threshold):

    respond, calls = _smooth_response(threshold)
    hi = _bisect.smallest_crossing(respond, 50.0, 2000.0, 1e-4, 10.0)
    assert (
        threshold <= hi <= threshold * (1 + 1e-4)
    )  # spikes, and within rel_tol
    n_bisect = [0]

    def fires(x):
        n_bisect[0] += 1
        return x >= threshold

    _bisect.smallest_firing(fires, 50.0, 2000.0, 1e-4)
    assert len(calls) < n_bisect[0]  # and in fewer runs than bisection


def test_smallest_crossing_gives_up_past_the_ceiling():

    respond, _ = _smooth_response(5000.0)
    assert (
        _bisect.smallest_crossing(respond, 50.0, 2000.0, 1e-4, 10.0) == np.inf
    )


def test_crossing_works_in_either_direction():
    # window() searches delay, where the spiking side is the smaller value

    def respond(
        d,
    ):  # the peak falls smoothly through 10 as the delay passes 0.1
        peak = 10.0 * np.exp(-(d - 0.1) * 8)
        return peak >= 10.0, min(peak, 10.0)

    far, near = _bisect.crossing(
        respond, 1.0, 0.0, lambda n, y: 1e-4, 10.0, [(1.0, respond(1.0)[1])]
    )
    assert near <= 0.1 <= far and far - near <= 1e-4
