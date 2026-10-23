"""Membrane overrides and the somatic-spike measurement.

The measurement follows Scott et al 2005's protocol.
"""

import numpy as np
import pytest

import msoaxon
from msoaxon import _solve, multi, somatic

DTX = {"soma_klt_scale": 0, "ais_klt_scale": 0, "rebalance_rest": True}
ARGS = ("EPSGpair", 5, 5.2, 60, 3, "active-full", 12, -68.0, 1)


def test_default_membrane_is_msoaxon_m():
    t1, y1 = msoaxon.mso_axon(*ARGS)
    t2, y2 = msoaxon.mso_axon(*ARGS, mem=multi.membrane(-68.0))
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


@pytest.mark.parametrize(
    "knobs",
    [
        {"rebalance_rest": True},
        DTX,
        {"soma_klt_scale": 3, "soma_na_vhalf": -77, "rebalance_rest": True},
    ],
)
def test_rebalanced_rest_holds(knobs):
    _, y = msoaxon.mso_axon(
        "step",
        5,
        10,
        0.0,
        3,
        "active-full",
        30,
        -68.0,
        1,
        mem=multi.membrane(-68.0, **knobs),
    )
    assert np.abs(y[:, :45] + 68).max() < 1e-6


def test_unbalanced_klt_block_fires_on_its_own():
    # why rebalancing matters: without it the soma depolarises and fires with no
    # input
    mem = multi.membrane(-68.0, soma_klt_scale=0, ais_klt_scale=0)
    _, y = msoaxon.mso_axon(
        "step", 5, 10, 0.0, 3, "active-full", 20, -68.0, 1, mem=mem
    )
    assert (y[:, 2] > -20).any()


def test_prestimulus_cache_separates_membranes():
    _solve._QUIET_CACHE.clear()
    msoaxon.mso_axon(*ARGS)
    t1, y1 = msoaxon.mso_axon(*ARGS, mem=multi.membrane(-68.0, **DTX))
    _solve._QUIET_CACHE.clear()
    t2, y2 = msoaxon.mso_axon(*ARGS, mem=multi.membrane(-68.0, **DTX))
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


def test_somatic_spike_is_mature_and_klt_sets_it():
    rb = somatic.rheobase(
        "multi", rel_tol=1e-2
    )  # the margins below are far wider than 1%
    control = somatic.spike_amplitude("multi", 2 * rb)["amplitude"]
    assert (
        10 < control < 19
    )  # Scott et al 2005: 17 +/- 2 mV mature, 5-15 near threshold

    mem = multi.membrane(-68.0, **DTX)
    rb_dtx = somatic.rheobase("multi", mem=mem, rel_tol=1e-2)
    blocked = somatic.spike_amplitude("multi", 2 * rb_dtx, mem=mem)["amplitude"]
    assert blocked > 1.8 * control  # dendrotoxin: 15 -> 37 mV
    assert rb_dtx < rb / 5


# rheobases at rel_tol 1e-4 (CHANGELOG 1.2.0), so the multiples are exact
RHEOBASE = {"multi": 3784.1, "two": 3846.0}
# PYTHON.md "Somatic spike" table: mV at 1.5x, 2x and 3x rheobase
TABLE = {
    "multi": {1.5: 7.6, 2.0: 13.4, 3.0: 19.4},
    "two": {1.5: 6.2, 2.0: 15.6, 3.0: 17.6},
}


@pytest.mark.parametrize(
    "model,k",
    [(m, k) for m in TABLE for k in TABLE[m]],
)
def test_somatic_spike_table_points_are_distinct(model, k):
    r = somatic.spike_amplitude(model, k * RHEOBASE[model])
    assert r["distinct"]
    assert r["amplitude"] == pytest.approx(TABLE[model][k], abs=0.06)
    assert r["t_inflection"] is not None and r["v_inflection"] is not None


@pytest.mark.parametrize("model", ["multi", "two"])
def test_somatic_spike_at_10x_has_no_distinct_inflection(model):
    # the axon crosses the criterion before the search window opens, and the
    # old picker returned ~0 mV from the upstroke next to the peak
    r = somatic.spike_amplitude(model, 10 * RHEOBASE[model])
    assert not r["distinct"]
    assert r["amplitude"] is None
    assert r["t_inflection"] is None and r["v_inflection"] is None
    assert r["t_spike"] <= somatic.START + somatic.CHARGING
    assert r["peak_above_rest"] > 60  # it did spike
    assert len(r["trace"][0]) > 0


@pytest.mark.parametrize("model", ["multi", "two"])
def test_somatic_spike_near_threshold_is_a_smooth_hump(model):
    # 1.05x used to return 8.8-9.3 mV, bigger than 1.5x
    r = somatic.spike_amplitude(model, 1.05 * RHEOBASE[model])
    assert not r["distinct"]
    assert r["amplitude"] is None
    assert r["peak_above_rest"] > 10


def test_somatic_spike_without_a_step_spike_is_none():
    assert somatic.spike_amplitude("multi", 0.5 * RHEOBASE["multi"]) is None


@pytest.mark.parametrize("k", [1.5, 3.0])
def test_somatic_spike_with_klt_removed_stays_distinct(k):
    # the sharp takeoff keeps its inflection at every step size, unlike the
    # control membrane. The amplitude is only bounded: this membrane's response
    # is steep in the current, and CI platforms differ by ~1 mV at 1.5x
    mem = multi.membrane(-68.0, **DTX)
    r = somatic.spike_amplitude("multi", k * 199.524, mem=mem)
    assert r["distinct"]
    assert 20 < r["amplitude"] < 40  # dendrotoxin: 15 -> 37 mV
