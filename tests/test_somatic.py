"""Membrane overrides and the somatic-spike measurement (Scott et al 2005 protocol)."""

import numpy as np
import pytest

import msoaxon
from msoaxon import _solve, multi, somatic

DTX = dict(soma_klt_scale=0, ais_klt_scale=0, rebalance_rest=True)
ARGS = ("EPSGpair", 5, 5.2, 60, 3, "active-full", 12, -68.0, 1)


def test_default_membrane_is_msoaxon_m():
    t1, y1 = msoaxon.mso_axon(*ARGS)
    t2, y2 = msoaxon.mso_axon(*ARGS, mem=multi.membrane(-68.0))
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


@pytest.mark.parametrize(
    "knobs",
    [
        dict(rebalance_rest=True),
        DTX,
        dict(soma_klt_scale=3, soma_na_vhalf=-77, rebalance_rest=True),
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
    # why rebalancing matters: without it the soma depolarises and fires with no input
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
