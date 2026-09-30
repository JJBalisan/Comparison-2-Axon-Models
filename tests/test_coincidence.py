"""EPSG kinetics, two-compartment recalibration, and the coincidence-window analysis."""

import numpy as np
import pytest

from msoaxon import SynParams, two_cpt
from msoaxon import _solve
from msoaxon.synaptic import epsg_unitary
from msoaxon.coincidence import half_width, threshold
from msoaxon.measure import passive_step
from msoaxon.two import GOLDWYN_2019

GOLDWYN_KW = {k: v for k, v in GOLDWYN_2019.items() if k != "v0"}


@pytest.mark.parametrize("tau", [(0.1, 0.18), (0.1, 0.3), (0.05, 0.5)])
def test_epsg_waveforms_peak_at_one(tau):
    t = np.linspace(0, 3, 300001)
    assert np.isclose(epsg_unitary(t, tau).max(), 1.0, atol=1e-4)


def test_default_epsg_is_the_matlab_expression():
    t = np.linspace(-1, 3, 4001)
    expected = (1 / 0.21317) * (t > 0) * (np.exp(-t / 0.18) - np.exp(-t / 0.1))
    np.testing.assert_array_equal(epsg_unitary(t), expected)


def _passive_response(v0=-68.0, **kw):
    r = passive_step("two", v0, model_type="passive", n_grid=399001, **kw)
    return r["rin_steady"], r["tau"]


def test_goldwyn_calibration_hits_its_passive_targets():
    r_in, tau = _passive_response(**GOLDWYN_2019)
    assert np.isclose(r_in, 8.5, rtol=1e-3)
    assert np.isclose(tau, 0.34, atol=0.005)
    r_in, tau = _passive_response()  # TwoCpt.m's own values
    assert np.isclose(r_in, 10.0, rtol=1e-3)
    assert np.isclose(tau, 0.71, atol=0.01)


def test_prestimulus_cache_separates_calibrations():
    args = ("EPSGpair", 5, 5.2, 30, 3, "active-full", 20, -68.0, 1)
    _solve._QUIET_CACHE.clear()
    two_cpt(*args)  # fills the cache for the default calibration
    t1, x1 = two_cpt(*args, **GOLDWYN_KW)
    _solve._QUIET_CACHE.clear()
    t2, x2 = two_cpt(*args, **GOLDWYN_KW)
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(x1, x2)


def test_slower_epsg_lowers_the_pair_threshold():
    fast = threshold("two", 0.0, rel_tol=1e-3)
    slow = threshold("two", 0.0, epsg_tau=(0.1, 0.3), rel_tol=1e-3)
    assert 46.9 < fast <= 48.05  # inside the ported search's bracket
    assert slow < fast


def test_threshold_is_the_spike_boundary():
    th = threshold("two", 0.1, rel_tol=1e-4)
    syn = SynParams(t_end=20)
    for amp, expect in ((th * 1.001, True), (th * 0.999, False)):
        t, x = two_cpt("EPSGpair", 5, 5.1, amp, 3, "active-full", 20, -68, 1, syn)
        assert ((x[:, 1] - x[:, 0]) > 10).any() == expect


def test_half_width_on_a_known_curve():
    delays = np.linspace(0, 0.5, 51)
    thresholds = 50 * (1 + delays ** 2)  # crosses 1.03x at sqrt(0.03)
    assert np.isclose(half_width(delays, thresholds, 0.03), 2 * np.sqrt(0.03), atol=2e-3)
    assert np.isnan(half_width(delays, thresholds, 5.0))
