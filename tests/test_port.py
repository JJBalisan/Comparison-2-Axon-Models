"""Checks that need no MATLAB: constants vs the stored .mat files, rest state, helpers."""

from pathlib import Path

import h5py
import numpy as np
import pytest

from msoaxon import constants as C
from msoaxon import mso_axon, two_cpt
from msoaxon.spiking import count_spikes, matlab_round, spiking

REPO = Path(__file__).resolve().parents[1]


def _mat(name):
    with h5py.File(REPO / name) as h:
        return {k: np.asarray(h[k]).ravel() for k in h if not k.startswith("#")}


def test_area_matches_constants_m_output():
    m = _mat("Area.mat")
    for key, ours in [("SA", C.SA), ("SAcm", C.SA_CM), ("areaRatio", C.AREA_RATIO),
                      ("L", C.L), ("Lcm", C.L_CM), ("Rcyl", C.R_CYL), ("XA", C.XA),
                      ("XAcm", C.XA_CM)]:
        np.testing.assert_allclose(ours, m[key], rtol=1e-12, err_msg=key)


def test_fractions_match_constants_m_output():
    m = _mat("Fractions.mat")
    np.testing.assert_allclose(C.KLT_FRAC, m["KLT_frac"], rtol=1e-12)
    np.testing.assert_allclose(C.H_FRAC, m["H_frac"], rtol=1e-12)
    np.testing.assert_allclose(C.NA_FRAC, m["Na_frac"], rtol=1e-12)


def test_coupling_matches_mat():
    m = _mat("Coupling.mat")
    np.testing.assert_array_equal(C.COUPLING1, m["coupling1"])
    np.testing.assert_array_equal(C.COUPLING2, m["coupling2"])


@pytest.mark.parametrize("model_type", ["passive", "active-KLT", "active-full"])
def test_two_cpt_stays_at_rest_without_input(model_type):
    _, x = two_cpt("step", 5, 10, 0.0, 3, model_type, 20, -68, 1)
    assert np.abs(x[:, :2] + 68).max() < 1e-5


def test_mso_axon_near_rest_without_input():
    # Unlike TwoCptODE.m, msoAxon.m does not subtract resting currents and starts
    # gating at rounded values (m0=.12, h0=.67, ...), so it sits ~0.9 mV off -68
    _, y = mso_axon("step", 5, 10, 0.0, 3, "active-full", 20, -68, 1)
    assert np.abs(y[:, :45] + 68).max() < 1.5


@pytest.mark.parametrize("model", [two_cpt, mso_axon])
def test_strong_epsg_pair_spikes(model):
    _, x = model("EPSGpair", 5, 5, 100, 3, "active-full", 20, -68, 1)
    axon = 1 if model is two_cpt else 2
    assert count_spikes(x[:, 0], x[:, axon], 10) == 1


def test_matlab_round_is_half_away_from_zero():
    assert matlab_round(2.5) == 3 and matlab_round(-2.5) == -3
    assert matlab_round(15, -1) == 20 and matlab_round(25, -1) == 30


def test_spiking_reset_carries_across_columns():
    # column 1 of a 'Two' input is soma - soma = 0, which re-arms before column 2
    x = np.array([[0, 0], [0, 50], [0, 50], [0, -1], [0, 50.0]])
    np.testing.assert_array_equal(spiking(x, 30, "Two")[:, 1], [0, 1, 0, 0, 1])
