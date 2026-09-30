"""Dendritic morphology (Lehnert et al 2014 Fig. 8 variant) and EPSGbilateral."""

import numpy as np
import pytest

from msoaxon import _solve, mso_axon
from msoaxon import constants as C
from msoaxon.coincidence import threshold
from msoaxon.measure import soma_on_grid
from msoaxon.multi import LUMPED, _tree_jac_sparsity, axial_current, membrane, with_dendrites

D = with_dendrites()
MID_L, MID_M = int(D.lateral[2]), int(D.medial[2])
ARGS = ("EPSGpair", 5, 5.2, 60, 3, "active-full", 12, -68.0, 1)


def test_default_morphology_is_msoaxon_m():
    t1, y1 = mso_axon(*ARGS)
    t2, y2 = mso_axon(*ARGS, morph=LUMPED)
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)
    # the general tree builder reproduces the chain's sparsity pattern
    assert (_tree_jac_sparsity(45, LUMPED.parent, LUMPED.child) != LUMPED.jac).nnz == 0


def test_dendritic_geometry_matches_lehnert():
    assert D.n == 55 and list(D.lateral) == [46, 47, 48, 49, 50]
    assert np.isclose(D.sa[0], 2467, atol=0.5)  # "reduced to 2467 um^2"
    assert np.isclose(D.sa[0] + D.sa[45:].sum(), C.SA[0])  # total stays 8750 um^2
    assert np.isclose(D.g_na[0] * D.sa[0], C.G_NA[0] * C.SA[0])  # total somatic Na kept
    assert not D.g_na[45:].any() and not D.g_kht[45:].any()  # no Na / KHT in dendrites
    klt = D.g_klt[45:50]
    assert np.all(np.diff(klt) < 0)
    assert np.allclose(klt[1:] / klt[:-1], np.exp(-40 / 74))  # length constant 74 um


def test_axial_current_conserves_charge():
    V = np.random.default_rng(0).normal(-60, 10, D.n)
    assert abs(np.sum(axial_current(V, D) * D.sa)) < 1e-9 * np.sum(np.abs(axial_current(V, D) * D.sa))


def test_dendritic_rest_holds_when_rebalanced():
    mem = membrane(-68.0, rebalance_rest=True, morph=D)
    _, y = mso_axon("step", 5, 10, 0.0, 3, "active-full", 20, -68.0, 1, morph=D, mem=mem)
    assert y.shape[1] == 7 * 55
    assert np.abs(y[:, :55] + 68).max() < 1e-6


def test_dendritic_epsp_is_attenuated_and_slowed():
    mem = membrane(-68.0, rebalance_rest=True, morph=D)
    peaks, widths = [], []
    for site in (1, MID_L):
        t, y = mso_axon("EPSG", 5, 10, 26.7, 3, "active-full", 10, -68.0, site, morph=D, mem=mem,
                        max_step=0.01)
        g, v = soma_on_grid(t, y, -68.0, 5, 10)
        peaks.append(v.max())
        widths.append(np.ptp(g[v >= v.max() / 2]))
    assert peaks[1] < peaks[0] and widths[1] > widths[0]


def test_bilateral_inputs_beat_unilateral():
    # Scott et al 2010: EPSGs on opposite dendrites reach threshold more easily
    # (58.5 vs 82.2, so a 1% tolerance is plenty)
    bil = threshold("multi", 0.0, rel_tol=1e-2, model_kw=dict(morph=D),
                    stim="EPSGbilateral", input_node=MID_L, input_node2=MID_M)
    uni = threshold("multi", 0.0, rel_tol=1e-2, model_kw=dict(morph=D),
                    stim="EPSGbilateral", input_node=MID_L, input_node2=MID_L)
    assert bil < uni


def test_prestimulus_cache_separates_morphologies():
    _solve._QUIET_CACHE.clear()
    mso_axon(*ARGS)
    t1, y1 = mso_axon(*ARGS, morph=D)
    _solve._QUIET_CACHE.clear()
    t2, y2 = mso_axon(*ARGS, morph=D)
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


def test_argument_checks():
    with pytest.raises(ValueError, match="input_node2"):
        mso_axon("EPSGbilateral", 5, 5, 60, 3, "active-full", 10, -68.0, MID_L, morph=D)
    with pytest.raises(ValueError, match="different morphology"):
        mso_axon(*ARGS, morph=D, mem=membrane(-68.0))
    with pytest.raises(ValueError, match="input_node"):
        mso_axon("EPSG", 5, 10, 26.7, 3, "active-full", 10, -68.0, MID_L)  # 48 needs dendrites
    with pytest.raises(ValueError, match="node must be 1..45"):  # spikes are read on the axon
        mso_axon("EPSG", 5, 10, 26.7, MID_L, "active-full", 10, -68.0, 1, morph=D)
    with pytest.raises(ValueError, match="stimType"):
        from msoaxon import two_cpt
        two_cpt("EPSGbilateral", 5, 5, 60, 3, "active-full", 10, -68.0, 1)


def test_total_klt_and_h_conserved_by_default():
    region = np.r_[0, np.arange(45, 55)]
    for g, lumped in ((D.g_klt, C.G_KLT[0]), (D.g_h, C.G_H[0])):
        assert np.isclose(np.sum(g[region] * D.sa[region]), lumped * C.SA[0])
    alt = with_dendrites(conserve_totals=False)  # gradient from the soma's Table 2 density
    assert alt.g_klt[0] == C.G_KLT[0]
    assert np.sum(alt.g_klt[region] * alt.sa[region]) < 0.6 * C.G_KLT[0] * C.SA[0]
    assert alt.key != D.key


def test_dendrite_ra_changes_only_dendritic_edges():
    hi = with_dendrites(dendrite_ra=200)
    assert hi.key != D.key
    into = D.child >= 45
    assert into.sum() == 10
    np.testing.assert_array_equal(hi.g_ax[~into], D.g_ax[~into])  # soma and axon untouched
    distal = into & (D.parent >= 45)  # dendrite-to-dendrite edges: both halves at 200
    np.testing.assert_allclose(hi.g_ax[distal], D.g_ax[distal] / 2, rtol=1e-12)
    root = into & (D.parent < 45)  # soma-to-dendrite: only the dendritic half doubles
    assert np.all((hi.g_ax[root] < D.g_ax[root]) & (hi.g_ax[root] > D.g_ax[root] / 2))


def test_prestimulus_cache_separates_membrane_v0():
    # the leak reversals come from membrane()'s v0, so it must be part of the cache key
    assert membrane(-60.0).key != membrane(-68.0).key
    args = ("EPSG", 5, 10, 26.7, 3, "active-full", 8, -68.0, 1)
    _solve._QUIET_CACHE.clear()
    mso_axon(*args, mem=membrane(-60.0))
    t1, y1 = mso_axon(*args, mem=membrane(-68.0))
    _solve._QUIET_CACHE.clear()
    t2, y2 = mso_axon(*args, mem=membrane(-68.0))
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


def test_shared_arrays_are_read_only_and_dendrite_scale_needs_dendrites():
    with pytest.raises(ValueError, match="read-only"):
        C.G_KLT[0] = 0.0
    with pytest.raises(ValueError, match="read-only"):
        LUMPED.g_ax[0] = 0.0
    with pytest.raises(ValueError, match="dendrite_klt_scale"):
        membrane(-68.0, dendrite_klt_scale=0.5)


def test_morph_defaults_to_mem_morph():
    mem = membrane(-68.0, morph=D)
    t1, y1 = mso_axon(*ARGS, mem=mem)
    t2, y2 = mso_axon(*ARGS, mem=mem, morph=D)
    np.testing.assert_array_equal(y1, y2)
    assert y1.shape[1] == 7 * 55
