"""Resonance protocol (#83) and the KLT-to-leak substitution of membrane()."""

import numpy as np
import pytest

from msoaxon import measure, multi, resonance as R

REBALANCED = {"rebalance_rest": True}
# a coarse grid: the two-compartment model peaks at ~350 Hz
COARSE = (20.0, 200.0, 350.0, 500.0, 700.0)


# --- the measurement ---------------------------------------------------------


def test_fundamental_recovers_amplitude_and_phase_through_a_drift():
    phase = np.linspace(0, 2 * np.pi * 10, 640, endpoint=False)
    x = np.linspace(-1, 1, len(phase))
    y = 3 * np.sin(phase + 0.4) + 1.5 + 0.7 * x - 0.4 * x**2
    c = R._fundamental(phase, y)
    assert abs(c) == pytest.approx(3.0, rel=1e-9)
    assert np.angle(c) == pytest.approx(0.4, abs=1e-9)


def test_half_wave_sine_has_a_fundamental_of_half_its_peak():
    phase = np.linspace(0, 2 * np.pi * 10, 640, endpoint=False)
    c = R._fundamental(phase, 20 * np.maximum(np.sin(phase), 0))
    assert c == pytest.approx(10.0, abs=0.02)  # the detrend leaks a little


def test_vertex_of_a_parabola():
    f = np.linspace(100, 200, 11)
    y = -((f - 137.0) ** 2) + 5
    assert R._vertex(f, y, int(np.argmax(y))) == pytest.approx(137.0)
    assert R._vertex(f, -y, int(np.argmin(-y))) == pytest.approx(137.0)
    # five points at the end see a vertex they do not contain: the grid point
    assert R._vertex(f, y, 0) == f[0]


def test_passive_impedance_is_low_pass_and_matches_the_step_resistance():
    kw = {"model_type": "passive"}
    r = R.subthreshold("two", COARSE, kw, workers=1, n_fine=3)
    assert not r["resonant"]
    assert np.isnan(r["f_res"])
    assert r["q"] == pytest.approx(1.0)
    assert np.all(np.diff(r["z"]) < 0)
    rin = measure.passive_step("two", model_type="passive")["rin_steady"]
    assert r["z_low"] == pytest.approx(rin, rel=0.03)


def test_impedance_is_linear_in_the_amplitude():
    z20 = R.impedance("two", 350.0, amp=20.0)
    z40 = R.impedance("two", 350.0, amp=40.0)
    assert abs(z40) == pytest.approx(abs(z20), rel=0.02)


def test_impedance_checks_its_arguments():
    with pytest.raises(ValueError, match="f"):
        R.impedance("two", 0.0)
    with pytest.raises(ValueError, match="cycles"):
        R.impedance("two", 100.0, cycles=0)


def test_two_compartment_model_resonates_in_the_subthreshold_range():
    r = R.subthreshold("two", COARSE, workers=1, n_fine=3)
    assert r["resonant"]
    assert 250 < r["f_res"] < 450
    assert r["q"] > 1.5
    assert r["z_peak"] == pytest.approx(max(r["z"]), rel=0.02)
    assert r["f_low"] == COARSE[0]


def test_multi_compartment_model_resonates_too():
    f = np.array([20.0, 350.0, 800.0])
    z = [abs(R.impedance("multi", x)) for x in f]
    assert z[1] > 1.5 * z[0]
    assert z[1] > 1.5 * z[2]


# --- thresholds --------------------------------------------------------------


def test_epsg_train_has_one_event_per_period():
    ev = R.epsg_train(250.0, dur=100.0)
    assert len(ev) == 25
    assert [e.t for e in ev[:3]] == pytest.approx(
        [R.START, R.START + 4.0, R.START + 8.0]
    )
    assert len(R.epsg_train(10.0, dur=20.0)) == 1  # at least one
    assert {e.node for e in R.epsg_train(100.0, node=2)} == {2}


def test_sine_threshold_is_a_minimum_in_the_middle_of_the_range():
    kw = {"rel_tol": 1e-2, "dur": 40.0}
    th = [R.sine_threshold("two", f, **kw) for f in (100.0, 350.0, 800.0)]
    assert np.all(np.isfinite(th))
    assert th[1] < th[0]
    assert th[1] < th[2]


def test_epsg_train_threshold_dips_where_the_period_matches_the_ring():
    kw = {"rel_tol": 1e-2}
    one = R.epsg_train_threshold("two", 400.0, dur=1.0, **kw)  # one event
    ring = R.epsg_train_threshold("two", 400.0, dur=20.0, **kw)
    off = R.epsg_train_threshold("two", 300.0, dur=20.0, **kw)
    assert ring < 0.95 * one
    assert off == pytest.approx(one, rel=0.01)  # no summation at 300 Hz


def test_threshold_functions_check_their_arguments():
    with pytest.raises(ValueError, match="dur"):
        R.sine_threshold("two", 300.0, dur=0.0)
    with pytest.raises(ValueError, match="rate"):
        R.epsg_train_threshold("two", -1.0)


def test_scan_reports_an_extremum_at_the_edge():
    r = R.spike_resonance(
        "two",
        (400.0, 600.0, 800.0),  # beyond the minimum, which is ~350
        workers=1,
        dur=40.0,
        rel_tol=2e-2,
    )
    assert r["at_edge"]
    assert r["f_best"] == 400.0
    assert r["freqs"] == [400.0, 600.0, 800.0]  # no second pass
    assert not r["resonant"]
    assert r["dip"] == 1.0  # the lowest frequency is the minimum


def test_spike_and_epsg_resonance_find_an_interior_minimum():
    kw = {"workers": 1, "n_fine": 3, "dur": 20.0, "rel_tol": 1e-2}
    sine = R.spike_resonance("two", (200.0, 500.0, 800.0), **kw)
    epsg = R.epsg_resonance("two", (300.0, 400.0, 600.0), **kw)
    for r in (sine, epsg):
        assert not r["at_edge"]
        assert r["resonant"]
        assert r["dip"] > R.MIN_RATIO
        assert r["threshold"] == pytest.approx(min(r["thresholds"]))
        assert r["freqs"] == sorted(r["freqs"])
    assert 350 < sine["f_best"] < 700
    assert 350 < epsg["f_best"] < 450


# --- membrane(klt_to_leak=...) ----------------------------------------------


def test_default_membrane_is_untouched_by_the_option():
    m = multi.membrane(-68.0)
    assert m.morph is multi.LUMPED
    assert m.key == multi.membrane(-68.0, klt_to_leak=None).key
    assert m.g_klt.sum() > 0


@pytest.mark.parametrize("kind", ["chord", "slope"])
def test_klt_to_leak_moves_the_conductance_into_the_leak(kind):
    m = multi.membrane(-68.0, klt_to_leak=kind, **REBALANCED)
    ref = multi.membrane(-68.0, **REBALANCED)
    assert not m.g_klt.any()
    # the leak grew exactly where KLT was, and nowhere else
    grew = m.morph.g_lk > multi.LUMPED.g_lk
    np.testing.assert_array_equal(grew, ref.g_klt > 0)
    np.testing.assert_array_equal(m.morph.g_lk[~grew], multi.LUMPED.g_lk[~grew])
    assert m.morph_key != multi.LUMPED.key
    assert m.morph.key != ref.morph.key
    assert np.isfinite(m.vlk).all()
    # everything else in the morphology is shared, not copied
    assert m.morph.g_klt is multi.LUMPED.g_klt
    assert m.morph.jac is multi.LUMPED.jac


def test_chord_leak_is_the_klt_conductance_at_rest():
    m = multi.membrane(-68.0, klt_to_leak="chord")
    v = np.full(multi.N, -68.0)
    g_rest = multi.LUMPED.g_klt * multi.C.winf(v) ** 4 * multi.C.zinf(v)
    np.testing.assert_allclose(
        m.morph.g_lk - multi.LUMPED.g_lk, g_rest, rtol=1e-12
    )


@pytest.mark.parametrize("kind", ["chord", "slope"])
def test_klt_to_leak_keeps_the_rest_potential(kind):
    m = multi.membrane(-68.0, klt_to_leak=kind, **REBALANCED)
    _, y = multi.mso_axon(
        "step", 5, 10, 0.0, 3, "active-full", 60, -68.0, 1, mem=m
    )
    assert np.abs(y[:, 0] + 68.0).max() < 1e-4


def test_the_tangent_leak_keeps_the_input_resistance_the_frozen_one_does_not():
    rin = {
        kind: measure.passive_step(
            "multi", mem=multi.membrane(-68.0, klt_to_leak=kind, **REBALANCED)
        )
        for kind in (None, "chord", "slope")
    }
    assert rin["chord"]["rin_steady"] > 3 * rin[None]["rin_steady"]
    assert rin["slope"]["rin_steady"] == pytest.approx(
        rin[None]["rin_steady"], rel=0.15
    )
    assert (
        rin["slope"]["rin_peak"] < 1.1 * rin["slope"]["rin_steady"]
    )  # no hump
    assert rin[None]["rin_peak"] > 1.5 * rin[None]["rin_steady"]


def test_klt_scales_carry_into_the_substitution():
    base = multi.membrane(-68.0, klt_to_leak="chord", **REBALANCED)
    soma = multi.membrane(
        -68.0, soma_klt_scale=2.0, klt_to_leak="chord", **REBALANCED
    )
    assert soma.morph.key != base.morph.key
    extra = base.morph.g_lk[0] - multi.LUMPED.g_lk[0]
    assert soma.morph.g_lk[0] - multi.LUMPED.g_lk[0] == pytest.approx(2 * extra)
    assert soma.morph.g_lk[1] == base.morph.g_lk[1]


def test_klt_to_leak_with_dendrites():
    d = multi.with_dendrites()
    m = multi.membrane(-68.0, morph=d, klt_to_leak="slope", **REBALANCED)
    assert m.morph.n == d.n
    assert not m.g_klt.any()
    assert (m.morph.g_lk[multi.N :] > d.g_lk[multi.N :]).all()


def test_klt_to_leak_rejects_unknown_kinds():
    with pytest.raises(ValueError, match="klt_to_leak"):
        multi.membrane(-68.0, klt_to_leak="frozen")
    with pytest.raises(ValueError, match="klt_to_leak"):
        multi.membrane(-68.0, klt_to_leak=True)


def test_mem_with_the_leak_replaces_the_morphology_of_the_run():
    m = multi.membrane(-68.0, klt_to_leak="chord", **REBALANCED)
    args = ("step", 5, 10, 0.0, 3, "active-full", 20, -68.0, 1)
    multi.mso_axon(*args, mem=m)  # morph comes from mem
    with pytest.raises(ValueError, match="different morphology"):
        multi.mso_axon(*args, mem=m, morph=multi.LUMPED)


def test_the_leak_substitution_removes_the_subthreshold_resonance():
    mem = multi.membrane(-68.0, klt_to_leak="slope", **REBALANCED)
    r = R.subthreshold("multi", COARSE, {"mem": mem}, workers=1, n_fine=3)
    assert not r["resonant"]  # flat to 1%, then falling: Q 1.001, no peak
    assert r["q"] < R.MIN_RATIO
    assert r["z"][-1] < 0.8 * r["z"][0]


def test_a_flat_response_is_not_called_resonant():
    flat = {"at_edge": False, "value": 1.0, "values": [1.02, 1.0, 1.02]}
    assert not R._dip({**flat, "freqs": [1, 2, 3], "f_best": 2.0})["resonant"]
