"""Tests for msoaxon.measure.passive_step (issue #68)."""

import numpy as np
import pytest

from msoaxon import measure, multi

N_GRID = 20001  # coarser than the default; plenty for these tolerances


def _step(model="multi", **kw):
    return measure.passive_step(model, n_grid=N_GRID, **kw)


def test_default_membrane_matches_the_rebalanced_values():
    # Before the fix: 61.4 / 68.7 MOhm and tau 0.0, because the baseline was
    # v0 and the membrane drifts by more than the step's response (microvolts).
    # pinned on the reference, as in test_analysis: native's steps move
    # rin_peak by ~1e-3
    reb = _step(mem=multi.membrane(-68.0, rebalance_rest=True), mode="exact")
    r = _step(mode="exact")
    assert 2.4 <= reb["rin_steady"] < 2.5 and 4.6 <= reb["rin_peak"] < 4.7
    # The default membrane rests ~0.13 mV below v0 where the conductances
    # differ a little, so it reads ~2% higher: 2.5 / 4.8 at one decimal.
    assert np.isclose(r["rin_steady"], reb["rin_steady"], rtol=0.03)
    assert np.isclose(r["rin_peak"], reb["rin_peak"], rtol=0.03)
    assert np.isclose(r["tau"], reb["tau"], atol=0.005)
    assert np.isclose(r["rin_steady"], 2.4826, rtol=1e-3)
    assert np.isclose(r["rin_peak"], 4.7619, rtol=1e-3)
    assert np.isclose(r["tau"], 0.1217, atol=1e-3)


def test_rebalanced_membrane_is_unchanged():
    # The dendrites table values (see test_analysis) on the coarse grid.
    r = _step(mem=multi.membrane(-68.0, rebalance_rest=True), mode="exact")
    assert np.isclose(r["rin_steady"], 2.4307, rtol=1e-4)
    assert np.isclose(r["rin_peak"], 4.6814, rtol=1e-4)
    assert np.isclose(r["tau"], 0.1197, atol=1e-4)


def test_rebalancing_at_the_settled_voltage_explains_the_default_offset():
    # The default membrane's ~2% excess is its operating point, not drift:
    # rebalancing at the voltage it settles to reproduces it.
    v = -68.1286
    r = _step(v0=v, mem=multi.membrane(v, rebalance_rest=True))
    assert np.isclose(r["rin_steady"], _step()["rin_steady"], rtol=2e-3)


@pytest.mark.parametrize(
    "kw",
    [
        {"model": "multi", "model_type": "passive"},
        {"model": "two", "model_type": "passive"},
        {"model": "two", "model_type": "passive", "v0": -60.0},
    ],
)
def test_positive_step_matches_negative_for_a_passive_cell(kw):
    neg = _step(I=-2.0, **kw)
    pos = _step(I=2.0, **kw)
    for k in ("rin_steady", "rin_peak"):
        assert neg[k] > 0
        assert np.isclose(pos[k], neg[k], rtol=1e-4)
    assert np.isclose(pos["tau"], neg["tau"], atol=1e-4)


def test_positive_step_matches_negative_with_active_channels():
    # Rectification makes the two signs differ slightly (~0.1%).
    neg = _step(I=-2.0)
    pos = _step(I=2.0)
    assert np.isclose(pos["rin_steady"], neg["rin_steady"], rtol=5e-3)
    assert np.isclose(pos["rin_peak"], neg["rin_peak"], rtol=5e-3)
    assert np.isclose(pos["tau"], neg["tau"], atol=0.005)


def test_passive_membrane_gives_a_positive_resistance():
    # Before the fix: -384.8 MOhm.
    r = _step(model_type="passive")
    assert np.isclose(r["rin_steady"], 9.9707, rtol=1e-3)
    assert np.isclose(r["tau"], 0.706, atol=0.005)


def test_two_compartment_positive_step_has_a_peak():
    # Before the fix: peak 0.09 MOhm and tau 0.0 for I = +2.
    r = _step("two", I=2.0)
    assert np.isclose(r["rin_steady"], 2.4256, rtol=1e-3)
    assert np.isclose(r["rin_peak"], 4.8758, rtol=1e-3)
    assert np.isclose(r["tau"], 0.1157, atol=1e-3)


def test_response_does_not_depend_on_the_initial_potential():
    # A passive cell started 8 mV off rest relaxes during the run; the
    # control subtraction removes that.
    a = _step("two", v0=-68.0, model_type="passive")
    b = _step("two", v0=-60.0, model_type="passive")
    assert np.isclose(a["rin_steady"], b["rin_steady"], rtol=1e-4)


def test_zero_current_is_rejected():
    with pytest.raises(ValueError, match="nonzero"):
        measure.passive_step("multi", I=0.0)
