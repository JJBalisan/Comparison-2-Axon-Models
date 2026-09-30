"""The analysis layer: parallel helpers, threshold sweeps, coincidence trials,
somatic measurements and passive measurements."""

import numpy as np
import pytest

from msoaxon import _solve
from msoaxon._parallel import map_tasks
from msoaxon.coincidence import threshold_curve


def test_map_tasks_handles_no_tasks():
    assert map_tasks(abs, []) == []
    assert threshold_curve("two", []).shape == (0,)


# --- threshold sweeps (threshold.py) --------------------------------------------------

def test_sweep_settings_follow_binarysearch_m():
    from msoaxon.threshold import START, sweep_setting, sweep_x

    assert sweep_setting("step", 3)[0] == 15.0 and sweep_setting("EPSG", 3)[0] == 10.0
    assert sweep_setting("ramp", 4)[0] == 5.4 and sweep_x("ramp", 4) == 5.4
    stop, syn, _ = sweep_setting("sine", 2)
    assert syn.f == 200 and stop == START + 1000 / 400  # the positive half-cycle
    assert sweep_x("sine", 2) == 200
    assert sweep_setting("EPSGpair", 3, 0.1)[0] == START + 0.2 and sweep_x("EPSGpair", 3, 0.1) == 0.2
    stop, syn, I = sweep_setting("Synaptic", 1)
    assert (stop, syn.random_in, I) == (15.0, 13986, 0.0)
    with pytest.raises(ValueError, match="not supported"):
        sweep_setting("SynapticPair", 1)
    with pytest.raises(ValueError, match="no x-axis"):
        sweep_x("step", 1)


def test_sine_sweep_point_regression():
    # in-process (workers=1) so the halving search itself runs under the test;
    # values from the current code, which reproduces the MATLAB-faithful baseline
    from msoaxon.threshold import binary_search

    assert binary_search("sine", 1, 3, 30, 15000, workers=1, rounded=False) == (
        [11054.0771484375], [11541.1376953125])


# --- coincidence trials (coincidence.py) ----------------------------------------------

def test_noise_free_trials_agree_with_the_threshold_curve():
    from msoaxon.coincidence import probability_trials

    delays = [0.0, 0.1, 0.6]
    th = threshold_curve("two", delays, workers=1)
    amp = 1.05 * th[0]
    assert np.all(np.abs(th / amp - 1) > 1e-3)  # no threshold within rel_tol of amp
    prob = probability_trials("two", delays, amp, n_trials=2, amp_cv=0.0, jitter=0.0, workers=1)
    np.testing.assert_array_equal(prob, (th <= amp).astype(float))
    assert prob[0] == 1.0 and prob[-1] == 0.0


def test_noisy_trials_are_seeded_probabilities():
    from msoaxon.coincidence import probability_trials

    kw = dict(n_trials=6, amp_cv=0.05, jitter=0.02, workers=1)
    p1 = probability_trials("two", [0.0, 0.5], 52.0, seed=3, **kw)
    np.testing.assert_array_equal(p1, probability_trials("two", [0.0, 0.5], 52.0, seed=3, **kw))
    assert np.all((p1 >= 0) & (p1 <= 1)) and p1[0] >= p1[1]


# --- model entry points and solver ----------------------------------------------------

def test_spikes_is_an_early_stop():
    from msoaxon._dispatch import run_model, spikes

    for I, expected in ((150.0, True), (20.0, False)):
        assert spikes("two", "EPSGpair", 5, 5, I, 3, 20, -68.0, 1, factor=10) is expected
        t, _ = run_model("two", "EPSGpair", 5, 5, I, 3, 20, -68.0, 1, stop_on_spike=10)
        assert bool(t[-1] < 20) is expected


def test_bilateral_at_one_site_is_an_epsg_pair():
    # the same two EPSGs, summed per input instead of per waveform: equal up to rounding
    from msoaxon import mso_axon

    t1, y1 = mso_axon("EPSGpair", 5, 5.2, 30, 3, "active-full", 8, -68.0, 1)
    t2, y2 = mso_axon("EPSGbilateral", 5, 5.2, 30, 3, "active-full", 8, -68.0, 1, input_node2=1)
    g = np.linspace(0, 8, 801)
    np.testing.assert_allclose(np.interp(g, t2, y2[:, 0]), np.interp(g, t1, y1[:, 0]), atol=1e-6)


def test_solver_failure_is_reported():
    with pytest.raises(RuntimeError, match="integration failed"):
        _solve.integrate(lambda t, y: y ** 2, [1.0], 2.0, [], 1e-6, 1e-6, 0.1)  # blows up at t = 1


def test_quiet_cache_is_bounded():
    _solve._QUIET_CACHE.clear()
    rhs = lambda t, y: -y  # noqa: E731
    n = _solve._QUIET_CACHE_SIZE + 5
    for i in range(n):
        _solve.integrate(rhs, [1.0], 2.0, [1.0], 1e-6, 1e-6, 0.5, quiet=(("test", i), rhs))
    assert len(_solve._QUIET_CACHE) == _solve._QUIET_CACHE_SIZE
    keys = [k[0] for k in _solve._QUIET_CACHE]
    assert ("test", 0) not in keys and ("test", n - 1) in keys  # oldest out first
    _solve._QUIET_CACHE.clear()


def test_remaining_argument_errors():
    from msoaxon import mso_axon
    from msoaxon.multi import with_dendrites
    from msoaxon.two import applied_current

    # TwoCpt.m's capitalised 'Active-sodium' exists only in the 2-compartment model
    with pytest.raises(ValueError, match="unknown model type"):
        mso_axon("step", 5, 10, 100.0, 3, "Active-sodium", 10, -68.0, 1)
    with pytest.raises(ValueError, match="larger than the lumped"):
        with_dendrites(length=2000.0)
    with pytest.raises(ValueError, match="unknown stimType"):
        applied_current(0.0, -68.0, "bogus", None)


# --- somatic and passive measurements -------------------------------------------------

def test_subthreshold_step_has_no_spike_amplitude():
    from msoaxon.somatic import fires, spike_amplitude

    assert not fires("two", 100.0)
    assert spike_amplitude("two", 100.0) is None


def test_spike_amplitude_trace_is_the_measured_run():
    from msoaxon.somatic import spike_amplitude

    r = spike_amplitude("two", 6000.0)
    t, x = r["trace"]
    assert np.isclose(t[-1], r["t_spike"] + 2.0)
    assert np.isclose(np.interp(r["t_inflection"], t, x[:, 0]), r["v_inflection"], atol=0.05)


def test_passive_step_matches_the_dendrites_table():
    # PYTHON.md / figures/dendrites: lumped model, rest-rebalanced, 2.4 / 4.7 MOhm
    from msoaxon.measure import passive_step
    from msoaxon.multi import membrane

    r = passive_step("multi", mem=membrane(-68.0, rebalance_rest=True))
    assert np.isclose(r["rin_steady"], 2.4307, rtol=1e-4)
    assert np.isclose(r["rin_peak"], 4.6815, rtol=1e-4)
    assert np.isclose(r["tau"], 0.1191, atol=1e-4)


def test_count_rising_edges():
    from msoaxon.spiking import count_rising_edges

    assert count_rising_edges([0, 1, 1, 0, 1, 0, 0, 1]) == 3
    assert count_rising_edges([1, 1, 0]) == 0  # a marker at the first sample isn't an edge
