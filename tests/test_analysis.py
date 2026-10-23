"""The analysis layer: parallel helpers, threshold sweeps, coincidence trials,
somatic measurements and passive measurements."""

import warnings
from concurrent import futures

import numpy as np
import pytest

import msoaxon
from msoaxon import (
    _common,
    _dispatch,
    _parallel,
    _solve,
    coincidence,
    measure,
    multi,
    somatic,
    spiking,
    synaptic,
    threshold,
    two,
)


def test_map_tasks_handles_no_tasks():
    assert _parallel.map_tasks(abs, []) == []
    assert coincidence.threshold_curve("two", []).shape == (0,)


# --- synaptic conductance memo (synaptic.py) ----------------------------------


def test_synaptic_memo_returns_equal_read_only_arrays():
    syn = msoaxon.SynParams(t_end=12.5, random_in=4242)
    t, g = synaptic.synaptic(syn)
    t2, g2 = synaptic.synaptic(syn)
    # the cached arrays are the same objects, equal to a fresh computation,
    # and a caller cannot change them
    assert t2 is t and g2 is g
    t_new, g_new = synaptic._table.__wrapped__(12.5, syn.freq, syn.gE, 4242)
    assert np.array_equal(t, t_new) and np.array_equal(g, g_new)
    for arr in (t, g):
        assert not arr.flags.writeable
        with pytest.raises(ValueError, match="read-only"):
            arr[0] = 1.0
    assert synaptic._table.cache_info().maxsize is not None  # bounded


def test_synaptic_memo_keys_on_the_fields_that_fix_the_waveform():
    base = msoaxon.SynParams(t_end=10.0, random_in=11)
    _, g = synaptic.synaptic(base)
    # a different seed, frequency, gE or length is a different table
    for other in (
        msoaxon.SynParams(t_end=10.0, random_in=12),
        msoaxon.SynParams(t_end=10.0, random_in=11, freq=300.0),
        msoaxon.SynParams(t_end=10.0, random_in=11, gE=2e7),
    ):
        _, g_other = synaptic.synaptic(other)
        assert g_other is not g and not np.array_equal(g_other, g)
    t_long, g_long = synaptic.synaptic(
        msoaxon.SynParams(t_end=12.0, random_in=11)
    )
    assert len(g_long) > len(g) and len(t_long) == len(g_long)
    # the fields the waveform ignores share the entry
    same = msoaxon.SynParams(
        t_end=10.0, random_in=11, VsynE=-10.0, diff=2.0, f=100.0
    )
    assert synaptic.synaptic(same)[1] is g


def test_stimulus_shares_the_cached_conductance():
    syn = msoaxon.SynParams(t_end=8.0, random_in=77)
    s1 = _common.stimulus("Synaptic", 5.0, 5.0, 1.0, 8.0, syn)
    s2 = _common.stimulus("SynapticPair", 5.0, 5.0, 1.0, 8.0, syn)
    assert np.array_equal(s1.g_syn, s2.g_syn) and not s1.g_syn.flags.writeable


# --- threshold sweeps (threshold.py) ------------------------------------------


def test_sweep_settings_follow_binarysearch_m():

    assert (
        threshold.sweep_setting("step", 3)[0] == 15.0
        and threshold.sweep_setting("EPSG", 3)[0] == 10.0
    )
    assert (
        threshold.sweep_setting("ramp", 4)[0] == 5.4
        and threshold.sweep_x("ramp", 4) == 5.4
    )
    stop, syn, _ = threshold.sweep_setting("sine", 2)
    assert (
        syn.f == 200 and stop == threshold.START + 1000 / 400
    )  # the positive half-cycle
    assert threshold.sweep_x("sine", 2) == 200
    assert (
        threshold.sweep_setting("EPSGpair", 3, 0.1)[0] == threshold.START + 0.2
        and threshold.sweep_x("EPSGpair", 3, 0.1) == 0.2
    )
    for stim in ("Synaptic", "SynapticPair"):
        with pytest.raises(ValueError, match="not supported"):
            threshold.sweep_setting(stim, 1)
        with pytest.raises(ValueError, match="not supported"):
            threshold.binary_search(stim, 1, 3, 10.0, 15000.0)
    with pytest.raises(ValueError, match="no x-axis"):
        threshold.sweep_x("step", 1)


def test_sine_sweep_point_regression():
    # in-process (workers=1) so the halving search itself runs under the test.
    # Its result is the last value tested; a flipped next-to-last decision
    # moves that by 15000 / 2**13 ~ 1.83 (1.7e-4 relative), which cross-
    # platform solver drift could do, so rel 2e-4 allows that one step and
    # nothing more

    multi, two = threshold.binary_search(
        "sine", 1, 3, 30, 15000, workers=1, rounded=False
    )
    assert multi == pytest.approx([11054.0771484375], rel=2e-4)
    assert two == pytest.approx([11541.1376953125], rel=2e-4)


# --- coincidence trials (coincidence.py) --------------------------------------


def test_noise_free_trials_agree_with_the_threshold_curve():

    delays = [0.0, 0.1, 0.6]
    th = coincidence.threshold_curve("two", delays, workers=1)
    amp = 1.05 * th[0]
    assert np.all(
        np.abs(th / amp - 1) > 1e-3
    )  # no threshold within rel_tol of amp
    prob = coincidence.probability_trials(
        "two", delays, amp, n_trials=2, amp_cv=0.0, jitter=0.0, workers=1
    )
    np.testing.assert_array_equal(prob, (th <= amp).astype(float))
    assert prob[0] == 1.0 and prob[-1] == 0.0


def test_noisy_trials_are_seeded_probabilities():

    kw = {"n_trials": 6, "amp_cv": 0.05, "jitter": 0.02, "workers": 1}
    p1 = coincidence.probability_trials("two", [0.0, 0.5], 52.0, seed=3, **kw)
    np.testing.assert_array_equal(
        p1,
        coincidence.probability_trials("two", [0.0, 0.5], 52.0, seed=3, **kw),
    )
    assert np.all((p1 >= 0) & (p1 <= 1)) and p1[0] >= p1[1]


# --- model entry points and solver --------------------------------------------


def test_spikes_is_an_early_stop():

    for I, expected in ((150.0, True), (20.0, False)):
        assert (
            _dispatch.spikes(
                "two", "EPSGpair", 5, 5, I, 3, 20, -68.0, 1, factor=10
            )
            is expected
        )
        t, _ = _dispatch.run_model(
            "two", "EPSGpair", 5, 5, I, 3, 20, -68.0, 1, stop_on_spike=10
        )
        assert bool(t[-1] < 20) is expected


def test_bilateral_at_one_site_is_an_epsg_pair():
    # the same two EPSGs, summed per input instead of per waveform: equal up to
    # rounding

    t1, y1 = msoaxon.mso_axon(
        "EPSGpair", 5, 5.2, 30, 3, "active-full", 8, -68.0, 1
    )
    t2, y2 = msoaxon.mso_axon(
        "EPSGbilateral",
        5,
        5.2,
        30,
        3,
        "active-full",
        8,
        -68.0,
        1,
        input_node2=1,
    )
    g = np.linspace(0, 8, 801)
    np.testing.assert_allclose(
        np.interp(g, t2, y2[:, 0]), np.interp(g, t1, y1[:, 0]), atol=1e-6
    )


def test_solver_failure_is_reported():
    with pytest.raises(RuntimeError, match="integration failed"):
        _solve.integrate(
            lambda t, y: y**2, [1.0], 2.0, [], 1e-6, 1e-6, 0.1
        )  # blows up at t = 1


def test_quiet_cache_is_bounded(fresh_quiet_cache):
    rhs = lambda t, y: -y  # noqa: E731
    n = _solve._QUIET_CACHE_SIZE + 5
    for i in range(n):
        _solve.integrate(
            rhs, [1.0], 2.0, [1.0], 1e-6, 1e-6, 0.5, quiet=(("test", i), rhs)
        )
    assert len(_solve._QUIET_CACHE) == _solve._QUIET_CACHE_SIZE
    keys = [k[0] for k in _solve._QUIET_CACHE]
    assert ("test", 0) not in keys and (
        "test",
        n - 1,
    ) in keys  # oldest out first


def test_remaining_argument_errors():

    # TwoCpt.m's capitalised 'Active-sodium' exists only in the 2-compartment
    # model
    with pytest.raises(ValueError, match="unknown model type"):
        msoaxon.mso_axon("step", 5, 10, 100.0, 3, "Active-sodium", 10, -68.0, 1)
    with pytest.raises(ValueError, match="larger than the lumped"):
        multi.with_dendrites(length=2000.0)
    with pytest.raises(ValueError, match="unknown stimType"):
        two.applied_current(0.0, -68.0, "bogus", None)


# --- somatic and passive measurements -----------------------------------------


def test_subthreshold_step_has_no_spike_amplitude():

    assert not somatic.fires("two", 100.0)
    assert somatic.spike_amplitude("two", 100.0) is None


def test_spike_amplitude_trace_is_the_measured_run():

    r = somatic.spike_amplitude("two", 6000.0)
    t, x = r["trace"]
    assert np.isclose(t[-1], r["t_spike"] + 2.0)
    assert np.isclose(
        np.interp(r["t_inflection"], t, x[:, 0]), r["v_inflection"], atol=0.05
    )


def test_passive_step_matches_the_dendrites_table():
    # PYTHON.md / figures/dendrites: lumped model, rest-rebalanced, 2.4 / 4.7
    # MOhm

    r = measure.passive_step(
        "multi", mem=multi.membrane(-68.0, rebalance_rest=True)
    )
    assert np.isclose(r["rin_steady"], 2.4307, rtol=1e-4)
    assert np.isclose(r["rin_peak"], 4.6815, rtol=1e-4)
    assert np.isclose(r["tau"], 0.1191, atol=1e-4)


def test_count_rising_edges():

    assert spiking.count_rising_edges([0, 1, 1, 0, 1, 0, 0, 1]) == 3
    assert (
        spiking.count_rising_edges([1, 1, 0]) == 0
    )  # a marker at the first sample isn't an edge


def test_stimulus_at_time_zero_skips_the_prestimulus_cache(fresh_quiet_cache):
    t, y = msoaxon.mso_axon("step", 0, 5, 500.0, 3, "active-full", 5, -68.0, 1)
    assert not _solve._QUIET_CACHE  # nothing stimulus-free to share
    assert y[-1, 0] > -68.0  # the step is on from the start


def test_pools_never_fork():

    assert _parallel.pool_context("fork").get_start_method() == "forkserver"
    assert _parallel.pool_context("spawn").get_start_method() == "spawn"
    assert (
        _parallel.pool_context().get_start_method() != "fork"
    )  # this platform's choice
    # the forkserver context really runs tasks (it's what Linux gets)
    with futures.ProcessPoolExecutor(
        2, mp_context=_parallel.pool_context("fork")
    ) as pool:
        assert _parallel.map_tasks(abs, [-1, -2, -3], executor=pool) == [
            1,
            2,
            3,
        ]


# --- settled stop (_solve.settled_event, search-speed A1) ---------------------


def test_input_end_per_stimulus():

    tau = (0.1, 0.18)
    assert _solve.input_end("EPSG", 5, 10, 20, tau) == pytest.approx(
        5 + 3.6
    )  # onset, not the cut-off
    assert _solve.input_end("EPSGpair", 5, 5.3, 20, tau) == pytest.approx(
        5.3 + 3.6
    )
    assert _solve.input_end(
        "EPSGbilateral", 5, 5.3, 20, (0.1, 0.3)
    ) == pytest.approx(5.3 + 6.0)
    assert _solve.input_end("step", 5, 105, 110, tau) == 105
    assert (
        _solve.input_end("ramp2", 5, 5.5, 20, tau) == 20
    )  # holds its current: never ends


@pytest.mark.parametrize("model", ["multi", "two"])
def test_settled_stop_ends_quiet_runs_early_with_the_same_answer(model):

    args = (
        model,
        "EPSGpair",
        5,
        5.3,
        20.0,
        3,
        20,
        -68.0,
        1,
    )  # far below threshold
    t_full, _ = _dispatch.run_model(*args, stop_on_spike=10)
    t, x = _dispatch.run_model(*args, stop_on_spike=10, stop_when_settled=True)
    t_off = _solve.input_end("EPSGpair", 5, 5.3, 20, (0.1, 0.18))
    # stopped after the input ended (the root-finder may land a few 1e-15 ms
    # early when the condition already holds at t_off, where the event switches)
    assert t_full[-1] == 20 and t_off - 1e-9 <= t[-1] < 20
    axon = 2 if model == "multi" else 1
    assert (
        x[-1, axon] - x[-1, 0] < 5
    )  # below factor / 2: a settled stop, not a spike
    # identical solver steps up to the stop (the last point is the event time
    # itself)
    np.testing.assert_array_equal(t[:-1], t_full[: len(t) - 1])
    assert _dispatch.spikes(*args, factor=10) is False


@pytest.mark.parametrize("model", ["multi", "two"])
def test_spike_stops_are_still_spikes(model):
    # a spike happens before the input ends, so the settled event can't touch it

    args = (model, "EPSGpair", 5, 5.3, 150.0, 3, 20, -68.0, 1)
    t1, y1 = _dispatch.run_model(*args, stop_on_spike=10)
    t2, y2 = _dispatch.run_model(
        *args, stop_on_spike=10, stop_when_settled=True
    )
    np.testing.assert_array_equal(y1, y2)
    axon = 2 if model == "multi" else 1
    assert y2[-1, axon] - y2[-1, 0] == pytest.approx(10, abs=1e-6)
    assert _dispatch.spikes(*args, factor=10) is True


def test_runs_that_never_settle_are_unchanged():

    args = (
        "multi",
        "ramp2",
        5,
        5.5,
        3000.0,
        3,
        20,
        -68.0,
        1,
    )  # ramp2 never ends
    t1, y1 = _dispatch.run_model(*args, stop_on_spike=30)
    t2, y2 = _dispatch.run_model(
        *args, stop_on_spike=30, stop_when_settled=True
    )
    np.testing.assert_array_equal(y1, y2)
    with pytest.raises(ValueError, match="needs stop_on_spike"):
        _dispatch.run_model(*args, stop_when_settled=True)


# --- coincidence windows found directly (coincidence.window, C1) -------------


def test_window_brackets_the_spike_boundary():

    w, th0 = coincidence.window("two", 0.03)
    amp, d = 1.03 * th0, w / 2
    site = ("EPSGpair", 1, None)
    assert coincidence._spikes(
        "two", amp, d - 2e-4, 3, -68.0, synaptic.EPSG_TAU, site, {}
    )  # still spikes just inside
    assert not coincidence._spikes(
        "two", amp, d + 2e-4, 3, -68.0, synaptic.EPSG_TAU, site, {}
    )  # not just outside
    assert 0.19 < w < 0.21  # PYTHON.md: ~199 us for TwoCpt.m at 3%


def test_window_agrees_with_the_threshold_curve():
    # the grid method interpolates linearly between 20 us points and so sits a
    # little low; on a fine grid the two agree closely

    w, th0 = coincidence.window("two", 0.03)
    delays = np.round(np.arange(0, 0.2001, 0.005), 4)
    grid = coincidence.half_width(
        delays, coincidence.threshold_curve("two", delays), 0.03
    )
    assert abs(w - grid) < 1e-3  # within 1 us


def test_window_reuses_threshold0_and_reports_no_window():

    w, th0 = coincidence.window("two", 0.03)
    assert coincidence.window("two", 0.03, threshold0=th0) == (w, th0)
    w_big, _ = coincidence.window(
        "two", 1.5, threshold0=th0
    )  # 2.5x th0 spikes even on the plateau
    assert np.isnan(w_big)


@pytest.mark.parametrize("th0", [np.inf, np.nan])
def test_window_without_a_finite_threshold_reports_no_window(th0):
    w, got = coincidence.window("two", 0.03, threshold0=th0)
    assert np.isnan(w) and (np.isnan(got) if np.isnan(th0) else got == th0)


def test_window_when_nothing_fires_returns_no_window():
    # a fixed 30 mV criterion is out of reach for Goldwyn's calibration: the
    # peak stays at 16-20 mV all the way up the search to its ceiling, and
    # the flat curve must not make the prediction warn about overflow
    kw = {k: v for k, v in two.GOLDWYN_2019.items() if k != "v0"}
    w, th0 = coincidence.window(
        "two", 0.03, v0=-58, model_kw=kw, criterion=30.0
    )
    assert np.isnan(w) and th0 == np.inf


# --- epsg_unitary far from onset (synaptic.py) --------------------------------

# epsg_unitary(t) for t > 0, captured from the code before the t <= -71 ms fix
# (which clamps t before the exponentials and must not move a bit here)
EPSG_UNITARY_AFTER_ONSET = {
    0.05: 0.7080474207643337,
    0.1: 0.9657736997044162,
    0.18: 0.9503239337141989,
    1.0: 0.01792241032842484,
    5.0: 4.050776728174052e-12,
}

EPSG_FAR_FROM_ONSET = pytest.mark.xfail(
    strict=True, reason="epsg_unitary NaN, fixed by the single-source PR"
)


def test_epsg_unitary_after_onset_is_bit_for_bit_unchanged():
    for t, expected in EPSG_UNITARY_AFTER_ONSET.items():
        assert synaptic.epsg_unitary(t) == expected
    got = synaptic.epsg_unitary(np.array(list(EPSG_UNITARY_AFTER_ONSET)))
    np.testing.assert_array_equal(got, list(EPSG_UNITARY_AFTER_ONSET.values()))


def test_epsg_unitary_just_before_onset_is_a_signed_zero():
    # (t > 0) * (negative difference of exponentials): -0.0 for t < 0, +0.0 at
    # the onset. Reachable today, so the fix must keep these bits.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = synaptic.epsg_unitary(np.array([-70.0, -1e-9, 0.0]))
    assert np.all(got == 0)
    assert list(np.signbit(got)) == [True, True, False]


@EPSG_FAR_FROM_ONSET
def test_epsg_unitary_far_before_onset_is_zero_not_nan():
    # exp(-t / 0.1) overflows for t < about -71 ms, and inf * 0 is NaN
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = synaptic.epsg_unitary(
            np.array([-200.0, -72.0, -70.0, -1e-9, 0.0])
        )
    assert not np.isnan(got).any()
    assert np.all(got == 0)
    assert list(np.signbit(got)) == [True, True, True, True, False]


@EPSG_FAR_FROM_ONSET
@pytest.mark.parametrize("model", ["multi", "two"])
def test_epsg_pair_with_the_second_onset_80_ms_later_runs(model):
    # at the second onset the first EPSG is 80 ms old; the second EPSG's
    # waveform is evaluated 80 ms before its own onset until then
    args = ("EPSGpair", 5, 85, 1.0, 3, "passive", 90, -68.0, 1)
    if model == "multi":
        t, y = msoaxon.mso_axon(*args, mode="exact")
    else:
        t, y = msoaxon.two_cpt(*args)
    assert t[-1] == 90
    assert np.isfinite(y).all()
