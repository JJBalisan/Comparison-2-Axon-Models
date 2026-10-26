"""The event list reaching the searches, and the protocols on it (issue #78).

The builders are checked on their statistics and shapes, the plumbing against
the EPSG-pair searches it must agree with, and the protocols on the
two-compartment model, which runs in tens of milliseconds.
"""

from typing import Any

import numpy as np
import pytest

from msoaxon import Event, coincidence, measure, somatic, synaptic

TWO: dict[str, Any] = {"workers": 1}  # serial: the tests are small


# --- builders -----------------------------------------------------------------


def test_fibre_events_are_seeded_and_sorted():
    a = synaptic.fibre_events(6, 3, jitter=0.3, g_sd=4.0)
    assert a == synaptic.fibre_events(6, 3, jitter=0.3, g_sd=4.0)
    assert a != synaptic.fibre_events(6, 4, jitter=0.3, g_sd=4.0)
    assert [e.t for e in a] == sorted(e.t for e in a)
    assert len(a) == 6 and all(e.e_rev == 0.0 and e.node == 1 for e in a)


def test_fibre_events_statistics():
    rng = np.random.default_rng(0)
    ev = synaptic.fibre_events(20000, rng, t0=5.0, g=37.0, g_sd=4.0, jitter=0.3)
    t = np.array([e.t for e in ev])
    g = np.array([e.g for e in ev])
    assert abs(t.mean() - 5.0) < 0.01 and abs(t.std() - 0.3) < 0.01
    assert abs(g.mean() - 37.0) < 0.2 and abs(g.std() - 4.0) < 0.2
    thinned = synaptic.fibre_events(20000, rng, release_p=0.45)
    assert abs(len(thinned) / 20000 - 0.45) < 0.015


def test_fibre_events_without_jitter_share_an_onset():
    ev = synaptic.fibre_events(4, 0, t0=7.0)
    assert {e.t for e in ev} == {7.0} and {e.g for e in ev} == {37.0}


def test_fibre_events_poisson_counts_and_window():
    ev = synaptic.fibre_events(2000, 1, t0=5.0, rate=300.0, duration=10.0)
    t = np.array([e.t for e in ev])
    assert abs(len(ev) / 2000 - 3.0) < 0.15  # 300 Hz over 10 ms
    assert t.min() >= 5.0 and t.max() < 15.0


def test_fibre_events_clip_onsets_at_zero_and_peaks_at_zero():
    ev = synaptic.fibre_events(200, 0, t0=0.0, g=1.0, g_sd=5.0, jitter=1.0)
    assert min(e.t for e in ev) == 0.0 and min(e.g for e in ev) == 0.0


def test_fibre_events_can_release_nothing():
    assert synaptic.fibre_events(10, 0, release_p=0.0) == ()


@pytest.mark.parametrize(
    "kw",
    [
        {"n": -1},
        {"n": 2, "release_p": 1.5},
        {"n": 2, "rate": 100.0},
        {"n": 2, "rate": 100.0, "duration": -1.0},
    ],
)
def test_fibre_events_reject_bad_arguments(kw):
    with pytest.raises(ValueError):
        synaptic.fibre_events(rng=0, **kw)


def test_bilateral_events_offset_the_second_side():
    ev = synaptic.bilateral_events(0.4, 5, n=3, nodes=(1, 2), t0=5.0)
    assert len(ev) == 6 and [e.t for e in ev] == sorted(e.t for e in ev)
    assert [(e.t, e.node) for e in ev] == [
        (5.0, 1),
        (5.0, 1),
        (5.0, 1),
        (5.4, 2),
        (5.4, 2),
        (5.4, 2),
    ]


def test_current_event_injects_its_amplitude_whatever_the_voltage():
    ev = synaptic.event_arrays((synaptic.current_event(5.0, 2000.0),))
    peak = 5.0 + 0.2  # near the peak of the default waveform
    for v in (-90.0, -68.0, -20.0):
        i = synaptic.event_currents(peak, np.array([v]), ev, 1.0)[0]
        w = synaptic.epsg_unitary(peak - 5.0)
        assert i == pytest.approx(2000.0 * w, rel=1e-3)


def test_epsc_train_is_regular_and_current_like():
    train = synaptic.epsc_train(4, 200.0, -50.0, t0=6.0, node=2)
    assert [e.t for e in train] == pytest.approx([6.0, 11.0, 16.0, 21.0])
    assert {e.node for e in train} == {2}
    assert all(e.e_rev == synaptic.CURRENT_E_REV for e in train)
    assert train[0].g * synaptic.CURRENT_E_REV == pytest.approx(-50.0)


def test_fixed_scales_every_peak_and_has_no_delay():
    src = synaptic.as_source(
        (Event(5.0, 1, 2.0), Event(6.0, 3, 1.0, e_rev=-80.0))
    )
    assert isinstance(src, synaptic.Fixed)
    out = src(123.0, 10.0)  # the delay is ignored
    assert [e.g for e in out] == [20.0, 10.0]
    assert [(e.t, e.node, e.e_rev) for e in out] == [
        (5.0, 1, 0.0),
        (6.0, 3, -80.0),
    ]


def test_as_source_passes_sources_through_and_rejects_empty():
    src = synaptic.PairSource()
    assert synaptic.as_source(src) is src
    with pytest.raises(ValueError, match="at least one Event"):
        synaptic.as_source(())


def test_pair_source_builds_the_pair_and_its_pre_pulses():
    src = synaptic.PairSource(
        start=5.0, node=1, node2=4, pre=((-1.0, 30.0), (-0.5, 10.0))
    )
    ev = src(0.25, 50.0)
    assert [(e.t, e.node, e.g) for e in ev] == [
        (4.0, 1, 30.0),
        (4.5, 1, 10.0),
        (5.0, 1, 50.0),
        (5.25, 4, 50.0),
    ]
    assert src(-0.25, 50.0)[-1].t == 4.75  # a negative delay puts it first


def test_event_syn_runs_the_source():
    syn = synaptic.event_syn(synaptic.PairSource(), 0.1, 40.0, 20.0)
    assert syn.t_end == 20.0 and [e.g for e in syn.events] == [40.0, 40.0]


# --- the searches take an event list ----------------------------------------


def test_threshold_of_an_event_pair_is_the_pair_threshold():
    ref = coincidence.threshold("two", 0.1, rel_tol=1e-3)
    via_source = coincidence.threshold(
        "two", 0.1, rel_tol=1e-3, events=synaptic.PairSource()
    )
    assert via_source == pytest.approx(ref, rel=2e-3)


def test_threshold_of_a_fixed_tuple_scales_it():
    ref = coincidence.threshold("two", 0.0, rel_tol=1e-3)
    pair = (Event(5.0, 1, 1.0), Event(5.0, 1, 1.0))
    th = coincidence.threshold("two", 0.0, rel_tol=1e-3, events=pair)
    assert th == pytest.approx(ref, rel=2e-3)
    # twice the peaks, half the scale
    doubled = tuple(Event(e.t, e.node, 2.0) for e in pair)
    half = coincidence.threshold("two", 0.0, rel_tol=1e-3, events=doubled)
    assert half == pytest.approx(ref / 2, rel=2e-3)


def test_a_fixed_tuple_has_no_delay():
    pair = (Event(5.0, 1, 1.0),)
    with pytest.raises(ValueError, match="no delay"):
        coincidence.threshold("two", 0.1, events=pair)
    with pytest.raises(ValueError, match="no delay"):
        coincidence.threshold_curve("two", [0.0, 0.1], events=pair, **TWO)
    with pytest.raises(ValueError, match="no delay"):
        coincidence.probability_trials(
            "two", [0.1], 50.0, n_trials=1, events=pair, **TWO
        )


def test_threshold_curve_takes_a_source():
    delays = [0.0, 0.2]
    ref = coincidence.threshold_curve("two", delays, rel_tol=1e-3, **TWO)
    got = coincidence.threshold_curve(
        "two", delays, rel_tol=1e-3, events=synaptic.PairSource(), **TWO
    )
    np.testing.assert_allclose(got, ref, rtol=2e-3)


def test_a_fixed_peak_in_a_source_is_not_scaled():
    # a 60 nS IPSG ahead of the pair raises the threshold; the search varies
    # only the excitation
    inhibited = synaptic.PairSource(pre=((-0.5, 60.0),))
    th = coincidence.threshold("two", 0.0, rel_tol=1e-3, events=inhibited)
    ref = coincidence.threshold("two", 0.0, rel_tol=1e-3)
    assert th != pytest.approx(ref, rel=0.01)


def test_window_of_an_event_pair_matches_the_symmetric_window():
    width, th0 = coincidence.window(
        "two", 0.03, epsg_tau=(0.1, 0.3), rel_tol=1e-3, delay_tol=1e-3
    )
    edges = coincidence.window_edges(
        "two",
        0.03,
        epsg_tau=(0.1, 0.3),
        rel_tol=1e-3,
        delay_tol=1e-3,
        events=synaptic.PairSource(tau=(0.1, 0.3)),
    )
    assert edges.width == pytest.approx(width, abs=4e-3)  # 4 us
    assert abs(edges.centre) < 3e-3 and abs(edges.best) < 5e-3
    assert edges.th0 == pytest.approx(th0, rel=5e-3)
    via_window = coincidence.window(
        "two",
        0.03,
        epsg_tau=(0.1, 0.3),
        rel_tol=1e-3,
        delay_tol=1e-3,
        events=synaptic.PairSource(tau=(0.1, 0.3)),
    )
    assert via_window == pytest.approx((edges.width, edges.th0))


def test_window_edges_follow_an_asymmetric_source():
    # a pre-pulse before the first EPSG only: the second, later one benefits
    # less, so the window is not centred on 0
    edges = coincidence.window_edges(
        "two",
        0.03,
        rel_tol=1e-3,
        delay_tol=1e-3,
        events=synaptic.PairSource(pre=((-0.4, 37.0),)),
    )
    assert edges.lo < edges.best < edges.hi
    assert edges.width > 0.05 and abs(edges.centre) > 5e-3


def test_window_edges_with_a_reference_too_low_to_spike_are_nan():
    edges = coincidence.window_edges(
        "two",
        0.03,
        rel_tol=1e-3,
        reference=(0.0, 30.0),  # 31 nS pairs do not fire the cell
        events=synaptic.PairSource(),
    )
    assert np.isnan([edges.lo, edges.hi, edges.width]).all()
    assert edges.best == 0.0 and edges.th0 == 30.0


def test_window_rejects_threshold0_with_events():
    with pytest.raises(ValueError, match="threshold0"):
        coincidence.window(
            "two", 0.03, threshold0=48.0, events=synaptic.PairSource()
        )


def test_best_delay_of_the_symmetric_pair_is_zero():
    best, th = coincidence.best_delay(
        "two", rel_tol=1e-3, events=synaptic.PairSource()
    )
    assert abs(best) < 5e-3
    assert th == pytest.approx(
        coincidence.threshold("two", 0.0, rel_tol=1e-3), rel=2e-3
    )


def test_probability_trials_of_a_source():
    th = coincidence.threshold(
        "two", 0.0, rel_tol=1e-3, events=synaptic.PairSource()
    )
    p = coincidence.probability_trials(
        "two",
        [0.0, 0.8],
        1.05 * th,
        n_trials=6,
        amp_cv=0.0,
        jitter=0.005,
        events=synaptic.PairSource(),
        **TWO,
    )
    assert p.tolist() == [1.0, 0.0]


def test_probability_trials_with_a_source_is_seeded():
    kw = dict(n_trials=8, jitter=0.1, events=synaptic.PairSource(), **TWO)
    a = coincidence.probability_trials("two", [0.3], 50.0, seed=2, **kw)
    assert (
        a.tolist()
        == coincidence.probability_trials(
            "two", [0.3], 50.0, seed=2, **kw
        ).tolist()
    )


def test_trial_spikes_runs_each_list():
    pair = (Event(5.0, 1, 60.0), Event(5.0, 1, 60.0))
    weak = (Event(5.0, 1, 20.0),)
    hits = coincidence.trial_spikes("two", [pair, weak, pair], **TWO)
    assert hits.tolist() == [True, False, True]
    assert coincidence.trial_spikes("two", [], **TWO).tolist() == []


def test_trial_spikes_rejects_an_empty_list_before_running():
    with pytest.raises(ValueError, match="at least one Event"):
        coincidence.trial_spikes("two", [()], **TWO)


# --- the protocols --------------------------------------------------------


def test_a_pre_pulse_raises_the_pair_threshold():
    # Franken et al 2015: a preceding EPSG raises the threshold, by less the
    # closer it comes to the pair's own depolarisation
    th0, th = coincidence.prepulse_thresholds(
        "two", pre_dts=(-1.0, -0.5), rel_tol=1e-3, **TWO
    )
    assert th0 == pytest.approx(47.7, abs=0.2)
    assert th[0] - th0 > 5.0 and th[0] > th[1] > th0


def test_bilateral_probability_is_seeded_and_shaped():
    kw = dict(n_trials=4, n=2, jitter=0.1, **TWO)
    a = coincidence.bilateral_probability("two", [0.0, 0.5], repeats=2, **kw)
    assert a.shape == (2, 2) and ((a == 0) | (a > 0)).all()
    b = coincidence.bilateral_probability("two", [0.0, 0.5], repeats=2, **kw)
    np.testing.assert_array_equal(a, b)
    c = coincidence.bilateral_probability(
        "two", [0.0, 0.5], repeats=2, seed=1, **kw
    )
    assert c.shape == a.shape


def test_fibre_count_sets_the_spike_probability():
    # Couchman et al 2010: 2 to 4 fibres of 37 nS fire the cell
    delays = [0.0]
    p = {
        n: coincidence.bilateral_probability(
            "two", delays, n_trials=20, n=n, jitter=0.0, **TWO
        )[0, 0]
        for n in (1, 2)
    }
    assert p[1] == 0.0 and p[2] == 1.0  # 2 x 37 nS alone is below threshold


def test_fibre_window_trials_summarises_each_repeat():
    delays = np.round(np.linspace(-0.4, 0.4, 5), 3)
    prob, summaries = coincidence.fibre_window_trials(
        "two", delays, n_fibres=2, n_trials=3, repeats=2, jitter=0.1, **TWO
    )
    assert prob.shape == (2, 5) and len(summaries) == 2
    assert set(summaries[0]) == {
        "peak",
        "width",
        "centre",
        "centroid",
        "lo",
        "hi",
    }
    assert summaries[0]["peak"] == prob[0].max()


def test_itd_sweep_is_tuned_to_zero_itd():
    itds = [-0.5, 0.0, 0.5]
    prob, summary = coincidence.itd_sweep(
        "two", itds, n_trials=4, n_syn=40, g=3.4, **TWO
    )
    assert prob[1] > prob[0] and prob[1] > prob[2]
    assert summary["peak"] == prob.max()
    # ITD 0 puts both sides' synapses together; at +-0.5 ms they are apart
    assert prob[1] >= 0.5 and prob[0] == prob[2] == 0.0


def test_itd_sweep_defaults_are_the_small_casarez_sweep():
    # the real protocol (125 synapses, p 0.45, 400 trials, +-0.5 ms in 10 us
    # steps) is the same call with these arguments
    assert (
        coincidence.CASAREZ_N,
        coincidence.CASAREZ_P,
    ) == (125, 0.45)
    ev = synaptic.bilateral_events(
        0.0,
        0,
        n=coincidence.CASAREZ_N,
        release_p=coincidence.CASAREZ_P,
        g=coincidence.CASAREZ_G,
        jitter=0.1,
    )
    assert 90 < len(ev) < 135  # 250 synapses released with p 0.45


# --- measurements -----------------------------------------------------------


def test_curve_summary_of_a_gaussian():
    x = np.linspace(-1, 1, 401)
    sigma = 0.2
    s = measure.curve_summary(
        x, 0.9 * np.exp(-((x - 0.1) ** 2) / (2 * sigma**2))
    )
    fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma
    assert s["peak"] == pytest.approx(0.9)
    assert s["width"] == pytest.approx(fwhm, abs=2e-3)
    assert s["centre"] == pytest.approx(0.1, abs=2e-3)
    assert s["centroid"] == pytest.approx(0.1, abs=5e-3)


def test_curve_summary_edge_cases():
    x = np.array([0.0, 1.0, 2.0])
    flat = measure.curve_summary(x, [0.8, 0.8, 0.8])
    assert flat["peak"] == 0.8 and np.isnan(flat["width"])
    assert np.isnan(measure.curve_summary(x, [0.0, 0.0, 0.0])["centroid"])
    one_sided = measure.curve_summary(x, [1.0, 0.2, 0.0])
    assert np.isnan(one_sided["lo"]) and one_sided["hi"] == pytest.approx(0.625)


def test_train_rin_rest_is_the_passive_value():
    r = measure.train_rin("two", rate=200.0, n=12, amp=0.001)
    assert r["change"] == pytest.approx(0.0, abs=1e-3)
    assert r["depolarisation"] == pytest.approx(0.0, abs=1e-3)
    assert 2.3 < r["rin_rest"] < 3.2  # MOhm, the lumped soma


def test_train_rin_falls_as_the_train_depolarises():
    # KLT opens in the train, so the input resistance falls (the models do
    # not reproduce Khurana et al 2011's rise; see PYTHON.md)
    slow = measure.train_rin("two", rate=200.0, n=12, amp=3000.0)
    fast = measure.train_rin("two", rate=400.0, n=12, amp=3000.0)
    assert fast["depolarisation"] > slow["depolarisation"] > 0
    assert fast["change"] < slow["change"] < 0.01


def test_train_rin_rejects_a_zero_probe_and_a_short_train():
    with pytest.raises(ValueError, match="nonzero probe"):
        measure.train_rin("two", probe=0.0)
    with pytest.raises(ValueError, match="too short"):
        measure.train_rin("two", n=2, probe_tau=(10.0, 30.0))


# --- somatic ----------------------------------------------------------------


def test_somatic_search_takes_an_event_list():
    train = tuple(Event(5.0 + 2.5 * k, 1, 1.0) for k in range(20))
    assert not somatic.fires("two", 20.0, events=train)
    assert somatic.fires("two", 100.0, events=train)
    spiked, peak = somatic.response("two", 20.0, events=train)
    assert not spiked and 0 < peak < 10
    rheo = somatic.rheobase("two", guess=50.0, rel_tol=1e-3, events=train)
    assert 60 < rheo < 100
    assert somatic.fires("two", rheo * 1.001, events=train)
    assert not somatic.fires("two", rheo * 0.99, events=train)


def test_somatic_step_is_unchanged_without_events():
    assert somatic.fires("two", 4000.0) and not somatic.fires("two", 3000.0)
