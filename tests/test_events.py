"""The "Events" stimulus (issue #78): a list of synaptic conductance events.

Each event is an EPSG-shaped conductance with its own onset, compartment,
peak, kinetics and reversal. The existing EPSG stimuli keep their own code
paths (their float order is part of the golden traces), so the tests here
compare right-hand sides rather than traces: a last-bit difference moves
the solver's steps.
"""

import platform
import sys
import warnings

import numpy as np
import pytest

from msoaxon import (
    Event,
    SynParams,
    _common,
    _fast,
    _solve,
    multi,
    synaptic,
    two,
)

T_END = 20.0
MYOGA = (0.1, 0.3)
IPSG = (0.28, 1.5)  # Roberts et al 2013 rise, Myoga et al 2014 decay


def stim(stim_type, start, stop, amp, **syn_kw):
    return _common.stimulus(
        stim_type, start, stop, amp, T_END, SynParams(t_end=T_END, **syn_kw)
    )


def events(*ev):
    return stim("Events", 0.0, 0.0, 1.0, events=ev)


ACTIVE = tuple(bool(g) for g in multi.ACTIVE_GATES["active-full"])
MEM = multi.membrane(-68.0)
P = two.get_params(-68.0, 3, 1, "active-full")


def _states(rng, n=45, count=3):
    yield np.concatenate([np.full(n, -68.0), np.full(6 * n, 0.3)])
    for _ in range(count):
        yield np.concatenate([rng.uniform(-90, 40, n), rng.random(6 * n)])


# --- the stimulus record ------------------------------------------------------


def test_stimulus_derives_its_range_from_the_onsets():
    s = events(Event(7.0, 1, 10.0), Event(5.0, 3, 20.0), Event(6.0, 1, 5.0))
    assert (s.start, s.stop) == (5.0, 7.0)
    assert s.events["k"].tolist() == [0, 2, 0]
    assert s.events["g"].tolist() == [10.0, 20.0, 5.0]
    assert not s.events["t"].flags.writeable
    assert s.I == 1.0


def test_event_defaults_and_inv_peak_per_event():
    e = Event(5.0)
    assert (e.node, e.g, e.tau, e.e_rev) == (1, 1.0, (0.1, 0.18), 0.0)
    s = events(Event(5.0), Event(5.0, tau=MYOGA))
    assert s.events["inv_peak"][0] == 1 / 0.21317
    assert s.events["inv_peak"][1] == synaptic.epsg_inv_peak(MYOGA)
    assert s.events["inv_peak"][1] != s.events["inv_peak"][0]


def test_other_stimuli_carry_no_events():
    assert stim("EPSG", 5.0, 10.0, 60.0).events is None
    assert stim("step", 5.0, 10.0, 60.0, events=(Event(5.0),)).events is None


def test_an_empty_list_is_rejected_at_construction():
    with pytest.raises(ValueError, match="at least one Event"):
        stim("Events", 0.0, 0.0, 1.0)


# --- validation -------------------------------------------------------------


@pytest.mark.parametrize(
    "bad,match",
    [
        (Event(5.0, 0), r"node must be 1\.\.45"),
        (Event(5.0, 46), r"node must be 1\.\.45"),
        (Event(-0.1), r"t must be in \[0, t_end"),
        (Event(T_END), r"t must be in \[0, t_end"),
        (Event(5.0, tau=(0.2, 0.2)), r"0 < rise < decay"),
        (Event(5.0, tau=(0.3, 0.1)), r"0 < rise < decay"),
        (Event(5.0, tau=(0.0, 0.1)), r"0 < rise < decay"),
        (Event(5.0, g=float("nan")), r"finite g and e_rev"),
        (Event(5.0, e_rev=float("inf")), r"finite g and e_rev"),
    ],
)
@pytest.mark.parametrize("model", ["multi", "two"])
def test_models_reject_bad_events_before_running(model, bad, match):
    run = multi.mso_axon if model == "multi" else two.two_cpt
    syn = SynParams(t_end=T_END, events=(Event(5.0), bad))
    with pytest.raises(ValueError, match=match):
        run("Events", 0.0, 0.0, 1.0, 3, "active-full", T_END, -68.0, 1, syn)


def test_models_reject_a_missing_list():
    with pytest.raises(ValueError, match="at least one Event"):
        two.two_cpt("Events", 0.0, 0.0, 1.0, 3, "active-full", T_END, -68, 1)


# --- the right-hand sides -----------------------------------------------------


def test_two_soma_event_is_the_epsg_bit_for_bit():
    # same expression and float order; two_cpt never cuts an EPSG off
    s_ev = events(Event(5.0, 1, 60.0))
    s_ep = stim("EPSG", 5.0, 10.0, 60.0)
    rng = np.random.default_rng(7)
    for t in (4.9, 5.0, 5.05, 5.3, 7.0, 15.0):
        for _ in range(4):
            x = np.concatenate([rng.uniform(-90, 40, 2), rng.random(9)])
            a = np.asarray(two._rhs(t, x, P, "Events", s_ev))
            b = np.asarray(two._rhs(t, x, P, "EPSG", s_ep))
            assert a.tobytes() == b.tobytes()


def test_two_soma_pair_matches_epsgpair_closely():
    s_ev = events(Event(5.0, 1, 60.0), Event(5.3, 1, 60.0))
    s_pr = stim("EPSGpair", 5.0, 5.3, 60.0)
    rng = np.random.default_rng(8)
    for t in (5.1, 5.3, 5.35, 6.0, 12.0):
        x = np.concatenate([rng.uniform(-90, 40, 2), rng.random(9)])
        a = np.asarray(two._rhs(t, x, P, "Events", s_ev))
        b = np.asarray(two._rhs(t, x, P, "EPSGpair", s_pr))
        np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)


def test_two_axon_event_drives_the_axon_at_its_own_voltage():
    s = events(Event(5.0, 3, 60.0))
    x = np.array([-60.0, -40.0, 0.3, 0.4, 0.2, 0.1, 0.1, 0.5, 0.2, 0.3, 0.3])
    i1, i2 = two.event_currents(5.1, x[0], x[1], s)
    assert i1 == 0.0
    te = 5.1 - 5.0
    wave = (1 / 0.21317) * (np.exp(-te / 0.18) - np.exp(-te / 0.1))
    assert i2 == pytest.approx(60.0 * (0.0 - x[1]) * wave)


def test_multi_single_event_is_the_epsgpair_half_bit_for_bit():
    # one event at the input node against EPSGpair with its second EPSG
    # beyond the run: the pair's second waveform is an exact 0.0 there
    s_ev = events(Event(5.0, 3, 60.0))
    s_pr = stim("EPSGpair", 5.0, 19.0, 60.0)
    rng = np.random.default_rng(9)
    for t in (4.9, 5.0, 5.05, 5.3, 7.0, 15.0):
        for x in _states(rng):
            a = multi._rhs(t, x, "Events", s_ev, 3, ACTIVE, MEM)
            b = multi._rhs(t, x, "EPSGpair", s_pr, 3, ACTIVE, MEM)
            assert a.tobytes() == b.tobytes()


def test_multi_pair_matches_epsgpair_closely():
    s_ev = events(Event(5.0, 1, 60.0), Event(5.2, 1, 60.0))
    s_pr = stim("EPSGpair", 5.0, 5.2, 60.0)
    rng = np.random.default_rng(10)
    for t in (5.1, 5.2, 5.25, 6.0, 12.0):
        for x in _states(rng):
            a = multi._rhs(t, x, "Events", s_ev, 1, ACTIVE, MEM)
            b = multi._rhs(t, x, "EPSGpair", s_pr, 1, ACTIVE, MEM)
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)


def test_multi_two_sites_match_epsgbilateral_closely():
    D = multi.with_dendrites()
    lat, med = int(D.lateral[2]), int(D.medial[2])
    s_ev = events(Event(5.0, lat, 60.0), Event(5.1, med, 60.0))
    s_bi = stim("EPSGbilateral", 5.0, 5.1, 60.0)
    mem = multi.membrane(-68.0, morph=D)
    rng = np.random.default_rng(11)
    for t in (5.05, 5.1, 5.2, 7.0):
        for x in _states(rng, n=D.n):
            a = multi._rhs(t, x, "Events", s_ev, 1, ACTIVE, mem, D)
            b = multi._rhs(
                t, x, "EPSGbilateral", s_bi, lat, ACTIVE, mem, D, med
            )
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)


def test_amplitude_scales_every_peak():
    one = events(Event(5.0, 1, 30.0), Event(5.2, 3, 40.0))
    two_x = stim(
        "Events",
        0.0,
        0.0,
        2.0,
        events=(Event(5.0, 1, 30.0), Event(5.2, 3, 40.0)),
    )
    x = np.concatenate([np.full(45, -60.0), np.full(270, 0.3)])
    a = multi.events_current(5.3, x[:45], one, multi.LUMPED.sa)
    b = multi.events_current(5.3, x[:45], two_x, multi.LUMPED.sa)
    assert np.count_nonzero(a) == 2
    np.testing.assert_allclose(b, 2 * a, rtol=1e-15)


def test_dv_is_the_voltage_block_of_rhs_bit_for_bit():
    s = events(
        Event(5.0, 1, 60.0),
        Event(5.2, 3, 40.0, MYOGA),
        Event(4.8, 1, 30, IPSG, -80),
    )
    rng = np.random.default_rng(12)
    for t in (4.79, 4.8, 5.0, 5.3, 9.0, 19.0):
        for x in _states(rng):
            want = multi._rhs(t, x, "Events", s, 1, ACTIVE, MEM)
            got = multi._dv(t, x, "Events", s, 1, MEM)
            assert got.tobytes() == want[:45].tobytes()
        for _ in range(3):
            x2 = np.concatenate([rng.uniform(-90, 40, 2), rng.random(9)])
            want = np.asarray(two._rhs(t, x2, P, "Events", s))
            got = np.asarray(two._dv(t, x2, P, "Events", s))
            assert got.tobytes() == want[:2].tobytes()


# --- the solver's contracts ---------------------------------------------------


def test_breakpoints_are_every_onset():
    s = events(Event(7.0), Event(5.0, 3), Event(5.0), Event(12.5, 5))
    cuts = _solve.breakpoints("Events", s.start, s.stop, T_END, s.events["t"])
    assert cuts == [5.0, 7.0, 12.5]
    assert _solve.breakpoints("Events", 5.0, 7.0, 6.0, [5.0, 7.0]) == [5.0]


def test_input_end_uses_the_slowest_decay():
    s = events(Event(5.0), Event(8.0, 1, 30.0, IPSG, -80.0))
    got = _solve.input_end(
        "Events", s.start, s.stop, T_END, (0.1, 0.18), s.events["decay"]
    )
    assert got == pytest.approx(8.0 + 20 * 1.5)
    assert _solve.input_end("Events", 5.0, 8.0, T_END, (0.1, 0.18)) == T_END


@pytest.mark.parametrize("model", ["multi", "two"])
def test_runs_share_the_quiet_prefix_and_the_settled_stop_works(
    fresh_quiet_cache, model
):
    run = multi.mso_axon if model == "multi" else two.two_cpt
    kw = {"mode": "exact"} if model == "multi" else {}
    # a subthreshold pair: the settled stop ends the run well before t_end
    syn = SynParams(
        t_end=T_END, events=(Event(6.0, 1, 10.0), Event(6.2, 1, 10.0))
    )
    args = ("Events", 0.0, 0.0, 1.0, 3, "active-full", T_END, -68.0, 1, syn)
    t, x = run(*args, stop_on_spike=10.0, stop_when_settled=True, **kw)
    assert 6.2 < t[-1] < T_END
    col = 2 if model == "multi" else 1
    assert x[-1, col] - x[-1, 0] < 10.0 * _solve.SETTLED_FRACTION
    # the stimulus-free prefix up to the first onset is cached once and
    # shared by the next run with the same settings
    assert len(_solve._QUIET_CACHE) == 1
    run(*args, **kw)
    assert len(_solve._QUIET_CACHE) == 1


BITWISE = sys.platform == "darwin" and platform.machine() == "arm64"
needs_ext = pytest.mark.skipif(
    not _fast.available(), reason="needs the fast extra (msoaxon_fast)"
)
# the golden traces' tolerance (tests/test_golden.py), rtol and atol [mV]
TOLERANCE = {"rtol": 1e-4, "atol": 5e-3}
SAMPLE_TIMES = (1.0, 5.1, 5.3, 5.6, 6.0, 7.0, 9.0, 12.0)
# no Events golden trace exists, so the modes are compared with each other
# as test_golden compares them with the traces: soma and axon sampled by
# linear interpolation. A pair at the EPSGpair golden's 100 nS, a stronger one
# that drives the axon to a partial spike, and a mixed list: a leading IPSG, a
# pair and an axonal event
MODE_CASES = {
    "pair": (Event(5.0, 1, 100.0), Event(5.2, 1, 100.0)),
    "strong": (Event(5.0, 1, 300.0), Event(5.2, 1, 300.0)),
    "mixed": (
        Event(4.7, 1, 30.0, IPSG, -80.0),
        Event(5.0, 1, 50.0),
        Event(5.2, 1, 50.0),
        Event(5.1, 3, 20.0, MYOGA),
    ),
}


def _mode_run(events, mode):
    syn = SynParams(t_end=T_END, events=events)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", _fast.DriftWarning)
        return multi.mso_axon(
            "Events",
            0.0,
            0.0,
            1.0,
            3,
            "active-full",
            T_END,
            -68.0,
            1,
            syn,
            mode=mode,
        )


@needs_ext
def test_fast_mode_runs_events_in_rust():
    events = MODE_CASES["pair"]
    # nothing but the DriftWarning is raised: no fallback to "exact"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.simplefilter("ignore", _fast.DriftWarning)
        t_f, x_f = _mode_run(events, "fast")
    t_e, x_e = _mode_run(events, "exact")
    if BITWISE:  # numpy and Rust share a libm here (tests/test_fast.py)
        assert t_f.tobytes() == t_e.tobytes()
        assert x_f.tobytes() == x_e.tobytes()
    else:
        assert t_f[-1] == t_e[-1] == T_END
        np.testing.assert_allclose(x_f[-1], x_e[-1], **TOLERANCE)


@needs_ext
@pytest.mark.parametrize("mode", ["fast", "native"])
@pytest.mark.parametrize("case", sorted(MODE_CASES))
def test_events_in_every_mode_agree_at_the_golden_tolerance(case, mode):
    events = MODE_CASES[case]
    t_e, x_e = _mode_run(events, "exact")
    t_m, x_m = _mode_run(events, mode)
    assert t_m[-1] == t_e[-1] == T_END
    for col in (0, 2):  # the soma and the third compartment
        want = np.interp(SAMPLE_TIMES, t_e, x_e[:, col])
        got = np.interp(SAMPLE_TIMES, t_m, x_m[:, col])
        np.testing.assert_allclose(got, want, **TOLERANCE)


# --- a demonstration ----------------------------------------------------------


def test_a_leading_glycinergic_event_lowers_the_pair_response():
    # Roberts et al 2013: IPSPs lead EPSPs by ~0.3 ms; Myoga et al 2014:
    # IPSG decay 1.5 ms, 30 nS. Reversal below rest makes it inhibitory.
    pair = (Event(5.0, 1, 50.0), Event(5.2, 1, 50.0))
    ipsg = Event(4.7, 1, 30.0, IPSG, -80.0)
    peaks = []
    for ev in (pair, (ipsg, *pair)):
        syn = SynParams(t_end=T_END, events=ev)
        t, x = two.two_cpt(
            "Events", 0, 0, 1.0, 3, "active-full", T_END, -68, 1, syn
        )
        peaks.append(float((x[:, 1] - x[:, 0]).max()))
    without, with_inh = peaks
    assert with_inh < without
