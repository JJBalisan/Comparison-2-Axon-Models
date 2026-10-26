"""The IPSG event and the dt_inh axis of the coincidence window (issue #77)."""

import importlib.util
import pathlib
from typing import Any

import numpy as np
import pytest

from msoaxon import coincidence, constants, multi, synaptic

TWO: dict[str, Any] = {"workers": 1}
MYOGA = (0.1, 0.3)


# --- the event --------------------------------------------------------------


def test_ipsg_defaults_are_the_in_vitro_glycinergic_waveform():
    e = synaptic.ipsg(5.0)
    assert (e.t, e.node, e.g, e.tau, e.e_rev) == (
        5.0,
        1,
        40.0,
        (0.28, 1.5),
        -80.0,
    )
    assert synaptic.ipsg(5.0, 60.0, 3, (0.3, 1.85), -75.0) == synaptic.Event(
        5.0, 3, 60.0, (0.3, 1.85), -75.0
    )


def test_ipsg_waveform_peaks_at_its_peak_conductance():
    ev = synaptic.event_arrays((synaptic.ipsg(0.0, 40.0),))
    t = np.linspace(0.0, 8.0, 8001)
    cur = [synaptic.event_currents(x, np.array([-68.0]), ev, 1.0)[0] for x in t]
    # 40 nS * (-80 - (-68)) mV, the peak of a unit waveform
    assert min(cur) == pytest.approx(40.0 * -12.0, rel=1e-3)


def test_pair_source_puts_the_ipsg_against_the_first_epsg():
    inh = synaptic.Inhibition(dt=-0.2, g=45.0, node=2)
    src = synaptic.PairSource(start=5.0, tau=MYOGA, inh=inh)
    ev = src(0.3, 50.0)
    assert [(e.t, e.node, e.g, e.e_rev) for e in ev] == [
        (pytest.approx(4.8), 2, 45.0, -80.0),
        (5.0, 1, 50.0, 0.0),
        (5.3, 1, 50.0, 0.0),
    ]
    assert ev[0].tau == synaptic.IPSG_TAU and ev[1].tau == MYOGA
    # the IPSG keeps its peak and onset when the delay and amplitude change
    other = src(-0.3, 80.0)
    assert (other[0].t, other[0].g) == (ev[0].t, 45.0)


# --- the dt_inh axis --------------------------------------------------------


def test_events_and_dt_inh_are_exclusive():
    with pytest.raises(ValueError, match="not both"):
        coincidence.window(
            "two", 0.03, events=synaptic.PairSource(), dt_inh=0.0
        )
    with pytest.raises(ValueError, match="not both"):
        coincidence.probability_trials(
            "two",
            [0.0],
            50.0,
            n_trials=1,
            events=synaptic.PairSource(),
            inh=synaptic.Inhibition(),
            **TWO,
        )


def test_a_leading_ipsg_raises_the_threshold_and_a_late_one_does_not():
    kw = {"rel_tol": 1e-3, "epsg_tau": MYOGA}
    plain = coincidence.threshold("two", 0.0, **kw)
    led = coincidence.threshold(
        "two",
        0.0,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(-0.2, 40.0)
        ),
        **kw,
    )
    late = coincidence.threshold(
        "two",
        0.0,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(1.5, 40.0)
        ),
        **kw,
    )
    assert led > 1.15 * plain  # the IPSG's size is not scaled by the search
    assert late == pytest.approx(plain, rel=0.05)


def test_window_with_dt_inh_returns_the_window_edges_width():
    kw = {
        "epsg_tau": MYOGA,
        "rel_tol": 1e-3,
        "delay_tol": 1e-3,
        "dt_inh": -0.2,
    }
    width, th0 = coincidence.window("two", 0.03, **kw)
    edges = coincidence.window_edges(
        "two",
        0.03,
        epsg_tau=MYOGA,
        rel_tol=1e-3,
        delay_tol=1e-3,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(-0.2, 40.0)
        ),
    )
    assert (width, th0) == pytest.approx((edges.width, edges.th0))
    assert 0.15 < width < 0.3


def test_a_late_ipsg_leaves_the_window_where_it_was():
    plain = coincidence.window(
        "two", 0.03, epsg_tau=MYOGA, rel_tol=1e-3, delay_tol=1e-3
    )[0]
    kw = {"epsg_tau": MYOGA, "rel_tol": 1e-3, "delay_tol": 1e-3}
    far = coincidence.window_edges(
        "two",
        0.03,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(1.5, 40.0)
        ),
        **kw,
    )
    assert far.width == pytest.approx(plain, abs=5e-3)
    assert abs(far.centre) < 5e-3


def test_the_window_centre_moves_with_the_ipsg_timing():
    # the optimum shifts towards the later-arriving inhibition's blind side:
    # the sign flips between an IPSG that leads by 0.6 ms and one that lags
    kw = {"epsg_tau": MYOGA, "rel_tol": 1e-3, "delay_tol": 1e-3}
    lead = coincidence.window_edges(
        "two",
        0.03,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(-0.6, 40.0)
        ),
        **kw,
    )
    lag = coincidence.window_edges(
        "two",
        0.03,
        events=synaptic.PairSource(
            tau=MYOGA, inh=synaptic.Inhibition(0.2, 40.0)
        ),
        **kw,
    )
    assert lead.centre > 0.01 > -0.02 > lag.centre


def test_probability_trials_with_dt_inh():
    th = coincidence.threshold("two", 0.0, epsg_tau=MYOGA, rel_tol=1e-3)
    kw = dict(n_trials=5, amp_cv=0.0, jitter=0.005, epsg_tau=MYOGA, **TWO)
    # excitation 3% above the threshold without inhibition: a 40 nS IPSG
    # leading by 0.2 ms stops it spiking, one 1.5 ms late does not
    led = coincidence.probability_trials(
        "two", [0.0], 1.03 * th, dt_inh=-0.2, **kw
    )
    late = coincidence.probability_trials(
        "two", [0.0], 1.03 * th, dt_inh=1.5, **kw
    )
    assert led.tolist() == [0.0] and late.tolist() == [1.0]


def test_an_ipsg_in_an_event_list_jitters_like_the_excitation(monkeypatch):
    seen = []
    monkeypatch.setattr(
        coincidence,
        "trial_spikes",
        lambda model, trials, *a, **k: (
            seen.extend(trials) or np.zeros(len(trials), bool)
        ),
    )
    events = (synaptic.ipsg(4.8), synaptic.Event(5.0, 1, 1.0))
    coincidence.probability_trials(
        "two", [0.0], 50.0, n_trials=5, jitter=0.1, events=events, **TWO
    )
    assert len({next(e for e in tr if e.e_rev < 0).t for tr in seen}) == 5


def test_inh_jitter_moves_only_the_inhibitory_onsets(monkeypatch):
    seen = []

    def fake(model, trials, *args, **kw):
        seen.extend(trials)
        return np.zeros(len(trials), bool)

    monkeypatch.setattr(coincidence, "trial_spikes", fake)
    kw = dict(n_trials=6, jitter=0.1, dt_inh=-0.2, **TWO)
    coincidence.probability_trials("two", [0.0], 50.0, inh_jitter=0.0, **kw)
    ipsg = {next(e for e in tr if e.e_rev < 0).t for tr in seen}
    epsg = {next(e for e in tr if e.e_rev == 0).t for tr in seen}
    assert ipsg == {5.0 - 0.2} and len(epsg) == 6
    seen.clear()
    coincidence.probability_trials("two", [0.0], 50.0, inh_jitter=0.1, **kw)
    assert len({next(e for e in tr if e.e_rev < 0).t for tr in seen}) == 6
    # left at None, the IPSG jitters like the excitation (as for an event
    # list that holds one before the change)
    seen.clear()
    coincidence.probability_trials("two", [0.0], 50.0, **kw)
    assert len({next(e for e in tr if e.e_rev < 0).t for tr in seen}) == 6


# --- the equal-leak membrane of scripts/inhibition_window.py ---------------


def _script():
    path = (
        pathlib.Path(__file__).parents[1] / "scripts" / "inhibition_window.py"
    )
    spec = importlib.util.spec_from_file_location("inhibition_window", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_klt_as_leak_keeps_the_chord_conductance_at_rest():
    leak = _script().klt_as_leak(-68.0)
    base = multi.LUMPED
    assert not leak.g_klt.any() and leak.key != base.key
    v = np.full(base.n, -68.0)
    chord = base.g_klt * constants.winf(v) ** 4 * constants.zinf(v)
    np.testing.assert_allclose(leak.g_lk, base.g_lk + chord, rtol=1e-12)
    assert leak.key != _script().klt_as_leak(-60.0).key  # v0 is in the key


def test_taking_klt_out_halves_the_pair_threshold():
    script = _script()
    leak = script.klt_as_leak(-68.0)
    kw = {
        "morph": leak,
        "mem": multi.membrane(-68.0, rebalance_rest=True, morph=leak),
    }
    keep = coincidence.threshold("multi", 0.0, rel_tol=1e-2)
    out = coincidence.threshold("multi", 0.0, model_kw=kw, rel_tol=1e-2)
    assert out < 0.7 * keep  # 22 against 48 nS: the Kv1 history is gone
