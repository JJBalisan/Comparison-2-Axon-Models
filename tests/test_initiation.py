"""The initiation-site readout (issue #80): where each spike starts.

Most tests run on synthetic traces whose answer is known; two run the model
(an EPSG pair, and a short run of sustained Synaptic input).
"""

import math

import numpy as np
import pytest
from scipy import special

from msoaxon import initiation, multi

LABELS = multi.LUMPED.labels
NODES = [i + 1 for i, name in enumerate(LABELS) if name == "node"]
AIS = [i + 1 for i, name in enumerate(LABELS) if name == "AIS"]


def synthetic(onsets, v0=-68.0, rise=100.0, width=0.02):
    """Builds a trace in which compartment c's steepest rise is at onsets[c].

    A compartment rises along a logistic centred on its onset; the soma
    (column 0) gets a gentle EPSP-like ramp 0.5 ms before each spike.
    Compartments with no onset stay at rest.

    Args:
        onsets: Per spike, a dict {1-indexed compartment: time [ms]}.
        v0: The resting potential [mV].
        rise: The height of each logistic [mV].
        width: The logistic's time constant [ms].

    Returns:
        (t, x) with x shaped (n_times, 45).
    """
    t = np.arange(0.0, 20.0, 0.002)
    x = np.full((len(t), multi.LUMPED.n), v0)
    for spike in onsets:
        for c, t0 in spike.items():
            x[:, c - 1] += rise * special.expit((t - t0) / width)
            x[:, c - 1] -= rise * (t > t0 + 0.5)  # and back down again
        epsp = min(spike.values()) - 0.5
        x[:, 0] += 15 * special.expit((t - epsp) / 0.3) - 15 * (t > epsp + 1)
    return t, x


def chain(first, t0):
    """Returns one spike's onsets: `first` leads, the rest follow by 10 us."""
    return {c: t0 + 0.01 * abs(c - first) for c in AIS + NODES}


def pair_run():
    """Runs a coincident EPSG pair, 60 nS (1.25x threshold), into the soma."""
    return multi.mso_axon(
        "EPSGpair", 5.0, 5.0, 60.0, 3, "active-full", 8.0, -68.0, 1
    )


@pytest.fixture(scope="module")
def pair():
    return pair_run()


def test_pair_initiates_in_the_ais(pair):
    (spike,) = initiation.find_spikes(*pair)
    assert spike.compartment in AIS
    assert 5.0 < spike.time < 6.0
    # the whole profile is returned, and it climbs the axon
    assert spike.times.shape == (multi.LUMPED.n,)
    assert spike.times[NODES[0] - 1] > spike.time
    assert spike.times[NODES[1] - 1] > spike.times[NODES[0] - 1]


def test_readout_is_deterministic(pair):
    a, b = initiation.find_spikes(*pair), initiation.find_spikes(*pair)
    assert [s.compartment for s in a] == [s.compartment for s in b]
    assert [s.time for s in a] == [s.time for s in b]
    np.testing.assert_array_equal(a[0].times, b[0].times)


@pytest.mark.parametrize("first", [AIS[1], NODES[0], NODES[1], NODES[5]])
def test_synthetic_first_site_is_recovered(first):
    t, x = synthetic([chain(first, 8.0)])
    (spike,) = initiation.find_spikes(t, x)
    assert spike.compartment == first
    assert spike.time == pytest.approx(8.0, abs=2e-3)
    # the profile matches the construction
    assert spike.times[NODES[2] - 1] == pytest.approx(
        8.0 + 0.01 * abs(NODES[2] - first), abs=2e-3
    )


def test_the_soma_and_internodes_never_win():
    # the soma's EPSP is the earliest event and the internodes' rises are
    # sandwiched by their nodes', but only AIS and nodes compete
    t, x = synthetic([chain(AIS[0], 8.0)])
    x[:, 3] = x[:, 4]  # compartment 4, an internode, rises with node 1
    (spike,) = initiation.find_spikes(t, x)
    assert spike.compartment == AIS[0]


def test_a_compartment_slower_than_the_soma_is_not_firing():
    # the earliest rise is in AIS compartment 2, but gentler than the EPSP
    t, x = synthetic([chain(AIS[1], 8.0)])
    x[:, AIS[0] - 1] = -68.0 + 4 * special.expit((t - 7.9) / 0.3)
    (spike,) = initiation.find_spikes(t, x)
    assert spike.compartment == AIS[1]
    assert spike.times[AIS[0] - 1] == pytest.approx(7.9, abs=0.05)


def test_each_spike_gets_its_own_window():
    t, x = synthetic([chain(AIS[1], 6.0), chain(NODES[1], 10.0)])
    first, second = initiation.find_spikes(t, x)
    assert (first.compartment, second.compartment) == (AIS[1], NODES[1])
    assert second.time == pytest.approx(10.0, abs=2e-3)


def test_no_spike_gives_no_spikes_and_nan_fractions():
    t, x = synthetic([])
    assert initiation.find_spikes(t, x) == []
    fractions = initiation.group_fractions([])
    assert fractions.n_spikes == 0
    assert math.isnan(fractions.ais)


def test_a_spike_without_an_axon_to_fire_is_skipped():
    # the reference is the soma, the only thing that rises
    t, x = synthetic([{1: 8.0}])
    assert initiation.find_spikes(t, x, reference=1) == []


def test_a_trace_too_sparse_to_slope_is_skipped():
    n = multi.LUMPED.n
    last = n - 1  # the most distal node, the default reference
    for t, ref in (
        (np.array([0.0, 1.0]), [-68.0, 100.0]),  # two samples
        (np.arange(0.0, 2.4, 0.4), [-68, -68, -68, 100, 100, 100]),  # sparse
    ):
        x = np.full((len(t), n), -68.0)
        x[:, last] = ref
        assert initiation.find_spikes(t, x) == []


def test_a_repeated_time_is_tolerated():
    # the solver repeats a time at a restart
    t, x = synthetic([chain(NODES[0], 8.0)])
    t = np.insert(t, 100, t[100])
    x = np.insert(x, 100, x[100], axis=0)
    assert initiation.find_spikes(t, x)[0].compartment == NODES[0]


def test_a_spike_at_the_end_of_a_window_still_reads():
    # the steepest point is the last interval: no neighbour to fit a parabola
    t = np.array([0.0, 1.0, 2.0])
    assert initiation._steepest(t, np.array([0.0, 1.0, 5.0])) == (1.5, 4.0)


def test_the_last_peak_skips_what_is_still_climbing():
    t = np.arange(0.0, 2.0, 0.01)
    v = 100 * special.expit((t - 1.0) / 0.02) + 300 * t  # a ramp and a rise
    assert initiation._steepest(t, v, last=True)[0] == pytest.approx(
        1.0, abs=2e-3
    )
    ramp = np.cumsum(t)  # steepening to the end, no interior peak
    time, slope = initiation._steepest(t, ramp, last=True)
    assert math.isnan(time)
    assert slope == -math.inf


def test_an_earlier_epsp_peak_does_not_hide_the_upstroke():
    # AIS compartment 2 has an EPSP-like rise, steeper than its spike, that
    # comes 0.3 ms before node 1 fires; node 1 still leads the spike
    t, x = synthetic([chain(NODES[0], 8.0)])
    x[:, AIS[0] - 1] += 60 * special.expit((t - 7.8) / 0.01)
    x[:, AIS[0] - 1] -= 60 * (t > 8.3)
    (spike,) = initiation.find_spikes(t, x)
    assert spike.compartment == NODES[0]
    assert spike.times[AIS[0] - 1] == pytest.approx(8.03, abs=2e-3)


def test_a_rise_at_the_start_is_reported_where_it_is():
    t = np.arange(5.0)
    steepest = initiation._steepest(t, np.array([0.0, 1, 2, 2.5, 2.5]))
    assert steepest == (0.5, 1.0)


def test_unknown_reference_or_labels_are_rejected():
    t, x = synthetic([])
    with pytest.raises(ValueError, match="reference"):
        initiation.find_spikes(t, x, reference=99)
    with pytest.raises(ValueError, match="labelled"):
        initiation.find_spikes(t, x, labels=("dendrite",))


def test_fractions_are_by_label():
    spikes = [
        initiation.Spike(c, 0.0, np.zeros(multi.LUMPED.n))
        for c in (2, 3, 3, 5, 7, 9, 4, 5)
    ]
    f = initiation.group_fractions(spikes)
    assert f.n_spikes == 8
    assert (f.ais, f.node1, f.node2, f.other) == (3 / 8, 2 / 8, 1 / 8, 2 / 8)
    assert f.ais + f.node1 + f.node2 + f.other == pytest.approx(1.0)
    assert f.distal == f.node1 + f.node2


def test_sustained_input_fractions_sum_to_one():
    f = initiation.sustained_input(freq=350.0, g_e=1.04e8, t_end=20.0)
    assert f.n_spikes >= 1
    assert f.ais + f.node1 + f.node2 + f.other == pytest.approx(1.0)
    assert all(0.0 <= v <= 1.0 for v in (f.ais, f.node1, f.node2, f.other))


def test_sustained_input_without_spikes():
    f = initiation.sustained_input(g_e=1.0e5, t_end=8.0)
    assert f.n_spikes == 0
    assert math.isnan(f.other)
