"""threshold._search, driven by a fake spiked_at.

No ODE solves: the fake spikes at I >= threshold and records what it was
asked, so each test pins the probe sequence, the returned value and the run
count of BinarySearch.m's halving search. Every value below is a dyadic
fraction of MAX_I, so the comparisons are exact.
"""

import itertools

import pytest

from msoaxon import threshold

MAX_I = 100.0
FACTOR = 10.0
EPSG_PAIR_DT = 1 / 25
MAX_CALLS = 50


class Fake:
    """A model stand-in: spikes at I >= threshold, and logs every run."""

    def __init__(self, spike_at):
        self.spike_at = spike_at
        self.calls = []

    def __call__(self, stim_type, start, stop, current, syn, factor):
        self.calls.append((stim_type, start, stop, current, syn, factor))
        if len(self.calls) > MAX_CALLS:
            raise RuntimeError(f"search did not end ({MAX_CALLS} runs)")
        return current >= self.spike_at

    @property
    def probes(self):
        return [call[3] for call in self.calls]


def search(fake, stim_type="step", i=1, ceiling=MAX_I, zoom=1.0):
    return threshold._search(
        fake, stim_type, i, FACTOR, ceiling, zoom, EPSG_PAIR_DT
    )


@pytest.mark.parametrize(
    "spike_at, probes",
    [
        # last probe (29.6875) is below the threshold: the search returns it
        # anyway, it does not report the lowest spiking value
        (30.0, [100, 50, 25, 37.5, 31.25, 28.125, 29.6875]),
        # last probe (60.9375) is above the threshold
        (60.0, [100, 50, 75, 62.5, 56.25, 59.375, 60.9375]),
        # fires at max_I and nowhere below: every probe after the first misses
        (100.0, [100, 50, 75, 87.5, 93.75, 96.875, 98.4375]),
        # fires everywhere: the probe halves each time
        (0.0, [100, 50, 25, 12.5, 6.25, 3.125, 1.5625]),
    ],
)
def test_probe_sequence_and_returned_value(spike_at, probes):
    fake = Fake(spike_at)
    assert search(fake) == probes[-1]  # the last tested value
    assert fake.probes == probes
    assert len(fake.calls) == len(probes)


def test_step_is_half_the_distance_to_the_previous_probe():
    # after probing p the search moves by |p - previous| / 2, where previous
    # starts at 0: 100 -> 50 (|100 - 0| / 2), 50 -> 75 (|50 - 100| / 2), ...
    fake = Fake(60.0)
    search(fake)
    probes = [0.0] + fake.probes
    steps = [abs(b - a) / 2 for a, b in itertools.pairwise(probes)]
    assert steps == [50, 25, 12.5, 6.25, 3.125, 1.5625, 0.78125]
    moves = [abs(b - a) for a, b in itertools.pairwise(fake.probes)]
    assert moves == steps[:-1]  # the last step ends the search instead


def test_never_firing_returns_max_i_after_one_run():
    fake = Fake(float("inf"))
    assert search(fake) == MAX_I
    assert fake.probes == [MAX_I]  # the miss steps to 150 > max_I and stops


def test_threshold_above_max_i_is_the_never_firing_case():
    fake = Fake(MAX_I + 0.5)
    assert search(fake) == MAX_I
    assert len(fake.calls) == 1


@pytest.mark.parametrize(
    "zoom, runs",
    [
        (50.0, 1),  # the first step is already 50, which is not above zoom
        (30.0, 2),
        (1.0, 7),
        (0.1, 10),
    ],
)
def test_search_stops_once_the_step_is_at_most_zoom(zoom, runs):
    fake = Fake(0.0)  # fires everywhere, so the steps are 50, 25, 12.5, ...
    search(fake, zoom=zoom)
    assert len(fake.calls) == runs


def test_search_scales_with_max_i():
    fake = Fake(0.0)
    assert search(fake, ceiling=8.0, zoom=1.0) == 2.0
    assert fake.probes == [8, 4, 2]  # the step reaches 1 = zoom after 2


def test_spiked_at_gets_the_sweep_point_settings():
    fake = Fake(30.0)
    search(fake, "ramp", i=3)
    stop, syn, _ = threshold.sweep_setting("ramp", 3)
    for stim_type, start, run_stop, _, run_syn, factor in fake.calls:
        assert (stim_type, start, run_stop) == ("ramp", threshold.START, stop)
        assert factor == FACTOR
        assert run_syn.t_end == syn.t_end


def test_fixed_current_stimulus_ignores_the_search_location():
    # Synaptic fixes I at 0, but the search still walks its locations, so the
    # run count and returned value follow the spikes, not the I passed in
    fake = Fake(0.0)
    assert search(fake, "Synaptic") == 1.5625
    assert fake.probes == [0.0] * 7
