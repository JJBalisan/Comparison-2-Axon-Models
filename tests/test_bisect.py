"""The searches in msoaxon/_bisect.py, driven by cheap fake predicates.

No ODE solves: each fake stands in for a model run, records what it was asked,
and refuses to run more than MAX_CALLS times, so a search that would never end
fails instead of hanging.
"""

import math

import numpy as np
import pytest

from msoaxon import _bisect, coincidence, somatic

LEVEL = 10.0  # the searches' spike criterion, axon - soma [mV]
MAX_CALLS = 200


class Fake:
    """A model stand-in: spikes at x >= threshold, peak rising smoothly below.

    The peak below threshold is LEVEL * (x / threshold) ** power, like the
    models' axon - soma rising convexly through 10 mV; `peak` overrides it.
    """

    def __init__(self, threshold, power=12, peak=None):
        self.threshold, self.power, self.peak = threshold, power, peak
        self.calls = []

    def respond(self, x):
        self.calls.append(x)
        if len(self.calls) > MAX_CALLS:
            raise RuntimeError(f"search did not end ({MAX_CALLS} runs)")
        spiked = x >= self.threshold
        if self.peak is not None:
            return spiked, LEVEL if spiked else self.peak(x)
        if spiked:
            return True, LEVEL
        return False, LEVEL * (x / self.threshold) ** self.power

    def fires(self, x):
        return self.respond(x)[0]


# --- smallest_firing: doubling, then bisection --------------------------------


@pytest.mark.parametrize("threshold", [3.0, 37.0, 50.0, 1500.0])
def test_smallest_firing_brackets_the_threshold(threshold):
    fake = Fake(threshold)
    hi = _bisect.smallest_firing(fake.fires, 50.0, 2000.0, 1e-4)
    assert threshold <= hi <= threshold * (1 + 1e-4)
    # tolerance termination: the last bracket [lo, hi] is within rel_tol
    lo = max((x for x in fake.calls if x < threshold), default=0.0)
    assert hi - lo <= 1e-4 * hi


def test_smallest_firing_gives_up_when_nothing_fires():
    fake = Fake(math.inf)
    assert _bisect.smallest_firing(fake.fires, 1.0, 100.0, 1e-3) == np.inf
    assert max(fake.calls) <= 100.0  # never tries past the ceiling


def test_smallest_firing_does_not_test_zero_when_the_guess_is_quiet():
    # the guess not firing already shows the cell is quiet at rest (firing
    # only gets easier as x grows), so the guard costs no run
    fake = Fake(400.0)
    _bisect.smallest_firing(fake.fires, 50.0, 2000.0, 1e-4)
    assert 0.0 not in fake.calls


def test_smallest_firing_tests_zero_once_when_the_guess_fires():
    fake = Fake(3.0)
    _bisect.smallest_firing(fake.fires, 50.0, 2000.0, 1e-4)
    assert fake.calls.count(0.0) == 1


def test_smallest_firing_fails_when_the_cell_fires_at_rest():
    fake = Fake(0.0)  # fires with no input
    with pytest.raises(_bisect.FiresWithoutInputError, match="no input"):
        _bisect.smallest_firing(fake.fires, 50.0, 2000.0, 1e-4)
    assert len(fake.calls) == 2  # the guess, then zero


def test_fires_without_input_error_is_a_value_error():
    assert issubclass(_bisect.FiresWithoutInputError, ValueError)


# --- smallest_crossing: the predictive (B4) search ----------------------------


@pytest.mark.parametrize("threshold", [3.0, 37.0, 49.9, 50.0, 1500.0])
def test_smallest_crossing_brackets_the_threshold(threshold):
    fake = Fake(threshold)
    hi = _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
    assert threshold <= hi <= threshold * (1 + 1e-4)


def test_smallest_crossing_gives_up_past_the_ceiling():
    fake = Fake(5000.0)
    assert (
        _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
        == np.inf
    )
    assert max(fake.calls) <= 2000.0


def test_smallest_crossing_does_not_test_zero_when_the_guess_is_quiet():
    fake = Fake(400.0)
    _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
    assert 0.0 not in fake.calls


def test_smallest_crossing_tests_zero_once_when_the_guess_fires():
    fake = Fake(3.0)
    _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
    assert fake.calls.count(0.0) == 1


def test_smallest_crossing_fails_when_the_cell_fires_at_rest():
    # today's descent would halve towards zero for ever
    fake = Fake(0.0)
    with pytest.raises(_bisect.FiresWithoutInputError, match="no input"):
        _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
    assert len(fake.calls) == 2


def test_smallest_crossing_falls_back_when_peaks_say_nothing():
    # a flat sub-threshold peak gives no prediction; the search must still
    # step up and bisect to the threshold
    fake = Fake(300.0, peak=lambda x: 2.0)
    hi = _bisect.smallest_crossing(fake.respond, 50.0, 2000.0, 1e-4, LEVEL)
    assert 300.0 <= hi <= 300.0 * (1 + 1e-4)


# --- rheobase: smallest_crossing on step responses (issue #19) ---------------

RHEOBASE = 3784.1  # the 45-compartment model's [pA]


def _step_peak(threshold):
    """Like the multi and two step responses: smooth up to 10 mV."""
    return lambda x: LEVEL * np.exp(-6 * (1 - x / threshold))


def _all_or_none(threshold):
    """Like KLT-blocked steps: ~0.3 mV just below, then a spike."""
    return lambda x: 0.3 * (x / threshold) ** 4


def test_smallest_crossing_descent_sets_the_first_step_down():
    fake = Fake(300.0, peak=_step_peak(300.0))
    _bisect.smallest_crossing(
        fake.respond, 2000.0, 1e6, 1e-4, LEVEL, descent=2.0
    )
    assert fake.calls[:3] == [2000.0, 0.0, 1000.0]
    default = Fake(300.0, peak=_step_peak(300.0))
    _bisect.smallest_crossing(default.respond, 2000.0, 1e6, 1e-4, LEVEL)
    assert default.calls[:3] == [2000.0, 0.0, 1600.0]  # 1.25, as before


def test_smallest_crossing_on_an_all_or_none_response_is_bisection():
    # the peaks below threshold are tiny and flat, so every prediction falls
    # outside the bracket: with descent 2 the search makes bisection's runs
    fake = Fake(199.5, peak=_all_or_none(199.5))
    hi = _bisect.smallest_crossing(
        fake.respond, 2000.0, 1e6, 1e-4, LEVEL, descent=2.0
    )
    bisection = Fake(199.5)
    assert _bisect.smallest_firing(bisection.fires, 2000.0, 1e6, 1e-4) == hi
    assert fake.calls == bisection.calls


def test_smallest_crossing_finds_a_smooth_rheobase_in_fewer_runs():
    fake = Fake(RHEOBASE, peak=_step_peak(RHEOBASE))
    hi = _bisect.smallest_crossing(
        fake.respond, 2000.0, 1e6, 1e-4, LEVEL, descent=2.0
    )
    assert RHEOBASE <= hi <= RHEOBASE * (1 + 1e-4)
    bisection = Fake(RHEOBASE)
    _bisect.smallest_firing(bisection.fires, 2000.0, 1e6, 1e-4)
    assert len(fake.calls) < len(bisection.calls)


def test_rheobase_predicts_from_the_step_responses(monkeypatch):
    calls = []
    peak = _step_peak(RHEOBASE)

    def response(model, stim, start, stop, amp, *args, factor, **kw):
        assert stim == "step" and factor == LEVEL
        calls.append(amp)
        assert len(calls) <= MAX_CALLS
        return amp >= RHEOBASE, LEVEL if amp >= RHEOBASE else peak(amp)

    monkeypatch.setattr(somatic._dispatch, "response", response)
    rb = somatic.rheobase("two")
    assert RHEOBASE <= rb <= RHEOBASE * (1 + 1e-4)
    bisection = Fake(RHEOBASE)
    _bisect.smallest_firing(bisection.fires, 2000.0, 1e6, 1e-4)
    assert len(calls) < len(bisection.calls)


# --- crossing: narrowing a known bracket --------------------------------------


def _midpoints(lo, hi, tol):
    """The probes plain bisection makes on [lo, hi] for a boundary above."""
    out = []
    while hi - lo > tol:
        mid = (lo + hi) / 2
        out.append(mid)
        lo, hi = (lo, mid) if mid >= 0.7 else (mid, hi)
    return out


def test_crossing_bisects_when_there_is_no_prediction():
    fake = Fake(0.7, peak=lambda x: 2.0)  # equal peaks: nothing to fit
    no, yes = _bisect.crossing(fake.respond, 0.0, 1.0, lambda n, y: 1e-3, LEVEL)
    assert fake.calls == _midpoints(0.0, 1.0, 1e-3)
    assert no < 0.7 <= yes and yes - no <= 1e-3


def test_crossing_bisects_when_the_prediction_leaves_the_bracket():
    # peaks that fall as x rises point the fit away from the boundary, below
    # the bracket, so each such step must halve instead
    fake = Fake(0.7, peak=lambda x: LEVEL * (1.0 - x))
    no, yes = _bisect.crossing(
        fake.respond, 0.0, 1.0, lambda n, y: 1e-3, LEVEL, [(0.0, LEVEL)]
    )
    assert no < 0.7 <= yes and yes - no <= 1e-3
    assert fake.calls[0] == 0.5
    assert len(fake.calls) <= len(_midpoints(0.0, 1.0, 1e-3)) + 2


def test_crossing_predicts_a_smooth_boundary_in_few_runs():
    fake = Fake(0.7)
    no, yes = _bisect.crossing(
        fake.respond,
        0.5,
        1.0,
        lambda n, y: 1e-6,
        LEVEL,
        [(0.5, Fake(0.7).respond(0.5)[1])],
    )
    assert no < 0.7 <= yes and yes - no <= 1e-6
    assert len(fake.calls) < len(_midpoints(0.5, 1.0, 1e-6))


def test_crossing_stops_at_once_inside_the_tolerance():
    fake = Fake(0.7)
    assert _bisect.crossing(
        fake.respond, 0.69, 0.71, lambda n, y: 0.05, LEVEL
    ) == (0.69, 0.71)
    assert fake.calls == []


# --- the callers name the configuration ---------------------------------------


def test_threshold_reports_a_cell_that_fires_at_rest(monkeypatch):
    calls = []

    def respond(model, amp, *args, **kw):
        calls.append(amp)
        assert len(calls) <= MAX_CALLS
        return True, LEVEL

    monkeypatch.setattr(coincidence, "_respond", respond)
    with pytest.raises(_bisect.FiresWithoutInputError) as err:
        coincidence.threshold("multi", 0.1, node=5)
    msg = str(err.value)
    assert "'multi'" in msg and "node=5" in msg and "delay=0.1" in msg
    assert "rebalance_rest=True" in msg


def test_rheobase_reports_a_cell_that_fires_at_rest(monkeypatch):
    calls = []

    def response(model, amp, *args, **kw):
        calls.append(amp)
        assert len(calls) <= MAX_CALLS
        return True, LEVEL

    monkeypatch.setattr(somatic, "response", response)
    with pytest.raises(_bisect.FiresWithoutInputError) as err:
        somatic.rheobase("multi", v0=-60.0)
    msg = str(err.value)
    assert "'multi'" in msg and "v0=-60.0" in msg
    assert "rebalance_rest=True" in msg
