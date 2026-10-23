"""Direct tests of the solver driver's contracts (msoaxon/_solve.py).

These pin what the driver does today, one rule per test, so that a new
integrator (issue #28) can be held to the same rules. Most are pure
functions of their arguments; the last two run the models to check which
runs share the stimulus-free prefix.
"""

import numpy as np
import pytest

import msoaxon
from msoaxon import _common, _solve, synaptic

# --- breakpoints ----------------------------------------------------------


@pytest.mark.parametrize(
    "stim_type,start,stop,t_end,expected",
    [
        # a set: equal start and stop collapse to one point
        ("step", 5, 5, 20, [5]),
        ("EPSGpair", 5, 5, 20, [5]),
        # only points strictly inside (0, t_end) are kept
        ("step", 5, 50, 20, [5]),
        ("step", 0, 10, 20, [10]),
        ("step", -3, 10, 20, [10]),
        ("step", 5, 20, 20, [5]),
        ("step", 0, 20, 20, []),
        # the sort puts a later start after an earlier stop
        ("step", 7, 5, 20, [5, 7]),
        # ramp and ramp2 add 5 and stop + 5, whatever start is
        ("ramp", 2, 4, 20, [2, 4, 5, 9]),
        ("ramp2", 2, 4, 20, [2, 4, 5, 9]),
        ("ramp", 8, 10, 30, [5, 8, 10, 15]),
        ("ramp", 5, 10, 30, [5, 10, 15]),
        # ... and are clipped like the others: stop + 5 beyond t_end drops
        ("ramp", 2, 4, 8, [2, 4, 5]),
        ("ramp2", 2, 4, 9, [2, 4, 5]),
        ("ramp", 0, 4, 20, [4, 5, 9]),
        # the other stimuli add nothing
        ("sine", 2, 4, 20, [2, 4]),
        ("Synaptic", 2, 4, 20, [2, 4]),
        ("SynapticPair", 2, 4, 20, [2, 4]),
        ("EPSG", 5, 10, 20, [5, 10]),
        ("EPSGpair", 5, 5.3, 20, [5, 5.3]),
        ("EPSGpair", 9, 4, 20, [4, 9]),
        ("EPSGbilateral", 5, 5.3, 20, [5, 5.3]),
    ],
)
def test_breakpoints(stim_type, start, stop, t_end, expected):
    got = _solve.breakpoints(stim_type, start, stop, t_end)
    assert isinstance(got, list)
    assert got == sorted(got)
    assert got == pytest.approx(expected, abs=0, rel=0)


# --- input_end ------------------------------------------------------------

STIM_NAMES = _common.STIM_TYPES + ("EPSGbilateral",)
TAU = (0.1, 0.3)  # 20 * max(TAU) = 6 ms


def _expected_input_end(stim_type, start, stop, t_end):
    if stim_type == "EPSG":
        return start + 6.0
    if stim_type in ("EPSGpair", "EPSGbilateral"):
        return max(start, stop) + 6.0
    if stim_type == "ramp2":
        return t_end
    return max(start, stop)


@pytest.mark.parametrize("stim_type", STIM_NAMES)
@pytest.mark.parametrize("start,stop", [(5, 8), (9, 4), (5, 5)])
def test_input_end_for_every_stimulus(stim_type, start, stop):
    got = _solve.input_end(stim_type, start, stop, 50, TAU)
    assert got == pytest.approx(_expected_input_end(stim_type, start, stop, 50))


def test_input_end_covers_every_stimulus_name():
    # a stimulus added to the models has to get an input_end rule on purpose:
    # the fall-through (t_end, "never ends") is the safe but slow answer
    assert set(STIM_NAMES) == {
        "step",
        "ramp",
        "ramp2",
        "sine",
        "Synaptic",
        "SynapticPair",
        "EPSG",
        "EPSGpair",
        "EPSGbilateral",
    }


@pytest.mark.parametrize("tau", [(0.1, 0.3), (0.5, 0.2)])
def test_input_end_uses_the_slower_time_constant(tau):
    # 20 * max(rise, decay), whichever of the two is larger
    assert _solve.input_end("EPSG", 5, 8, 50, tau) == pytest.approx(
        5 + 20 * max(tau)
    )


@pytest.mark.parametrize("stim_type", ["step", "ramp", "sine", "Synaptic"])
def test_input_end_ignores_the_time_constants_and_t_end(stim_type):
    a = _solve.input_end(stim_type, 5, 8, 50, (0.1, 0.18))
    b = _solve.input_end(stim_type, 5, 8, 7, (1.0, 2.0))
    assert a == b == 8


# --- pre_stimulus_is_quiet --------------------------------------------------


@pytest.mark.parametrize(
    "cuts,start,stop,expected",
    [
        ([], 5, 10, False),  # no breakpoint, no prefix to share
        ([], 0, 0, False),
        ([3.0], 5, 10, True),  # a cut before the stimulus
        ([5.0], 5, 10, True),  # a cut at the stimulus start counts
        ([8.0], 5, 10, False),  # a cut after the start
        ([5.0, 10.0], 8, 5, True),  # min(start, stop) is the earlier onset
        ([6.0, 8.0], 8, 5, False),
        ([2.0], 5, 1, False),  # stop before start: min is stop
        ([1.0], 5, 1, True),
    ],
)
def test_pre_stimulus_is_quiet(cuts, start, stop, expected):
    got = _solve.pre_stimulus_is_quiet(cuts, start, stop)
    assert got is expected


# --- spike_event ------------------------------------------------------------


def test_spike_event_is_the_axon_minus_soma_minus_factor():
    event = _solve.spike_event(2, 10.0)
    assert vars(event) == {"terminal": True, "direction": 1}
    x = np.array([-60.0, 99.0, -45.0])
    assert event(0.0, x) == pytest.approx(-45.0 - -60.0 - 10.0)
    assert event(123.4, x) == event(0.0, x)  # time plays no part
    x[2] = -50.0  # exactly factor above the soma: the crossing
    assert event(0.0, x) == 0.0
    x[2] = -40.0
    assert event(0.0, x) == pytest.approx(10.0)


# --- settled_event ------------------------------------------------------------


def _dv_of(values):
    return lambda t, x: np.array(values, dtype=float)


def test_settled_event_is_terminal_and_rising():
    event = _solve.settled_event(_dv_of([-1.0]), 5.0, 1, 10.0)
    assert vars(event) == {"terminal": True, "direction": 1}


def test_settled_event_waits_for_the_input_to_end():
    event = _solve.settled_event(_dv_of([-1.0, -1.0]), 5.0, 1, 10.0)
    x = np.array([-60.0, -60.0])  # quiet and falling
    assert event(4.999, x) == -1.0
    assert event(5.0, x) == 1.0  # t_off itself counts


def test_settled_event_needs_the_axon_gap_below_half_the_factor():
    event = _solve.settled_event(_dv_of([-1.0, -1.0]), 5.0, 1, 10.0)
    assert event(6.0, np.array([-60.0, -50.0])) == -1.0  # gap 10
    assert event(6.0, np.array([-60.0, -55.0])) == -1.0  # gap 5 = factor / 2
    assert event(6.0, np.array([-60.0, -55.1])) == 1.0  # just under
    assert event(6.0, np.array([-60.0, -70.0])) == 1.0  # axon below the soma


@pytest.mark.parametrize(
    "dv,expected",
    [
        ([-1.0, -2.0, -0.1], 1.0),
        ([-1.0, 0.0, -1.0], -1.0),  # a flat compartment is not falling
        ([-1.0, 1e-9, -1.0], -1.0),  # one rising compartment
        ([1.0, 1.0, 1.0], -1.0),
    ],
)
def test_settled_event_needs_every_voltage_falling(dv, expected):
    event = _solve.settled_event(_dv_of(dv), 5.0, 1, 10.0)
    assert event(6.0, np.array([-60.0, -62.0, -61.0])) == expected


def test_settled_event_does_not_call_dv_before_t_off_or_while_spiking():
    calls = []

    def dv(t, x):
        calls.append(t)
        return np.array([-1.0, -1.0])

    event = _solve.settled_event(dv, 5.0, 1, 10.0)
    event(1.0, np.array([-60.0, -60.0]))
    event(6.0, np.array([-60.0, -40.0]))
    assert calls == []
    event(6.0, np.array([-60.0, -60.0]))
    assert calls == [6.0]


def test_settled_event_stops_a_real_integration_once_it_has_settled():
    # x' = -x from 1: every voltage falls from the start, so t_off = 2 alone
    # sets the stop, and the root is located where the event switches sign
    event = _solve.settled_event(lambda t, x: -x, 2.0, 0, 10.0)
    t, y = _solve.integrate(
        lambda t, x: -x,
        np.array([1.0, 1.0]),
        10.0,
        [],
        rtol=1e-8,
        atol=1e-8,
        max_step=0.1,
        stop_event=event,
    )
    assert t[-1] == pytest.approx(2.0, abs=1e-6)
    assert y.shape == (len(t), 2)


# --- synaptic.interp_g --------------------------------------------------------


def test_interp_g_is_zero_outside_the_table_and_linear_inside():
    t_syn = np.array([0.0, 1.0, 2.0])
    g_syn = np.array([0.0, 10.0, 4.0])
    for t in (-1e6, -1e-9, 2.0 + 1e-9, 3.0, 1e6):  # left = right = 0
        assert synaptic.interp_g(t_syn, g_syn, t) == 0.0
    # exact at the knots, the ends included
    for t, g in zip(t_syn, g_syn, strict=True):
        assert synaptic.interp_g(t_syn, g_syn, float(t)) == g
    assert synaptic.interp_g(t_syn, g_syn, 0.5) == pytest.approx(5.0)
    assert synaptic.interp_g(t_syn, g_syn, 1.5) == pytest.approx(7.0)
    assert synaptic.interp_g(t_syn, g_syn, 0.25) == pytest.approx(2.5)


def test_interp_g_returns_a_python_float():
    got = synaptic.interp_g(np.array([0.0, 1.0]), np.array([1.0, 3.0]), 0.5)
    assert type(got) is float and got == 2.0


def test_interp_g_does_not_extrapolate_a_nonzero_end_value():
    # unlike np.interp's default, which holds the end values; the table's last
    # conductance is not carried past its end
    t_syn = np.array([0.0, 1.0])
    g_syn = np.array([5.0, 5.0])
    assert synaptic.interp_g(t_syn, g_syn, 1.5) == 0.0
    assert synaptic.interp_g(t_syn, g_syn, -0.5) == 0.0


# --- the two/multi quiet-prefix asymmetry -------------------------------------
#
# mso_axon switches stimuli on for t > start, so the equations at t == start
# are still stimulus-free and every stimulus can share the quiet prefix.
# two_cpt switches step and the synaptic inputs on for t >= start, so sharing
# would change their results and they skip the cache.

# (stim_type, start, stop, I, node, model_type, t_end, v0, input_node)
QUIET_ARGS = {
    "step": ("step", 5, 8, 50.0, 3, "passive", 12, -68.0, 1),
    "Synaptic": ("Synaptic", 5, 8, 0.0, 3, "passive", 12, -68.0, 1),
    "EPSG": ("EPSG", 5, 8, 1.0, 3, "passive", 12, -68.0, 1),
}


@pytest.mark.parametrize("stim_type", ["step", "Synaptic"])
def test_two_cpt_does_not_share_the_quiet_prefix_for_step_and_synaptic(
    fresh_quiet_cache, stim_type
):
    msoaxon.two_cpt(*QUIET_ARGS[stim_type])
    assert len(_solve._QUIET_CACHE) == 0


@pytest.mark.parametrize("stim_type", ["step", "Synaptic", "EPSG"])
def test_mso_axon_shares_the_quiet_prefix_when_the_stimulus_starts_late(
    fresh_quiet_cache, stim_type
):
    msoaxon.mso_axon(*QUIET_ARGS[stim_type], mode="exact")
    assert len(_solve._QUIET_CACHE) == 1


def test_two_cpt_shares_the_quiet_prefix_for_the_other_stimuli(
    fresh_quiet_cache,
):
    # the control: the asymmetry is specific to step and the synaptic inputs
    msoaxon.two_cpt(*QUIET_ARGS["EPSG"])
    assert len(_solve._QUIET_CACHE) == 1


@pytest.mark.parametrize("model", ["multi", "two"])
def test_no_quiet_prefix_when_the_stimulus_starts_at_zero(
    fresh_quiet_cache, model
):
    # no breakpoint before the stimulus (start == 0 is not a cut), so there is
    # no prefix to share in either model
    args = ("EPSG", 0, 8, 1.0, 3, "passive", 12, -68.0, 1)
    if model == "multi":
        msoaxon.mso_axon(*args, mode="exact")
    else:
        msoaxon.two_cpt(*args)
    assert len(_solve._QUIET_CACHE) == 0
