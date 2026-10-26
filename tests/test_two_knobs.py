"""The free two-compartment knobs: coupling, fractions, gNa2, Na V-half.

Issues #84 and #79. Every existing call must keep its bits (the node lookup
supplies whatever a knob leaves unset), a bad value must raise before any
run, and the V-half must hold rest at v0 without a rebalance.
"""

import dataclasses

import numpy as np
import pytest

from msoaxon import _solve, coincidence, constants as C, two

RUN = ("EPSG", 5, 5, 1.0, 3, "active-full", 20, -68.0, 1)


def run(*args, **kw):
    return two.two_cpt(*(args or RUN), **kw)


def _bits(a):
    return np.ascontiguousarray(a).tobytes()


def same(a, b):
    return _bits(a[0]) == _bits(b[0]) and _bits(a[1]) == _bits(b[1])


@pytest.fixture
def no_run(monkeypatch):
    """Fails the test if a model run is reached."""

    def reached(*args, **kw):
        raise AssertionError("a run happened before the check")

    monkeypatch.setattr(_solve, "integrate", reached)


# --- every existing result holds ---


@pytest.mark.parametrize("node", [2, 3, 5, 9])
@pytest.mark.parametrize("model_type", ["active-full", "Active-sodium"])
def test_the_lookups_own_values_given_explicitly_are_the_same_bits(
    node, model_type
):
    P = two.get_params(-68.0, node, 1, model_type)
    given = {
        "gC": P.gC,
        "couple12": P.couple12,
        "couple21": P.couple21,
        "area_ratio": C.AREA_RATIO[node - 1],
        "gNa2": P.gNa2,
        "soma_na_vhalf": -62.5,
        "axon_na_vhalf": -62.5,
    }
    args = ("EPSG", 5, 5, 1.0, node, model_type, 20, -68.0, 1)
    assert same(run(*args), run(*args, **given))


def test_the_default_parameters_have_the_default_vhalf():
    P = two.get_params(-68.0, 3, 1, "active-full")
    assert (P.na_vhalf1, P.na_vhalf2) == (62.5, 62.5)


@pytest.mark.parametrize("v0", [-68.0, -58.0])
def test_a_shifted_vhalf_leaves_the_calibration_alone(v0):
    a = two.get_params(v0, 3, 1, "active-full")
    b = two.get_params(v0, 3, 1, "active-full", soma_na_vhalf=-77.0)
    for f in dataclasses.fields(two.TwoCptParams):
        if f.name not in ("INa0_1", "na_vhalf1"):
            assert getattr(a, f.name) == getattr(b, f.name), f.name
    assert b.na_vhalf1 == 77.0
    assert b.INa0_1 != a.INa0_1


# --- validation, before any run ---

BAD = [0.0, -1.0, float("nan")]


@pytest.mark.parametrize("bad", BAD)
@pytest.mark.usefixtures("no_run")
def test_gc_must_be_positive(bad):
    with pytest.raises(ValueError, match=r"gC must be > 0"):
        run(gC=bad)


@pytest.mark.parametrize("bad", BAD + [1.0, 1.5])
@pytest.mark.parametrize("name", ["couple12", "couple21"])
@pytest.mark.usefixtures("no_run")
def test_coupling_constants_must_be_in_0_1(name, bad):
    with pytest.raises(ValueError, match=rf"{name} must be in \(0, 1\)"):
        run(**{name: bad})


@pytest.mark.parametrize("bad", BAD)
@pytest.mark.usefixtures("no_run")
def test_area_ratio_must_be_positive(bad):
    with pytest.raises(ValueError, match=r"area_ratio must be > 0"):
        run(area_ratio=bad)


@pytest.mark.parametrize("bad", [-0.1, 1.1, float("nan")])
@pytest.mark.parametrize("name", ["klt1", "klt2", "h1", "h2", "na1"])
@pytest.mark.usefixtures("no_run")
def test_fractions_must_be_in_0_1(name, bad):
    with pytest.raises(ValueError, match=rf"{name} must be in \[0, 1\]"):
        run(**{name: bad})


@pytest.mark.usefixtures("no_run")
def test_a_fraction_is_checked_even_where_the_model_type_switches_it_off():
    with pytest.raises(ValueError, match=r"klt1 must be in"):
        run("step", 5, 6, 100.0, 3, "passive", 20, -68.0, 1, klt1=2.0)


@pytest.mark.parametrize(
    "kw, who",
    [
        ({"klt1": 0.7, "h1": 0.2, "na1": 0.2}, "soma"),
        ({"klt2": 0.6, "h2": 0.5}, "axon"),
    ],
)
@pytest.mark.usefixtures("no_run")
def test_fractions_may_not_leave_a_compartment_without_leak(kw, who):
    with pytest.raises(ValueError, match=rf"the {who}'s KLT, h and Na"):
        run(**kw)


def test_fractions_of_channels_that_are_off_do_not_count():
    # h1 = 1 would leave no leak, but 'active-KLT' has no h
    P = two.get_params(-68.0, 3, 1, "active-KLT", h1=1.0, h2=1.0)
    assert (P.gh1, P.gh2) == (0.0, 0.0)


@pytest.mark.parametrize("bad", [-1.0, float("nan")])
@pytest.mark.usefixtures("no_run")
def test_gna2_must_not_be_negative(bad):
    with pytest.raises(ValueError, match=r"gNa2 must be >= 0"):
        run(gNa2=bad)


def test_gna2_may_be_zero():
    P = two.get_params(-68.0, 3, 1, "active-full", gNa2=0.0)
    assert P.gNa2 == 0.0 and P.INa0_2 == 0.0


@pytest.mark.parametrize("name", ["soma_na_vhalf", "axon_na_vhalf"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
@pytest.mark.usefixtures("no_run")
def test_vhalf_must_be_finite(name, bad):
    with pytest.raises(ValueError, match=rf"{name} must be finite"):
        run(**{name: bad})


@pytest.mark.usefixtures("no_run")
def test_gc_and_r1_are_one_or_the_other():
    with pytest.raises(ValueError, match=r"gC or r1, not both"):
        run(gC=50.0, r1=8.5)


@pytest.mark.usefixtures("no_run")
def test_a_search_refuses_a_bad_knob_before_it_runs():
    with pytest.raises(ValueError, match=r"gC must be > 0"):
        coincidence.threshold("two", 0.1, model_kw={"gC": -1.0})


# --- what the knobs do ---


def test_gc_replaces_the_one_r1_gives_and_leaves_the_node_coupling():
    P = two.get_params(-68.0, 3, 1, "active-full")
    Q = two.get_params(-68.0, 3, 1, "active-full", gC=56.4)  # Jing et al 2025
    assert Q.gC == 56.4
    assert (Q.couple12, Q.couple21) == (P.couple12, P.couple21)
    # the constants fix the ratios, gC the scale
    assert Q.gTot1 / Q.gC == pytest.approx(P.gTot1 / P.gC)
    assert Q.gTot2 / Q.gC == pytest.approx(P.gTot2 / P.gC)
    assert Q.gTot1 == pytest.approx(P.gTot1 * 56.4 / P.gC)


def test_the_docstring_example_runs():
    t, x = two.two_cpt("EPSG", 5, 5, 1.0, 3, "active-full", 20, -68, 1, gC=56.4)
    assert np.isfinite(x).all() and x[:, 0].max() > -68


def test_couplings_and_a_fraction_move_the_response():
    base = run()
    for kw in ({"couple12": 0.5}, {"couple21": 0.1}, {"klt1": 0.3}):
        t, x = run(**kw)
        assert not same(base, (t, x)), kw


def test_free_coupling_gives_the_goldwyn_plane_a_valid_cell():
    for c12, c21 in [(0.05, 0.05), (0.95, 0.05), (0.05, 0.95), (0.9, 0.9)]:
        P = two.get_params(
            -68.0, 3, 1, "active-full", couple12=c12, couple21=c21
        )
        assert P.gC > 0 and P.gTot1 > 0 and P.gTot2 > 0
        assert P.cap1 > 0 and P.cap2 > 0


def test_the_quiet_prefix_key_follows_every_knob(monkeypatch):
    keys = []
    real = _solve.integrate

    def spy(*args, quiet=None, **kw):
        keys.append(quiet[0])
        return real(*args, quiet=quiet, **kw)

    monkeypatch.setattr(_solve, "integrate", spy)
    run()
    for kw in (
        {"gC": 56.4},
        {"couple12": 0.5},
        {"klt1": 0.3},
        {"gNa2": 50.0},
        {"soma_na_vhalf": -77.0},
        {"axon_na_vhalf": -77.0},
    ):
        run(**kw)
    assert len(keys) == 7 and len(set(keys)) == 7


@pytest.mark.parametrize("v0", [-68.0, -58.0])
@pytest.mark.parametrize("knob", ["soma_na_vhalf", "axon_na_vhalf"])
def test_rest_holds_at_v0_with_a_shifted_vhalf(v0, knob):
    t, x = two.two_cpt(
        "step", 5, 6, 0.0, 3, "active-full", 40, v0, 1, **{knob: -77.0}
    )
    # the solver's own residual (6e-7 mV at -68, 5e-6 at -58), which the
    # unshifted model has too
    assert np.abs(x[:, :2] - v0).max() < 1e-4


def test_the_gates_start_at_their_shifted_steady_state():
    P = two.get_params(-58.0, 3, 1, "active-full", soma_na_vhalf=-77.0)
    x0 = two.initial_state(P)
    assert x0[3] == C.hinf(-58.0, 77.0)  # h1
    assert x0[7] == C.hinf(-58.0, 62.5)  # h2: the axon keeps the default
    assert x0[3] < x0[7]


def test_the_somatic_vhalf_shifts_the_soma_and_the_axon_only_the_axon():
    base = run()
    soma = run(soma_na_vhalf=-77.0)
    axon = run(axon_na_vhalf=-77.0)
    assert not same(base, soma) and not same(base, axon)
    # less available somatic Na: a smaller soma depolarization
    assert soma[1][:, 0].max() < base[1][:, 0].max()


def test_the_inactivated_fractions_of_the_issue_table():
    # the issue's table: 1 - hinf at rest, the models' own slope
    got = {
        (rest, vh): 1 - C.hinf(rest, -vh)
        for rest in (-68.0, -58.0)
        for vh in (-62.5, -77.0)
    }
    want = {
        (-68.0, -62.5): 0.33,
        (-68.0, -77.0): 0.76,
        (-58.0, -62.5): 0.64,
        (-58.0, -77.0): 0.92,
    }
    for k, v in want.items():
        assert got[k] == pytest.approx(v, abs=0.005)
