"""The settled stop's voltage-only derivatives (multi._dv, two._dv).

_solve.settled_event only needs the sign of every dV/dt, so the models give
it a function that computes the voltage derivatives alone. It must return
exactly the leading entries of the full right-hand side, bit for bit (same
code, same float order), so the event and the stop time do not move.
"""

import numpy as np
import pytest

from msoaxon import (
    _common,
    _dispatch,
    _fast,
    _solve,
    constants,
    multi,
    synaptic,
    two,
)

D = multi.with_dendrites()
D_RA = multi.with_dendrites(conserve_totals=False, dendrite_ra=200.0)
LAT, MED = int(D.lateral[2]), int(D.medial[2])
MYOGA = (0.1, 0.3)

# (stim, start, stop, I, input_node, input_node2, SynParams keywords)
STIMULI = (
    ("none", 5.0, 10.0, 0.0, 1, None, {}),
    ("step", 5.0, 10.0, 8500.0, 3, None, {}),
    ("ramp", 5.0, 8.0, 3000.0, 1, None, {}),
    ("ramp2", 5.0, 5.5, 8500.0, 1, None, {}),
    ("sine", 5.0, 7.5, 8500.0, 1, None, {}),
    ("Synaptic", 5.0, 15.0, 0.0, 1, None, {}),
    ("SynapticPair", 5.0, 15.0, 0.0, 3, None, {"diff": 1.0}),
    ("EPSG", 5.0, 10.0, 108.0, 1, None, {}),
    ("EPSG", 5.0, 10.0, 60.0, 3, None, {"epsg_tau": MYOGA}),
    ("EPSGpair", 5.0, 5.2, 100.0, 1, None, {}),
    ("EPSGpair", 5.0, 5.3, 49.0, 3, None, {"epsg_tau": MYOGA}),
    ("EPSGbilateral", 5.0, 5.2, 60.0, 1, 3, {}),
    (
        "Events",
        5.0,
        5.2,
        1.0,
        1,
        None,
        {
            "events": (
                synaptic.Event(5.0, 1, 60.0),
                synaptic.Event(5.2, 3, 40.0, MYOGA),
                synaptic.Event(4.8, 1, 30.0, (0.28, 1.5), -80.0),
            )
        },
    ),
)


def _multi_configs():
    """Yields (label, stimulus tuple, mem, morph)."""
    for stim in STIMULI:
        yield f"lumped/{stim[0]}", stim, multi.membrane(-68.0), multi.LUMPED
    epsg_pair = STIMULI[9]
    for label, mem in (
        ("rebalance", multi.membrane(-68.0, rebalance_rest=True)),
        (
            "dendrotoxin",
            multi.membrane(
                -68.0, soma_klt_scale=0, ais_klt_scale=0, rebalance_rest=True
            ),
        ),
        (
            "soma-na-vhalf",
            multi.membrane(
                -60.0, soma_klt_scale=3, soma_na_vhalf=-77, rebalance_rest=True
            ),
        ),
    ):
        yield f"mem/{label}", epsg_pair, mem, multi.LUMPED
    bilateral = ("EPSGbilateral", 5.0, 5.1, 60.0, LAT, MED, {})
    for label, morph in (("dendrites", D), ("dendrites-ra200", D_RA)):
        for stim in (STIMULI[1], STIMULI[5], STIMULI[7], STIMULI[9]):
            name, start, stop, amp, _, k2, kw = stim
            yield (
                f"{label}/{name}",
                (name, start, stop, amp, LAT, k2, kw),
                multi.membrane(-68.0, morph=morph),
                morph,
            )
        yield f"{label}/EPSGbilateral", bilateral, None, morph
    yield (
        "dendrites/mem",
        bilateral,
        multi.membrane(
            -68.0, morph=D, dendrite_klt_scale=2.0, rebalance_rest=True
        ),
        D,
    )


MULTI = {c[0]: c[1:] for c in _multi_configs()}


def _times(stim, start, stop, epsg_tau, rng):
    """The window edges, both sides, times after the input, and random."""
    t_off = _solve.input_end(stim, start, stop, 20.0, epsg_tau)
    edges = [start, stop, 5.0, stop + 5.0, t_off]
    times = [e + d for e in edges for d in (0.0, -1e-9, 1e-9)]
    times += list(rng.uniform(t_off, 20.0, 6))  # when the event looks
    times += list(rng.uniform(0.0, 20.0, 6))
    return [float(t) for t in times]


def _same_bits(got, want):
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    assert got.shape == want.shape
    # tobytes, not array_equal: -0.0 and 0.0 must not compare equal
    assert got.tobytes() == want.tobytes(), np.max(np.abs(got - want))


@pytest.mark.parametrize("name", sorted(MULTI))
def test_multi_dv_is_the_voltage_block_of_rhs_bit_for_bit(name):
    (stim, start, stop, amp, k, k2, syn_kw), mem, morph = MULTI[name]
    syn = synaptic.SynParams(t_end=20.0, **syn_kw)
    s = _common.stimulus(stim, start, stop, amp, 20.0, syn)
    mem = mem or multi.membrane(-68.0, morph=morph)
    n = morph.n
    rng = np.random.default_rng(20261013)
    for model_type in ("active-full", "passive", "active-KLT+H"):
        active = tuple(bool(g) for g in multi.ACTIVE_GATES[model_type])
        for t in _times(stim, start, stop, syn.epsg_tau, rng):
            xs = [np.concatenate([np.full(n, -68.0), np.full(6 * n, 0.3)])]
            xs += [
                np.concatenate([rng.uniform(-90, 40, n), rng.random(6 * n)])
                for _ in range(3)
            ]
            for x in xs:
                want = multi._rhs(t, x, stim, s, k, active, mem, morph, k2)
                got = multi._dv(t, x, stim, s, k, mem, morph, k2)
                _same_bits(got, want[:n])


TWO_STIMULI = [s for s in STIMULI if s[0] != "EPSGbilateral"]


@pytest.mark.parametrize("node", [3, 5, 7])
@pytest.mark.parametrize("input_node", [1, 3])
@pytest.mark.parametrize(
    "model_type", ["active-full", "passive", "Active-sodium", "active-KLT+H"]
)
@pytest.mark.parametrize("goldwyn", [False, True])
def test_two_dv_is_the_voltage_block_of_rhs_bit_for_bit(
    node, input_node, model_type, goldwyn
):
    kw = {"r1": 8.5, "tau_est": 0.34} if goldwyn else {}
    v0 = two.GOLDWYN_2019["v0"] if goldwyn else -68.0
    P = two.get_params(v0, node, input_node, model_type, **kw)
    rng = np.random.default_rng(node * 10 + input_node)
    for stim, start, stop, amp, _, _, syn_kw in TWO_STIMULI:
        syn = synaptic.SynParams(t_end=20.0, **syn_kw)
        s = _common.stimulus(stim, start, stop, amp, 20.0, syn)
        for t in _times(stim, start, stop, syn.epsg_tau, rng):
            xs = [
                np.array(
                    [v0, v0]
                    + [constants.winf(v0), constants.hinf(v0)] * 4
                    + [constants.ainf(v0)]
                )
            ]
            xs += [
                np.concatenate([rng.uniform(-90, 40, 2), rng.random(9)])
                for _ in range(3)
            ]
            for x in xs:
                want = two._rhs(t, x, P, stim, s)[:2]
                _same_bits(two._dv(t, x, P, stim, s), want)


# --- the stop itself: identical to the full right-hand side's -------------

# (model, stim, start, stop, I, node, t_end, v0, input_node), factor, kw
RUNS = {
    "multi EPSGpair sub": (
        ("multi", "EPSGpair", 5, 5.3, 20.0, 3, 20, -68.0, 1),
        10,
        {},
    ),
    "multi EPSGpair near": (
        ("multi", "EPSGpair", 5, 5.1, 48.0, 3, 20, -68.0, 1),
        10,
        {},
    ),
    "multi ramp factor 30": (
        ("multi", "ramp", 5, 8, 3000.0, 3, 30, -68.0, 1),
        30,
        {},
    ),
    "multi dendrites bilateral": (
        ("multi", "EPSGbilateral", 5, 5.1, 30.0, 3, 20, -68.0, LAT),
        10,
        {"morph": D, "input_node2": MED},
    ),
    "two EPSGpair near": (
        ("two", "EPSGpair", 5, 5.1, 47.0, 3, 20, -68.0, 1),
        10,
        {},
    ),
    "two ramp factor 30": (
        ("two", "ramp", 5, 8, 3000.0, 3, 30, -68.0, 1),
        30,
        {},
    ),
}


def _old_settled_event(monkeypatch, n_v, calls):
    """Restores the settled check that sliced the full right-hand side.

    settled_event as it was before _dv: the event ignores the dv it is
    given and calls the rhs that _solve.integrate receives, keeping its
    first n_v entries. Counts those calls in calls["n"].
    """
    holder = {}
    integrate = _solve.integrate

    def settled_event(dv, t_off, axon_col, factor):  # noqa: ARG001
        def event(t, x):
            if t < t_off or x[axon_col] - x[0] >= factor / 2:
                return -1.0
            calls["n"] += 1
            d = np.asarray(holder["rhs"](t, x)[:n_v])
            return 1.0 if np.all(d < 0) else -1.0

        event.terminal = True
        event.direction = 1
        return event

    def keep_rhs(rhs, *args, **kw):
        holder["rhs"] = rhs
        return integrate(rhs, *args, **kw)

    monkeypatch.setattr(_solve, "settled_event", settled_event)
    monkeypatch.setattr(_solve, "integrate", keep_rhs)


CASES = [(name, "exact") for name in sorted(RUNS)] + [
    pytest.param(
        name,
        "fast",
        marks=pytest.mark.skipif(
            not _fast.available(), reason="needs the fast extra"
        ),
    )
    for name in sorted(RUNS)
    if name.startswith("multi")
]


@pytest.mark.parametrize(("name", "mode"), CASES)
def test_settled_stop_is_unchanged_by_the_voltage_only_check(
    name, mode, monkeypatch
):
    args, factor, kw = RUNS[name]
    if args[0] == "multi":
        kw = {**kw, "mode": mode}
    t, y = _dispatch.run_model(
        *args, stop_on_spike=factor, stop_when_settled=True, **kw
    )
    # a settled stop: before t_end, axon - soma at or below factor / 2
    axon = args[5] - 1 if args[0] == "multi" else 1
    assert t[-1] < args[6]
    assert y[-1, axon] - y[-1, 0] <= factor / 2
    calls = {"n": 0}
    n_v = kw.get("morph", multi.LUMPED).n if args[0] == "multi" else 2
    with monkeypatch.context() as m:
        _old_settled_event(m, n_v, calls)
        t_full, y_full = _dispatch.run_model(
            *args, stop_on_spike=factor, stop_when_settled=True, **kw
        )
    assert calls["n"] > 0  # the old check did run
    assert t.tobytes() == t_full.tobytes()
    assert y.tobytes() == y_full.tobytes()
