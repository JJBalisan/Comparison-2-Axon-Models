"""The parameter dataclasses: frozen, identity-hashed, keys unchanged."""

import dataclasses

import numpy as np
import pytest

from msoaxon import (
    _common,
    _dispatch,
    coincidence,
    constants as C,
    multi,
    synaptic,
    two,
)


def stimulus(stim_type="step"):
    return _common.stimulus(
        stim_type, 5.0, 6.0, 60.0, 20.0, synaptic.SynParams()
    )


def objects():
    """One instance of each bundle, with a field to try assigning."""
    return {
        "Stimulus": (stimulus(), "start"),
        "Morphology": (multi.with_dendrites(), "n"),
        "Membrane": (multi.membrane(-68.0), "vlk"),
        "TwoCptParams": (two.get_params(-68.0, 3, 1, "active-full"), "gC"),
    }


@pytest.mark.parametrize("name", list(objects()))
def test_frozen(name):
    obj, field = objects()[name]
    assert dataclasses.is_dataclass(obj)
    with pytest.raises(dataclasses.FrozenInstanceError):
        setattr(obj, field, 0)


@pytest.mark.parametrize("name", list(objects()))
def test_identity_equality_and_hash(name):
    # eq=False: no generated __eq__ (it would compare ndarray fields), so
    # equal contents are still different objects and hash by identity
    a, _ = objects()[name]
    b, _ = objects()[name]
    assert a is not b
    assert a != b
    assert a == a
    assert hash(a) == object.__hash__(a)
    assert len({a, b, a}) == 2


def test_equal_membranes_are_not_equal():
    a, b = multi.membrane(-68.0), multi.membrane(-68.0)
    assert a.key == b.key
    assert a != b
    assert hash(a) != hash(b)


def test_dataclasses_declare_every_field():
    for cls in (
        _common.Stimulus,
        multi.Morphology,
        multi.Membrane,
        two.TwoCptParams,
    ):
        for field in dataclasses.fields(cls):
            assert field.type  # annotated
    # a missing field fails at construction, not inside a right-hand side
    with pytest.raises(TypeError):
        _common.Stimulus(**{"start": 1.0})


def test_stimulus_optional_fields():
    plain = stimulus("EPSG")
    assert (plain.f, plain.t_syn, plain.g_syn) == (None, None, None)
    assert (plain.VsynE, plain.diff, plain.events) == (None, None, None)
    assert plain.epsg_tau == synaptic.EPSG_TAU
    assert stimulus("sine").f == synaptic.SynParams().f
    syn = stimulus("SynapticPair")
    assert syn.t_syn is not None and syn.g_syn is not None
    assert (syn.VsynE, syn.diff) == (0.0, 1.0)


def _bits(x):
    return np.float64(x).tobytes()


@pytest.mark.parametrize("stim_type", ["ramp", "ramp2"])
def test_ramp_fields_equal_the_old_inline_expressions(stim_type):
    # what multi.external_current and two.applied_current recomputed per call
    start, stop = 5.3, 8.1
    s = _common.stimulus(
        stim_type, start, stop, 60.0, 20.0, synaptic.SynParams()
    )
    slope = 1 / (stop - start)
    assert _bits(s.slope) == _bits(slope)
    if stim_type == "ramp":
        assert _bits(s.t_top) == _bits(start + 1 / slope)
        assert _bits(s.t_end_local) == _bits(stop + 5)
    else:
        assert (s.t_top, s.t_end_local) == (None, None)
    assert s.two_pi_f is None


def test_sine_field_equals_the_old_inline_expression():
    syn = synaptic.SynParams(f=317.0)
    s = _common.stimulus("sine", 5.0, 9.0, 60.0, 20.0, syn)
    assert _bits(s.two_pi_f) == _bits(2 * np.pi * syn.f)
    for t in (5.0, 5.37, 8.123):
        old = np.sin(2 * np.pi * s.f * (t - s.start) / 1000)
        assert _bits(np.sin(s.two_pi_f * (t - s.start) / 1000)) == _bits(old)
    assert (s.slope, s.t_top, s.t_end_local) == (None, None, None)


def test_inv_peak_equals_the_old_inline_expressions():
    default = stimulus("EPSG")
    assert _bits(default.inv_peak) == _bits(1 / 0.21317)
    tau = (0.1, 0.3)  # Myoga et al 2014: decay 0.3
    s = _common.stimulus(
        "EPSG", 5.0, 9.0, 60.0, 20.0, synaptic.SynParams(epsg_tau=tau)
    )
    rise, decay = tau
    t_peak = np.log(decay / rise) * rise * decay / (decay - rise)
    peak = np.exp(-t_peak / decay) - np.exp(-t_peak / rise)
    assert _bits(s.inv_peak) == _bits(1 / peak)
    # epsg_unitary gives the same bits whether it is handed the field or not
    t = np.linspace(-2.0, 6.0, 801)
    expected = (1 / peak) * (t > 0) * (np.exp(-t / decay) - np.exp(-t / rise))
    for given in (s.inv_peak, None):
        got = synaptic.epsg_unitary(t, s.epsg_tau, given)
        assert got.tobytes() == expected.tobytes()


def test_morphology_keys():
    assert multi.LUMPED.key == "lumped"
    assert multi.with_dendrites().key == (
        "dendrites",
        200.0,
        5.0,
        5,
        74.0,
        True,
        100.0,
    )
    assert multi.with_dendrites(
        length=150, n_seg=3, klt_lambda=50, conserve_totals=False
    ).key == ("dendrites", 150.0, 5.0, 3, 50.0, False, 100.0)


def test_morphology_dendrite_fields():
    assert multi.LUMPED.lateral is None and multi.LUMPED.medial is None
    d = multi.with_dendrites()
    assert d.lateral is not None and d.medial is not None
    assert list(d.lateral) == [46, 47, 48, 49, 50]
    assert list(d.medial) == [51, 52, 53, 54, 55]


def test_membrane_keys():
    m = multi.membrane(-68.0)
    assert m.key == (-68.0, 1.0, 1.0, -62.5, False, 1.0, -68.0)
    assert m.morph_key == "lumped"
    assert m.morph is multi.LUMPED
    assert m.y0_gates is None
    d = multi.with_dendrites()
    m = multi.membrane(
        -60,
        soma_klt_scale=3,
        ais_klt_scale=2,
        soma_na_vhalf=-77,
        rebalance_rest=True,
        morph=d,
        dendrite_klt_scale=0.5,
    )
    assert m.key == (-60.0, 3.0, 2.0, -77.0, True, 0.5, -68.0)
    assert m.morph_key == d.key
    assert m.y0_gates is not None and len(m.y0_gates) == 6


def test_keys_feed_the_caches():
    m = multi.membrane(-68.0)
    assert coincidence._frozen("mem", m) == ("key", m.key, "lumped")
    assert coincidence._frozen("morph", multi.LUMPED) == (
        "key",
        "lumped",
        None,
    )
    assert len(_dispatch._MEMBRANE_KEY) == len(m.key)
    assert _dispatch.describe("multi", 3, -68.0, mem=m).endswith(
        "mem=membrane(v0=-68.0, soma_klt_scale=1.0, ais_klt_scale=1.0, "
        "soma_na_vhalf=-62.5, rebalance_rest=False, dendrite_klt_scale=1.0, "
        "axon_rest=-68.0)"
    )


def _old_multi_y0(morph, mem, v0):
    """mso_axon's inline initial state before initial_state existed."""
    n = morph.n
    if mem.y0_gates is None:
        return np.concatenate(
            [np.full(n, float(v0))] + [np.full(n, g) for g in multi._GATE0]
        )
    return np.concatenate([np.full(n, float(v0)), *mem.y0_gates])


@pytest.mark.parametrize("rebalance_rest", [False, True])
@pytest.mark.parametrize("dendrites", [False, True])
@pytest.mark.parametrize("v0", [-68.0, -60.0, -65])
def test_multi_initial_state_equals_the_old_inline_code(
    rebalance_rest, dendrites, v0
):
    morph = multi.with_dendrites() if dendrites else multi.LUMPED
    mem = multi.membrane(v0, morph=morph, rebalance_rest=rebalance_rest)
    if mem.v_rest is not None:
        # a rest rebalanced away from the axon's -68 mV is root-found, not
        # flat (#128): V starts from the solved profile, gates at its steady
        # state; tests/test_rest.py checks that it is one
        new = multi.initial_state(morph, mem, v0)
        V = mem.v_rest
        # the old inline code's gate blocks, at the solved V
        gates = (
            C.minf(V),
            C.hinf(V, mem.na_vhalf),
            C.pinf(V),
            C.winf(V),
            C.zinf(V),
            C.ainf(V),
        )
        want = np.concatenate([V, *gates])
        assert new.dtype == want.dtype
        assert new.shape == (7 * morph.n,)
        assert new.tobytes() == want.tobytes()
        assert new[0] == v0
        return
    new, old = (
        multi.initial_state(morph, mem, v0),
        _old_multi_y0(morph, mem, v0),
    )
    assert new.dtype == old.dtype
    assert new.shape == (7 * morph.n,)
    assert np.array_equal(new, old)
    assert new.tobytes() == old.tobytes()


def _old_two_x0(params):
    """two_cpt's inline initial state before initial_state existed."""
    Vr = params.Vrest
    return [
        Vr,
        Vr,
        C.winf(Vr),
        C.hinf(Vr),
        C.winf(Vr),
        C.minf(Vr),
        C.minf(Vr),
        C.hinf(Vr),
        C.pinf(Vr),
        C.ainf(Vr),
        C.ainf(Vr),
    ]


@pytest.mark.parametrize("v0", [-68.0, -58.0, -65])
@pytest.mark.parametrize("model_type", ["passive", "active-full"])
def test_two_initial_state_equals_the_old_inline_code(v0, model_type):
    P = two.get_params(v0, 3, 1, model_type)
    new = two.initial_state(P)
    old = np.asarray(_old_two_x0(P), dtype=float)  # as integrate() takes it
    assert new.dtype == old.dtype
    assert np.array_equal(new, old)
    assert new.tobytes() == old.tobytes()


@pytest.mark.parametrize(
    "bundle",
    [
        multi.LUMPED,
        multi.with_dendrites(),
        multi.membrane(-68.0),
        multi.membrane(-60.0, soma_klt_scale=0.5, rebalance_rest=True),
        multi.membrane(-68.0, morph=multi.with_dendrites()),
    ],
    ids=["lumped", "dendrites", "membrane", "rebalanced", "dendritic-mem"],
)
def test_bundle_arrays_are_read_only(bundle):
    # a mutated bundle would poison the quiet-prefix cache and the spike
    # heights, which are keyed by the bundle's settings
    arrays = []
    for field in dataclasses.fields(bundle):
        value = getattr(bundle, field.name)
        for a in value if isinstance(value, tuple) else (value,):
            if isinstance(a, np.ndarray):
                arrays.append(a)
    assert arrays
    for a in arrays:
        assert not a.flags.writeable
        with pytest.raises(ValueError, match="read-only"):
            a[...] = 0
