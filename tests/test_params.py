"""The parameter dataclasses: frozen, identity-hashed, keys unchanged."""

import dataclasses

import pytest

from msoaxon import _common, _dispatch, coincidence, multi, synaptic, two


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
    assert (plain.VsynE, plain.diff) == (None, None)
    assert plain.epsg_tau == synaptic.EPSG_TAU
    assert stimulus("sine").f == synaptic.SynParams().f
    syn = stimulus("SynapticPair")
    assert syn.t_syn is not None and syn.g_syn is not None
    assert (syn.VsynE, syn.diff) == (0.0, 1.0)


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
    assert m.key == (-68.0, 1.0, 1.0, -62.5, False, 1.0)
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
    assert m.key == (-60.0, 3.0, 2.0, -77.0, True, 0.5)
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
        "soma_na_vhalf=-62.5, rebalance_rest=False, dendrite_klt_scale=1.0)"
    )
