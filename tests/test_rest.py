"""The rebalanced rest (#128): soma-only rebalance, steady state, stability.

membrane(rebalance_rest=True) rebalances the soma, AIS and dendrites to v0
and the axon to axon_rest (default -68 mV), root-finds the steady state when
the two differ, and warns when the rest is unstable. A whole cable rebalanced
to -58 mV fires from its nodes at 461 Hz.
"""

import importlib.util
import pathlib
import warnings

import numpy as np
import pytest

from msoaxon import _fast, constants as C, multi

REPO = pathlib.Path(__file__).resolve().parents[1]
D = multi.with_dendrites()
needs_ext = pytest.mark.skipif(
    not _fast.available(), reason="needs the fast extra (msoaxon_fast)"
)


def _uniform_vlk(v0, morph=multi.LUMPED, **knobs):
    """membrane()'s leak reversals before #128: every compartment at v0."""
    mem = multi.membrane(v0, morph=morph, **knobs)  # for g_klt, na_vhalf
    V = np.full(morph.n, float(v0))
    mi, hi, pi, wi, zi = (
        C.minf(V),
        C.hinf(V, mem.na_vhalf),
        C.pinf(V),
        C.winf(V),
        C.zinf(V),
    )
    I_ion = (
        C.na_open(mi, hi, morph.g_na) * (V - multi.V_NA)
        + morph.g_kht * pi * (V - multi.V_K)
        + mem.g_klt * wi**4 * zi * (V - multi.V_K)
        + morph.g_h * C.ainf(V) * (V - multi.V_H)
    )
    return V + I_ion / morph.g_lk


def _quiet(mem):
    """Builds mem with UnstableRestWarning raised as an error."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", multi.UnstableRestWarning)
        return mem()


def _residual(mem):
    x = multi.initial_state(mem.morph, mem, mem.key[0])
    return np.abs(
        multi._rhs(0.0, x, "none", None, 1, (True,) * 6, mem, mem.morph)
    ).max()


@pytest.mark.parametrize(
    "v0, morph, knobs",
    [
        (-68.0, multi.LUMPED, {}),
        (-68.0, multi.LUMPED, {"soma_klt_scale": 0, "ais_klt_scale": 0}),
        (-68.0, multi.LUMPED, {"soma_klt_scale": 3, "soma_na_vhalf": -77}),
        (-68.0, D, {"dendrite_klt_scale": 2.0}),
        (-58.0, multi.LUMPED, {"axon_rest": -58.0}),
    ],
)
def test_uniform_rebalance_is_unchanged_bit_for_bit(v0, morph, knobs):
    # at v0 = -68 the default axon_rest is v0, so the cable is rebalanced as
    # before #128 (the rebalanced golden traces depend on it); so is a whole
    # cable rebalanced to v0 by asking for it
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", multi.UnstableRestWarning)
        mem = multi.membrane(v0, morph=morph, rebalance_rest=True, **knobs)
    rest = {k: v for k, v in knobs.items() if k != "axon_rest"}
    assert mem.vlk.tobytes() == _uniform_vlk(v0, morph, **rest).tobytes()
    assert mem.v_rest is None
    x = multi.initial_state(morph, mem, v0)
    assert np.array_equal(x[: morph.n], np.full(morph.n, v0))


def test_axon_rest_needs_a_rebalance():
    with pytest.raises(ValueError, match="axon_rest needs rebalance_rest"):
        multi.membrane(-58.0, axon_rest=-58.0)


def test_axon_rest_is_part_of_the_key():
    a = multi.membrane(-60.0, rebalance_rest=True)
    b = _quiet(
        lambda: multi.membrane(-60.0, rebalance_rest=True, axon_rest=-64)
    )
    assert a.key[-1] == multi.AXON_REST == -68.0
    assert b.key[-1] == -64.0
    assert a.key != b.key


@pytest.mark.parametrize("vhalf", [-62.5, -77.0])
def test_soma_only_rebalance_holds_the_soma_and_rests_the_axon(vhalf):
    mem = _quiet(
        lambda: multi.membrane(-58.0, soma_na_vhalf=vhalf, rebalance_rest=True)
    )
    V = mem.v_rest
    assert V[0] == -58.0  # held exactly: vlk[0] is the unknown
    assert np.allclose(V[1:3], -58.0, atol=0.5)  # the AIS follows the soma
    nodes = np.array(mem.morph.labels) == "node"
    # the nodes stay below their Hopf point (-60.5 mV), the distal ones
    # near the axon's own -68
    assert V[nodes].max() < -60.0
    assert V[-1] == pytest.approx(-67.87, abs=0.05)
    # the axon's leak reversals are the -68 mV balance, untouched
    axon = np.isin(mem.morph.labels, ("internode", "node"))
    assert np.array_equal(mem.vlk[axon], _uniform_vlk(-68.0)[axon])


@pytest.mark.parametrize(
    "build",
    [
        lambda: multi.membrane(-58.0, rebalance_rest=True),
        lambda: multi.membrane(-58.0, soma_na_vhalf=-77, rebalance_rest=True),
        lambda: multi.membrane(
            -60.0, soma_klt_scale=3, soma_na_vhalf=-77, rebalance_rest=True
        ),
        lambda: multi.membrane(-58.0, rebalance_rest=True, morph=D),
        lambda: multi.membrane(-58.0, rebalance_rest=True, klt_to_leak="slope"),
        lambda: multi.membrane(-58.0, rebalance_rest=True, klt_to_leak="chord"),
    ],
    ids=["-58", "-58-vhalf", "-60-klt3", "-58-dendrites", "slope", "chord"],
)
def test_rebalanced_rest_is_a_steady_state(build):
    mem = _quiet(build)
    assert mem.v_rest is not None
    assert _residual(mem) < 1e-9


@pytest.mark.parametrize(
    "mode", ["exact", pytest.param("native", marks=needs_ext)]
)
@pytest.mark.parametrize(
    "morph", [multi.LUMPED, D], ids=["lumped", "dendrites"]
)
def test_rebalanced_cell_stays_at_rest(mode, morph):
    mem = _quiet(
        lambda: multi.membrane(-58.0, rebalance_rest=True, morph=morph)
    )
    n = morph.n
    x0 = multi.initial_state(morph, mem, -58.0)
    t, y = multi.mso_axon(
        "step", 5, 10, 0.0, 3, "active-full", 20, -58.0, 1, mem=mem, mode=mode
    )
    assert t[-1] == 20
    assert np.abs(y[:, :n] - x0[:n]).max() < 1e-6


def test_whole_cable_rebalance_to_minus_58_warns():
    with pytest.warns(multi.UnstableRestWarning, match="node compartments"):
        multi.membrane(-58.0, rebalance_rest=True, axon_rest=-58.0)


@pytest.mark.parametrize("vhalf", [-62.5, -77.0])
def test_soma_only_rebalance_to_minus_58_does_not_warn(vhalf):
    _quiet(
        lambda: multi.membrane(-58.0, soma_na_vhalf=vhalf, rebalance_rest=True)
    )


def test_rest_jacobian_finds_the_nodal_hopf_mode():
    # #128's linearisation: +1.088 +- 2.250j /ms at a uniform -58, and the
    # soma-only cell's slowest mode decays (-0.003 /ms)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", multi.UnstableRestWarning)
        uniform = multi.membrane(-58.0, rebalance_rest=True, axon_rest=-58.0)
    lam = np.linalg.eigvals(multi.rest_jacobian(uniform)[0])
    lead = lam[np.argmax(lam.real)]
    assert lead.real == pytest.approx(1.088, abs=2e-3)
    assert abs(lead.imag) == pytest.approx(2.250, abs=2e-3)
    soma_only = multi.membrane(-58.0, rebalance_rest=True)
    lam = np.linalg.eigvals(multi.rest_jacobian(soma_only)[0])
    assert -0.01 < lam.real.max() < 0


def _na_availability():
    spec = importlib.util.spec_from_file_location(
        "na_availability", REPO / "scripts" / "na_availability.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@needs_ext
def test_ripple_watches_every_compartment():
    script = _na_availability()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", multi.UnstableRestWarning)
        firing = multi.membrane(-58.0, rebalance_rest=True, axon_rest=-58.0)
    # the nodes fire 110 mV spikes; the soma sees 0.18 mV of them
    rip = script.ripple("multi", -58.0, {"mem": firing, "mode": "native"})
    assert rip > 50
    # each compartment is compared with its own rest: the soma-only cell's
    # axon rests 10 mV below the soma and is steady
    steady = multi.membrane(-58.0, rebalance_rest=True)
    rip = script.ripple("multi", -58.0, {"mem": steady, "mode": "native"})
    assert rip < script.RIPPLE_MAX
