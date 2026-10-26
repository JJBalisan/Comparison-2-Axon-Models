"""The Rust right-hand side (mode="fast") against multi._rhs, and the dispatch.

The native mode (the Rust integrator over the same right-hand side) is
tested at the end: its plumbing, its segment contract and its events.

Equivalence: the Rust port must match multi._rhs to rtol 1e-12 on every
configuration it covers, and bit for bit on darwin-arm64, where it was
written and where numpy and Rust call the same libm. Those tests need the
`fast` extra and skip without it; the dispatch tests run either way (they
fake or hide the extension).
"""

import multiprocessing
import os
import platform
import sys
import warnings

import numpy as np
import pytest

import msoaxon
from msoaxon import (
    _common,
    _dispatch,
    _fast,
    _parallel,
    _solve,
    coincidence,
    constants,
    multi,
    synaptic,
    two,
)

BITWISE = sys.platform == "darwin" and platform.machine() == "arm64"
needs_ext = pytest.mark.skipif(
    not _fast.available(), reason="needs the fast extra (msoaxon_fast)"
)

D = multi.with_dendrites()
D_RA = multi.with_dendrites(conserve_totals=False, dendrite_ra=200.0)
LAT, MED = int(D.lateral[2]), int(D.medial[2])
DTX = {"soma_klt_scale": 0, "ais_klt_scale": 0, "rebalance_rest": True}
MYOGA = (0.1, 0.3)
IPSG = (0.28, 1.5)  # Roberts et al 2013 rise, Myoga et al 2014 decay
# the "Events" lists: mixed kinetics, reversals and compartments, then
# events that share a compartment and an onset (np.bincount adds them in
# list order, so the sum's last bits depend on it)
EVENTS = (
    synaptic.Event(5.0, 1, 60.0),
    synaptic.Event(5.2, 3, 40.0, MYOGA),
    synaptic.Event(4.8, 1, 30.0, IPSG, -80.0),
)
EVENTS_SHARED = (
    synaptic.Event(5.0, 1, 30.0),
    synaptic.Event(5.0, 1, 20.0, MYOGA, -20.0),
    synaptic.Event(5.0, 2, 10.0),
    synaptic.Event(5.1, 1, 25.0),
    synaptic.Event(5.0, 1, 0.1),
)

# (stim, start, stop, I, input_node, input_node2, SynParams keywords); ints
# where the scripts and goldens pass ints
STIMULI = {
    "none": ("none", 5, 10, 0, 1, None, {}),
    "step": ("step", 5, 10, 8500, 3, None, {}),
    "ramp": ("ramp", 5, 5.5, 8500, 1, None, {}),
    "ramp-long": ("ramp", 5.0, 5.3, 12000.0, 3, None, {}),
    "ramp2": ("ramp2", 5, 5.5, 8500, 1, None, {}),
    "sine": ("sine", 5, 7.5, 8500, 1, None, {}),
    "sine-int-f": ("sine", 5.0, 6.25, 9000.0, 3, None, {"f": 400}),
    "Synaptic": ("Synaptic", 5, 15, 0, 1, None, {"random_in": 1396}),
    "SynapticPair": ("SynapticPair", 5, 15, 0, 3, None, {"diff": 1.0}),
    "EPSG": ("EPSG", 5, 10, 108, 1, None, {}),
    "EPSG-myoga": ("EPSG", 5.0, 10.0, 60.0, 3, None, {"epsg_tau": MYOGA}),
    "EPSGpair": ("EPSGpair", 5, 5.2, 100, 1, None, {}),
    "EPSGpair-myoga": (
        "EPSGpair",
        5.0,
        5.3,
        49.0,
        1,
        None,
        {"epsg_tau": MYOGA},
    ),
    "EPSGbilateral": ("EPSGbilateral", 5.0, 5.2, 60.0, 1, 3, {}),
    "EPSGbilateral-same": ("EPSGbilateral", 5.0, 5.1, 60.0, 2, 2, {}),
    # start and stop are derived from the onsets, input_node is ignored; the
    # run's I (the amplitude) scales every peak
    "Events": ("Events", 0.0, 0.0, 1.0, 1, None, {"events": EVENTS}),
    "Events-scaled": ("Events", 0.0, 0.0, 1.7, 1, None, {"events": EVENTS}),
    "Events-shared": (
        "Events",
        0.0,
        0.0,
        1.0,
        1,
        None,
        {"events": EVENTS_SHARED},
    ),
}


def configs():
    """name -> (stim settings, model_type, mem, morph)."""
    out = {}
    for name, stim in STIMULI.items():
        out[name] = (stim, "active-full", None, multi.LUMPED)
    for mt in multi.ACTIVE_GATES:
        out[f"type/{mt}"] = (STIMULI["EPSGpair"], mt, None, multi.LUMPED)
    out["mem/dendrotoxin"] = (
        STIMULI["step"],
        "active-full",
        multi.membrane(-68.0, **DTX),
        multi.LUMPED,
    )
    out["mem/soma-na-vhalf"] = (
        STIMULI["EPSGpair"],
        "active-full",
        multi.membrane(
            -60.0, soma_klt_scale=3, soma_na_vhalf=-77, rebalance_rest=True
        ),
        multi.LUMPED,
    )
    bilateral = ("EPSGbilateral", 5.0, 5.1, 60.0, LAT, MED, {})
    for label, morph in (("dendrites", D), ("dendrites-ra200", D_RA)):
        for name in ("step", "sine", "Synaptic", "EPSG", "EPSGpair"):
            stim, start, stop, amp, _, k2, syn_kw = STIMULI[name]
            out[f"{label}/{name}"] = (
                (stim, start, stop, amp, LAT, k2, syn_kw),
                "active-full",
                None,
                morph,
            )
        out[f"{label}/EPSGbilateral"] = (bilateral, "active-full", None, morph)
        # one event on each dendrite, one on the soma and one in the axon
        events = (
            synaptic.Event(5.0, LAT, 60.0),
            synaptic.Event(5.1, MED, 60.0),
            synaptic.Event(5.1, 1, 20.0, IPSG, -80.0),
            synaptic.Event(5.3, 3, 10.0),
        )
        out[f"{label}/Events"] = (
            ("Events", 0.0, 0.0, 1.0, 1, None, {"events": events}),
            "active-full",
            None,
            morph,
        )
    out["dendrites/mem"] = (
        bilateral,
        "active-KLT+H",
        multi.membrane(
            -68.0, morph=D, dendrite_klt_scale=2.0, rebalance_rest=True
        ),
        D,
    )
    return out


CONFIGS = configs()


def setup(name):
    """Returns (rhs args after t, x; the times worth sampling)."""
    (stim, start, stop, amp, k, k2, syn_kw), mt, mem, morph = CONFIGS[name]
    syn = synaptic.SynParams(t_end=20.0, **syn_kw)
    s = _common.stimulus(stim, start, stop, amp, 20.0, syn)
    mem = mem or multi.membrane(-68.0, morph=morph)
    active = tuple(bool(g) for g in multi.ACTIVE_GATES[mt])
    # the window edges, both sides, and the ramp's own breakpoints
    edges = [start, stop, 5.0, stop + 5.0, 0.0]
    if stim in ("ramp", "ramp2"):
        slope = 1 / (stop - start)
        edges.append(start + 1 / slope)
    times = [e + d for e in edges for d in (0.0, -1e-9, 1e-9)]
    if s.events is not None:  # every onset, both sides
        times += [t + d for t in s.events["t"] for d in (0.0, -1e-9, 1e-9)]
    times += list(np.random.default_rng(len(name)).uniform(0, 20, 12))
    if stim in ("Synaptic", "SynapticPair"):
        times += list(s.t_syn[500:1500:97])  # exactly on table points
    return (stim, s, k, active, mem, morph, k2), times


def states(n, rng, count):
    """Random states: V in -90..40 mV, gates in 0..1, plus rest."""
    out = [np.concatenate([np.full(n, -68.0), np.full(6 * n, 0.3)])]
    for _ in range(count):
        out.append(np.concatenate([rng.uniform(-90, 40, n), rng.random(6 * n)]))
    return out


def assert_same(got, want, n):
    if BITWISE:
        np.testing.assert_array_equal(got, want)
        return
    for block in range(7):
        g, w = (
            got[block * n : (block + 1) * n],
            want[block * n : (block + 1) * n],
        )
        np.testing.assert_allclose(
            g, w, rtol=1e-12, atol=1e-12 * float(np.max(np.abs(w)))
        )


@needs_ext
@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_rust_rhs_matches_numpy(name):
    args, times = setup(name)
    fast = _fast._build(*args)
    morph = args[5]
    rng = np.random.default_rng(12)
    for t in times:
        for x in states(morph.n, rng, 3):
            want = multi._rhs(float(t), x, *args)
            assert_same(fast(float(t), x), want, morph.n)


@needs_ext
def test_rust_rhs_reads_strided_states():
    # scipy's num_jac passes columns of a 2-D array
    args, _ = setup("EPSGpair")
    fast = _fast._build(*args)
    X = np.stack(states(multi.N, np.random.default_rng(3), 4), axis=1)
    for j in range(X.shape[1]):
        col = X[:, j]
        assert not col.flags.c_contiguous
        np.testing.assert_array_equal(
            fast(5.1, col), fast(5.1, np.ascontiguousarray(col))
        )


@needs_ext
def test_rust_epsg_pair_is_finite_far_before_the_second_onset():
    # the second EPSG starts 80 ms after the first: the exponentials of its
    # waveform overflow at te - td < -71 unless clamped (as in the Python)
    syn = synaptic.SynParams(t_end=20.0)
    s = _common.stimulus("EPSGpair", 5.0, 85.0, 100.0, 20.0, syn)
    active = (True,) * 6
    mem = multi.membrane(-68.0)
    fast = _fast._build("EPSGpair", s, 1, active, mem, multi.LUMPED, None)
    x = states(multi.N, np.random.default_rng(8), 1)[1]
    assert np.all(np.isfinite(fast(5.1, x)))


@needs_ext
def test_rust_rhs_rejects_a_wrong_state_length():
    args, _ = setup("EPSGpair")
    with pytest.raises(ValueError, match="315"):
        _fast._build(*args)(5.0, np.zeros(10))


# the constants the Rust side takes from Python: (module, name, new value)
PASSED_CONSTANTS = [
    (multi, "V_NA", multi.V_NA + 1.0),
    (multi, "V_K", multi.V_K + 1.0),
    (multi, "V_H", multi.V_H + 1.0),
    (constants, "_P_TEMP", constants._P_TEMP * 1.5),
    (constants, "_A_TEMP", constants._A_TEMP * 1.5),
]


@needs_ext
@pytest.mark.parametrize(
    ("module", "name", "value"),
    PASSED_CONSTANTS,
    ids=[f"{m.__name__}.{n}" for m, n, _ in PASSED_CONSTANTS],
)
def test_rust_rhs_follows_a_changed_python_constant(
    monkeypatch, module, name, value
):
    # the Rust side has no numbers of its own: move one in Python and its
    # right-hand side must move exactly as multi._rhs does, so a change to
    # the model can't be mistaken for hardware drift
    args, _ = setup("EPSGpair")
    n = args[5].n
    x = states(n, np.random.default_rng(5), 1)[1]
    before = multi._rhs(5.1, x, *args)
    with monkeypatch.context() as mp:
        mp.setattr(module, name, value)
        want = multi._rhs(5.1, x, *args)
        assert not np.array_equal(want, before)
        assert_same(_fast._build(*args)(5.1, x), want, n)
    # and back to the original once the patch is undone
    assert_same(_fast._build(*args)(5.1, x), before, n)


@needs_ext
@pytest.mark.parametrize("change", ["unknown", "missing"])
def test_rust_rhs_rejects_a_mismatched_gating_dict(monkeypatch, change):
    args, _ = setup("EPSGpair")
    gating = dict(_fast._GATING)
    if change == "unknown":
        gating["tauq_c"] = 1.0
        match = 'unknown gating constant "tauq_c"'
    else:
        del gating["tauz_c"]
        match = 'missing gating constant "tauz_c"'
    monkeypatch.setattr(_fast, "_GATING", gating)
    with pytest.raises(ValueError, match=match):
        _fast._build(*args)


@needs_ext
@pytest.mark.skipif(not BITWISE, reason="bit-identical on darwin-arm64 only")
@pytest.mark.parametrize(
    "name", ["EPSGpair", "Synaptic", "Events", "dendrites/mem"]
)
def test_runs_are_bit_identical_here(name):
    (stim, start, stop, amp, k, k2, syn_kw), mt, mem, morph = CONFIGS[name]
    syn = synaptic.SynParams(t_end=15.0, **syn_kw)
    kw = {"syn": syn, "mem": mem, "morph": morph, "input_node2": k2}
    args = (stim, start, stop, amp, 3, mt, 15.0, -68.0, k)
    t1, y1 = multi.mso_axon(*args, **kw, mode="exact")
    t2, y2 = multi.mso_axon(*args, **kw, mode="fast")
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_array_equal(y1, y2)


@needs_ext
def test_self_check_passes():
    differ, total, worst = _fast.self_check()
    assert total == 24
    if BITWISE:
        assert (differ, worst) == (0, 0.0)
    else:
        # elementwise, so cancellation in (x_inf - x) / tau shows; GitHub's
        # ubuntu x86_64 runner measured 2.7e-13 (every state differed)
        assert worst < 1e-10


@needs_ext
def test_self_check_is_cheap():
    # a wall-clock bound, so it is skipped (not loosened) when the machine
    # is busier than its CPUs: it read 0.45 s at a load average near 180
    import time  # noqa: PLC0415

    getloadavg = getattr(os, "getloadavg", None)
    if getloadavg is not None:
        load, cpus = getloadavg()[0], os.cpu_count() or 1
        if load > cpus:
            pytest.skip(f"load average {load:.0f} on {cpus} CPUs")
    best = np.inf
    for _ in range(3):
        t0 = time.perf_counter()
        _fast.self_check()
        best = min(best, time.perf_counter() - t0)
    assert best < 0.05


# --- Table::interp against np.interp -----------------------------------------


def interp_g_like(xs, ys, x):
    """synaptic.interp_g's call: zero outside the table."""
    return float(np.interp(x, xs, ys, left=0.0, right=0.0))


def assert_interp_same(got, want):
    if np.isnan(want):
        assert np.isnan(got)
    elif BITWISE:
        # np.interp's fused multiply-add on arm64 is reproduced in Rust
        assert np.float64(got).tobytes() == np.float64(want).tobytes()
    else:
        # Rust fuses slope * dx + y on every aarch64 target; numpy's build
        # need not (linux-aarch64 differs): one rounding against two, so
        # they can be an ulp of the operands apart. A relative tolerance
        # alone fails when the sum cancels (got -0.04247429203836899, want
        # -0.042474292038368944); the largest gap over these tables is
        # 4.4e-16 absolute.
        np.testing.assert_allclose(got, want, rtol=1e-15, atol=1e-15)


@needs_ext
def test_rust_interp_matches_numpy_on_the_edges():
    xs = np.array([0.0, 0.3, 1.0, 2.5, 2.5000000000000004])
    ys = np.array([0.0, 2.0, 1.5, -1.0, 4.0])
    ulp = np.nextafter
    points = [*xs, -1e-300, -0.1, 2.6, 9e99, -9e99, np.nan]
    points += [ulp(x, np.inf) for x in xs] + [ulp(x, -np.inf) for x in xs]
    points += [-np.inf, np.inf, 0.1, 0.65, 1.7, 2.0]
    for x in points:
        want = interp_g_like(xs, ys, float(x))
        assert_interp_same(_ext_interp(xs, ys, x), want)


@needs_ext
@pytest.mark.parametrize("x", [1.4, 1.5, 1.6, np.nan, -np.inf, np.inf])
def test_rust_interp_on_a_one_point_table(x):
    xs, ys = np.array([1.5]), np.array([7.0])
    assert_interp_same(_ext_interp(xs, ys, x), interp_g_like(xs, ys, x))


@needs_ext
def test_rust_interp_matches_numpy_on_random_tables():
    rng = np.random.default_rng(2610)
    for _ in range(40):
        n = int(rng.integers(2, 60))
        xs = np.cumsum(rng.uniform(1e-3, 1.0, n))
        ys = rng.normal(size=n)
        ys[rng.integers(0, n, 2)] = ys[0]  # some flat spans
        for x in rng.uniform(xs[0] - 0.2, xs[-1] + 0.2, 60):
            assert_interp_same(
                _ext_interp(xs, ys, float(x)), interp_g_like(xs, ys, float(x))
            )


@needs_ext
def test_rust_interp_matches_a_synaptic_table():
    t_syn, g_syn = synaptic.synaptic(synaptic.SynParams(t_end=20.0))
    rng = np.random.default_rng(14)
    for x in [-0.5, 0.0, 20.0, 20.5, *rng.uniform(0.0, 20.0, 300)]:
        want = synaptic.interp_g(t_syn, g_syn, float(x))
        assert_interp_same(_ext_interp(t_syn, g_syn, float(x)), want)


@needs_ext
def test_rust_interp_rejects_empty_and_ragged_tables():
    with pytest.raises(ValueError, match="non-empty"):
        _ext_interp(np.array([]), np.array([]), 0.0)
    with pytest.raises(ValueError, match="same length"):
        _ext_interp(np.array([0.0, 1.0]), np.array([0.0]), 0.5)


def _ext_interp(xs, ys, x):
    return _fast._ext.interp(xs, ys, float(x))


# --- the drift warning -------------------------------------------------------


@needs_ext
def test_drift_warning_fires_once_per_process(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned", False)
    monkeypatch.delenv(_fast.MODE_ENV, raising=False)
    with pytest.warns(msoaxon.DriftWarning) as record:
        for _ in range(2):
            msoaxon.mso_axon(
                "EPSGpair",
                5,
                5.1,
                30,
                3,
                "active-full",
                8,
                -68.0,
                1,
                mode="fast",
            )
    drift = [w for w in record if issubclass(w.category, _fast.DriftWarning)]
    assert len(drift) == 1
    text = str(drift[0].message)
    assert 'mode="exact"' in text and "filterwarnings" in text
    assert "other hardware" in text
    assert 'the default is mode="native"' in text
    if BITWISE:
        assert "bit-identical" in text
    assert drift[0].filename == __file__  # points at the caller


@needs_ext
def test_drift_warning_reports_a_mismatch(monkeypatch):
    monkeypatch.setattr(_fast, "self_check", lambda: (2, 24, 3e-16))
    with pytest.warns(_fast.DriftWarning, match="2 of 24") as record:
        _fast._warn_drift()
    assert "3.0e-16" in str(record[0].message)


def test_exact_mode_never_warns(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned", False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode="exact"
        )
    assert not _fast._STATE["warned"]


# --- mode resolution ---------------------------------------------------------


@pytest.fixture
def no_env(monkeypatch):
    monkeypatch.delenv(_fast.MODE_ENV, raising=False)


@pytest.fixture
def fake_ext(monkeypatch, no_env):
    """Pretends the extension is installed; its RHS is multi._rhs.

    Records (mode, stim_type) per build; a fake native RHS has no
    `segment`, so scipy integrates it.
    """
    built = []

    def builder(mode):
        def build(
            stim_type, s, input_node, active, mem, morph, input_node2=None
        ):
            built.append((mode, stim_type))
            return lambda t, x: multi._rhs(
                t, x, stim_type, s, input_node, active, mem, morph, input_node2
            )

        return build

    monkeypatch.setattr(_fast, "_ext", object())
    monkeypatch.setattr(_fast, "build", builder("fast"))
    monkeypatch.setattr(_fast, "build_native", builder("native"))
    return built


def test_fallback_without_the_extension_is_silent(monkeypatch, no_env):
    monkeypatch.setattr(_fast, "_ext", None)
    assert not _fast.available()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _fast.resolve(None) == "exact"
        assert _fast.resolve("exact") == "exact"
        t, _ = multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1
        )
    assert t[-1] == 8


def test_asking_for_fast_without_the_extension_warns(monkeypatch, no_env):
    monkeypatch.setattr(_fast, "_ext", None)
    with pytest.warns(UserWarning, match="not installed") as record:
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode="fast"
        )
    assert record[0].filename == __file__
    monkeypatch.setenv(_fast.MODE_ENV, "fast")
    with pytest.warns(UserWarning, match="not installed"):
        assert _fast.resolve(None) == "exact"


def test_default_is_native_when_installed(fake_ext):
    assert _fast.effective(None) == _fast.resolve(None) == "native"
    multi.mso_axon("EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1)
    # run and quiet prefix
    assert fake_ext == [("native", "EPSGpair"), ("native", "none")]


def test_env_overrides_the_default_not_the_argument(fake_ext, monkeypatch):
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    assert _fast.resolve(None) == "exact"
    assert _fast.resolve("fast") == "fast"  # explicit wins
    multi.mso_axon("EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1)
    assert fake_ext == []
    monkeypatch.setenv(_fast.MODE_ENV, "")  # empty means unset
    assert _fast.resolve(None) == "native"


def test_unknown_modes_are_rejected(monkeypatch, no_env):
    with pytest.raises(ValueError, match="mode must be"):
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode="rust"
        )
    monkeypatch.setenv(_fast.MODE_ENV, "turbo")
    with pytest.raises(ValueError, match="turbo"):
        _fast.resolve(None)


@pytest.mark.parametrize("stim", ["ramp", "ramp2"])
@pytest.mark.parametrize("mode", [None, "exact", "fast"])
@pytest.mark.parametrize("stop", [5.0, 4.0, float("nan")])
def test_ramp_needs_stop_after_start_before_any_dispatch(stim, mode, stop):
    with pytest.raises(ValueError, match="ramp needs stop > start"):
        multi.mso_axon(
            stim, 5.0, stop, 30, 3, "active-full", 8, -68.0, 1, mode=mode
        )
    with pytest.raises(ValueError, match="ramp needs stop > start"):
        two.two_cpt(stim, 5.0, stop, 30, 3, "active-full", 8, -68.0, 1)


def test_mode_is_part_of_the_quiet_cache_key(fake_ext, fresh_quiet_cache):
    for mode in ("exact", "fast"):
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode=mode
        )
    modes = sorted(k[0][-1] for k in _solve._QUIET_CACHE)
    assert modes == ["exact", "fast"]


def test_run_model_passes_mode_to_multi_only(fake_ext):
    args = ("EPSGpair", 5, 5.1, 30, 3, 8, -68.0, 1)
    _dispatch.run_model("multi", *args, mode="fast")
    assert fake_ext == [("fast", "EPSGpair"), ("fast", "none")]
    _dispatch.run_model("multi", *args, mode="exact")
    _dispatch.run_model("two", *args, mode="fast")  # no fast two_cpt: ignored
    _dispatch.run_model("two", *args)  # nor a native one
    assert fake_ext == [("fast", "EPSGpair"), ("fast", "none")]


def test_pool_workers_inherit_the_mode_env(monkeypatch, no_env):
    # a forkserver worker forks from a server that keeps the environment it
    # started with; the pool initializer hands it the parent's
    ctx = multiprocessing.get_context("forkserver")
    monkeypatch.setattr(_parallel, "pool_context", lambda: ctx)
    monkeypatch.setitem(_fast._STATE, "warned", False)
    monkeypatch.setitem(_fast._STATE, "warned_native", False)
    with _parallel.process_pool(1) as pool:
        assert list(pool.map(_fast.requested, [None])) == [None]
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    monkeypatch.setitem(_fast._STATE, "warned", True)
    monkeypatch.setitem(_fast._STATE, "warned_native", False)
    with _parallel.process_pool(1) as pool:
        state = pool.submit(_fast.parent_state).result()
    assert state == ("exact", (True, False))


@pytest.mark.skipif(
    "forkserver" not in multiprocessing.get_all_start_methods(),
    reason="no forkserver start method (Windows); the spawn path is covered",
)
def test_pool_workers_inherit_the_settled_stop(monkeypatch):
    ctx = multiprocessing.get_context("forkserver")
    monkeypatch.setattr(_parallel, "pool_context", lambda: ctx)
    for value in (False, True):
        monkeypatch.setattr(_dispatch, "STOP_WHEN_SETTLED", value)
        with _parallel.process_pool(1) as pool:
            assert pool.submit(_dispatch.settled_stop).result() is value


def test_a_missing_extension_imports_silently():
    def missing(name):
        raise ModuleNotFoundError(f"No module named {name!r}", name=name)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _fast._import_ext(missing) is None


@pytest.mark.parametrize(
    "error",
    [
        ImportError("dlopen(msoaxon_fast.abi3.so): symbol not found"),
        ModuleNotFoundError("No module named 'numpy.core'", name="numpy.core"),
    ],
)
def test_a_broken_extension_warns_and_falls_back(error):
    def broken(name):
        raise error

    with pytest.warns(UserWarning, match="installed but failed to import"):
        assert _fast._import_ext(broken) is None


def test_inherit_unsets_a_stale_env(monkeypatch):
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    monkeypatch.setitem(_fast._STATE, "warned", False)
    _fast.inherit(None, False)
    assert _fast.requested() is None and not _fast._STATE["warned"]
    _fast.inherit("fast", True)
    assert _fast.requested() == "fast" and _fast._STATE["warned"]


def test_inherit_takes_both_warned_flags(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned", False)
    monkeypatch.setitem(_fast._STATE, "warned_native", False)
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    _fast.inherit("native", (False, True))
    assert _fast.requested() == "native"
    assert _fast._STATE == {"warned": False, "warned_native": True}
    _fast.inherit("native", True)  # a bool is the fast flag alone
    assert _fast._STATE == {"warned": True, "warned_native": True}


def test_epsg_constants_follow_synaptic():
    assert _fast._epsg(synaptic.EPSG_TAU) == (1 / 0.21317, 0.18, 0.1)
    inv_peak, decay, rise = _fast._epsg(MYOGA)
    assert (decay, rise) == (0.3, 0.1)
    t = np.linspace(0.01, 3, 50)
    np.testing.assert_array_equal(
        inv_peak * (np.exp(-t / decay) - np.exp(-t / rise)),
        synaptic.epsg_unitary(t, MYOGA),
    )


# --- coincidence.spike_height's cache ----------------------------------------


@pytest.fixture
def fake_runs(monkeypatch):
    """Fakes _dispatch: threshold 40, height 30; counts the height runs."""
    runs = []

    def response(model, stim, start, stop, amp, *args, factor, **kw):
        if amp >= 40.0:
            return True, factor
        return False, factor * min(amp / 40.0, 0.999) ** 6

    def run_model(*args, **kw):
        runs.append(kw.get("mode"))
        return np.array([0.0, 20.0]), np.array([[0.0, 0.0], [0.0, 30.0]])

    monkeypatch.setattr(_dispatch, "response", response)
    monkeypatch.setattr(_dispatch, "run_model", run_model)
    monkeypatch.setattr(coincidence, "_HEIGHTS", {})
    return runs


def test_spike_height_cache_keys_the_resolved_mode(
    fake_ext, fake_runs, monkeypatch
):
    # unset, the mode resolves to "native" here (fake extension); once
    # MSOAXON_MODE=exact it resolves to "exact", so the cached native
    # height must not serve it
    assert coincidence.spike_height("multi", node=2) == 30.0
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    assert coincidence.spike_height("multi", node=2) == 30.0
    assert len(fake_runs) == 2
    # an explicit mode="exact" is the same configuration: a hit
    coincidence.spike_height("multi", node=2, model_kw={"mode": "exact"})
    assert len(fake_runs) == 2
    monkeypatch.delenv(_fast.MODE_ENV)
    coincidence.spike_height("multi", node=2, model_kw={"mode": "native"})
    assert len(fake_runs) == 2  # the first, unset-mode entry
    coincidence.spike_height("multi", node=2, model_kw={"mode": "fast"})
    assert len(fake_runs) == 3  # fast is not the default: a miss


# --- the native mode ----------------------------------------------------------

EPSG8 = ("EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1)
SPIKING = ("EPSGpair", 5, 5.1, 100, 3, "active-full", 20, -68.0, 1)


@pytest.fixture
def quiet_native(monkeypatch):
    """Native runs without the first-use warning."""
    monkeypatch.setitem(_fast._STATE, "warned_native", True)


def test_modes_add_native_between_fast_and_exact():
    assert _fast.MODES == ("fast", "native", "exact")


def test_native_is_the_default_fast_and_exact_stay(fake_ext, monkeypatch):
    assert _fast.effective(None) == "native"
    assert _fast.effective("fast") == "fast"
    assert _fast.effective("exact") == "exact"
    monkeypatch.setenv(_fast.MODE_ENV, "fast")
    assert _fast.effective(None) == "fast"
    assert _fast.resolve(None) == "fast"
    assert _fast.effective("native") == "native"  # explicit wins


def test_asking_for_native_without_the_extension_warns(monkeypatch, no_env):
    monkeypatch.setattr(_fast, "_ext", None)
    assert _fast.effective("native") == "exact"
    with pytest.warns(UserWarning, match='mode="native" asked for') as rec:
        t, _ = multi.mso_axon(*EPSG8, mode="native")
    assert rec[0].filename == __file__
    assert t[-1] == 8
    monkeypatch.setenv(_fast.MODE_ENV, "native")
    with pytest.warns(UserWarning, match="not installed"):
        assert _fast.resolve(None) == "exact"


def test_rust_modes_build_events_for_the_run_and_none_for_the_prefix(
    fake_ext,
):
    syn = synaptic.SynParams(
        t_end=8, events=(synaptic.Event(t=5.0, node=1, g=30.0),)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no fallback warning any more
        multi.mso_axon(
            "Events", 0, 0, 1.0, 3, "active-full", 8, -68.0, 1, syn=syn
        )
        multi.mso_axon(
            "Events",
            0,
            0,
            1.0,
            3,
            "active-full",
            8,
            -68.0,
            1,
            syn=syn,
            mode="fast",
        )
    assert fake_ext == [
        ("native", "Events"),
        ("native", "none"),
        ("fast", "Events"),
        ("fast", "none"),
    ]


def test_segment_without_a_segment_method_is_scipy():
    calls = []

    def rhs(t, x):
        calls.append(t)
        return -x

    t, y, stopped = _solve._segment(
        rhs, 0.0, 1.0, np.ones(2), 1e-8, 1e-8, 0.1, None
    )
    assert calls and not stopped and t[-1] == 1.0
    assert y.shape == (2, len(t))


def test_segment_hands_off_to_a_segment_method():
    class Fake:
        def segment(self, *args):
            self.args = args
            return "t", "y", True

    fake = Fake()
    out = _solve._segment(fake, 0.0, 1.0, "y0", 1e-8, 1e-9, 0.1, "jac", "ev")
    assert out == ("t", "y", True)
    assert fake.args == (0.0, 1.0, "y0", 1e-8, 1e-9, 0.1, "ev")


def test_events_carry_their_spec():
    assert _solve.spike_event(2, 20.0).spec == ("spike", 2, 20.0)
    dv = lambda t, x: x  # noqa: E731
    spec = _solve.settled_event(dv, 7.5, 2, 20.0).spec
    assert spec == ("settled", 7.5, 2, 20.0 * _solve.SETTLED_FRACTION)


def _native_rhs(stim="EPSGpair", amp=30.0, t_end=8.0):
    syn = synaptic.SynParams(t_end=t_end)
    s = _common.stimulus(stim, 5.0, 5.1, amp, t_end, syn)
    mem = multi.membrane(-68.0)
    active = (True,) * 6
    args = (stim, s, 1, active, mem, multi.LUMPED)
    return _fast._build(*args), _fast.Native(_fast._build(*args)), mem


@needs_ext
def test_native_calls_the_rust_rhs():
    fast, native, mem = _native_rhs()
    x = multi.initial_state(multi.LUMPED, mem, -60.0)
    np.testing.assert_array_equal(native(5.05, x), fast(5.05, x))


@needs_ext
def test_native_segment_contract():
    _, native, mem = _native_rhs()
    y0 = multi.initial_state(multi.LUMPED, mem, -68.0)
    t, y, stopped = native.segment(0.25, 5.0, y0, 1e-8, 1e-8, 0.8)
    assert not stopped
    assert t[0] == 0.25 and t[-1] == 5.0  # exactly
    assert np.all(np.diff(t) > 0) and np.all(np.diff(t) <= 0.8)
    assert y.shape == (7 * multi.LUMPED.n, len(t))  # solve_ivp's (n, n_t)
    np.testing.assert_array_equal(y[:, 0], y0)
    # the same integrator serves the next segment, from where this ended
    t2, y2, _ = native.segment(5.0, 8.0, y[:, -1], 1e-8, 1e-8, 0.8)
    assert t2[0] == 5.0 and t2[-1] == 8.0
    np.testing.assert_array_equal(y2[:, 0], y[:, -1])


@needs_ext
def test_native_rejects_an_event_it_cannot_run():
    _, native, mem = _native_rhs()
    y0 = multi.initial_state(multi.LUMPED, mem, -68.0)
    with pytest.raises(ValueError, match="spike_event"):
        native.segment(0.0, 1.0, y0, 1e-8, 1e-8, 0.8, [lambda t, x: 1.0])


@pytest.mark.parametrize(
    "mode",
    [
        "exact",
        pytest.param("fast", marks=needs_ext),
        pytest.param("native", marks=needs_ext),
    ],
)
def test_a_solver_failure_is_a_runtime_error_in_every_mode(mode, quiet_native):
    # an EPSG pair of 1e30 nS: scipy's step size underflows, diffsol's too
    with pytest.raises(RuntimeError, match=r"integration failed on \[5") as e:
        multi.mso_axon(*EPSG8[:3], 1e30, *EPSG8[4:], mode=mode)
    if mode == "native":
        assert "diffsol: " in str(e.value)  # diffsol's own message, kept
        assert isinstance(e.value.__cause__, ValueError)


@needs_ext
def test_native_argument_errors_stay_value_errors():
    _, native, mem = _native_rhs()
    y0 = multi.initial_state(multi.LUMPED, mem, -68.0)
    with pytest.raises(ValueError, match="need t0 < t1"):
        native.segment(1.0, 1.0, y0, 1e-8, 1e-8, 0.8)


@pytest.mark.parametrize(
    ("message", "raised"),
    [
        ("diffsol: Sparsity error", RuntimeError),
        ("the pattern indexes 3 of its 4 entries", ValueError),
    ],
)
def test_building_the_integrator_maps_diffsol_errors_too(message, raised):
    # the first segment builds the integrator; diffsol failing there is a
    # solver failure like one during the run, a check of ours is not
    class Rhs:
        def integrator(self):
            raise ValueError(message)

    with pytest.raises(raised, match=message) as e:
        _fast.Native(Rhs()).segment(0.0, 1.0, np.zeros(7), 1e-8, 1e-8, 0.8)
    if raised is RuntimeError:
        assert str(e.value).startswith("integration failed on [0.0, 1.0]")
        assert isinstance(e.value.__cause__, ValueError)


@needs_ext
def test_native_spike_stop_agrees_with_scipys_event(quiet_native):
    kw = {"stop_on_spike": 20.0}
    t_f, x_f = multi.mso_axon(*SPIKING, mode="fast", **kw)
    t_n, x_n = multi.mso_axon(*SPIKING, mode="native", **kw)
    assert t_f[-1] < 20 and t_n[-1] < 20
    assert abs(t_n[-1] - t_f[-1]) < 1e-6
    np.testing.assert_allclose(x_n[-1], x_f[-1], rtol=0, atol=1e-4)
    assert abs(x_n[-1, 2] - x_n[-1, 0] - 20.0) < 1e-6


@needs_ext
def test_native_settled_stop_ends_a_subthreshold_run(quiet_native):
    kw = {"stop_on_spike": 20.0, "stop_when_settled": True}
    t, x = multi.mso_axon(*EPSG8[:6], 40, *EPSG8[7:], mode="native", **kw)
    t_off = 5.1 + _solve.EPSG_TAIL * max(synaptic.EPSG_TAU)
    assert t_off <= t[-1] < 40
    assert x[-1, 2] - x[-1, 0] < 20.0 * _solve.SETTLED_FRACTION


@needs_ext
def test_native_max_step_bounds_every_step(quiet_native):
    t, _ = multi.mso_axon(
        *EPSG8[:6], 2, *EPSG8[7:], max_step=0.005, mode="native"
    )
    assert t[-1] == 2
    assert np.max(np.diff(t)) <= 0.005 + 1e-12


@needs_ext
def test_native_quiet_cache_key_differs_from_fast(
    fresh_quiet_cache, quiet_native, monkeypatch
):
    monkeypatch.setitem(_fast._STATE, "warned", True)
    for mode in ("fast", "native"):
        multi.mso_axon(*EPSG8, mode=mode)
    modes = sorted(k[0][-1] for k in _solve._QUIET_CACHE)
    assert modes == ["fast", "native"]


@needs_ext
def test_native_matches_exact_with_a_synaptic_input_off_the_soma(
    quiet_native,
):
    # the synaptic input drives its own compartment, off the soma
    kw = {"morph": D, "syn": synaptic.SynParams(t_end=15.0)}
    args = ("Synaptic", 5, 15, 0, 3, "active-full", 15, -68.0, LAT)
    t_e, x_e = multi.mso_axon(*args, mode="exact", **kw)
    t_n, x_n = multi.mso_axon(*args, mode="native", **kw)
    assert t_n[-1] == t_e[-1] == 15
    times = np.linspace(0, 15, 151)
    for col in (0, 2, LAT - 1):
        np.testing.assert_allclose(
            np.interp(times, t_n, x_n[:, col]),
            np.interp(times, t_e, x_e[:, col]),
            rtol=1e-4,
            atol=5e-3,
        )


@needs_ext
def test_native_drift_warning_fires_once_and_says_tolerance(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned_native", False)
    monkeypatch.setitem(_fast._STATE, "warned", True)
    monkeypatch.delenv(_fast.MODE_ENV, raising=False)
    with pytest.warns(msoaxon.DriftWarning) as record:
        for _ in range(2):
            msoaxon.mso_axon(*EPSG8, mode="native")
    drift = [w for w in record if issubclass(w.category, _fast.DriftWarning)]
    assert len(drift) == 1
    text = str(drift[0].message)
    assert 'mode="native"' in text and "within the solver" in text
    assert "not bit for bit" in text and "filterwarnings" in text
    assert "the default when the msoaxon-fast extension is installed" in text
    assert 'Pass mode="fast" (the same right-hand side) or' in text
    assert 'mode="exact" (the numpy reference)' in text
    assert drift[0].filename == __file__


@needs_ext
def test_a_bare_default_run_issues_the_native_warning(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned_native", False)
    monkeypatch.setitem(_fast._STATE, "warned", False)
    monkeypatch.delenv(_fast.MODE_ENV, raising=False)
    with pytest.warns(msoaxon.DriftWarning) as record:
        msoaxon.mso_axon(*EPSG8)  # no mode, no MSOAXON_MODE
    drift = [w for w in record if issubclass(w.category, _fast.DriftWarning)]
    assert len(drift) == 1
    assert "Rust BDF integrator" in str(drift[0].message)
    assert _fast._STATE == {"warned": False, "warned_native": True}
