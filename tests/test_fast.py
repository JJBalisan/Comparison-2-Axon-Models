"""The Rust right-hand side (mode="fast") against multi._rhs, and the dispatch.

Equivalence: the Rust port must match multi._rhs to rtol 1e-12 on every
configuration it covers, and bit for bit on darwin-arm64, where it was
written and where numpy and Rust call the same libm. Those tests need the
`fast` extra and skip without it; the dispatch tests run either way (they
fake or hide the extension).
"""

import multiprocessing
import platform
import sys
import warnings

import numpy as np
import pytest

import msoaxon
from msoaxon import _common, _fast, _parallel, _solve, multi, synaptic

BITWISE = sys.platform == "darwin" and platform.machine() == "arm64"
needs_ext = pytest.mark.skipif(
    not _fast.available(), reason="needs the fast extra (msoaxon_fast)"
)

D = multi.with_dendrites()
D_RA = multi.with_dendrites(conserve_totals=False, dendrite_ra=200.0)
LAT, MED = int(D.lateral[2]), int(D.medial[2])
DTX = {"soma_klt_scale": 0, "ais_klt_scale": 0, "rebalance_rest": True}
MYOGA = (0.1, 0.3)

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
def test_rust_rhs_rejects_a_wrong_state_length():
    args, _ = setup("EPSGpair")
    with pytest.raises(ValueError, match="315"):
        _fast._build(*args)(5.0, np.zeros(10))


@needs_ext
@pytest.mark.skipif(not BITWISE, reason="bit-identical on darwin-arm64 only")
@pytest.mark.parametrize("name", ["EPSGpair", "Synaptic", "dendrites/mem"])
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
def test_self_check_passes_and_is_cheap():
    import time  # noqa: PLC0415

    best = np.inf
    for _ in range(3):
        t0 = time.perf_counter()
        differ, total, worst = _fast.self_check()
        best = min(best, time.perf_counter() - t0)
    assert total == 24
    assert best < 0.05
    if BITWISE:
        assert (differ, worst) == (0, 0.0)
    else:
        # elementwise, so cancellation in (x_inf - x) / tau shows; GitHub's
        # ubuntu x86_64 runner measured 2.7e-13 (every state differed)
        assert worst < 1e-10


# --- the drift warning -------------------------------------------------------


@needs_ext
def test_drift_warning_fires_once_per_process(monkeypatch):
    monkeypatch.setitem(_fast._STATE, "warned", False)
    monkeypatch.delenv(_fast.MODE_ENV, raising=False)
    with pytest.warns(msoaxon.DriftWarning) as record:
        for _ in range(2):
            msoaxon.mso_axon(
                "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1
            )
    drift = [w for w in record if issubclass(w.category, _fast.DriftWarning)]
    assert len(drift) == 1
    text = str(drift[0].message)
    assert 'mode="exact"' in text and "filterwarnings" in text
    assert "other hardware" in text
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
    """Pretends the extension is installed; its RHS is multi._rhs."""
    built = []

    def build(stim_type, s, input_node, active, mem, morph, input_node2=None):
        built.append(stim_type)
        return lambda t, x: multi._rhs(
            t, x, stim_type, s, input_node, active, mem, morph, input_node2
        )

    monkeypatch.setattr(_fast, "_ext", object())
    monkeypatch.setattr(_fast, "build", build)
    return built


def test_fallback_without_the_extension_is_silent(monkeypatch, no_env):
    monkeypatch.setattr(_fast, "_ext", None)
    assert not _fast.available()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _fast.resolve(None, True) == "exact"
        assert _fast.resolve("exact", True) == "exact"
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
        assert _fast.resolve(None, True) == "exact"


def test_default_is_fast_when_installed_and_covered(fake_ext):
    assert _fast.resolve(None, True) == "fast"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _fast.resolve(None, False) == "exact"  # silently
    with pytest.warns(UserWarning, match="does not cover"):
        assert _fast.resolve("fast", False) == "exact"
    multi.mso_axon("EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1)
    assert fake_ext == ["EPSGpair", "none"]  # run and quiet prefix


def test_env_overrides_the_default_not_the_argument(fake_ext, monkeypatch):
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    assert _fast.resolve(None, True) == "exact"
    assert _fast.resolve("fast", True) == "fast"  # explicit wins
    multi.mso_axon("EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1)
    assert fake_ext == []
    monkeypatch.setenv(_fast.MODE_ENV, "")  # empty means unset
    assert _fast.resolve(None, True) == "fast"


def test_unknown_modes_are_rejected(monkeypatch, no_env):
    with pytest.raises(ValueError, match="mode must be"):
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode="rust"
        )
    monkeypatch.setenv(_fast.MODE_ENV, "turbo")
    with pytest.raises(ValueError, match="turbo"):
        _fast.resolve(None, True)


def test_coverage_excludes_only_degenerate_ramps():
    for stim in multi.MSO_STIM_TYPES:
        assert _fast.covers(stim, 5.0, 5.5)
    assert not _fast.covers("ramp", 5.0, 5.0)
    assert not _fast.covers("ramp2", 5, 5)
    assert _fast.covers("EPSGpair", 5.0, 5.0)


def test_mode_is_part_of_the_quiet_cache_key(fake_ext):
    _solve._QUIET_CACHE.clear()
    for mode in ("exact", "fast"):
        multi.mso_axon(
            "EPSGpair", 5, 5.1, 30, 3, "active-full", 8, -68.0, 1, mode=mode
        )
    modes = sorted(k[0][-1] for k in _solve._QUIET_CACHE)
    _solve._QUIET_CACHE.clear()
    assert modes == ["exact", "fast"]


def test_run_model_passes_mode_to_multi_only(fake_ext):
    from msoaxon import _dispatch  # noqa: PLC0415

    args = ("EPSGpair", 5, 5.1, 30, 3, 8, -68.0, 1)
    _dispatch.run_model("multi", *args, mode="fast")
    assert fake_ext == ["EPSGpair", "none"]
    _dispatch.run_model("multi", *args, mode="exact")
    _dispatch.run_model("two", *args, mode="fast")  # no fast two_cpt: ignored
    assert fake_ext == ["EPSGpair", "none"]


def test_pool_workers_inherit_the_mode_env(monkeypatch, no_env):
    # a forkserver worker forks from a server that keeps the environment it
    # started with; the pool initializer hands it the parent's
    ctx = multiprocessing.get_context("forkserver")
    monkeypatch.setattr(_parallel, "pool_context", lambda: ctx)
    monkeypatch.setitem(_fast._STATE, "warned", False)
    with _parallel.process_pool(1) as pool:
        assert list(pool.map(_fast.requested, [None])) == [None]
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    monkeypatch.setitem(_fast._STATE, "warned", True)
    with _parallel.process_pool(1) as pool:
        state = pool.submit(_fast.parent_state).result()
    assert state == ("exact", True)


def test_inherit_unsets_a_stale_env(monkeypatch):
    monkeypatch.setenv(_fast.MODE_ENV, "exact")
    monkeypatch.setitem(_fast._STATE, "warned", False)
    _fast.inherit(None, False)
    assert _fast.requested() is None and not _fast._STATE["warned"]
    _fast.inherit("fast", True)
    assert _fast.requested() == "fast" and _fast._STATE["warned"]


def test_epsg_constants_follow_synaptic():
    assert _fast._epsg(synaptic.EPSG_TAU) == (1 / 0.21317, 0.18, 0.1)
    inv_peak, decay, rise = _fast._epsg(MYOGA)
    assert (decay, rise) == (0.3, 0.1)
    t = np.linspace(0.01, 3, 50)
    np.testing.assert_array_equal(
        inv_peak * (np.exp(-t / decay) - np.exp(-t / rise)),
        synaptic.epsg_unitary(t, MYOGA),
    )
