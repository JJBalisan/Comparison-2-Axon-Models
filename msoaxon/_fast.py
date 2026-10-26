"""The optional Rust right-hand side of the 45-compartment model.

`msoaxon_fast` (rust/, the `fast` extra) ports multi._rhs to Rust,
expression for expression in numpy's floating-point order. multi._rhs stays
the reference: mso_axon(mode="exact") runs it, mode="fast" runs the port
under scipy's BDF, and mode="native" runs the port under a Rust BDF
integrator (diffsol), so no Python is left in the time loop (Native).

Mode, per mso_axon call: an explicit mode= wins, then the MSOAXON_MODE
environment variable, then the default, which is "native" when the
extension imports, else "exact" (silently when the extra is not installed;
with a warning when it is installed but fails to load). Every configuration
is ported. Asking for "fast" or "native" (by argument or environment) when
it can't run warns and runs "exact".

The first fast run in a process checks the port against multi._rhs on a few
fixed states and reports the result in a DriftWarning: bit-identical on the
machine the port was checked on (darwin-arm64), it may differ in the last
bits on other hardware, where numpy and Rust can use different exp/pow.
The first native run issues its own DriftWarning: its steps are diffsol's,
not scipy's, so it agrees with the reference within the solver tolerance,
not bit for bit.
"""

import importlib
import os
import warnings
from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np

from msoaxon import _common, constants, synaptic

if TYPE_CHECKING:  # multi imports this module
    from msoaxon import multi


def _import_ext(
    import_module: Callable[[str], object] = importlib.import_module,
) -> object | None:
    """Imports msoaxon_fast, or returns None.

    Not installed (ModuleNotFoundError for the package itself) is the
    numpy-only configuration and stays silent. Anything else, an
    extension that is present but will not load (a wheel built for
    another Python or architecture, a missing symbol), is a broken
    install: it warns with the loader's message, because falling back
    silently would hide a 2.7x slowdown and last-bit differences behind
    results that look fine.

    Args:
        import_module: The importer (a parameter for the tests).

    Returns:
        The module, or None.
    """
    try:
        return import_module("msoaxon_fast")
    except ImportError as e:
        if isinstance(e, ModuleNotFoundError) and e.name == "msoaxon_fast":
            return None
        warnings.warn(
            "msoaxon-fast is installed but failed to import "
            f'({e}); running mode="exact"',
            stacklevel=2,
        )
        return None


_ext = _import_ext()

MODES = ("fast", "native", "exact")
MODE_ENV = "MSOAXON_MODE"

# per process: whether the fast and the native DriftWarning have been issued
# (pool workers inherit the parent's flags, see inherit())
_STATE = {"warned": False, "warned_native": False}

Rhs = Callable[[float, np.ndarray], np.ndarray]


class DriftWarning(UserWarning):
    """The Rust right-hand side is in use (mode="fast" or "native").

    Fast results can differ from the numpy reference (mode="exact") in the
    last bits on hardware where the two use different maths libraries;
    native results, from another integrator, within the solver tolerance.
    Silence it with warnings.filterwarnings("ignore",
    category=msoaxon.DriftWarning).
    """


def available() -> bool:
    """Returns whether the msoaxon_fast extension imported."""
    return _ext is not None


def requested(mode: str | None = None) -> str | None:
    """Returns the mode asked for: `mode`, else MSOAXON_MODE, else None.

    Args:
        mode: The mode passed to mso_axon, or None.

    Returns:
        "fast", "native", "exact" or None (no preference).

    Raises:
        ValueError: For any other value.
    """
    if mode is None:
        mode = os.environ.get(MODE_ENV) or None
    if mode is not None and mode not in MODES:
        raise ValueError(f"mode must be one of {MODES} or None, got {mode!r}")
    return mode


def effective(mode: str | None) -> str:
    """Returns the mode a run will use, without warning.

    For cache keys: what mso_axon will run, not what was asked for, so a
    result cached under one MSOAXON_MODE never serves another.

    Args:
        mode: The mode passed to mso_axon, or None.

    Returns:
        "fast" when asked for and the extension imports; "native" when
        nothing or "native" is asked for and it imports; else "exact".
    """
    want = requested(mode)
    if want == "exact" or _ext is None:
        return "exact"
    return "fast" if want == "fast" else "native"


def resolve(mode: str | None) -> str:
    """Returns the mode one run uses, warning if "fast"/"native" can't run.

    Args:
        mode: The mode passed to mso_axon, or None.

    Returns:
        "fast", "native" or "exact".
    """
    got = effective(mode)
    want = requested(mode)
    if got == "exact" and want in ("fast", "native"):
        warnings.warn(
            f'mode="{want}" asked for, but the msoaxon-fast extension is '
            'not installed (uv sync --extra fast); running mode="exact"',
            stacklevel=3,
        )
    return got


def _epsg(tau: tuple[float, float]) -> tuple[float, float, float]:
    """Returns synaptic.epsg_unitary's (1 / peak, decay, rise) for `tau`."""
    if tau == synaptic.EPSG_TAU:
        return 1 / 0.21317, 0.18, 0.1
    rise, decay = tau
    t_peak = np.log(decay / rise) * rise * decay / (decay - rise)
    peak = np.exp(-t_peak / decay) - np.exp(-t_peak / rise)
    return float(1 / peak), float(decay), float(rise)


def _f64(a: np.ndarray) -> np.ndarray:
    return np.asarray(a, dtype=np.float64)


# The numbers inside constants.py's gating functions are literals there, so
# they are repeated here once, grouped by the function each belongs to; the
# Rust side (rust/src/lib.rs, `Gating`) has none of its own. A name is the
# function's and the number's role in it, the signs are as the function
# writes them. Changing a literal in constants.py means changing it here;
# the bitwise tests (tests/test_fast.py) fail if the two disagree. The
# values the module does name (the reversal potentials in multi, the Q10
# factors in constants) are read when a run is built, see _gating().
_GATING = {
    # na_open: g * m**4 * (0.993 h + 0.007)
    "na_open_h": 0.993,
    "na_open_floor": 0.007,
    # minf: 1 / (1 + exp((V + 46) / -11))
    "minf_v": 46.0,
    "minf_k": -11.0,
    # taum: (0.141 + -0.0826 / (1 + exp((-20.5 - V) / 10.8))) / 3
    "taum_a": 0.141,
    "taum_b": -0.0826,
    "taum_v": -20.5,
    "taum_k": 10.8,
    "taum_div": 3.0,
    # hinf: 1 / (1 + exp((V + vhalf) / 7.77)); vhalf is Membrane.na_vhalf
    "hinf_k": 7.77,
    # tauh: (4 + -3.74 / (1 + exp((-40.6 - V) / 5.05))) / 3
    "tauh_a": 4.0,
    "tauh_b": -3.74,
    "tauh_v": -40.6,
    "tauh_k": 5.05,
    "tauh_div": 3.0,
    # pinf: 1 / (1 + exp(-(V + 23) / 6))
    "pinf_v": 23.0,
    "pinf_k": 6.0,
    # taup: _P_TEMP * (100 / (4 exp((V + 60) / 32) + 5 exp(-(V + 60) / 22))
    # + 5); both 60s are one value here
    "taup_num": 100.0,
    "taup_a": 4.0,
    "taup_ka": 32.0,
    "taup_b": 5.0,
    "taup_kb": 22.0,
    "taup_v": 60.0,
    "taup_c": 5.0,
    # winf: 1 / (1 + exp((V + 57.34) / -11.7))
    "winf_v": 57.34,
    "winf_k": -11.7,
    # tauw: 21.5 / (6 exp((V + 60) / 7) + 24 exp(-(V + 60) / 50.6)) + 0.35;
    # both 60s are one value here
    "tauw_num": 21.5,
    "tauw_a": 6.0,
    "tauw_ka": 7.0,
    "tauw_b": 24.0,
    "tauw_kb": 50.6,
    "tauw_v": 60.0,
    "tauw_c": 0.35,
    # zinf: (1 - 0.27) / (1 + exp((V + 67) / 6.16)) + 0.27
    "zinf_floor": 0.27,
    "zinf_v": 67.0,
    "zinf_k": 6.16,
    # tauz: 170 / (5 exp((V + 60) / 10) + exp(-(V + 70) / 8)) + 10.7
    "tauz_num": 170.0,
    "tauz_a": 5.0,
    "tauz_va": 60.0,
    "tauz_ka": 10.0,
    "tauz_vb": 70.0,
    "tauz_kb": 8.0,
    "tauz_c": 10.7,
    # ainf: 1 / (1 + exp(0.1 (V + 80.4)))
    "ainf_k": 0.1,
    "ainf_v": 80.4,
    # taua: _A_TEMP * (79 + 417 exp(-(V + 61.5)**2 / 800))
    "taua_c": 79.0,
    "taua_a": 417.0,
    "taua_v": 61.5,
    "taua_k": 800.0,
}


def _gating() -> dict[str, float]:
    """Returns every constant the Rust right-hand side takes from Python.

    Read when a run is built, not at import, so a change to multi's
    reversal potentials or constants' Q10 factors reaches the Rust side
    exactly as it reaches multi._rhs.

    Returns:
        name -> value, the names of `Gating` in rust/src/lib.rs.
    """
    multi = importlib.import_module("msoaxon.multi")
    return {
        "v_na": float(multi.V_NA),
        "v_k": float(multi.V_K),
        "v_h": float(multi.V_H),
        "taup_temp": float(constants._P_TEMP),
        "taua_temp": float(constants._A_TEMP),
        **_GATING,
    }


def build(
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    active: tuple[bool, ...],
    mem: 'multi.Membrane',
    morph: 'multi.Morphology',
    input_node2: int | None = None,
) -> Rhs:
    """Builds the Rust right-hand side for one run.

    The arguments are multi._rhs's; every per-run constant is copied into
    the Rust object once, here.

    Args:
        stim_type: The stimulus name, or "none".
        s: The stimulus settings (_common.stimulus); for "Events" its
            event table (s.events) goes to Rust whole.
        input_node: Where the stimulus enters, 1-indexed.
        active: Which gates (m, h, p, w, z, a) evolve.
        mem: Channel overrides (multi.membrane()).
        morph: The morphology.
        input_node2: The second input site, for "EPSGbilateral".

    Returns:
        f(t, x) -> dx/dt, bit for bit multi._rhs on darwin-arm64.
    """
    if not _STATE["warned"]:
        _STATE["warned"] = True
        _warn_drift()
    return _build(stim_type, s, input_node, active, mem, morph, input_node2)


def build_native(
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    active: tuple[bool, ...],
    mem: 'multi.Membrane',
    morph: 'multi.Morphology',
    input_node2: int | None = None,
) -> 'Native':
    """Builds the Rust right-hand side and integrator for one run.

    As build(), for mode="native"; the first call in a process issues the
    native DriftWarning.

    Returns:
        A Native: callable as f(t, x), and solving whole segments in Rust.
    """
    if not _STATE["warned_native"]:
        _STATE["warned_native"] = True
        _warn_native()
    return Native(
        _build(stim_type, s, input_node, active, mem, morph, input_node2)
    )


class Native:
    """A Rust right-hand side whose segments its Rust integrator solves.

    _solve._segment hands a segment to segment() whenever the right-hand
    side has one, instead of calling scipy; called as f(t, x) it is the
    plain Rust right-hand side (the settled stop's dv, direct callers).
    The integrator (diffsol's BDF with the run's Jacobian pattern) is
    built on the first segment and reused for every later one.
    """

    def __init__(self, rhs: Rhs) -> None:
        """Wraps `rhs`, an msoaxon_fast.MultiRhs."""
        self._rhs = rhs
        self._integrator = None

    def __call__(self, t: float, x: np.ndarray) -> np.ndarray:
        """Returns dx/dt at (t, x), the Rust right-hand side."""
        return self._rhs(t, x)

    def segment(
        self,
        a: float,
        b: float,
        y: np.ndarray,
        rtol: float,
        atol: float,
        max_step: float,
        events: object = None,
    ) -> tuple[np.ndarray, np.ndarray, bool]:
        """Solves [a, b] from y in Rust; _solve._segment's contract.

        Args:
            a: Segment start [ms].
            b: Segment end [ms].
            y: The state at a.
            rtol: Relative tolerance.
            atol: Absolute tolerance.
            max_step: Largest step [ms].
            events: None, or one or a list of _solve.spike_event /
                settled_event closures (read through their `spec`).

        Returns:
            (t, y, stopped): t from a (y included) to b or the event, y
            shaped (n_states, len(t)) like solve_ivp's, stopped whether
            an event fired.

        Raises:
            ValueError: For an event without a `spec`.
            RuntimeError: If diffsol fails on the segment (or building
                its problem, on the first), as _solve._segment raises
                for scipy's BDF; the message keeps diffsol's.
        """
        spike = settled = None
        if events is not None and not isinstance(events, list | tuple):
            events = [events]
        for event in events or ():
            spec = getattr(event, "spec", None)
            if spec is None:
                raise ValueError(
                    'mode="native" runs only _solve.spike_event and '
                    f"settled_event, got {event!r}"
                )
            if spec[0] == "spike":
                spike = (int(spec[1]), float(spec[2]))
            else:
                settled = (float(spec[1]), int(spec[2]), float(spec[3]))
        try:
            if self._integrator is None:
                self._integrator = self._rhs.integrator()
            t, ys, stopped = self._integrator.run(
                np.ascontiguousarray(y, dtype=np.float64),
                float(a),
                float(b),
                float(rtol),
                float(atol),
                float(max_step),
                spike,
                settled,
            )
        except ValueError as e:
            # diffsol's own failures, building the problem or solving it
            # (rust/src/integrate.rs prefixes them); the argument checks
            # stay ValueError
            if not str(e).startswith("diffsol:"):
                raise
            raise RuntimeError(f"integration failed on [{a}, {b}]: {e}") from e
        return t, ys.T, stopped


def _build(
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    active: tuple[bool, ...],
    mem: 'multi.Membrane',
    morph: 'multi.Morphology',
    input_node2: int | None = None,
) -> Rhs:
    """build() without the first-use DriftWarning."""
    stim = {
        "stim_type": stim_type,
        "start": float(s.start),
        "stop": float(s.stop),
        "amp": float(s.I),
        "k": input_node - 1,
        "k2": (input_node2 or 1) - 1,
        "epsg": _epsg(tuple(s.epsg_tau)),
    }
    if stim_type == "sine":
        stim["f"] = float(s.f)
    if stim_type in ("Synaptic", "SynapticPair"):
        stim.update(
            t_syn=_f64(s.t_syn),
            g_syn=_f64(s.g_syn),
            vsyn=float(s.VsynE),
            diff=float(s.diff),
        )
    if stim_type == "Events":
        # "amp" is the run's I, which scales every event's peak
        stim.update(
            {
                f"ev_{name}": (
                    np.asarray(a, dtype=np.int64) if name == "k" else _f64(a)
                )
                for name, a in s.events.items()
            }
        )
    return _ext.MultiRhs(
        g_na=_f64(mem.g_na),
        g_kht=_f64(morph.g_kht),
        g_klt=_f64(mem.g_klt),
        g_h=_f64(morph.g_h),
        g_lk=_f64(morph.g_lk),
        vlk=_f64(mem.vlk),
        na_vhalf=_f64(mem.na_vhalf),
        cap=_f64(morph.cap),
        sa=_f64(morph.sa),
        g_ax=_f64(morph.g_ax),
        parent=np.asarray(morph.parent, dtype=np.int64),
        child=np.asarray(morph.child, dtype=np.int64),
        active=tuple(bool(a) for a in active),
        gating=_gating(),
        **stim,
    )


def self_check(n_states: int = 8) -> tuple[int, int, float]:
    """Compares the Rust and numpy right-hand sides on fixed states.

    Three configurations (the chain with an EPSG pair, active-full; a
    rebalanced membrane with a step; the dendrites with bilateral EPSGs),
    each on `n_states` seeded random states while the input is on.

    Args:
        n_states: States per configuration.

    Returns:
        (states that differ in any bit, states compared, largest relative
        difference).
    """
    multi = importlib.import_module("msoaxon.multi")
    common = importlib.import_module("msoaxon._common")
    rng = np.random.default_rng(20261012)
    syn = synaptic.SynParams()
    dend = multi.with_dendrites()
    configs = (
        ("EPSGpair", multi.membrane(-68.0), multi.LUMPED, 1, None),
        (
            "step",
            multi.membrane(-68.0, soma_klt_scale=3, rebalance_rest=True),
            multi.LUMPED,
            1,
            None,
        ),
        (
            "EPSGbilateral",
            multi.membrane(-68.0, morph=dend),
            dend,
            int(dend.lateral[2]),
            int(dend.medial[2]),
        ),
    )
    active = (True,) * 6
    differ, total, worst = 0, 0, 0.0
    for stim, mem, morph, k, k2 in configs:
        s = common.stimulus(stim, 5.0, 5.2, 60.0, 20.0, syn)
        fast = _build(stim, s, k, active, mem, morph, k2)
        for _ in range(n_states):
            t = float(rng.uniform(5.0, 6.0))
            x = np.concatenate(
                [rng.uniform(-90.0, 30.0, morph.n), rng.random(6 * morph.n)]
            )
            want = multi._rhs(t, x, stim, s, k, active, mem, morph, k2)
            got = fast(t, x)
            total += 1
            if not np.array_equal(got, want):
                differ += 1
                scale = np.maximum(np.abs(want), np.finfo(float).tiny)
                worst = max(worst, float(np.max(np.abs(got - want) / scale)))
    return differ, total, worst


def _warn_drift() -> None:
    """Runs self_check() and issues the DriftWarning with its result."""
    differ, total, worst = self_check()
    if differ == 0:
        verdict = (
            f"on this machine it is bit-identical to the numpy reference "
            f"(self-check: {total} states)"
        )
    else:
        verdict = (
            f"on this machine it differs from the numpy reference in "
            f"{differ} of {total} self-check states (largest relative "
            f"difference {worst:.1e})"
        )
    warnings.warn(
        f'msoaxon is running its Rust right-hand side (mode="fast", under '
        f"scipy's BDF; the default is mode=\"native\"); "
        f"{verdict}. Results may differ from the numpy reference in the "
        f'last bits on other hardware. Pass mode="exact" (or set '
        f"{MODE_ENV}=exact) to run the reference; silence this with "
        f'warnings.filterwarnings("ignore", category=msoaxon.DriftWarning).',
        DriftWarning,
        stacklevel=4,
    )


def _warn_native() -> None:
    """Issues the native DriftWarning."""
    warnings.warn(
        'msoaxon is running its Rust BDF integrator (mode="native", the '
        "default when the msoaxon-fast extension is installed): diffsol's "
        "BDF over the Rust right-hand side, in place of scipy's. Results "
        "agree with the numpy reference within the solver tolerance, not "
        "bit for bit: the steps, and so the last digits and event times, "
        'are its own. Pass mode="fast" (the same right-hand side) or '
        f'mode="exact" (the numpy reference), or set {MODE_ENV}, for '
        "scipy's BDF; silence this with "
        'warnings.filterwarnings("ignore", category=msoaxon.DriftWarning).',
        DriftWarning,
        stacklevel=4,
    )


def parent_state() -> tuple[str | None, tuple[bool, bool]]:
    """Returns what pool workers inherit: MSOAXON_MODE and the warned flags.

    Returns:
        (MSOAXON_MODE or None, (fast warned, native warned)).
    """
    return os.environ.get(MODE_ENV), (
        _STATE["warned"],
        _STATE["warned_native"],
    )


def inherit(env_mode: str | None, warned: bool | tuple[bool, bool]) -> None:
    """Pool-worker initializer: adopts the parent's parent_state().

    A forkserver worker forks from a server started with the environment
    of the pool's first use, so a later MSOAXON_MODE would not reach it
    otherwise. A parent that already issued the DriftWarning spares its
    workers from repeating it.

    Args:
        env_mode: The parent's MSOAXON_MODE (None if unset).
        warned: Whether the parent issued the fast and the native
            DriftWarning, as parent_state() pairs them; a bool is the
            fast one alone.
    """
    if env_mode is None:
        os.environ.pop(MODE_ENV, None)
    else:
        os.environ[MODE_ENV] = env_mode
    fast, native = warned if isinstance(warned, tuple) else (warned, False)
    _STATE["warned"] = _STATE["warned"] or fast
    _STATE["warned_native"] = _STATE["warned_native"] or native
