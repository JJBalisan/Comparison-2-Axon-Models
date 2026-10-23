"""The optional Rust right-hand side of the 45-compartment model.

`msoaxon_fast` (rust/, the `fast` extra) ports multi._rhs to Rust,
expression for expression in numpy's floating-point order. multi._rhs stays
the reference: mso_axon(mode="exact") runs it, mode="fast" runs the port.

Mode, per mso_axon call: an explicit mode= wins, then the MSOAXON_MODE
environment variable, then the default, which is "fast" when the extension
imports, else "exact" (silently). Every configuration is ported. Asking
for "fast" (by argument or environment) when it can't run warns and runs
"exact".

The first fast run in a process checks the port against multi._rhs on a few
fixed states and reports the result in a DriftWarning: bit-identical on the
machine the port was checked on (darwin-arm64), it may differ in the last
bits on other hardware, where numpy and Rust can use different exp/pow.
"""

import importlib
import os
import types
import warnings
from collections.abc import Callable

import numpy as np

from msoaxon import constants, synaptic

try:
    import msoaxon_fast as _ext
except ImportError:  # the fast extra isn't installed: numpy only
    _ext = None

MODES = ("fast", "exact")
MODE_ENV = "MSOAXON_MODE"

# per process: whether the DriftWarning has been issued (pool workers inherit
# the parent's flag, see inherit())
_STATE = {"warned": False}

Rhs = Callable[[float, np.ndarray], np.ndarray]


class DriftWarning(UserWarning):
    """The Rust right-hand side is in use (mode="fast").

    Its results can differ from the numpy reference (mode="exact") in the
    last bits on hardware where the two use different maths libraries.
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
        "fast", "exact" or None (no preference).

    Raises:
        ValueError: For any other value.
    """
    if mode is None:
        mode = os.environ.get(MODE_ENV) or None
    if mode is not None and mode not in MODES:
        raise ValueError(f"mode must be one of {MODES} or None, got {mode!r}")
    return mode


def effective(mode: str | None, covered: bool = True) -> str:
    """Returns the mode a run will use, without warning.

    For cache keys: what mso_axon will run, not what was asked for, so a
    result cached under one MSOAXON_MODE never serves another.

    Args:
        mode: The mode passed to mso_axon, or None.
        covered: Whether the configuration is ported (always, today).

    Returns:
        "fast" or "exact".
    """
    want = requested(mode)
    return (
        "fast" if want != "exact" and _ext is not None and covered else "exact"
    )


def resolve(mode: str | None, covered: bool = True) -> str:
    """Returns the mode one run uses, warning if "fast" can't run.

    Args:
        mode: The mode passed to mso_axon, or None.
        covered: Whether the configuration is ported (always, today).

    Returns:
        "fast" or "exact".
    """
    if effective(mode, covered) == "fast":
        return "fast"
    if requested(mode) == "fast":
        why = (
            "the msoaxon-fast extension is not installed (uv sync --extra fast)"
            if _ext is None
            else "the Rust port does not cover this configuration"
        )
        warnings.warn(
            f'mode="fast" asked for, but {why}; running mode="exact"',
            stacklevel=3,
        )
    return "exact"


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


def build(
    stim_type: str,
    s: types.SimpleNamespace,
    input_node: int,
    active: tuple[bool, ...],
    mem: types.SimpleNamespace,
    morph: types.SimpleNamespace,
    input_node2: int | None = None,
) -> Rhs:
    """Builds the Rust right-hand side for one run.

    The arguments are multi._rhs's; every per-run constant is copied into
    the Rust object once, here.

    Args:
        stim_type: The stimulus name, or "none".
        s: The stimulus settings (_common.stimulus).
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


def _build(
    stim_type: str,
    s: types.SimpleNamespace,
    input_node: int,
    active: tuple[bool, ...],
    mem: types.SimpleNamespace,
    morph: types.SimpleNamespace,
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
        p_temp=constants._P_TEMP,
        a_temp=constants._A_TEMP,
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
        f'msoaxon is running its Rust right-hand side (mode="fast"); '
        f"{verdict}. Results may differ from the numpy reference in the "
        f'last bits on other hardware. Pass mode="exact" (or set '
        f"{MODE_ENV}=exact) to run the reference; silence this with "
        f'warnings.filterwarnings("ignore", category=msoaxon.DriftWarning).',
        DriftWarning,
        stacklevel=4,
    )


def parent_state() -> tuple[str | None, bool]:
    """Returns what pool workers inherit: MSOAXON_MODE and the warned flag."""
    return os.environ.get(MODE_ENV), _STATE["warned"]


def inherit(env_mode: str | None, warned: bool) -> None:
    """Pool-worker initializer: adopts the parent's parent_state().

    A forkserver worker forks from a server started with the environment
    of the pool's first use, so a later MSOAXON_MODE would not reach it
    otherwise. A parent that already issued the DriftWarning spares its
    workers from repeating it.

    Args:
        env_mode: The parent's MSOAXON_MODE (None if unset).
        warned: Whether the parent issued the DriftWarning.
    """
    if env_mode is None:
        os.environ.pop(MODE_ENV, None)
    else:
        os.environ[MODE_ENV] = env_mode
    _STATE["warned"] = _STATE["warned"] or warned
