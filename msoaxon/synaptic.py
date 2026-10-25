"""Synaptic input conductance (port of Synaptic.m)."""

import dataclasses
import functools

import numpy as np

EPSG_TAU = (
    0.1,
    0.18,
)  # (rise, decay) [ms], the Franken et al 2015 waveform used by the MATLAB

# The unitary EPSG's peak for EPSG_TAU, from the MATLAB code. Kept as the
# literal (fast mode's _fast._epsg and rust/src/lib.rs hold their own copies).
EPSG_NORM = 0.21317
# epsg_unitary clamps t to at least this [ms] inside its exponentials (see
# there); below onset the waveform is masked to zero whatever they hold.
EPSG_T_FLOOR = -1.0


@dataclasses.dataclass(frozen=True)
class Event:
    """One synaptic conductance event of the "Events" stimulus.

    The conductance is g times the unitary EPSG waveform (epsg_unitary,
    peak 1) from t, driving the compartment's voltage towards e_rev; an
    event never ends (it decays). Events add; any number may share a
    compartment or an onset.

    Attributes:
        t: Onset [ms], at or after 0 and before the run's t_end.
        node: The compartment it enters, 1-indexed. In two_cpt, 1 is the
            soma and anything else the axon compartment.
        g: Peak conductance [nS]; mso_axon's and two_cpt's `I` scales
            every event's g, so pass I=1.0 to use the peaks as given.
        tau: The (rise, decay) time constants [ms], 0 < rise < decay.
        e_rev: Reversal potential [mV]: 0 for an EPSG, below rest for a
            glycinergic IPSG.
    """

    t: float
    node: int = 1
    g: float = 1.0
    tau: tuple[float, float] = EPSG_TAU
    e_rev: float = 0.0


def check_events(
    events: tuple[Event, ...], n_max: int | None, t_end: float
) -> None:
    """Rejects an "Events" list the right-hand sides could not evaluate.

    Args:
        events: The SynParams.events tuple.
        n_max: The number of compartments an event may enter, or None
            to skip that check.
        t_end: End of the run [ms].

    Raises:
        ValueError: If there is no event, a node is out of range, an
            onset is before 0 or not before t_end, rise is not below
            decay (epsg_inv_peak divides by their difference) or either
            is not positive, or g or e_rev is not finite.
    """
    if not events:
        raise ValueError('"Events" needs at least one Event in syn.events')
    for i, e in enumerate(events):
        rise, decay = e.tau
        if n_max is not None and not 1 <= e.node <= n_max:
            raise ValueError(
                f"events[{i}].node must be 1..{n_max} (1-indexed), got {e.node}"
            )
        if not 0 <= e.t < t_end:
            raise ValueError(
                f"events[{i}].t must be in [0, t_end={t_end}), got {e.t}"
            )
        if not 0 < rise < decay:
            raise ValueError(
                f"events[{i}].tau needs 0 < rise < decay, got {e.tau}"
            )
        if not (np.isfinite(e.g) and np.isfinite(e.e_rev)):
            raise ValueError(
                f"events[{i}] needs finite g and e_rev, got {e.g}, {e.e_rev}"
            )


def event_arrays(events: tuple[Event, ...]) -> dict[str, np.ndarray]:
    """Lays the events out as read-only arrays for the right-hand sides.

    Args:
        events: The SynParams.events tuple (see check_events).

    Returns:
        dict with t, k (0-based compartment), g, rise, decay, inv_peak
        (epsg_inv_peak per event) and e_rev, one entry per event.
    """
    out = {
        "t": np.array([e.t for e in events], dtype=float),
        "k": np.array([e.node - 1 for e in events], dtype=np.intp),
        "g": np.array([e.g for e in events], dtype=float),
        "rise": np.array([e.tau[0] for e in events], dtype=float),
        "decay": np.array([e.tau[1] for e in events], dtype=float),
        "inv_peak": np.array(
            [epsg_inv_peak(tuple(e.tau)) for e in events], dtype=float
        ),
        "e_rev": np.array([e.e_rev for e in events], dtype=float),
    }
    for a in out.values():
        a.flags.writeable = False
    return out


def event_currents(
    t: float, v: np.ndarray, ev: dict[str, np.ndarray], scale: float
) -> np.ndarray:
    """Returns each event's current [pA] at time t; positive depolarises.

    scale * g * (e_rev - v) * waveform per event, the waveform being
    epsg_unitary's expression with the same float order (the mask
    (t > 0) and the clamp at EPSG_T_FLOOR included).

    Args:
        t: Time [ms].
        v: The voltage [mV] of each event's compartment, one per event.
        ev: event_arrays(...).
        scale: The run's I, multiplying every peak.

    Returns:
        The currents, one per event.
    """
    te = t - ev["t"]
    tc = np.maximum(te, EPSG_T_FLOOR)
    wave = (
        ev["inv_peak"]
        * (te > 0)
        * (np.exp(-tc / ev["decay"]) - np.exp(-tc / ev["rise"]))
    )
    return scale * ev["g"] * (ev["e_rev"] - v) * wave


def epsg_inv_peak(tau: tuple[float, float] = EPSG_TAU) -> float:
    """Returns the factor that makes epsg_unitary's peak 1.

    The default kinetics keep the MATLAB expression, 1 / EPSG_NORM, so
    results are unchanged; other (rise, decay) pairs are normalised
    exactly.

    Args:
        tau: The (rise, decay) time constants [ms].

    Returns:
        1 / peak of the difference of exponentials.
    """
    if tau == EPSG_TAU:
        return 1 / EPSG_NORM
    rise, decay = tau
    t_peak = np.log(decay / rise) * rise * decay / (decay - rise)
    peak = np.exp(-t_peak / decay) - np.exp(-t_peak / rise)
    return float(1 / peak)


def epsg_unitary(
    t: float | np.ndarray,
    tau: tuple[float, float] = EPSG_TAU,
    inv_peak: float | None = None,
) -> float | np.ndarray:
    """Returns the unitary EPSG waveform (difference of exponentials).

    The peak is 1 (see epsg_inv_peak).

    Args:
        t: Time since onset [ms]; zero for t <= 0.
        tau: The (rise, decay) time constants [ms].
        inv_peak: epsg_inv_peak(tau) if the caller has it already (the
            models read it from Stimulus.inv_peak, computed once per
            run); None computes it here.

    Returns:
        The waveform at t. Negative t gives -0.0, not 0.0, as it always
        has: the mask is 0 times a negative difference of exponentials.
    """
    if inv_peak is None:
        inv_peak = epsg_inv_peak(tau)
    rise, decay = tau
    # The exponentials are evaluated at t clamped to >= EPSG_T_FLOOR: for t
    # below about -71 ms exp(-t / rise) overflows and 0 * inf is nan. The mask
    # zeroes those values anyway, and the clamped difference is still negative,
    # so negative t keeps its -0.0.
    tc = np.maximum(t, EPSG_T_FLOOR)
    return inv_peak * (t > 0) * (np.exp(-tc / decay) - np.exp(-tc / rise))


@dataclasses.dataclass
class SynParams:
    """Stand-in for the MATLAB `Syn` struct.

    Only the fields the chosen stimulus uses need to be meaningful.

    Attributes:
        t_end: End of the conductance table [ms].
        freq: Input frequency [Hz].
        gE: EPSG conductance (non-meaningful units).
        VsynE: Excitatory reversal potential [mV].
        random_in: RNG seed.
        diff: Time offset of the second input for SynapticPair [ms].
        f: Sine frequency [Hz].
        epsg_tau: EPSG/EPSGpair (rise, decay) [ms]; Myoga et al 2014
            used decay 0.3.
        events: The "Events" stimulus: any number of Event records (see
            Event); ignored by the other stimuli.
    """

    t_end: float = 20.0
    freq: float = 350.0  # input frequency [Hz]
    gE: float = 1.3e7  # epsg conductance (non-meaningful units)
    VsynE: float = 0.0  # excitatory reversal potential [mV]
    random_in: int = 1396  # RNG seed
    diff: float = 1.0  # time offset of the second input for SynapticPair [ms]
    f: float = 200.0  # sine frequency [Hz]
    # EPSG/EPSGpair (rise, decay) [ms]; Myoga et al 2014 used decay 0.3
    epsg_tau: tuple = EPSG_TAU
    events: tuple[Event, ...] = ()  # the "Events" stimulus (see Event)


def synaptic(syn: SynParams) -> tuple[np.ndarray, np.ndarray]:
    """Makes the synaptic conductance: filtered noise convolved with an EPSG.

    The noise is gammatone-filtered and rectified. numpy's RNG differs
    from MATLAB's rng/randn, so the same seed gives a different
    (statistically equivalent) noise realisation.

    The table depends only on t_end, freq, gE and random_in, and every
    Synaptic or SynapticPair run asks for it again, so recent tables are
    kept (see _table). The arrays are shared between calls and read-only.

    Args:
        syn: The settings (t_end, freq, gE, random_in).

    Returns:
        (tSyn, gSyn): the time grid [ms] and the conductance on it, both
        read-only.
    """
    return _table(syn.t_end, syn.freq, syn.gE, syn.random_in)


@functools.lru_cache(maxsize=16)
def _table(
    t_end: float, freq: float, g_e: float, random_in: int
) -> tuple[np.ndarray, np.ndarray]:
    """Computes synaptic()'s tables, cached on the fields that fix them."""
    dt = 0.01
    t = np.arange(0, t_end + dt / 2, dt)
    nt = len(t)
    n = np.random.default_rng(random_in).standard_normal(nt)

    wc = freq * 2 * np.pi / 1000  # center frequency
    gw = 24.7 * (4.37 * wc / (2 * np.pi) + 1)
    gam_filter = (
        (t / 1000) ** 4 * np.exp(-2 * np.pi * t * gw / 1000) * np.cos(t * wc)
    )
    rect_gam = np.maximum(np.convolve(gam_filter, n)[:nt], 0)

    # 17.57 normalises so max is G, 37 is the unitary epsg amplitude (Lehnert et
    # al)
    epsg = (
        g_e * 37 * 17.57 * ((1 - np.exp(-t / 1.0)) ** 1.3 * np.exp(-t / 0.27))
    )
    g = np.convolve(rect_gam, epsg)[:nt]
    t.flags.writeable = g.flags.writeable = False  # shared by every caller
    return t, g


def interp_g(t_syn: np.ndarray, g_syn: np.ndarray, t: float) -> float:
    """Interpolates the conductance table like interp1q, but 0 outside it.

    MATLAB's interp1q gives NaN outside the table; TwoCptODE's guard
    against that is applied everywhere.

    Args:
        t_syn: The table's times [ms].
        g_syn: The table's conductances.
        t: The time [ms].

    Returns:
        The conductance at t as a float.
    """
    return float(np.interp(t, t_syn, g_syn, left=0.0, right=0.0))
