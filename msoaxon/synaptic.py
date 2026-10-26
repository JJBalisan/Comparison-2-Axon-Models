"""Synaptic input conductance (port of Synaptic.m)."""

import dataclasses
import functools
from collections.abc import Callable
from typing import Any

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


# The unitary inhibitory waveform of the in vitro MSO papers: rise 0.28 ms
# (Roberts et al 2013) and decay 1.5 ms (Myoga et al 2014 report 1.5 to 1.85
# ms), glycinergic, reversing below rest. In vivo IPSG kinetics are not in the
# literature: these are slice values, and the dt_inh results inherit that.
IPSG_TAU = (0.28, 1.5)
IPSG_E_REV = -80.0  # [mV]
IPSG_G = 40.0  # [nS]; the papers span 30 to 60
# A conductance event driven from this far away is a current source to within
# 1e-4 of its driving force: current_event uses it for EPSC-shaped currents.
CURRENT_E_REV = 1e5  # [mV]
# The default onset [ms] of the protocol builders (_common.START, which this
# module cannot import).
T0 = 5.0


def ipsg(
    t: float,
    g: float = IPSG_G,
    node: int = 1,
    tau: tuple[float, float] = IPSG_TAU,
    e_rev: float = IPSG_E_REV,
) -> Event:
    """Returns one inhibitory conductance event.

    Args:
        t: Onset [ms].
        g: Peak conductance [nS].
        node: The compartment it enters, 1-indexed.
        tau: The (rise, decay) time constants [ms].
        e_rev: Reversal potential [mV].

    Returns:
        An Event with the glycinergic defaults.
    """
    return Event(t, node, g, tuple(tau), e_rev)


def current_event(
    t: float,
    amp: float,
    node: int = 1,
    tau: tuple[float, float] = EPSG_TAU,
) -> Event:
    """Returns an event that injects an EPSC-shaped current, not a conductance.

    The event is a conductance of amp / CURRENT_E_REV nS driven from a
    reversal of CURRENT_E_REV mV, so its current, g * (e_rev - v), stays
    amp * waveform to within the swing of v (1e-3 of it). A run's I
    scales it like any event.

    Args:
        t: Onset [ms].
        amp: Peak current [pA]; negative hyperpolarises.
        node: The compartment it enters, 1-indexed.
        tau: The (rise, decay) time constants [ms].

    Returns:
        The Event.
    """
    return Event(t, node, amp / CURRENT_E_REV, tuple(tau), CURRENT_E_REV)


def epsc_train(
    n: int,
    rate: float,
    amp: float,
    t0: float = T0,
    node: int = 1,
    tau: tuple[float, float] = EPSG_TAU,
) -> tuple[Event, ...]:
    """Returns a regular train of EPSC-shaped current events.

    Args:
        n: The number of pulses.
        rate: The pulse rate [Hz].
        amp: Each pulse's peak current [pA].
        t0: The first onset [ms].
        node: The compartment it enters, 1-indexed.
        tau: The (rise, decay) time constants of each pulse [ms].

    Returns:
        n current events, 1000 / rate ms apart.
    """
    return tuple(
        current_event(t0 + k * 1000.0 / rate, amp, node, tau) for k in range(n)
    )


@dataclasses.dataclass(frozen=True)
class Inhibition:
    """An inhibitory event timed against the first EPSG of a pair.

    Attributes:
        dt: Onset relative to the pair's first EPSG [ms]; negative leads.
        g: Peak conductance [nS]. It does not scale with the searched
            amplitude: the search varies the excitation only.
        tau: The (rise, decay) time constants [ms].
        e_rev: Reversal potential [mV].
        node: The compartment it enters, 1-indexed.
    """

    dt: float = 0.0
    g: float = IPSG_G
    tau: tuple[float, float] = IPSG_TAU
    e_rev: float = IPSG_E_REV
    node: int = 1


@dataclasses.dataclass(frozen=True)
class Fixed:
    """A fixed event list whose peaks scale with the searched amplitude.

    The simplest event source (see Source): the list is the shape of
    the input, and the search's amplitude multiplies every g, as `I`
    does in mso_axon("Events", ...). It has no delay.

    Attributes:
        events: The events; their g is the peak at amplitude 1.
    """

    events: tuple[Event, ...]

    def __call__(self, _delay: float, amp: float) -> tuple[Event, ...]:
        """Returns the events with every peak multiplied by amp."""
        return tuple(dataclasses.replace(e, g=e.g * amp) for e in self.events)


@dataclasses.dataclass(frozen=True)
class PairSource:
    """Two EPSGs `delay` apart, with optional pre-pulses and inhibition.

    The first EPSG is at `start` and the second at start + delay, so a
    negative delay puts the second first; both have the searched amplitude.
    Pre-pulses and the IPSG are timed against the first EPSG and keep
    their own peaks. With a node2 the two EPSGs enter different
    compartments, as EPSGbilateral does.

    Attributes:
        start: The first EPSG's onset [ms].
        tau: The EPSGs' (rise, decay) time constants [ms].
        node: Where the first EPSG enters, 1-indexed.
        node2: Where the second enters, or None for node.
        pre: (dt, g) pairs: an EPSG of peak g [nS] at start + dt, with
            the pair's kinetics and node.
        inh: An Inhibition, or None.
    """

    start: float = T0
    tau: tuple[float, float] = EPSG_TAU
    node: int = 1
    node2: int | None = None
    pre: tuple[tuple[float, float], ...] = ()
    inh: Inhibition | None = None

    def __call__(self, delay: float, amp: float) -> tuple[Event, ...]:
        """Returns the events of one trial at this delay and amplitude."""
        tau = tuple(self.tau)
        out = [Event(self.start + dt, self.node, g, tau) for dt, g in self.pre]
        if self.inh is not None:
            i = self.inh
            out.append(
                ipsg(self.start + i.dt, i.g, i.node, i.tau, i.e_rev),
            )
        out.append(Event(self.start, self.node, amp, tau))
        node2 = self.node if self.node2 is None else self.node2
        out.append(Event(self.start + delay, node2, amp, tau))
        return tuple(out)


# An event source: (delay, amp) -> the events of one run. The coincidence and
# somatic searches run it with I = 1 and vary amp, so a source can hold some
# events fixed (an IPSG, a pre-pulse) while it scales the rest.
Source = Callable[[float, float], tuple[Event, ...]]
# what the searches' `events=` takes: a source, or a tuple of Events (see Fixed)
EventsLike = tuple[Event, ...] | Source


def as_source(events: EventsLike) -> Source:
    """Returns an event source for a list of events or a source.

    Args:
        events: A source (returned as given), or a non-empty sequence of
            Events, wrapped in Fixed.

    Returns:
        The source.

    Raises:
        ValueError: If events is an empty sequence.
    """
    if callable(events):
        return events
    events = tuple(events)
    if not events:
        raise ValueError('"Events" needs at least one Event')
    return Fixed(events)


def fibre_events(
    n: int,
    rng: np.random.Generator | int | None,
    t0: float = T0,
    g: float = 37.0,
    g_sd: float = 0.0,
    jitter: float = 0.0,
    node: int = 1,
    tau: tuple[float, float] = EPSG_TAU,
    release_p: float = 1.0,
    rate: float | None = None,
    duration: float = 0.0,
) -> tuple[Event, ...]:
    """Draws the excitatory events of n fibres (or synapses).

    Couchman et al 2010: MSO neurons get 4 to 8 fibres of 37 +/- 4 nS.
    With rate None every fibre fires once, at t0 plus a Gaussian jitter
    (the first-spike jitter of an AVCN bushy cell); with a rate each
    fibre fires a Poisson process of that rate over [t0, t0 + duration),
    its spikes uniform in the window. Each spike releases with
    probability release_p, and each release is an EPSG of peak
    g * N(1, g_sd / g) (never below 0).

    Args:
        n: The number of fibres.
        rng: A numpy Generator, or a seed for one.
        t0: The nominal onset [ms] (the start of the window with a rate).
        g: The mean peak [nS].
        g_sd: The standard deviation of the peak [nS] (Couchman: 4).
        jitter: The standard deviation of each onset [ms], with rate
            None.
        node: The compartment the events enter, 1-indexed.
        tau: The (rise, decay) time constants [ms].
        release_p: The probability a spike releases.
        rate: The Poisson rate per fibre [Hz], or None for one jittered
            spike each.
        duration: The window of the Poisson process [ms].

    Returns:
        The events, sorted by onset (none if nothing released), onsets
        clipped to be at least 0.

    Raises:
        ValueError: If n is negative, release_p is outside [0, 1], or
            rate is given with a duration that is not positive.
    """
    if n < 0 or not 0 <= release_p <= 1:
        raise ValueError(
            f"fibre_events needs n >= 0 and 0 <= release_p <= 1, got "
            f"{n}, {release_p}"
        )
    rng = np.random.default_rng(rng)
    if rate is None:
        times = (
            t0 + jitter * rng.standard_normal(n) if jitter else np.full(n, t0)
        )
    else:
        if duration <= 0:
            raise ValueError(
                f"a Poisson drive needs duration > 0, got {duration}"
            )
        counts = rng.poisson(rate * duration / 1000.0, n)
        times = t0 + duration * rng.random(counts.sum())
    released = rng.random(len(times)) < release_p
    times = times[released]
    peaks = g + g_sd * rng.standard_normal(len(times))
    order = np.argsort(times, kind="stable")
    tau = tuple(tau)
    return tuple(
        Event(float(max(times[i], 0.0)), node, float(max(peaks[i], 0.0)), tau)
        for i in order
    )


def bilateral_events(
    delay: float,
    rng: np.random.Generator | int | None,
    n: int = 4,
    nodes: tuple[int, int] = (1, 1),
    t0: float = T0,
    **kw: Any,
) -> tuple[Event, ...]:
    """Draws the events of two sides, the second `delay` after the first.

    Each side is fibre_events(n, ...): the sides differ in the onset
    window (t0 and t0 + delay) and node. With the default 4 fibres of
    37 nS this is the Couchman et al 2010 input; with n = 125 and a
    release_p of 0.45 it is the Casarez et al 2026 one.

    Args:
        delay: Onset of the second side relative to the first [ms];
            negative puts it first.
        rng: A numpy Generator, or a seed for one; both sides draw from it.
        n: The number of fibres (synapses) per side.
        nodes: The compartment of each side, 1-indexed.
        t0: The first side's nominal onset [ms].
        **kw: Further fibre_events arguments (g, g_sd, jitter, tau,
            release_p, rate, duration), the same for both sides.

    Returns:
        The events of both sides, sorted by onset.
    """
    rng = np.random.default_rng(rng)
    first = fibre_events(n, rng, t0, node=nodes[0], **kw)
    second = fibre_events(n, rng, t0 + delay, node=nodes[1], **kw)
    return tuple(sorted(first + second, key=lambda e: e.t))


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


def event_syn(
    source: Source, delay: float, amp: float, t_end: float
) -> SynParams:
    """Returns the SynParams of one run of an event source.

    Args:
        source: See as_source.
        delay: The delay passed to the source [ms].
        amp: The amplitude passed to the source.
        t_end: End of the run [ms].

    Returns:
        SynParams(t_end, events=source(delay, amp)).
    """
    return SynParams(t_end=t_end, events=tuple(source(delay, amp)))


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
