"""Passive and synaptic measurements on the soma.

Shared by the scripts and the tests.
"""

from typing import Any

import numpy as np

from msoaxon import _dispatch, synaptic


def soma_on_grid(
    t: np.ndarray,
    x: np.ndarray,
    v0: float,
    t0: float,
    t1: float,
    dt: float = 0.001,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolates the soma voltage onto a uniform grid.

    Args:
        t: Times of the run [ms].
        x: The run's state, soma voltage in column 0.
        v0: The resting potential to subtract [mV].
        t0: Grid start [ms].
        t1: Grid end (excluded) [ms].
        dt: Grid spacing [ms].

    Returns:
        (grid, soma voltage relative to v0 on it).
    """
    g = np.arange(t0, t1, dt)
    return g, np.interp(g, t, x[:, 0]) - v0


def passive_step(
    model: str,
    v0: float = -68.0,
    model_type: str = "active-full",
    I: float = -2.0,
    n_grid: int = 200001,
    **model_kw: Any,
) -> dict[str, float]:
    """Measures soma input resistance and time constant from a current step.

    The step lasts 40 ms (5 to 45 ms). The response is the soma voltage
    minus that of an identical run without the step, so any drift of the
    membrane away from v0 cancels. The responses are microvolts (2 pA
    through a few MOhm), while a membrane that is not rest-balanced
    drifts by about 0.1 mV before the step and 0.01 mV during it; a
    constant baseline (v0, or even the voltage at the step onset) is
    swamped by that. The grid stops at 44.9 ms: the implicit solver's
    step that ends exactly at switch-off is evaluated with the current
    already off, so that one point has decayed a little.

    Assumes the step is small enough to stay in the linear range, and
    that the response is a monotone relaxation to a steady state: the
    last grid point is taken as steady, and tau is the first crossing of
    63% of it. Either sign of I gives positive resistances; the peak is
    the largest deflection in the direction of the step. Resistances
    depend on the operating point, so the default membrane (resting a
    little below v0) reads about 2% above the rest-balanced one.

    Args:
        model: "multi" or "two".
        v0: Initial potential [mV].
        model_type: Which channels are active.
        I: Step amplitude [pA], nonzero.
        n_grid: Points of the analysis grid.
        **model_kw: Further run_model keywords (mem, morph, r1, ...).

    Returns:
        dict(rin_steady, rin_peak [MOhm], tau [ms]); tau is the time to
        63% of the steady response.

    Raises:
        ValueError: If I is zero.
    """
    if I == 0:
        raise ValueError("passive_step needs a nonzero step current I")
    start = 5

    def soma(current: float) -> np.ndarray:
        t, x = _dispatch.run_model(
            model,
            "step",
            start,
            45,
            current,
            3,
            50,
            v0,
            1,
            model_type=model_type,
            **model_kw,
        )
        return np.interp(g, t, x[:, 0])

    g = np.linspace(start, 44.9, n_grid)
    dv = soma(I) - soma(0.0)
    w = dv * np.sign(I)  # the response in the direction of the step
    tau = g[np.argmax(w >= w[-1] * (1 - np.exp(-1)))] - start
    return {
        "rin_steady": dv[-1] / I * 1e3,
        "rin_peak": dv[np.argmax(w)] / I * 1e3,
        "tau": tau,
    }


def curve_summary(
    x: np.ndarray | list[float], p: np.ndarray | list[float]
) -> dict[str, float]:
    """Summarises a tuning curve (spike probability against delay or ITD).

    The width is the full width at half maximum: the distance between the
    two points where the curve, walking out from its highest sample,
    first falls below half of it, each found by linear interpolation. The
    centre is the midpoint of those two points; the centroid is the
    probability-weighted mean of x over the whole curve.

    Args:
        x: The sample positions, increasing [ms].
        p: The curve at them.

    Returns:
        dict(peak, width, centre, centroid, lo, hi), lo and hi being the
        half-maximum points in x. width, centre, lo and hi are nan when
        the curve never falls below half its peak on that side (so all
        four are nan for a flat curve), and everything but peak is nan
        if the curve is zero throughout.
    """
    x, p = np.asarray(x, float), np.asarray(p, float)
    i = int(np.argmax(p))
    peak = float(p[i])
    nan = float("nan")
    if peak <= 0:
        return {
            "peak": peak,
            "width": nan,
            "centre": nan,
            "centroid": nan,
            "lo": nan,
            "hi": nan,
        }
    half = peak / 2

    def edge(step: int) -> float:
        j = i
        while 0 <= j + step < len(p):
            j += step
            if p[j] < half:
                # interpolate between the last sample at or above half and j
                return float(
                    np.interp(half, [p[j], p[j - step]], [x[j], x[j - step]])
                )
        return nan

    lo, hi = edge(-1), edge(1)
    return {
        "peak": peak,
        "width": hi - lo,
        "centre": (lo + hi) / 2,
        "centroid": float(np.sum(x * p) / np.sum(p)),
        "lo": lo,
        "hi": hi,
    }


def train_rin(
    model: str,
    rate: float = 300.0,
    n: int = 10,
    amp: float = 30.0,
    probe: float = -100.0,
    probe_tau: tuple[float, float] = (2.0, 4.0),
    v0: float = -68.0,
    model_type: str = "active-full",
    **model_kw: Any,
) -> dict[str, float]:
    """Measures the soma input resistance during an EPSC train.

    Khurana et al 2011 saw it rise by up to 60% during trains while the
    summed Ih and KLT conductance stayed nearly constant. The train is n
    EPSC-shaped currents (synaptic.epsc_train: current events, not
    conductances, so it adds no shunt). The probe is a slow hyperpolarising
    current event, slow enough (probe_tau against a membrane time constant
    under a millisecond) that its response at its own peak is quasi-steady,
    and timed so that peak falls in the middle of the train. The input
    resistance there is the probe's response in the train, the difference
    of the voltage with and without the probe, over its peak current.
    Rest is the same without the train.

    Assumes the train stays subthreshold (the model does not spike during
    the probe's peak) and the probe is small enough to stay linear.

    Args:
        model: "multi" or "two".
        rate: The train rate [Hz].
        n: The number of pulses.
        amp: Each pulse's peak current [pA].
        probe: The probe's peak current [pA], nonzero.
        probe_tau: The probe's (rise, decay) time constants [ms].
        v0: Initial potential [mV].
        model_type: Which channels are active.
        **model_kw: Further run_model keywords (mem, morph, r1, ...).

    Returns:
        dict(rin_rest, rin_train [MOhm], change = rin_train / rin_rest -
        1, depolarisation [mV], the train's own soma depolarisation at
        the probe's peak).

    Raises:
        ValueError: If probe is zero or the train is too short to hold
            the probe's peak.
    """
    if probe == 0:
        raise ValueError("train_rin needs a nonzero probe current")
    rise, decay = probe_tau
    t_peak = np.log(decay / rise) * rise * decay / (decay - rise)
    dt = 1000.0 / rate
    t0 = synaptic.T0
    t_probe = t0 + (n - 1) * dt / 2 - t_peak
    if t_probe < 0.5:
        raise ValueError(
            f"a train of {n} pulses at {rate} Hz is too short for a "
            f"probe peaking {t_peak:.1f} ms after its onset"
        )
    t_end = t0 + n * dt + 10.0
    model_kw = {"max_step": 0.01, **model_kw}
    when = t_probe + t_peak
    train = synaptic.epsc_train(n, rate, amp, t0)
    pulse = synaptic.current_event(t_probe, probe, tau=probe_tau)
    null = synaptic.current_event(t0, 0.0)

    # one train period around the probe's peak: the train's own transients
    # ride on the response, and over a whole period they average out
    grid = np.linspace(when - dt / 2, when + dt / 2, 201)
    i_probe = probe * synaptic.epsg_unitary(grid - t_probe, probe_tau)  # pA

    def soma(events: tuple[synaptic.Event, ...]) -> np.ndarray:
        t, x = _dispatch.run_model(
            model,
            "Events",
            0.0,
            0.0,
            1.0,
            3,
            t_end,
            v0,
            1,
            synaptic.SynParams(t_end=t_end, events=events),
            model_type=model_type,
            **model_kw,
        )
        return np.interp(grid, t, x[:, 0])

    def resistance(dv: np.ndarray) -> float:
        return float(np.mean(dv) / np.mean(i_probe) * 1e3)

    v_rest = soma((null,))
    rin_rest = resistance(soma((pulse,)) - v_rest)
    v_train = soma(train)
    rin_train = resistance(soma((*train, pulse)) - v_train)
    return {
        "rin_rest": rin_rest,
        "rin_train": rin_train,
        "change": rin_train / rin_rest - 1,
        "depolarisation": float(np.mean(v_train - v_rest)),
    }
