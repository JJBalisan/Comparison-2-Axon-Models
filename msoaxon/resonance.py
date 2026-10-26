"""Resonance of the MSO soma: impedance, sine-train and EPSG-train thresholds.

Mikiel-Hunter, Kotak & Rinzel 2016 (PLoS Comput Biol 12:e1005166) report three
frequencies for their gerbil MSO model at 35 C, all tied to KLT: subthreshold
resonance 242 Hz (the peak of the soma's impedance), spike resonance 285 Hz
(the sine frequency with the lowest firing threshold) and the EPSC-train
threshold minimum 317 Hz. This module measures the same three in both models:

- impedance(), subthreshold(): the soma's impedance from small sine currents,
  one frequency per run;
- sine_threshold(), spike_resonance(): the smallest sine train that fires, over
  frequency;
- epsg_train_threshold(), epsg_resonance(): the smallest EPSG train (an
  "Events" list at one rate) that fires, over rate.

The stimulus is the models' "sine": a half-wave rectified sine, I * max(sin, 0),
on from start to stop. There is no chirp (a new stim_type would need the Rust
side), so the impedance comes from a series of runs, and the paper's own
discrete-frequency protocol is a series too. Rectification does not bias the
estimate: the half-wave sine is I/pi + (I/2) sin(wt) plus even harmonics, so
over whole cycles the response at the fundamental divides by I/2, which is read
off the stimulus itself, and the mean and harmonics are orthogonal to it. They
do add a mean current of I/pi, 6 pA at the default 20 pA, about 20 uV across
the soma. The paper's sine is symmetric apart from a hyperpolarising half
scaled by 0.5; the rectified one is the nearest the models have.

Spike resonance differs from the paper's in four ways: the criterion is the
repo's (an axon compartment 10 mV above the soma, not dV/dt); the train lasts
`dur` ms (50, the paper's 500: the threshold keeps falling slowly with the
train length, by 4 to 9% from 50 to 500 ms); the threshold is the smallest peak
current with any spike, not the lowest amplitude that evokes spikes per cycle;
and the sine is rectified. The EPSG train has the repo's EPSG (rise and decay
0.1 and 0.18 ms), not the paper's 0.3 ms alpha.
"""

from collections.abc import Callable
from concurrent import futures
from typing import Any

import numpy as np

from msoaxon import _bisect, _common, _dispatch, _parallel, synaptic

START = _common.START  # stimulus onset [ms]
V0 = -68.0
NODE = 3  # the compartment compared with the soma, the AIS
FACTOR = _common.FACTOR  # spike = NODE rises this far above the soma [mV]

# impedance: the small-signal protocol
AMP = 20.0  # sine amplitude [pA]; the response is linear to ~1% up to 80 pA
SETTLE = 30.0  # discarded after onset [ms]
CYCLES = 10  # whole cycles analysed
STEPS_PER_CYCLE = 64  # solver max_step = period / this
# coarse frequencies [Hz]; 20 is the reference for Q (the lowest frequency)
IMPEDANCE_FREQS = (20.0, 50.0, *np.arange(100.0, 801.0, 50.0))

# thresholds
DUR = 50.0  # train length [ms]; thresholds fall 4-9% by 500 ms
TAIL = 5.0  # run time after the last input [ms]
SINE_FREQS = tuple(np.arange(100.0, 801.0, 50.0))
# the EPSG-train threshold dips only near the period of the membrane's ring
EPSG_RATES = tuple(np.arange(100.0, 801.0, 25.0))
SINE_GUESS = 2000.0  # first sine amplitude tried [pA]
SINE_CEILING = 2e4
EPSG_GUESS = 50.0  # first EPSG peak conductance tried [nS]
EPSG_CEILING = 2000.0
REL_TOL = 1e-3
N_FINE = 9  # points of the second, local pass
# A peak of |Z| (or dip of a threshold) is resonance only if it differs from
# the lowest frequency's value by more than this ratio: a membrane that is flat
# to 0.1% would otherwise have a "resonance" at a random frequency.
MIN_RATIO = 1.05


def _fundamental(phase: np.ndarray, y: np.ndarray, detrend: int = 2) -> complex:
    """Returns the complex amplitude c of y's component at the sine's phase.

    y ~ Im(c * exp(i * phase)) = a sin(phase) + b cos(phase), c = a + ib,
    fitted by least squares together with a polynomial in time, which
    removes the slow drift (the default membrane relaxes for tens of ms
    after the start, and the half-wave sine's mean current shifts it).
    Over whole cycles the polynomial is nearly orthogonal to the sine.

    Args:
        phase: The sine's phase at each sample [rad].
        y: The samples.
        detrend: Degree of the polynomial.

    Returns:
        c; abs(c) is the amplitude and angle(c) the phase lead of y over
        sin(phase).
    """
    x = np.linspace(-1.0, 1.0, len(y))
    cols = [np.sin(phase), np.cos(phase)] + [x**k for k in range(detrend + 1)]
    a, b = np.linalg.lstsq(np.stack(cols, axis=1), y, rcond=None)[0][:2]
    return complex(a, b)


def impedance(
    model: str,
    f: float,
    amp: float = AMP,
    settle: float = SETTLE,
    cycles: int = CYCLES,
    v0: float = V0,
    model_kw: dict[str, Any] | None = None,
) -> complex:
    """Returns the soma's impedance at f from a small rectified-sine run.

    The current enters the soma; the voltage is the soma's. The response
    over the last `cycles` whole cycles is projected on the sine (see
    _fundamental) and divided by the stimulus' own projection, which for
    the half-wave sine is amp / 2.

    Args:
        model: "multi" or "two".
        f: Frequency [Hz].
        amp: Peak current [pA]; small enough to stay linear (the result
            moves by 1% from 20 to 80 pA).
        settle: Time from the onset to the analysed window [ms].
        cycles: Whole cycles in the window.
        v0: Initial potential [mV].
        model_kw: Further model keywords, e.g. mem, mode, model_type,
            or r1/tau_est for "two".

    Returns:
        The impedance [MOhm]; angle() is the voltage's lead over the
        current.

    Raises:
        ValueError: If f, amp, settle or cycles is not positive.
    """
    for name, value in (
        ("f", f),
        ("amp", amp),
        ("settle", settle),
        ("cycles", cycles),
    ):
        _bisect.check_positive(name, value)
    period = 1000.0 / f
    t0 = START + settle
    t_end = t0 + cycles * period
    kw: dict[str, Any] = {
        "max_step": period / STEPS_PER_CYCLE,
        **(model_kw or {}),
    }
    t, x = _dispatch.run_model(
        model,
        "sine",
        START,
        t_end,
        amp,
        NODE,
        t_end,
        v0,
        1,
        synaptic.SynParams(t_end=t_end, f=f),
        **kw,
    )
    g = np.arange(t0, t_end, period / STEPS_PER_CYCLE)
    phase = 2 * np.pi * f * (g - START) / 1000.0
    current = amp * np.maximum(np.sin(phase), 0.0)
    # mV / pA = GOhm
    return (
        _fundamental(phase, np.interp(g, t, x[:, 0]))
        / _fundamental(phase, current)
        * 1e3
    )


def sine_threshold(
    model: str,
    f: float,
    dur: float = DUR,
    node: int = NODE,
    v0: float = V0,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = REL_TOL,
    guess: float = SINE_GUESS,
    ceiling: float = SINE_CEILING,
) -> float:
    """Finds the smallest rectified-sine train that makes `node` spike.

    The sine is on for `dur` ms from START and the run lasts TAIL ms
    longer; a spike is `node` rising FACTOR mV above the soma.

    Args:
        model: "multi" or "two".
        f: Frequency [Hz].
        dur: Train length [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        model_kw: Further model keywords, e.g. mem or mode.
        rel_tol: The relative tolerance of the threshold.
        guess: The first amplitude tried [pA].
        ceiling: The largest amplitude worth trying [pA].

    Returns:
        The threshold peak current [pA], or inf if nothing up to
        `ceiling` spikes.

    Raises:
        ValueError: If f or dur is not positive.
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
    _bisect.check_positive("f", f)
    _bisect.check_positive("dur", dur)
    stop = START + dur
    t_end = stop + TAIL
    syn = synaptic.SynParams(t_end=t_end, f=f)
    kw = model_kw or {}

    def respond(amp: float) -> tuple[bool, float]:
        return _dispatch.response(
            model,
            "sine",
            START,
            stop,
            amp,
            node,
            t_end,
            v0,
            1,
            syn,
            factor=FACTOR,
            **kw,
        )

    return _bisect.smallest_crossing(respond, guess, ceiling, rel_tol, FACTOR)


def epsg_train(rate: float, dur: float = DUR, node: int = 1) -> tuple:
    """Returns an EPSG train: events every 1/rate from START for dur ms.

    Args:
        rate: Events per second.
        dur: Train length [ms]; the train has round(rate * dur / 1000)
            events, at least one.
        node: The compartment every event enters, 1-indexed.

    Returns:
        A tuple of synaptic.Event, each of peak 1 nS (the `I` of the
        run scales them).
    """
    n = max(1, round(rate * dur / 1000.0))
    return tuple(
        synaptic.Event(START + k * 1000.0 / rate, node, 1.0) for k in range(n)
    )


def epsg_train_threshold(
    model: str,
    rate: float,
    dur: float = DUR,
    node: int = NODE,
    v0: float = V0,
    model_kw: dict[str, Any] | None = None,
    rel_tol: float = REL_TOL,
    guess: float = EPSG_GUESS,
    ceiling: float = EPSG_CEILING,
) -> float:
    """Finds the smallest EPSG train, at one rate, that makes `node` spike.

    The events are epsg_train(rate, dur) at the soma, all of one peak
    conductance (the waveform of the Events stimulus, rise and decay
    0.1 and 0.18 ms, reversal 0 mV).

    Args:
        model: "multi" or "two".
        rate: Events per second.
        dur: Train length [ms].
        node: The compartment compared with the soma, 1-indexed.
        v0: Initial potential [mV].
        model_kw: Further model keywords, e.g. mem or mode.
        rel_tol: The relative tolerance of the threshold.
        guess: The first peak conductance tried [nS].
        ceiling: The largest peak conductance worth trying [nS].

    Returns:
        The threshold peak conductance per event [nS], or inf if nothing
        up to `ceiling` spikes.

    Raises:
        ValueError: If rate or dur is not positive.
        _bisect.FiresWithoutInputError: If the cell fires with no input.
    """
    _bisect.check_positive("rate", rate)
    _bisect.check_positive("dur", dur)
    events = epsg_train(rate, dur)
    t_end = events[-1].t + TAIL
    syn = synaptic.SynParams(t_end=t_end, events=events)
    kw = model_kw or {}

    def respond(peak: float) -> tuple[bool, float]:
        return _dispatch.response(
            model,
            "Events",
            START,
            t_end,
            peak,
            node,
            t_end,
            v0,
            1,
            syn,
            factor=FACTOR,
            **kw,
        )

    return _bisect.smallest_crossing(respond, guess, ceiling, rel_tol, FACTOR)


def _call(task: tuple[Callable[..., Any], str, float, dict[str, Any]]) -> Any:
    """Runs fn(model, f, **kw); module-level so workers can import it."""
    fn, model, f, kw = task
    return fn(model, f, **kw)


def _vertex(f: np.ndarray, y: np.ndarray, i: int) -> float:
    """Returns the vertex of the parabola through y around index i.

    A least-squares parabola over up to five points centred on i. If it
    is not convex towards the extremum, or its vertex is outside the
    points, f[i] is returned.

    Args:
        f: The abscissae, increasing.
        y: The values (finite near i).
        i: The index of the extremum, with a neighbour on both sides.

    Returns:
        The vertex abscissa.
    """
    lo, hi = max(i - 2, 0), min(i + 3, len(f))
    a, b, _ = np.polyfit(f[lo:hi], y[lo:hi], 2)
    if a == 0 or (a < 0) != (y[i] >= np.median(y[lo:hi])):
        return float(f[i])
    v = -b / (2 * a)
    return float(v) if f[lo] <= v <= f[hi - 1] else float(f[i])


def scan(
    fn: Callable[..., float | complex],
    model: str,
    freqs: tuple[float, ...] | np.ndarray,
    maximum: bool,
    kw: dict[str, Any] | None = None,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    key: Callable[[Any], float] = float,
    n_fine: int = N_FINE,
) -> dict[str, Any]:
    """Finds the frequency that extremises fn(model, f), in two passes.

    A coarse pass over `freqs`, then `n_fine` points between the
    extremum's two neighbours, and the vertex of a parabola over the five
    points around the best.

    Args:
        fn: impedance, sine_threshold or epsg_train_threshold.
        model: "multi" or "two".
        freqs: The coarse frequencies, increasing [Hz].
        maximum: True to find the largest key(fn(...)), False the
            smallest.
        kw: Keywords for fn.
        workers: Processes for the independent runs (default: all CPUs;
            1 runs in-process).
        executor: An existing process pool to use instead.
        key: Maps fn's result to the number being extremised (abs for
            the impedance).
        n_fine: Points of the second pass (at least 3).

    Returns:
        dict(f_best, value, at_edge, freqs, values): the extremum's
        frequency (the vertex) and value there (the best run's), whether
        it sits at an end of `freqs` (then f_best is that end and no
        second pass ran), and every run, sorted by frequency (values
        are fn's results).
    """
    kw = kw or {}
    freqs = np.asarray(freqs, dtype=float)

    def run(fs: np.ndarray) -> list[Any]:
        tasks = [(fn, model, float(f), kw) for f in fs]
        return _parallel.map_tasks(_call, tasks, workers, executor)

    def numbers(vals: list[Any]) -> np.ndarray:
        y = np.array([key(v) for v in vals], dtype=float)
        return np.where(np.isfinite(y), y, np.nan)

    vals = run(freqs)
    y = numbers(vals)
    pick = np.nanargmax if maximum else np.nanargmin
    i = int(pick(y))
    at_edge = i in (0, len(freqs) - 1)
    f_best = float(freqs[i])
    if not at_edge:
        lo, hi = freqs[i - 1], freqs[i + 1]
        fine = np.linspace(lo, hi, n_fine + 2)[1:-1]
        fine = fine[~np.isin(fine, freqs)]
        vals += run(fine)
        allf = np.concatenate([freqs, fine])
        order = np.argsort(allf)
        allf, vals = allf[order], [vals[k] for k in order]
        y = numbers(vals)
        near = np.flatnonzero((allf >= lo) & (allf <= hi))
        j = int(near[pick(y[near])])
        f_best = _vertex(allf, y, j)
        i, freqs = j, allf
    return {
        "f_best": f_best,
        "value": vals[i],
        "at_edge": at_edge,
        "freqs": [float(f) for f in freqs],
        "values": vals,
    }


def _dip(r: dict[str, Any]) -> dict[str, Any]:
    """Returns scan's result for a threshold curve, with its dip."""
    th = [float(v) for v in r["values"]]
    dip = th[0] / float(r["value"])
    return {
        "f_best": r["f_best"],
        "threshold": float(r["value"]),
        "dip": dip,
        "resonant": not r["at_edge"] and dip > MIN_RATIO,
        "at_edge": r["at_edge"],
        "freqs": r["freqs"],
        "thresholds": th,
    }


def subthreshold(
    model: str,
    freqs: tuple[float, ...] | np.ndarray = IMPEDANCE_FREQS,
    model_kw: dict[str, Any] | None = None,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    n_fine: int = N_FINE,
    **impedance_kw: Any,
) -> dict[str, Any]:
    """Measures the soma's resonance: where |Z| peaks, and how sharply.

    Q is the peak |Z| over |Z| at the lowest frequency of `freqs` (20 Hz
    by default, so it is not over the DC input resistance the paper uses;
    the slow gates, Ih and KLT inactivation, lower |Z| there a little
    below it). Resonance needs an interior peak with Q > MIN_RATIO; a
    passive low-pass membrane peaks at the lowest frequency, and one
    that is flat has a Q of 1.

    Args:
        model: "multi" or "two".
        freqs: The coarse frequencies, increasing [Hz].
        model_kw: Further model keywords (see impedance).
        workers: Processes (see scan).
        executor: An existing process pool.
        n_fine: Points of the second pass (see scan).
        **impedance_kw: amp, settle, cycles or v0 for impedance().

    Returns:
        dict(f_res, z_peak, q, z_low, f_low, resonant, freqs, z, phase):
        f_res and z_peak [MOhm] (nan unless resonant), Q, |Z| at f_low,
        and the curve (|Z| in MOhm, phase in degrees).
    """
    r = scan(
        impedance,
        model,
        freqs,
        True,
        {"model_kw": model_kw, **impedance_kw},
        workers,
        executor,
        key=abs,
        n_fine=n_fine,
    )
    z = np.array([abs(v) for v in r["values"]])
    q = float(np.nanmax(z) / z[0])
    resonant = not r["at_edge"] and q > MIN_RATIO
    nan = float("nan")
    return {
        "f_res": r["f_best"] if resonant else nan,
        "z_peak": float(abs(r["value"])) if resonant else nan,
        "q": q,
        "z_low": float(z[0]),
        "f_low": r["freqs"][0],
        "resonant": resonant,
        "freqs": r["freqs"],
        "z": [float(v) for v in z],
        "phase": [float(np.degrees(np.angle(v))) for v in r["values"]],
    }


def spike_resonance(
    model: str,
    freqs: tuple[float, ...] | np.ndarray = SINE_FREQS,
    model_kw: dict[str, Any] | None = None,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    n_fine: int = N_FINE,
    **threshold_kw: Any,
) -> dict[str, Any]:
    """Finds the sine frequency with the lowest firing threshold.

    Args:
        model: "multi" or "two".
        freqs: The coarse frequencies, increasing [Hz].
        model_kw: Further model keywords (see sine_threshold).
        workers: Processes (see scan).
        executor: An existing process pool.
        n_fine: Points of the second pass (see scan).
        **threshold_kw: dur, node, v0, rel_tol, guess or ceiling.

    Returns:
        dict(f_best, threshold, dip, resonant, at_edge, freqs,
        thresholds): the frequency of the lowest threshold (a parabola's
        vertex), that threshold [pA], the threshold at the lowest
        frequency over it (the dip), whether the minimum is interior
        with a dip above MIN_RATIO, and the curve.
    """
    r = scan(
        sine_threshold,
        model,
        freqs,
        False,
        {"model_kw": model_kw, **threshold_kw},
        workers,
        executor,
        n_fine=n_fine,
    )
    return _dip(r)


def epsg_resonance(
    model: str,
    rates: tuple[float, ...] | np.ndarray = EPSG_RATES,
    model_kw: dict[str, Any] | None = None,
    workers: int | None = None,
    executor: futures.Executor | None = None,
    n_fine: int = N_FINE,
    **threshold_kw: Any,
) -> dict[str, Any]:
    """Finds the EPSG-train rate with the lowest firing threshold.

    Args:
        model: "multi" or "two".
        rates: The coarse rates, increasing [Hz].
        model_kw: Further model keywords (see epsg_train_threshold).
        workers: Processes (see scan).
        executor: An existing process pool.
        n_fine: Points of the second pass (see scan).
        **threshold_kw: dur, node, v0, rel_tol, guess or ceiling.

    Returns:
        dict(f_best, threshold, dip, resonant, at_edge, freqs,
        thresholds): the rate of the lowest threshold (a parabola's
        vertex), that threshold [nS per event], the threshold at the
        lowest rate over it (the dip), whether the minimum is interior
        with a dip above MIN_RATIO, and the curve.
    """
    r = scan(
        epsg_train_threshold,
        model,
        rates,
        False,
        {"model_kw": model_kw, **threshold_kw},
        workers,
        executor,
        n_fine=n_fine,
    )
    return _dip(r)
