"""Event-list protocols (issue #78): fibres, pre-pulses, EPSC trains, ITD sweep.

    uv run scripts/event_protocols.py --out-dir figures/event_protocols

Four protocols on the "Events" stimulus, each for the 45-compartment model and
the two-compartment model of TwoCpt.m, printed as tables and saved as JSON:

1. Fibre drive (Couchman et al 2010, Myoga et al 2014): n fibres of 37 nS per
   side with Gaussian onset jitter, spike probability per delay, and the
   half-width of that curve over repeated experiments.
2. Pre-pulse (Franken et al 2015): the EPSG-pair threshold after an EPSG at
   -1.0 / -0.75 / -0.5 ms; the target is +7.2 / +6.6 nS at 1.0 / 0.75 ms.
3. EPSC trains (Khurana et al 2011): the input resistance during a train; the
   target is a rise of up to 60%.
4. A Casarez-style trial-based ITD sweep (125 synapses per side, release
   probability 0.45), with a small default number of trials.

Runs in the default mode (fast, bit-identical to the reference here, when the
Rust extension is installed). --mode native uses the Rust BDF integrator for
the 45-compartment model instead (within solver tolerance, about 8x faster),
and says so in the tables.
"""  # noqa: E501 (printed verbatim as --help)

import argparse
import json
import pathlib
from concurrent import futures
from typing import Any

import numpy as np

import msoaxon
from msoaxon import _fast, coincidence, measure

MODELS = {"45-compartment": "multi", "2-cpt, TwoCpt.m": "two"}
PRE_DTS = (-1.0, -0.75, -0.5)
FRANKEN_NS = {-1.0: 7.2, -0.75: 6.6}  # Franken et al 2015: the rise [nS]


def main() -> None:
    """Parses the command line and runs the protocols."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", default="figures/event_protocols")
    ap.add_argument(
        "--mode",
        choices=("fast", "native", "exact"),
        help="mso_axon's mode for the 45-compartment model (default: fast "
        "when the extension is installed)",
    )
    ap.add_argument("--repeats", type=int, default=8)
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument(
        "--quick",
        action="store_true",
        help="coarse grids, few trials and 1%% threshold tolerance: a smoke "
        "test, not results",
    )
    a = ap.parse_args()
    with msoaxon.process_pool() as pool:
        run(a, pool)


def _kw(a: argparse.Namespace, model: str) -> dict[str, Any]:
    """Returns the model keywords: the mode, for the 45-compartment model."""
    return {"mode": a.mode} if a.mode and model == "multi" else {}


def fibre_window(
    a: argparse.Namespace, pool: futures.Executor
) -> dict[str, Any]:
    """Spike probability per delay for fibre drive, and its half-width."""
    delays = np.round(np.arange(-0.8, 0.8001, 0.2 if a.quick else 0.1), 3)
    n_trials = 4 if a.quick else a.trials
    repeats = 2 if a.quick else a.repeats
    out: dict[str, Any] = {"delays_ms": delays.tolist()}
    print(
        f"\n1. Fibre drive: n fibres x 37 nS per side, Gaussian jitter "
        f"0.32 ms; {n_trials} trials per delay, {repeats} repeats"
    )
    print(
        "model | fibres/side | P(ITD 0) | P(|ITD| 0.8 ms, monaural) | "
        "half-width [us], mean +- sd (min-max over repeats)"
    )
    for name, model in MODELS.items():
        for n in (2, 3):
            prob, summ = coincidence.fibre_window_trials(
                model,
                delays,
                n_fibres=n,
                jitter=0.32,
                n_trials=n_trials,
                repeats=repeats,
                executor=pool,
                model_kw=_kw(a, model),
            )
            width = np.array([s["width"] for s in summ]) * 1e3
            ok = np.isfinite(width)
            m = prob.mean(axis=0)
            wtxt = (
                f"{width[ok].mean():.0f} +- {width[ok].std():.0f} "
                f"({width[ok].min():.0f}-{width[ok].max():.0f}), "
                f"{ok.sum()}/{repeats} repeats fall below half maximum"
                if ok.any()
                else f"never falls below half maximum (0/{repeats})"
            )
            print(
                f"{name} | {n} | {m[len(m) // 2]:.2f} | "
                f"{(m[0] + m[-1]) / 2:.2f} | {wtxt}",
                flush=True,
            )
            out[f"{model}_{n}"] = {
                "probability": prob.tolist(),
                "half_widths_us": width.tolist(),
            }
    return out


def prepulse(a: argparse.Namespace, pool: futures.Executor) -> dict[str, Any]:
    """The EPSG-pair threshold after a preceding EPSG."""
    tol = 1e-2 if a.quick else 1e-4
    gs = (37.0,) if a.quick else (10.0, 15.0, 20.0, 37.0)
    print(
        "\n2. Pre-pulse: pair threshold rise [nS] after an EPSG at -1.0 / "
        "-0.75 / -0.5 ms (Franken 2015: +7.2 / +6.6 at 1.0 / 0.75 ms)"
    )
    print(
        "model | pre-pulse peak [nS] | threshold, no pre-pulse | "
        + " | ".join(f"rise at {dt} ms" for dt in PRE_DTS)
    )
    out: dict[str, Any] = {}
    for name, model in MODELS.items():
        for g in gs:
            th0, th = coincidence.prepulse_thresholds(
                model,
                PRE_DTS,
                g,
                rel_tol=tol,
                executor=pool,
                model_kw=_kw(a, model),
            )
            print(
                f"{name} | {g:g} | {th0:.2f} | "
                + " | ".join(f"{x:+.1f}" for x in th - th0),
                flush=True,
            )
            out[f"{model}_{g:g}"] = {"th0": th0, "rise": (th - th0).tolist()}
    return out


def trains(a: argparse.Namespace, pool: futures.Executor) -> dict[str, Any]:
    """The input resistance during EPSC trains."""
    rates = (200.0, 400.0) if a.quick else (100.0, 200.0, 300.0, 400.0)
    n = 12 if a.quick else 20
    print(
        f"\n3. EPSC trains: {n} current pulses of 3 nA, input resistance "
        "during the train (Khurana 2011: up to +60%)"
    )
    print(
        "model | rate [Hz] | Rin rest [MOhm] | Rin train | change | depol [mV]"
    )
    out: dict[str, Any] = {}
    jobs = {
        (name, rate): pool.submit(
            measure.train_rin,
            model,
            rate,
            n,
            3000.0,
            **_kw(a, model),
        )
        for name, model in MODELS.items()
        for rate in rates
    }
    for (name, rate), job in jobs.items():
        r = job.result()
        print(
            f"{name} | {rate:g} | {r['rin_rest']:.2f} | {r['rin_train']:.2f} | "
            f"{r['change']:+.0%} | {r['depolarisation']:+.2f}",
            flush=True,
        )
        out[f"{MODELS[name]}_{rate:g}"] = r
    return out


def itd(a: argparse.Namespace, pool: futures.Executor) -> dict[str, Any]:
    """A Casarez-style trial-based ITD sweep."""
    itds = np.round(np.linspace(-0.5, 0.5, 6 if a.quick else 11), 4)
    n_trials = 4 if a.quick else a.trials
    print(
        f"\n4. ITD sweep: {coincidence.CASAREZ_N} synapses per side of "
        f"{coincidence.CASAREZ_G} nS, release probability "
        f"{coincidence.CASAREZ_P}, jitter 0.1 ms, {n_trials} trials per ITD"
    )
    print("model | peak | half-width [us] | centre [us] | centroid [us]")
    out: dict[str, Any] = {"itds_ms": itds.tolist()}
    for name, model in MODELS.items():
        prob, s = coincidence.itd_sweep(
            model,
            itds,
            n_trials=n_trials,
            executor=pool,
            model_kw=_kw(a, model),
        )
        print(
            f"{name} | {s['peak']:.2f} | {s['width'] * 1e3:.0f} | "
            f"{s['centre'] * 1e3:+.0f} | {s['centroid'] * 1e3:+.0f}",
            flush=True,
        )
        out[model] = {"probability": prob.tolist(), **s}
    return out


def run(a: argparse.Namespace, pool: futures.Executor) -> None:
    """Runs the four protocols, prints their tables and saves the JSON.

    Args:
        a: The parsed command-line arguments.
        pool: The process pool every parallel call shares.
    """
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mode = _fast.effective(a.mode)
    print(
        f"45-compartment model run in mode {mode!r}"
        + (
            " (within solver tolerance of the reference)"
            if mode == "native"
            else ""
        )
        + "; the two-compartment model has only its numpy path"
    )
    res = {
        "mode": mode,
        "fibres": fibre_window(a, pool),
        "prepulse": prepulse(a, pool),
        "trains": trains(a, pool),
        "itd": itd(a, pool),
    }
    (out / "event_protocols.json").write_text(json.dumps(res, indent=1))
    print(f"\nsaved {out}/event_protocols.json")


if __name__ == "__main__":
    main()
