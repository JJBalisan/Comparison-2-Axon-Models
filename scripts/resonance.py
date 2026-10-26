"""Resonance of the soma against Mikiel-Hunter et al 2016 (issue #83).

    uv run scripts/resonance.py --out-dir figures/resonance

Their gerbil MSO model at 35 C resonates at 242 Hz subthreshold (the peak of
the soma's impedance), 285 Hz for spikes (the sine frequency with the lowest
threshold) and 317 Hz for EPSC trains (the rate with the lowest threshold), all
KLT dependent. This measures the same three in the 45-compartment model, with
KLT replaced by a leak (multi.membrane's klt_to_leak) and in the two-compartment
model (which has no membrane(), so no replacement). See msoaxon/resonance.py for
the protocols and how they differ from the paper's.
"""

import argparse
import json
import math
import pathlib
from typing import Any

import msoaxon
from msoaxon import measure, multi, resonance

TARGETS = {"subthreshold": 242.0, "spike": 285.0, "epsg": 317.0}  # [Hz]
REBALANCED = {"rebalance_rest": True}
ARMS = {
    "45-compartment (msoAxon.m)": ("multi", None),
    "45-compartment, rest rebalanced": ("multi", REBALANCED),
    "45-compartment, KLT frozen at its resting conductance": (
        "multi",
        {**REBALANCED, "klt_to_leak": "chord"},
    ),
    "45-compartment, KLT replaced by its slope conductance": (
        "multi",
        {**REBALANCED, "klt_to_leak": "slope"},
    ),
    "2-compartment (TwoCpt.m)": ("two", None),
}
# --quick: a smoke test with coarse grids, short trains and loose tolerances
QUICK_ARMS = (
    "45-compartment (msoAxon.m)",
    "45-compartment, KLT replaced by its slope conductance",
    "2-compartment (TwoCpt.m)",
)
QUICK: dict[str, dict[str, Any]] = {
    "zap": {"freqs": (20.0, 200.0, 350.0, 500.0), "n_fine": 3},
    "sine": {"freqs": (250.0, 450.0, 650.0), "n_fine": 3},
    "epsg": {"rates": (300.0, 400.0, 500.0), "n_fine": 3},
    "thresholds": {"dur": 20.0, "rel_tol": 2e-2},
}


def _clean(x: Any) -> Any:
    """Returns x with nan and inf as None, so the JSON is strict."""
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_clean(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def _best(d: dict[str, Any]) -> str:
    """Names a threshold curve's resonance: its frequency, or that it has none.

    A minimum at the end of the swept range, or one less than
    resonance.MIN_RATIO below the lowest frequency's threshold, is no
    resonance; the text then gives the dip and where the minimum sat.
    """
    if d["resonant"]:
        return f"{d['f_best']:.0f} Hz"
    return f"none (dip {d['dip']:.2f}, minimum at {d['f_best']:.0f} Hz)"


def measure_arm(
    model: str,
    knobs: dict[str, Any] | None,
    mode: str | None,
    quick: bool,
    pool: Any,
) -> dict[str, Any]:
    """Measures the three resonances, and the input resistance, of one arm.

    Args:
        model: "multi" or "two".
        knobs: multi.membrane keywords, or None for the default membrane.
        mode: The model mode (None: the default).
        quick: Use QUICK's grids and tolerances.
        pool: A process pool.

    Returns:
        dict(rin, subthreshold, spike, epsg): passive_step's input
        resistance and time constant, and the three measurements of
        msoaxon.resonance.
    """
    kw: dict[str, Any] = {"mode": mode}
    if knobs:
        kw["mem"] = multi.membrane(-68.0, **knobs)
    q: dict[str, dict[str, Any]] = (
        QUICK
        if quick
        else {"zap": {}, "sine": {}, "epsg": {}, "thresholds": {}}
    )
    return {
        "rin": {
            k: float(v) for k, v in measure.passive_step(model, **kw).items()
        },
        "subthreshold": resonance.subthreshold(
            model, model_kw=kw, executor=pool, **q["zap"]
        ),
        "spike": resonance.spike_resonance(
            model, model_kw=kw, executor=pool, **q["sine"], **q["thresholds"]
        ),
        "epsg": resonance.epsg_resonance(
            model, model_kw=kw, executor=pool, **q["epsg"], **q["thresholds"]
        ),
    }


def summary(name: str, r: dict[str, Any]) -> str:
    """Returns the arm's one-line summary."""
    z, s, e = r["subthreshold"], r["spike"], r["epsg"]
    return (
        f"{name}: Rin {r['rin']['rin_steady']:.1f} MOhm; "
        + (
            f"subthreshold {z['f_res']:.0f} Hz (Q {z['q']:.2f})"
            if z["resonant"]
            else f"no subthreshold resonance (Q {z['q']:.2f})"
        )
        + f"; spike {_best(s)}; EPSG train {_best(e)}"
    )


def markdown(results: dict[str, Any]) -> str:
    """Returns the results as a Markdown table."""
    rows = [
        "| Model | R_in (MOhm) | Subthreshold resonance, Q "
        "| Spike resonance | EPSG-train minimum |",
        "|---|---|---|---|---|",
        f"| Mikiel-Hunter 2016 | | {TARGETS['subthreshold']:.0f} Hz "
        f"| {TARGETS['spike']:.0f} Hz | {TARGETS['epsg']:.0f} Hz |",
    ]
    for name, r in results.items():
        z = r["subthreshold"]
        sub = (
            f"{z['f_res']:.0f} Hz, Q {z['q']:.2f}"
            if z["resonant"]
            else f"none (Q {z['q']:.2f})"
        )
        rows.append(
            f"| {name} | {r['rin']['rin_steady']:.1f} | {sub} "
            f"| {_best(r['spike'])} | {_best(r['epsg'])} |"
        )
    return "\n".join(rows)


def plot(results: dict[str, Any], path: pathlib.Path) -> None:
    """Draws the three curves of every arm."""
    # imported here rather than at the top, like the pool scripts, so importing
    # this module stays cheap
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    panels = (
        ("subthreshold", "z", "Soma impedance |Z| (MOhm)", "Frequency (Hz)"),
        ("spike", "thresholds", "Sine-train threshold (pA)", "Frequency (Hz)"),
        (
            "epsg",
            "thresholds",
            "EPSG-train threshold (nS per event)",
            "Rate (Hz)",
        ),
    )
    for ax, (key, ykey, ylabel, xlabel) in zip(axes, panels, strict=True):
        for name, r in results.items():
            d = r[key]
            ax.plot(d["freqs"], d[ykey], ".-", ms=4, lw=1.2, label=name)
        ax.axvline(TARGETS[key], color="k", ls=":", lw=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if key != "subthreshold":
            ax.set_yscale("log")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(
        "Dotted lines: Mikiel-Hunter et al 2016 (242, 285, 317 Hz)", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main() -> None:
    """Runs the resonance protocols and saves the figure and JSON."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", default="figures/resonance")
    ap.add_argument(
        "--mode",
        choices=("exact", "fast", "native"),
        help="model mode (default: fast); native is about 8x faster and "
        "agrees within solver tolerance",
    )
    ap.add_argument("--workers", type=int, help="processes (default: all)")
    ap.add_argument(
        "--quick",
        action="store_true",
        help="coarse grids, 20 ms trains and 2%% tolerance, three models: "
        "a smoke test, not results",
    )
    a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    with msoaxon.process_pool(a.workers) as pool:
        for name, (model, knobs) in ARMS.items():
            if a.quick and name not in QUICK_ARMS:
                continue
            results[name] = measure_arm(model, knobs, a.mode, a.quick, pool)
            print(summary(name, results[name]), flush=True)
    print()
    print(markdown(results))
    plot(results, out / "resonance.png")
    (out / "resonance.json").write_text(
        json.dumps(
            _clean(
                {
                    "mode": a.mode or "default",
                    "quick": a.quick,
                    "targets_hz": TARGETS,
                    "results": results,
                }
            ),
            indent=1,
        )
    )
    print(f"saved {out}/resonance.png and .json")


if __name__ == "__main__":
    main()
