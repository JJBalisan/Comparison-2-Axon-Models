"""Resting Na availability: rest x Na inactivation V-half, both models (#79).

    uv run scripts/na_availability.py --out-dir figures/na_availability

Scott, Mathews & Golding 2010 measured somatic Na inactivation with
V-half = -77 mV and ~92% of channels inactivated at rest (-58 mV); the
models use -62.5 mV and rest at -68 mV. The 2 x 2 (rest {-68, -58} x V-half
{-62.5, -77}) shows what the somatic spike, rheobase, input resistance,
EPSP half-width and 3% coincidence window do under each combination.

The V-half moves the soma only, in both models:
multi.membrane(soma_na_vhalf=..., rebalance_rest=True) and
two_cpt(soma_na_vhalf=...) (which holds rest at v0 by subtracting the
resting currents). Nothing here changes a default; the table is a
proposal for the PYTHON.md "Somatic spike" section.
"""

import argparse
import json
import math
import os
import pathlib
from typing import Any

import numpy as np

import msoaxon
from msoaxon import (
    _dispatch,
    coincidence,
    constants,
    measure,
    multi,
    somatic,
    synaptic,
)

RESTS = (-68.0, -58.0)
VHALFS = (-62.5, -77.0)
MULTS = (1.5, 2.0, 3.0)
MODELS = {"45-compartment": "multi", "2-compartment": "two"}
CELLS = len(MODELS) * len(RESTS) * len(VHALFS)
# a rest this unsteady (largest deviation from v0 between 100 and 200 ms after
# a kick, [mV]) has no input resistance to measure: a 2 pA step moves a steady
# cell by 5 uV
RIPPLE_MAX = 0.01

# a single EPSG this fraction of the coincident threshold: deeply
# subthreshold, so the half-width is the EPSP's and not a spike's
EPSP_FRACTION = 0.25


def inactivated(rest: float, vhalf: float) -> float:
    """Returns the fraction of Na channels inactivated at rest.

    Args:
        rest: The resting potential [mV].
        vhalf: The signed inactivation midpoint [mV].

    Returns:
        1 - constants.hinf(rest) with the model's own slope.
    """
    return 1.0 - float(constants.hinf(rest, -vhalf))


def model_kw(
    model: str, rest: float, vhalf: float, mode: str | None
) -> dict[str, Any]:
    """Returns the model keywords of one cell of the table.

    Args:
        model: "multi" or "two".
        rest: The resting potential [mV].
        vhalf: The somatic Na inactivation midpoint [mV].
        mode: The 45-compartment model's solver mode; None for
            MSOAXON_MODE or the default.

    Returns:
        The keywords for run_model; the two-compartment model ignores
        mode, so it is left out.
    """
    if model == "two":
        return {"soma_na_vhalf": vhalf}
    mem = multi.membrane(rest, soma_na_vhalf=vhalf, rebalance_rest=True)
    return {"mem": mem, "mode": mode}


def ripple(model: str, rest: float, kw: dict[str, Any]) -> float:
    """Measures whether rest is steady: a kick, then how far V strays.

    The run starts at the steady state, which an unstable one never
    leaves without a nudge, so a 2 pA, 5 ms step perturbs it. A stable
    rest settles within milliseconds and what remains is the solver's
    1e-6 mV; the 45-compartment model at -58 mV does not settle: it
    rings with 0.18 mV at a 1.4 ms period, whatever the kick, solver
    mode or step size.

    Args:
        model: "multi" or "two".
        rest: The resting potential [mV].
        kw: The model keywords (model_kw).

    Returns:
        The largest |V - rest| of the soma between 100 and 200 ms [mV].
    """
    t, x = _dispatch.run_model(
        model,
        "step",
        5,
        10,
        2.0,
        3,
        200,
        rest,
        1,
        **{**kw, "max_step": 0.05},
    )
    return float(np.abs(x[t >= 100, 0] - rest).max())


def epsp(
    model: str, rest: float, kw: dict[str, Any], amplitude: float
) -> dict[str, float]:
    """Measures a single somatic EPSG's soma depolarization.

    Args:
        model: "multi" or "two".
        rest: The resting potential [mV].
        kw: The model keywords (model_kw).
        amplitude: The EPSG amplitude, in the units of the "EPSG"
            stimulus.

    Returns:
        dict(peak [mV above rest], half_width [us]).
    """
    start = coincidence.START
    syn = synaptic.SynParams(t_end=coincidence.T_END)
    t, x = _dispatch.run_model(
        model,
        "EPSG",
        start,
        coincidence.T_END,  # mso_axon cuts EPSG off at stop
        amplitude,
        3,
        coincidence.T_END,
        rest,
        1,
        syn,
        **{**kw, "max_step": 0.005},
    )
    grid = np.arange(start - 0.05, start + 10.0, 1e-3)
    v = np.interp(grid, t, x[:, 0])
    dv = v - np.interp(start - 0.05, t, x[:, 0])
    peak = dv.max()
    above = np.nonzero(dv >= peak / 2)[0]
    return {
        "peak": float(peak),
        "half_width": float((grid[above[-1]] - grid[above[0]]) * 1e3),
    }


def measure_cell(
    task: tuple[str, float, float, str | None, bool],
) -> dict[str, Any]:
    """Runs the measurements of one cell of the table.

    With quick, a rest that is not steady gets only its ripple: the
    searches cost minutes there (the 45-compartment model at -58 mV
    takes tiny steps through its oscillation), and none of them means
    much without a steady rest.

    Args:
        task: (model, rest, vhalf, mode, quick): "multi" or "two"; the
            resting potential and the somatic Na inactivation midpoint
            [mV]; the 45-compartment model's solver mode (None:
            MSOAXON_MODE or the default); and whether to use coarse
            tolerances and two multiples (a smoke test). One tuple so a
            process pool can map it.

    Returns:
        A dict of the cell's numbers.
    """
    model, rest, vhalf, mode, quick = task
    kw = model_kw(model, rest, vhalf, mode)
    rip = ripple(model, rest, kw)
    if quick and rip >= RIPPLE_MAX:
        return {
            "rest_mV": rest,
            "vhalf_mV": vhalf,
            "inactivated": inactivated(rest, vhalf),
            "rest_ripple_mV": rip,
            "rheobase_pA": None,
            "rin_steady_MOhm": None,
        }
    tol = 1e-2 if quick else 1e-4
    rb = somatic.rheobase(model, v0=rest, rel_tol=tol, model_kw=kw)
    amps, peaks = {}, {}
    for k in (1.5, 3.0) if quick else MULTS:
        r = somatic.spike_amplitude(model, rb * k, v0=rest, model_kw=kw)
        # None: no spike, or no distinct inflection (peak is still valid)
        amps[str(k)] = (
            r["amplitude"] if r is not None and r["distinct"] else None
        )
        peaks[str(k)] = None if r is None else float(r["peak_above_rest"])
    rin = (
        measure.passive_step(model, v0=rest, **kw) if rip < RIPPLE_MAX else None
    )
    width, th0 = coincidence.window(
        model, 0.03, v0=rest, model_kw=kw, rel_tol=tol, delay_tol=tol
    )
    e = epsp(model, rest, kw, EPSP_FRACTION * th0)
    return {
        "rest_mV": rest,
        "vhalf_mV": vhalf,
        "inactivated": inactivated(rest, vhalf),
        "rheobase_pA": float(rb),
        "spike_amplitude_mV": amps,
        "spike_peak_above_rest_mV": peaks,
        "rest_ripple_mV": rip,
        "rin_steady_MOhm": None if rin is None else float(rin["rin_steady"]),
        "tau_ms": None if rin is None else float(rin["tau"]),
        "epsp_peak_mV": e["peak"],
        "epsp_half_width_us": e["half_width"],
        "window_3pct_us": float(width) * 1e3,
        "threshold0": float(th0),
    }


def fmt(x: float | None, digits: int = 1) -> str:
    """Formats a number, or "n/a" for None and nan."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "n/a"
    return f"{x:.{digits}f}"


def table(results: dict) -> str:
    """Renders the results as one Markdown table per measurement."""
    lines = []
    for label, key, digits in (
        ("Na inactivated at rest", "inactivated", 2),
        ("Rheobase [pA]", "rheobase_pA", 0),
        ("Input resistance [MOhm]", "rin_steady_MOhm", 2),
        ("Rest ripple, max |V - rest| with no input [mV]", "rest_ripple_mV", 3),
        ("EPSP half-width [us]", "epsp_half_width_us", 0),
        ("3% coincidence window [us]", "window_3pct_us", 0),
    ):
        lines += [
            f"{label}",
            "| model | rest | V-half | value |",
            "|---|---|---|---|",
        ]
        for name, cells in results.items():
            for c in cells:
                lines.append(
                    f"| {name} | {c['rest_mV']:.0f} | {c['vhalf_mV']:g} "
                    f"| {fmt(c.get(key), digits)} |"
                )
        lines.append("")
    lines += [
        "Somatic spike amplitude [mV]",
        "| model | rest | V-half | 1.5x | 2x | 3x |",
        "|---|---|---|---|---|---|",
    ]
    for name, cells in results.items():
        for c in cells:
            a = c.get("spike_amplitude_mV", {})
            p = c.get("spike_peak_above_rest_mV", {})
            lines.append(
                f"| {name} | {c['rest_mV']:.0f} | {c['vhalf_mV']:g} | "
                + " | ".join(
                    fmt(a[str(k)])
                    if a.get(str(k)) is not None
                    else (
                        f"no inflection (peak {fmt(p[str(k)])})"
                        if str(k) in a
                        else "-"
                    )
                    for k in MULTS
                )
                + " |"
            )
    return "\n".join(lines)


def main() -> None:
    """Runs the 2 x 2 for both models and prints and saves the table."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", default="figures/na_availability")
    ap.add_argument(
        "--mode",
        default=None,
        choices=("exact", "fast", "native"),
        help="the 45-compartment model's solver mode (default: "
        "MSOAXON_MODE, else native when the extension is built; the "
        "2-compartment model runs scipy whatever this is)",
    )
    ap.add_argument(
        "--quick",
        action="store_true",
        help="1%% tolerances, two multiples of rheobase, and no searches "
        "where the rest is not steady: a smoke test, not results",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=min(CELLS, os.cpu_count() or 1),
        help="processes (default: one per cell, or the cores if fewer)",
    )
    a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    tasks = [
        (model, rest, vhalf, a.mode, a.quick)
        for model in MODELS.values()
        for rest in RESTS
        for vhalf in VHALFS
    ]
    cells = msoaxon.map_tasks(measure_cell, tasks, a.workers)
    results = {
        name: [c for c, t in zip(cells, tasks, strict=True) if t[0] == model]
        for name, model in MODELS.items()
    }
    for name, rows in results.items():
        for c in rows:
            print(
                f"{name}, rest {c['rest_mV']:.0f}, V-half {c['vhalf_mV']:g}: "
                f"rheobase {fmt(c['rheobase_pA'], 0)} pA"
            )
    print()
    print(table(results))
    (out / "na_availability.json").write_text(json.dumps(results, indent=1))
    print(f"\nsaved {out}/na_availability.json (mode={a.mode or 'default'})")


if __name__ == "__main__":
    main()
