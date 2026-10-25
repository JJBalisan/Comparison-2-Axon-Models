"""Compares diffsol-real traces with the golden traces and with scipy.

Each golden case is sampled at tests/test_golden.TIMES with np.interp,
exactly as test_golden.sample does, and scored with the margin consumed of
scripts/phase0_harness.py: max |got - want| / (atol + rtol |want|) over the
soma and axon samples (1.0 = the tolerance edge). Prints a markdown table.

Usage (from the repository root):

    uv run python prototype/diffsol-real/compare.py CASE_DIR OUT_DIR

CASE_DIR holds dump_cases.py's JSON files, OUT_DIR the driver's outputs
under the same file names.
"""

import argparse
import importlib
import json
import pathlib
import sys
from typing import Any

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
TESTS = REPO / "tests"


def _golden_module() -> Any:
    sys.path.insert(0, str(TESTS))
    try:
        return importlib.import_module("test_golden")
    finally:
        sys.path.remove(str(TESTS))


GOLDEN_MOD = _golden_module()
RTOL = GOLDEN_MOD.TOLERANCE["rtol"]
ATOL = GOLDEN_MOD.TOLERANCE["atol"]


def sample(t: list[float], soma: list[float], axon: list[float]) -> dict:
    """test_golden.sample on a stored trace."""
    times = [s for s in GOLDEN_MOD.TIMES if s <= t[-1]]
    return {
        "t_last": float(t[-1]),
        "soma": list(np.interp(times, t, soma)),
        "axon": list(np.interp(times, t, axon)),
    }


def margin(got: dict, want: dict) -> float:
    """Margin consumed over the soma and axon samples."""
    out = 0.0
    for k in ("soma", "axon"):
        g, w = np.asarray(got[k]), np.asarray(want[k])
        if g.shape != w.shape:
            return float("inf")
        out = max(out, float(np.max(np.abs(g - w) / (ATOL + RTOL * abs(w)))))
    return out


def rows(case_dir: pathlib.Path, out_dir: pathlib.Path) -> list[dict]:
    """Returns one result row per golden case found in both directories."""
    golden = json.loads((TESTS / "golden_traces.json").read_text())
    result = []
    for path in sorted(case_dir.glob("multi-*.json")):
        out_path = out_dir / path.name
        if not out_path.exists():
            continue
        case = json.loads(path.read_text())
        ds = json.loads(out_path.read_text())
        sp = case["scipy"]
        want = golden[case["name"]]
        got = sample(ds["t"], ds["soma"], ds["axon"])
        ref = sample(sp["t"], sp["soma"], sp["axon"])
        result.append(
            {
                "case": case["name"],
                "scipy_steps": sp["steps"],
                "diffsol_steps": ds["steps"],
                "diffsol_ms": ds["wall_us"] / 1e3,
                "python_ms": sp["wall_s_cold"] * 1e3,
                "margin": margin(got, want),
                "scipy_margin": margin(ref, want),
                "vs_scipy": margin(got, ref),
                "t_last_ok": got["t_last"] == want["t_last"],
                "t_last_raw_ok": ds["t_last_raw"] == want["t_last"],
            }
        )
    return result


def table(result: list[dict]) -> str:
    """Formats rows() as markdown."""
    lines = [
        "| case | scipy steps | diffsol steps | diffsol wall ms "
        "| python fast wall ms | ratio | margin | scipy margin "
        "| margin vs scipy | t_last ok (raw) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in result:
        lines.append(
            f"| {r['case']} | {r['scipy_steps']} | {r['diffsol_steps']} "
            f"| {r['diffsol_ms']:.1f} | {r['python_ms']:.1f} "
            f"| {r['python_ms'] / r['diffsol_ms']:.1f} "
            f"| {r['margin']:.3g} | {r['scipy_margin']:.3g} "
            f"| {r['vs_scipy']:.3g} "
            f"| {r['t_last_ok']} ({r['t_last_raw_ok']}) |"
        )
    worst = max(r["margin"] for r in result)
    lines.append("")
    lines.append(
        f"{len(result)} cases; max margin {worst:.3g}; t_last exact in "
        f"{sum(r['t_last_ok'] for r in result)} (raw "
        f"{sum(r['t_last_raw_ok'] for r in result)})."
    )
    return "\n".join(lines)


def main() -> None:
    """Prints the comparison table."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", type=pathlib.Path, help="dump_cases.py output")
    ap.add_argument("outs", type=pathlib.Path, help="driver outputs")
    a = ap.parse_args()
    sys.stdout.write(table(rows(a.cases, a.outs)) + "\n")


if __name__ == "__main__":
    main()
