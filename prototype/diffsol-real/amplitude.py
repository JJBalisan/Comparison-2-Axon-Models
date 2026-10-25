"""Runs somatic.spike_amplitude with its capped second run on diffsol.

spike_amplitude runs the step twice: once to find the spike time ts, then
from 0 to ts + 2 ms with max_step 0.005 ms, and takes the inflection from a
cubic spline through that run's own steps. This script replaces the second
run with the diffsol-real driver (the cap emulated by set_stop_time hops,
or no cap at all) and compares the result with scipy's.

Usage (from the repository root, after cargo build --release):

    uv run python prototype/diffsol-real/amplitude.py SCRATCH_DIR
"""

import argparse
import json
import pathlib
import statistics
import subprocess
import sys
import time
import warnings
from typing import Any

import dump_cases
import numpy as np

import msoaxon
from msoaxon import _dispatch, _fast, _solve, somatic

HERE = pathlib.Path(__file__).resolve().parent
BINARY = HERE / "target" / "release" / "diffsol-real"
STEP_PA = 8500.0
CAP = 0.005


def _with_diffsol(scratch: pathlib.Path, max_step: str) -> dict[str, Any]:
    """spike_amplitude with the capped run replaced by diffsol's."""
    real = _dispatch.run_model
    info: dict[str, Any] = {}

    def run_model(*a: Any, **kw: Any) -> Any:
        if kw.get("max_step") != CAP:
            return real(*a, **kw)
        dump_cases.REC.reset()
        _solve._QUIET_CACHE.clear()
        t, y = real(*a, **kw)
        case = {
            "name": f"amplitude/second-run/{max_step}",
            "rhs": next(
                r for r in dump_cases.REC.rhs if r["stim_type"] != "none"
            ),
            "quiet_rhs": next(
                r for r in dump_cases.REC.rhs if r["stim_type"] == "none"
            ),
            **dump_cases.REC.solve,
            "stops": None,
            "axon_col": 2,
            "scipy": {"steps": len(t), "t_last": float(t[-1])},
        }
        src = scratch / "amplitude-case.json"
        out = scratch / f"amplitude-out-{max_step}.json"
        src.write_text(json.dumps(case))
        subprocess.run(  # noqa: S603 (our own binary, our own files)
            [str(BINARY), str(src), str(out), "--max-step", max_step]
            + ["--repeat", "5"],
            check=True,
            capture_output=True,
        )
        ds = json.loads(out.read_text())
        info.update(
            scipy_steps=len(t),
            diffsol_steps=ds["steps"],
            diffsol_ms=ds["wall_us"] / 1e3,
            lu=ds["stats"]["lu_setups"],
            t_end=float(t[-1]),
        )
        x = np.zeros((len(ds["t"]), y.shape[1]))
        x[:, 0] = ds["soma"]
        x[:, 2] = ds["axon"]
        return np.asarray(ds["t"]), x

    _dispatch.run_model = run_model
    try:
        res = somatic.spike_amplitude(
            "multi", STEP_PA, model_kw={"mode": "fast"}
        )
    finally:
        _dispatch.run_model = real
    return {**info, **_summary(res)}


def _summary(res: dict[str, Any] | None) -> dict[str, Any]:
    if res is None:
        return {"fires": False}
    keys = ("amplitude", "distinct", "peak_above_rest", "t_inflection")
    return {k: res[k] for k in keys} | {"t_spike": res["t_spike"]}


def _scipy() -> dict[str, Any]:
    """spike_amplitude as it is, timed (median of 3, cache cleared)."""
    walls, second = [], {}
    real = _dispatch.run_model

    def run_model(*a: Any, **kw: Any) -> Any:
        t0 = time.perf_counter()
        out = real(*a, **kw)
        if kw.get("max_step") == CAP:
            second["steps"] = len(out[0])
            second["ms"] = (time.perf_counter() - t0) * 1e3
        return out

    _dispatch.run_model = run_model
    try:
        for _ in range(3):
            _solve._QUIET_CACHE.clear()
            t0 = time.perf_counter()
            res = somatic.spike_amplitude(
                "multi", STEP_PA, model_kw={"mode": "fast"}
            )
            walls.append(time.perf_counter() - t0)
    finally:
        _dispatch.run_model = real
    return {
        "wall_ms": statistics.median(walls) * 1e3,
        "second_run_steps": second["steps"],
        "second_run_ms_last": second["ms"],
        **_summary(res),
    }


def main() -> None:
    """Prints scipy's and diffsol's spike_amplitude results as JSON."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("scratch", type=pathlib.Path, help="scratch directory")
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=msoaxon.DriftWarning)
    _fast._STATE["warned"] = True
    dump_cases._install()
    out = {
        "scipy": _scipy(),
        "diffsol_hops": _with_diffsol(a.scratch, str(CAP)),
        "diffsol_uncapped": _with_diffsol(a.scratch, "none"),
    }
    sys.stdout.write(json.dumps(out, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
