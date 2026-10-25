"""Dumps msoaxon runs as JSON cases for the diffsol-real driver.

For every lumped active-full golden case of tests/test_golden.py (the 15 of
scripts/phase0_harness.py), plus the 110 ms somatic step (with and without
max_step=0.005), the golden EPSGpair run stopped by spike_event, and every
run of one coincidence.threshold search, this runs msoaxon in fast mode and
captures exactly what reaches the solver: the keyword arguments of both
msoaxon_fast.MultiRhs objects (the stimulus and the stimulus-free "none"
right-hand side of the quiet prefix), y0, the breakpoints, t_end, the
tolerances, max_step and the stop events. It also stores scipy's own result
(every accepted step of the soma and the compared axon column, the final
state) and its wall time. One JSON file per case goes into OUT_DIR.

Usage (from the repository root):

    uv run python prototype/diffsol-real/dump_cases.py OUT_DIR [--fixture F]

--fixture also writes the small test fixture of the Rust crate (the
EPSGpair case without its trace).
"""

import argparse
import importlib
import json
import pathlib
import statistics
import sys
import time
import warnings
from collections.abc import Callable
from typing import Any

import numpy as np

import msoaxon
from msoaxon import _common, _fast, _solve, coincidence, somatic

REPO = pathlib.Path(__file__).resolve().parents[2]
TESTS = REPO / "tests"
REPEATS = 3
SOMATIC_STEP_PA = 8500.0
SPIKE_FACTOR = coincidence.FACTOR  # the threshold search's criterion [mV]
FIXTURE_CASE = "multi/EPSGpair"


def _plain(v: Any) -> Any:
    """Returns v as JSON-ready Python values (arrays become lists)."""
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [_plain(x) for x in v]
    return v


class Recorder:
    """Captures what one mso_axon call hands to the Rust RHS and solver."""

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        """Forgets the previous run."""
        self.rhs: list[dict[str, Any]] = []
        self.solve: dict[str, Any] = {}
        self.stops: dict[str, Any] | None = None


REC = Recorder()


def _install() -> None:
    """Wraps msoaxon_fast.MultiRhs, _solve.integrate and _common._stops."""
    real_cls = _fast._ext.MultiRhs

    def multi_rhs(**kw: Any) -> Any:
        REC.rhs.append(_plain(kw))
        return real_cls(**kw)

    _fast._ext.MultiRhs = multi_rhs

    real_integrate = _solve.integrate

    def integrate(
        rhs: Any, y0: Any, t_end: float, cuts: list, **kw: Any
    ) -> Any:
        REC.solve = {
            "y0": _plain(np.asarray(y0, dtype=float)),
            "t_end": float(t_end),
            "cuts": [float(c) for c in cuts],
            "rtol": float(kw["rtol"]),
            "atol": float(kw["atol"]),
            "max_step": float(kw["max_step"]),
            "use_quiet": kw.get("quiet") is not None,
        }
        return real_integrate(rhs, y0, t_end, cuts, **kw)

    _solve.integrate = integrate
    # multi.py calls _solve.integrate through the module, so this suffices.

    real_stops = _common._stops

    def stops(
        dv: Any,
        axon_col: int,
        stop_on_spike: float | None,
        stop_when_settled: bool,
        stim_type: str,
        s: Any,
    ) -> Any:
        if stop_on_spike is None:
            REC.stops = None
        else:
            t_off = _solve.input_end(
                stim_type, s.start, s.stop, s.t_end, s.epsg_tau
            )
            REC.stops = {
                "axon_col": int(axon_col),
                "factor": float(stop_on_spike),
                "settled": bool(stop_when_settled and t_off < s.t_end),
                "t_off": float(t_off),
                "settled_fraction": _solve.SETTLED_FRACTION,
            }
        return real_stops(
            dv, axon_col, stop_on_spike, stop_when_settled, stim_type, s
        )

    _common._stops = stops


def _golden() -> Any:
    """Imports tests/test_golden.py as phase0_harness does."""
    sys.path.insert(0, str(TESTS))
    try:
        return importlib.import_module("test_golden")
    finally:
        sys.path.remove(str(TESTS))


def _lumped_active_full(kw: dict[str, Any]) -> bool:
    return (
        "morph" not in kw
        and "mem" not in kw
        and kw.get("model_type", "active-full") == "active-full"
    )


def _timed(run: Callable[[], Any], cold: bool) -> tuple[Any, float]:
    """Returns (result, median wall [s]) of REPEATS runs after a warm-up.

    cold clears the quiet-prefix cache before every run.
    """
    _solve._QUIET_CACHE.clear()
    out = run()
    walls = []
    for _ in range(REPEATS):
        if cold:
            _solve._QUIET_CACHE.clear()
        t0 = time.perf_counter()
        run()
        walls.append(time.perf_counter() - t0)
    return out, statistics.median(walls)


def capture(
    name: str, run: Callable[[], Any], axon_col: int, trace: bool = True
) -> dict[str, Any]:
    """Runs one case, timed, and returns its JSON record."""
    _, warm = _timed(run, cold=False)
    (t, y), cold = _timed(run, cold=True)
    REC.reset()
    _solve._QUIET_CACHE.clear()
    t, y = run()
    stim = [r for r in REC.rhs if r["stim_type"] != "none"]
    quiet = [r for r in REC.rhs if r["stim_type"] == "none"]
    case = {
        "name": name,
        "rhs": stim[0] if stim else quiet[0],
        "quiet_rhs": quiet[0] if quiet else None,
        **REC.solve,
        "stops": REC.stops,
        "axon_col": axon_col,
        "scipy": {
            "steps": len(t),
            "t_last": float(t[-1]),
            "final": _plain(y[-1]),
            "wall_s_cold": cold,
            "wall_s_warm": warm,
        },
    }
    if trace:
        case["scipy"].update(
            t=_plain(t), soma=_plain(y[:, 0]), axon=_plain(y[:, axon_col])
        )
    return case


def golden_runs() -> list[tuple[str, Callable[[], Any], int]]:
    """Returns (name, run, axon column) for the 15 golden cases."""
    mod = _golden()
    out = []
    for name, (model, kw, axon) in sorted(mod.CASES.items()):
        if model != "multi" or not _lumped_active_full(kw):
            continue

        def run(kw: dict[str, Any] = kw) -> Any:
            return mod._run("multi", **kw, mode="fast")

        out.append((name, run, axon))
    return out


def _somatic(**kw: Any) -> Any:
    return msoaxon.mso_axon(
        "step",
        somatic.START,
        somatic.STOP,
        SOMATIC_STEP_PA,
        3,
        "active-full",
        somatic.T_END,
        -68.0,
        1,
        msoaxon.SynParams(t_end=somatic.T_END),
        mode="fast",
        **kw,
    )


def _epsgpair_spike() -> Any:
    return msoaxon.mso_axon(
        "EPSGpair",
        5,
        5.2,
        100,
        3,
        "active-full",
        20,
        -68.0,
        1,
        msoaxon.SynParams(t_end=20),
        stop_on_spike=SPIKE_FACTOR,
        mode="fast",
    )


def search_runs(out: pathlib.Path) -> dict[str, Any]:
    """Dumps every run of one threshold search; returns its summary."""
    real = _solve.integrate
    runs: list[dict[str, Any]] = []

    def counting(*a: Any, **kw: Any) -> Any:
        # REC.solve/stops/rhs are filled by the wrappers inside real()
        res = real(*a, **kw)
        stim = [r for r in REC.rhs if r["stim_type"] != "none"]
        quiet = [r for r in REC.rhs if r["stim_type"] == "none"]
        runs.append(
            {
                "rhs": stim[-1],
                "quiet_rhs": quiet[-1] if quiet else None,
                **REC.solve,
                "stops": REC.stops,
                "axon_col": 2,
                "scipy": {
                    "steps": len(res[0]),
                    "t_last": float(res[0][-1]),
                    "final": _plain(res[1][-1]),
                },
            }
        )
        REC.rhs.clear()
        return res

    def search() -> float:
        return coincidence.threshold("multi", 0.1, model_kw={"mode": "fast"})

    _solve._QUIET_CACHE.clear()
    search()  # warm-up
    walls = []
    for _ in range(REPEATS):
        _solve._QUIET_CACHE.clear()
        t0 = time.perf_counter()
        search()
        walls.append(time.perf_counter() - t0)
    _solve._QUIET_CACHE.clear()
    REC.reset()
    _solve.integrate = counting
    try:
        amp = search()
    finally:
        _solve.integrate = real
    for i, r in enumerate(runs):
        r["name"] = f"search/{i:02d}"
        (out / f"search-{i:02d}.json").write_text(json.dumps(r))
    return {
        "threshold": amp,
        "runs": len(runs),
        "wall_s_best": min(walls),
        "wall_s_median": statistics.median(walls),
        "scipy_steps": [r["scipy"]["steps"] for r in runs],
    }


def _fname(name: str) -> str:
    return name.replace("/", "-") + ".json"


def main() -> None:
    """Writes every case into the output directory."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=pathlib.Path, help="output directory")
    ap.add_argument("--fixture", type=pathlib.Path, help="test fixture path")
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=msoaxon.DriftWarning)
    _fast._STATE["warned"] = True  # no self_check() objects in the capture
    _install()
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [
        *golden_runs(),
        ("extra/somatic-step-110ms", _somatic, 2),
        ("extra/EPSGpair-spike-event", _epsgpair_spike, 2),
    ]
    for name, run, axon in jobs:
        case = capture(name, run, axon)
        (a.out / _fname(name)).write_text(json.dumps(case))
        sys.stdout.write(f"{name}: {case['scipy']['steps']} steps\n")
        if a.fixture and name == FIXTURE_CASE:
            small = {k: v for k, v in case.items() if k != "scipy"}
            small["scipy"] = {
                k: case["scipy"][k] for k in ("steps", "t_last", "final")
            }
            a.fixture.write_text(json.dumps(small, indent=0) + "\n")
    name = "extra/somatic-step-110ms-maxstep0.005"
    case = capture(name, lambda: _somatic(max_step=0.005), 2, trace=False)
    (a.out / _fname(name)).write_text(json.dumps(case))
    sys.stdout.write(f"{name}: {case['scipy']['steps']} steps\n")
    summary = search_runs(a.out)
    (a.out / "search-summary.json").write_text(json.dumps(summary, indent=1))
    sys.stdout.write(f"search: {json.dumps(summary)}\n")


if __name__ == "__main__":
    main()
