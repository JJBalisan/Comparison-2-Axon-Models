"""Phase 0 measurement harness for the Rust integrator (issue #28).

    uv run scripts/phase0_harness.py                      # every mode
    uv run scripts/phase0_harness.py --mode fast --json out.json

Per mode it reports three things, as a markdown table and optionally JSON:

1. Margin consumed: for each lumped active-full case of tests/test_golden.py,
   max over the golden samples (soma and axon, at the golden sample grid) of
   |got - want| / (atol + rtol * |want|), rtol 1e-4 and atol 5e-3 mV, with
   `want` read from tests/golden_traces.json exactly as test_golden does.
   1.0 is the edge of the tolerance; 0.0 is a bit-identical trace.
2. Wall time (best of 3 after one warm-up) of the multi-compartment
   coincidence.threshold at delay 0.1 ms and coincidence.window at margin
   0.03. The quiet-prefix cache is emptied before every timed call, so each
   call pays for its own quiet prefix and none is served from an earlier
   call, case or mode.
3. Solver step count of one 20 ms EPSGpair run and one 110 ms somatic step
   run: a different BDF that takes more steps gives back part of the gain
   per step.
"""

import argparse
import importlib
import json
import pathlib
import platform
import sys
import time
import types
import warnings
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

import msoaxon
from msoaxon import _fast, _solve, coincidence, somatic

REPO = pathlib.Path(__file__).resolve().parents[1]
TESTS = REPO / "tests"
RTOL, ATOL = 1e-4, 5e-3  # tests/test_golden.py TOLERANCE [mV]
REPEATS = 3
SOMATIC_STEP_PA = 8500.0  # a 100 ms somatic step that spikes (golden step)
EPSG_PAIR = ("EPSGpair", 5, 5.2, 100, 20)  # the golden multi/EPSGpair case


def _golden_module() -> types.ModuleType | None:
    """Imports tests/test_golden.py, or returns None if that fails."""
    sys.path.insert(0, str(TESTS))
    try:
        return importlib.import_module("test_golden")
    except ImportError:
        return None
    finally:
        sys.path.remove(str(TESTS))


def _is_lumped_active_full(kw: dict[str, Any]) -> bool:
    return (
        "morph" not in kw
        and "mem" not in kw
        and kw.get("model_type", "active-full") == "active-full"
    )


def golden_cases() -> tuple[list[str], Callable[..., dict], dict]:
    """Returns (case names, sample function, golden dict).

    The names are the lumped active-full multi cases of
    tests/test_golden.py and the sampler is its own sample(), so the
    compared quantity and sample grid are the test's. If the test module
    can't be imported the 11 STIMULI cases are rebuilt here (the same
    parameters, sampled the same way).
    """
    golden = json.loads((TESTS / "golden_traces.json").read_text())
    mod = _golden_module()
    if mod is not None:
        names = [
            n
            for n, (model, kw, _) in sorted(mod.CASES.items())
            if model == "multi" and _is_lumped_active_full(kw)
        ]
        return names, mod.sample, golden

    stimuli = {  # copy of test_golden.STIMULI; keep in step with it
        "step": ("step", 5, 10, 8500, 20, {}),
        "step-subthreshold": ("step", 5, 10, 500, 20, {}),
        "ramp": ("ramp", 5, 5.5, 8500, 20, {}),
        "ramp2": ("ramp2", 5, 5.5, 8500, 20, {}),
        "sine": ("sine", 5, 7.5, 8500, 20, {}),
        "Synaptic": ("Synaptic", 5, 15, 0, 20, {"random_in": 1396}),
        "SynapticPair": ("SynapticPair", 5, 15, 0, 20, {"diff": 1.0}),
        "EPSG": ("EPSG", 5, 10, 108, 20, {}),
        "EPSG-subthreshold": ("EPSG", 5, 10, 26.7, 20, {}),
        "EPSGpair": ("EPSGpair", 5, 5.2, 100, 20, {}),
        "EPSGpair-subthreshold": ("EPSGpair", 5, 5.3, 30, 20, {}),
    }
    times = (1.0, 5.1, 5.3, 5.6, 6.0, 7.0, 9.0, 12.0)

    def sample(name: str, mode: str = "exact") -> dict[str, Any]:
        stim, start, stop, amp, t_end, syn_kw = stimuli[name.split("/")[1]]
        syn = msoaxon.SynParams(t_end=t_end, **syn_kw)
        t, x = msoaxon.mso_axon(
            stim, start, stop, amp, 3, "active-full", t_end, -68.0, 1, syn,
            mode=mode,
        )  # fmt: skip
        ts = [s for s in times if s <= t[-1]]
        return {
            "t_last": float(t[-1]),
            "soma": list(np.interp(ts, t, x[:, 0])),
            "axon": list(np.interp(ts, t, x[:, 2])),
        }

    return [f"multi/{n}" for n in sorted(stimuli)], sample, golden


def margin_consumed(got: Sequence[float], want: Sequence[float]) -> float:
    """Returns max |got - want| / (ATOL + RTOL * |want|) over a trace pair."""
    got, want = np.asarray(got), np.asarray(want)
    return float(np.max(np.abs(got - want) / (ATOL + RTOL * np.abs(want))))


def margins(
    mode: str, names: list[str], sample: Callable[..., dict], golden: dict
) -> dict[str, dict[str, Any]]:
    """Returns {case: {"margin": float, "t_last_equal": bool}} for a mode."""
    out = {}
    for name in names:
        _solve._QUIET_CACHE.clear()
        got, want = sample(name, mode), golden[name]
        out[name] = {
            "margin": max(
                margin_consumed(got[k], want[k]) for k in ("soma", "axon")
            ),
            "t_last_equal": got["t_last"] == want["t_last"],
        }
    return out


def best_of(fn: Callable[[], object], repeats: int = REPEATS) -> float:
    """Returns the best wall time [s] of fn() after one warm-up call."""
    _solve._QUIET_CACHE.clear()
    fn()
    best = float("inf")
    for _ in range(repeats):
        _solve._QUIET_CACHE.clear()
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def timings(mode: str) -> dict[str, float]:
    """Returns the best-of-3 seconds of one threshold and one window search."""
    kw = {"mode": mode}
    return {
        "threshold_s": best_of(
            lambda: coincidence.threshold("multi", 0.1, model_kw=kw)
        ),
        "window_s": best_of(
            lambda: coincidence.window("multi", 0.03, model_kw=kw)
        ),
    }


def step_counts(mode: str) -> dict[str, int]:
    """Returns len(t) of the full EPSGpair run and the somatic step run."""
    stim, start, stop, amp, t_end = EPSG_PAIR
    _solve._QUIET_CACHE.clear()
    t_epsg, _ = msoaxon.mso_axon(
        stim, start, stop, amp, 3, "active-full", t_end, -68.0, 1,
        msoaxon.SynParams(t_end=t_end), mode=mode,
    )  # fmt: skip
    _solve._QUIET_CACHE.clear()
    t_som, _ = msoaxon.mso_axon(
        "step", somatic.START, somatic.STOP, SOMATIC_STEP_PA, 3,
        "active-full", somatic.T_END, -68.0, 1,
        msoaxon.SynParams(t_end=somatic.T_END), mode=mode,
    )  # fmt: skip
    return {"steps_epsgpair_20ms": len(t_epsg), "steps_step_110ms": len(t_som)}


def measure(
    mode: str, names: list[str], sample: Callable[..., dict], golden: dict
) -> dict[str, Any]:
    """Returns every number for one mode."""
    m = margins(mode, names, sample, golden)
    return {
        "margins": m,
        "max_margin": max(v["margin"] for v in m.values()),
        "t_last_all_equal": all(v["t_last_equal"] for v in m.values()),
        **timings(mode),
        **step_counts(mode),
    }


def table(results: dict[str, dict[str, Any]]) -> str:
    """Returns the markdown summary and per-case margin tables."""
    modes = list(results)
    lines = [
        "| mode | max margin | t_last equal | threshold 0.1 ms [s] "
        "| window 0.03 [s] | steps EPSGpair 20 ms | steps step 110 ms |",
        "|---|---|---|---|---|---|---|",
    ]
    for mode, r in results.items():
        lines.append(
            f"| {mode} | {r['max_margin']:.3g} | {r['t_last_all_equal']} "
            f"| {r['threshold_s']:.3f} | {r['window_s']:.3f} "
            f"| {r['steps_epsgpair_20ms']} | {r['steps_step_110ms']} |"
        )
    lines += [
        "",
        "Margin consumed per case (1.0 = tolerance edge):",
        "",
        "| case | " + " | ".join(modes) + " |",
        "|---|" + "---|" * len(modes),
    ]
    for case in results[modes[0]]["margins"]:
        cells = (f"{results[m]['margins'][case]['margin']:.3g}" for m in modes)
        lines.append(f"| {case} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main() -> None:
    """Measures the requested modes and prints (and saves) the results."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--mode",
        action="append",
        choices=_fast.MODES,
        help="mode to measure (repeatable; default: every mode)",
    )
    ap.add_argument("--json", metavar="PATH", help="also write the numbers")
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=msoaxon.DriftWarning)
    modes = a.mode or list(_fast.MODES)
    if "fast" in modes and not _fast.available():
        ap.error('mode "fast" needs the fast extra (uv sync --extra fast)')

    names, sample, golden = golden_cases()
    results = {m: measure(m, names, sample, golden) for m in modes}
    print(table(results))
    if a.json:
        meta = {
            "machine": platform.node(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "repeats": REPEATS,
            "tolerance": {"rtol": RTOL, "atol": ATOL},
        }
        pathlib.Path(a.json).write_text(
            json.dumps({"meta": meta, "modes": results}, indent=1) + "\n"
        )


if __name__ == "__main__":
    main()
