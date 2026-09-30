"""Benchmark the threshold searches: simulation runs, simulated time, wall time, error.

    uv run scripts/benchmark_search.py                    # current settings
    uv run scripts/benchmark_search.py --no-settled-stop  # every non-spiking run to t_end
    uv run scripts/benchmark_search.py --check            # also run each yes/no both ways

Runs a fixed set of searches serially in this process, counting every simulation
through _dispatch.run_model. Each value is compared with a reference computed at
rel_tol 1e-7 without the settled stop (cached in figures/benchmark/reference.json;
--reference recomputes it), so methods that change results can be compared with
methods that don't. --check makes every spikes() call run with and without the
settled stop and counts the decisions that differ (it must be 0).
"""

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path

from msoaxon import _dispatch, coincidence, somatic
from msoaxon import threshold as threshold_mod
from msoaxon.coincidence import threshold
from msoaxon.multi import membrane, with_dendrites
from msoaxon.somatic import rheobase
from msoaxon.threshold import binary_search

D = with_dendrites()
MID_L, MID_M = int(D.lateral[2]), int(D.medial[2])


def cases():
    """(name, stimulus, search(rel_tol) -> value). rel_tol is ignored by the
    BinarySearch.m port, which has its own fixed grid."""
    out = []
    for model in ("multi", "two"):
        for d in (0.0, 0.1, 0.3, 0.6, 1.0):
            out.append((f"EPSG pair {model}, delay {d}", "EPSGpair",
                        lambda tol, m=model, d=d: threshold(m, d, rel_tol=tol)))
    for d in (0.0, 0.3):
        out.append((f"EPSG bilateral dendrites, delay {d}", "EPSGbilateral",
                    lambda tol, d=d: threshold("multi", d, rel_tol=tol, model_kw=dict(morph=D),
                                               stim="EPSGbilateral", input_node=MID_L,
                                               input_node2=MID_M)))
    dtx = membrane(-68.0, soma_klt_scale=0, ais_klt_scale=0, rebalance_rest=True)
    for label, model, mem in (("multi", "multi", None), ("two", "two", None),
                              ("multi, KLT blocked", "multi", dtx)):
        out.append((f"rheobase {label}", "step (100 ms)",
                    lambda tol, m=model, mem=mem: rheobase(m, mem=mem, rel_tol=tol)))
    for stim, n, factor, max_I, kw in (("EPSGpair", 3, 10, 150, dict(epsg_pair_dt=0.1)),
                                       ("step", 1, 30, 15000, {}), ("sine", 1, 30, 15000, {}),
                                       ("ramp", 2, 30, 15000, {}), ("EPSG", 1, 30, 15000, {})):
        out.append((f"BinarySearch.m {stim}", stim,
                    lambda tol, s=stim, n=n, f=factor, M=max_I, kw=kw: binary_search(
                        s, n, 3, f, M, workers=1, rounded=False, **kw)))
    return out


def instrument(stats, check):
    """Count runs and simulated time; with check, compare both stop settings."""
    run_model, spikes = _dispatch.run_model, _dispatch.spikes

    def counted(*a, **kw):
        t, x = run_model(*a, **kw)
        stats["runs"] += 1
        stats["sim_ms"] += t[-1]
        return t, x

    def checked(*a, **kw):
        answer = spikes(*a, **kw)
        if check:
            saved, _dispatch.STOP_WHEN_SETTLED = _dispatch.STOP_WHEN_SETTLED, False
            n, ms = stats["runs"], stats["sim_ms"]
            if spikes(*a, **kw) != answer:
                stats["mismatches"] += 1
            stats["runs"], stats["sim_ms"] = n, ms  # the check's runs don't count
            _dispatch.STOP_WHEN_SETTLED = saved
        return answer

    # the analysis modules imported spikes by name, so replace it there too
    _dispatch.run_model = counted
    for mod in (_dispatch, coincidence, somatic, threshold_mod):
        mod.spikes = checked


def as_list(v):
    return [float(x) for x in (v[0] + v[1] if isinstance(v, tuple) else [v])]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-settled-stop", action="store_true", help="run non-spiking runs to t_end")
    ap.add_argument("--check", action="store_true", help="run every yes/no both ways")
    ap.add_argument("--rel-tol", type=float, default=1e-4)
    ap.add_argument("--reference", action="store_true", help="recompute the reference values")
    ap.add_argument("--out", default="figures/benchmark")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ref_path = out / "reference.json"
    todo = cases()

    if a.reference or not ref_path.exists():
        _dispatch.STOP_WHEN_SETTLED = False
        ref = {}
        for name, _, search in todo:
            ref[name] = as_list(search(1e-7))
            print(f"reference {name}: {ref[name]}", flush=True)
        ref_path.write_text(json.dumps(ref, indent=1))
    ref = json.loads(ref_path.read_text())

    _dispatch.STOP_WHEN_SETTLED = not a.no_settled_stop
    stats = defaultdict(float)
    instrument(stats, a.check)
    rows, by_stim = [], defaultdict(lambda: defaultdict(float))
    for name, stim, search in todo:
        before = dict(stats)
        t0 = time.perf_counter()
        value = as_list(search(a.rel_tol))
        wall = time.perf_counter() - t0
        err = max(abs(v - r) / abs(r) for v, r in zip(value, ref[name]) if r not in (0, float("inf")))
        row = dict(case=name, stimulus=stim, runs=int(stats["runs"] - before.get("runs", 0)),
                   sim_ms=stats["sim_ms"] - before.get("sim_ms", 0), wall_s=wall, max_rel_err=err,
                   mismatches=int(stats["mismatches"] - before.get("mismatches", 0)))
        rows.append(row)
        for k in ("runs", "sim_ms", "wall_s", "mismatches"):
            by_stim[stim][k] += row[k]
        print(f"{name:38s} runs {row['runs']:4d}  sim {row['sim_ms']:8.1f} ms  wall {wall:6.2f} s  "
              f"err {err:.1e}" + (f"  mismatches {row['mismatches']}" if a.check else ""), flush=True)

    print("\nby stimulus:")
    for stim, s in by_stim.items():
        print(f"  {stim:15s} runs {int(s['runs']):4d}  sim {s['sim_ms']:9.1f} ms  wall {s['wall_s']:7.2f} s"
              + (f"  mismatches {int(s['mismatches'])}" if a.check else ""))
    total = {k: sum(r[k] for r in rows) for k in ("runs", "sim_ms", "wall_s", "mismatches")}
    print(f"  {'total':15s} runs {total['runs']:4d}  sim {total['sim_ms']:9.1f} ms  "
          f"wall {total['wall_s']:7.2f} s" + (f"  mismatches {total['mismatches']}" if a.check else ""))
    method = "full runs" if a.no_settled_stop else "settled stop"
    tag = method.replace(" ", "-") + f"_tol{a.rel_tol:g}"
    (out / f"{tag}.json").write_text(json.dumps(dict(method=method, rel_tol=a.rel_tol, rows=rows,
                                                     total=total), indent=1))
    print(f"saved {out / (tag + '.json')}")


if __name__ == "__main__":
    main()
