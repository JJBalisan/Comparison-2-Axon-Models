"""Where spikes start in the 45-compartment model (issue #80).

    uv run scripts/initiation_site.py --json initiation.json

Sweeps the sustained Synaptic input over input frequency and EPSG
conductance, pooling several noise seeds per point, and prints the share of
spikes starting in the AIS, at node 1 and at node 2 against Lehnert et al
2014's 5-60% at nodes 1 and 2. Then reads the EPSG-pair protocol at
multiples of the coincident threshold. See msoaxon.initiation for the
criterion.
"""

import argparse
import itertools
import json
from typing import Any

import msoaxon
from msoaxon import coincidence, initiation, multi

FREQS = (200.0, 350.0, 500.0, 800.0)  # input frequency [Hz]
GES = (0.65e7, 1.3e7, 2.6e7, 5.2e7, 1.04e8)  # 1.3e7 is SynParams' default
SEEDS = (1396, 1, 2)  # noise seeds; 1396 is SynParams' default
T_END = 200.0  # [ms] per run
MULTS = (1.05, 1.5, 2.0, 3.0)  # pair amplitude / coincident threshold


def sustained(task: tuple[float, float, int, float]) -> dict[str, Any]:
    """Runs one point of the sweep (a pool worker).

    Args:
        task: (freq, gE, seed, t_end).

    Returns:
        initiation.Fractions as a dict.
    """
    freq, g_e, seed, t_end = task
    f = initiation.sustained_input(freq, g_e, t_end, seed)
    return {
        "n_spikes": f.n_spikes,
        "ais": f.ais,
        "node1": f.node1,
        "node2": f.node2,
        "other": f.other,
    }


def pooled(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Pools the seeds of one point: counts add, shares follow.

    Args:
        rows: One sustained() result per seed.

    Returns:
        The same keys as a row, over all their spikes.
    """
    n = sum(r["n_spikes"] for r in rows)
    out = {"n_spikes": n}
    for key in ("ais", "node1", "node2", "other"):
        count = sum(r["n_spikes"] * r[key] for r in rows if r["n_spikes"])
        out[key] = count / n if n else float("nan")
    return out


def pair_sites() -> list[dict[str, Any]]:
    """Reads the EPSG-pair protocol at multiples of the coincident threshold.

    Returns:
        One dict per multiple: the threshold, the multiple, the site and
        its 1-indexed compartment.
    """
    th = coincidence.threshold("multi", 0.0)
    rows = []
    for mult in MULTS:
        t, x = multi.mso_axon(
            "EPSGpair", 5.0, 5.0, mult * th, 3, "active-full", 9.0, -68.0, 1
        )
        (spike,) = initiation.find_spikes(t, x)
        rows.append(
            {
                "threshold_nS": th,
                "mult": mult,
                "compartment": spike.compartment,
                "label": multi.LUMPED.labels[spike.compartment - 1],
            }
        )
    return rows


def main() -> None:
    """Runs the sweep and prints (and optionally saves) the results."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--json", help="save the results here")
    ap.add_argument("--freqs", type=float, nargs="+", default=FREQS, help="Hz")
    ap.add_argument("--ges", type=float, nargs="+", default=GES, help="EPSG")
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    ap.add_argument("--t-end", type=float, default=T_END, help="ms per run")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument(
        "--quick",
        action="store_true",
        help="350 Hz, 2.6e7, one seed, 50 ms: a smoke test, not results",
    )
    a = ap.parse_args()
    freqs, ges, seeds, t_end = a.freqs, a.ges, a.seeds, a.t_end
    if a.quick:
        freqs, ges, seeds, t_end = (350.0,), (2.6e7,), seeds[:1], 50.0

    grid = list(itertools.product(freqs, ges, seeds))
    tasks = [(f, g, s, t_end) for f, g, s in grid]
    with msoaxon.process_pool(a.workers) as pool:
        rows = msoaxon.map_tasks(sustained, tasks, executor=pool)
    points = {}
    for (f, g, _), row in zip(grid, rows, strict=True):
        points.setdefault((f, g), []).append(row)

    print(f"sustained Synaptic input, {t_end:g} ms x {len(seeds)} seeds")
    print("freq Hz      gE  spikes   AIS  node1  node2  other  node1+2")
    results = []
    for (f, g), group in points.items():
        p = pooled(group)
        print(
            f"{f:7g} {g:7.2e} {p['n_spikes']:6d} {p['ais']:5.2f} "
            f"{p['node1']:6.2f} {p['node2']:6.2f} {p['other']:6.2f} "
            f"{p['node1'] + p['node2']:8.2f}"
        )
        results.append({"freq": f, "gE": g, **p})
    pair = pair_sites() if not a.quick else []
    for r in pair:
        print(
            f"pair {r['mult']:g}x threshold ({r['threshold_nS']:.1f} nS): "
            f"compartment {r['compartment']} ({r['label']})"
        )
    if a.json:
        with open(a.json, "w") as f:
            json.dump({"sustained": results, "pair": pair}, f, indent=1)


if __name__ == "__main__":
    main()
