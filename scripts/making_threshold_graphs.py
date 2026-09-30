"""Spiking-threshold curves for both models (port of Making_Threshold_graphs.m).

    uv run scripts/making_threshold_graphs.py --node 3 --out-dir figures
    uv run scripts/making_threshold_graphs.py --only EPSGpair --jpg-grid   # the repo's jpgs

Each binary-search point runs ~8 simulations per model; points run in parallel
across all CPUs unless --workers says otherwise.
"""

import argparse
from pathlib import Path

import numpy as np

from msoaxon._parallel import process_pool
from msoaxon.threshold import binary_search, sweep_x

# stimType: (n_points, factor, max, title, xlabel); x-values come from threshold.sweep_x
SWEEPS = {
    "ramp": (15, 30, 15000, "Ramp", "Stop Value"),
    "sine": (12, 30, 15000, "Sine", "Frequency"),
    "ramp2": (15, 30, 15000, "Ramp2", "Stop Value"),
    "EPSGpair": (25, 10, 150, "EPSGpair", "Time Difference"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--node", type=int, default=3)
    ap.add_argument("--only", choices=SWEEPS, action="append", help="run just these sweeps")
    ap.add_argument("--jpg-grid", action="store_true",
                    help="EPSGpair on the grid the repo's jpgs used: 11 points 0.1 ms apart, unrounded")
    ap.add_argument("--out-dir", help="save PNGs here instead of showing them")
    ap.add_argument("--workers", type=int, help="processes for sweep points (default: all CPUs)")
    ap.add_argument("--points", type=int, help="run only the first N points of each sweep")
    a = ap.parse_args()
    if a.workers == 1:  # in-process, as binary_search(workers=1) does
        run(a, None)
    else:
        with process_pool(a.workers) as pool:  # shared by every sweep
            run(a, pool)


def run(a, pool):
    # imported here rather than at the top: pool workers re-import this script,
    # and only the main process draws (msoaxon.plotting imports matplotlib)
    import matplotlib.pyplot as plt

    from msoaxon.plotting import threshold_figure

    for stim in a.only or SWEEPS:
        n, factor, max_I, title, xlabel = SWEEPS[stim]
        kw, dt = {}, 1 / 25
        if stim == "EPSGpair":
            dt = 0.1 if a.jpg_grid else dt
            n = 11 if a.jpg_grid else n
            kw = dict(epsg_pair_dt=dt, rounded=not a.jpg_grid)
        n = min(n, a.points) if a.points else n
        x_values = np.array([sweep_x(stim, i, dt) for i in range(1, n + 1)])
        multi, two = binary_search(stim, n, a.node, factor, max_I, workers=a.workers,
                                   executor=pool, **kw)
        print(f"{stim}: multi={multi}\n{' ' * len(stim)}  two  ={two}")
        fig = threshold_figure(x_values, multi, two, f"{title} Thresholds Compartment {a.node}", xlabel)
        if a.out_dir:
            Path(a.out_dir).mkdir(parents=True, exist_ok=True)
            path = Path(a.out_dir) / f"{stim}_thresholds_cpt{a.node}.png"
            fig.savefig(path, dpi=150)
            print(f"saved {path}")
    if not a.out_dir:
        plt.show()


if __name__ == "__main__":
    main()
