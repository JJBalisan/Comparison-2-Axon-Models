"""Spiking-threshold curves for both models (port of Making_Threshold_graphs.m).

    uv run scripts/making_threshold_graphs.py --node 3 --out-dir figures
    uv run scripts/making_threshold_graphs.py --only EPSGpair --jpg-grid   # the repo's jpgs

Each binary-search point runs ~8 simulations per model; points run in parallel
across all CPUs unless --workers says otherwise.
"""  # noqa: E501 (printed verbatim as --help)

import argparse
import pathlib
from concurrent import futures

import numpy as np

from msoaxon import _parallel, threshold

# stimType: (n_points, factor, max, title, xlabel); x-values come from
# threshold.sweep_x
SWEEPS = {
    "ramp": (15, 30, 15000, "Ramp", "Stop Value"),
    "sine": (12, 30, 15000, "Sine", "Frequency"),
    "ramp2": (15, 30, 15000, "Ramp2", "Stop Value"),
    "EPSGpair": (25, 10, 150, "EPSGpair", "Time Difference"),
}


def epsgpair_grid(jpg_grid: bool) -> tuple[int, dict[str, float | bool]]:
    """Returns the EPSGpair sweep's grid.

    tests/test_coverage.py pins the jpg grid's thresholds through this.

    Args:
        jpg_grid: Use the grid the repo's EPSGpair_Thresholds*.jpg were
            made on (11 points 0.1 ms apart, unrounded) instead of
            BinarySearch.m's.

    Returns:
        (n_points, keywords), the keywords being binary_search's
        epsg_pair_dt and rounded.
    """
    if jpg_grid:
        return 11, {"epsg_pair_dt": 0.1, "rounded": False}
    return SWEEPS["EPSGpair"][0], {"epsg_pair_dt": 1 / 25, "rounded": True}


def main() -> None:
    """Parses the command line and runs the threshold sweeps."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--node", type=int, default=3)
    ap.add_argument(
        "--only", choices=SWEEPS, action="append", help="run just these sweeps"
    )
    ap.add_argument(
        "--jpg-grid",
        action="store_true",
        help="EPSGpair on the grid the repo's jpgs used: 11 points 0.1 ms "
        "apart, unrounded",
    )
    ap.add_argument("--out-dir", help="save PNGs here instead of showing them")
    ap.add_argument(
        "--workers",
        type=int,
        help="processes for sweep points (default: all CPUs)",
    )
    ap.add_argument(
        "--points", type=int, help="run only the first N points of each sweep"
    )
    a = ap.parse_args()
    if a.workers == 1:  # in-process, as binary_search(workers=1) does
        run(a, None)
    else:
        with _parallel.process_pool(a.workers) as pool:  # shared by every sweep
            run(a, pool)


def run(
    a: argparse.Namespace, pool: futures.ProcessPoolExecutor | None
) -> None:
    # imported here rather than at the top: pool workers re-import this script,
    # and only the main process draws (msoaxon.plotting imports matplotlib)
    """Runs the threshold sweeps and draws or saves one figure per sweep.

    Args:
        a: The parsed command-line arguments.
        pool: The process pool every sweep shares, or None to run
            in-process.
    """
    import matplotlib.pyplot as plt

    from msoaxon import plotting

    stims = a.only or list(SWEEPS)

    def sweep(stim: str) -> tuple[np.ndarray, list[float], list[float]]:
        n, factor, max_I, _, _ = SWEEPS[stim]
        kw, dt = {}, 1 / 25
        if stim == "EPSGpair":
            n, kw = epsgpair_grid(a.jpg_grid)
            dt = kw["epsg_pair_dt"]
        n = min(n, a.points) if a.points else n
        x_values = np.array(
            [threshold.sweep_x(stim, i, dt) for i in range(1, n + 1)]
        )
        multi, two = threshold.binary_search(
            stim,
            n,
            a.node,
            factor,
            max_I,
            workers=a.workers,
            executor=pool,
            **kw,
        )
        return x_values, multi, two

    # binary_search maps its points onto the pool and waits for them, so called
    # one sweep at a time it would drain the pool between sweeps. From one
    # thread per sweep, every sweep's points are queued at once (submitting to
    # the pool is thread-safe) and the pool stays busy until the last point;
    # the figures are still drawn here, in order. Without a pool the sweeps run
    # in-process, one after another.
    if pool is None:
        sweeps = {stim: sweep(stim) for stim in stims}
    else:
        with futures.ThreadPoolExecutor(len(stims)) as threads:
            jobs = {stim: threads.submit(sweep, stim) for stim in stims}
            sweeps = {stim: job.result() for stim, job in jobs.items()}
    for stim in stims:
        _, _, _, title, xlabel = SWEEPS[stim]
        x_values, multi, two = sweeps[stim]
        print(f"{stim}: multi={multi}\n{' ' * len(stim)}  two  ={two}")
        fig = plotting.threshold_figure(
            x_values,
            multi,
            two,
            f"{title} Thresholds Compartment {a.node}",
            xlabel,
        )
        if a.out_dir:
            pathlib.Path(a.out_dir).mkdir(parents=True, exist_ok=True)
            path = (
                pathlib.Path(a.out_dir) / f"{stim}_thresholds_cpt{a.node}.png"
            )
            fig.savefig(path, dpi=150)
            print(f"saved {path}")
    if not a.out_dir:
        plt.show()


if __name__ == "__main__":
    main()
