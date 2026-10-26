"""Window and threshold over the two-compartment coupling plane (issue #84).

    uv run scripts/coupling_plane.py --out-dir figures/coupling_plane

Goldwyn, Remme & Rinzel 2019 describe a two-compartment neuron by its forward
(soma to axon, couple12) and backward (axon to soma, couple21) coupling
constants. Their Fig 3 plane has forward coupling no weaker than backward
("the upper left half of the parameter space is empty") and, in Methods, a
ratio couple12 / couple21 of at most ten. This sweeps that region on a regular
grid and, at each point, finds the EPSG-pair threshold (coincident inputs) and
the 3% coincidence window (coincidence.threshold and coincidence.window, as in
scripts/coincidence_window.py), and says whether the cell fires at all. The
area ratio, the channel fractions and the axonal Na conductance stay at node
3's values, so unlike Goldwyn et al, who re-reference g_Na for every
configuration, a weakly coupled cell may simply not fire.

It marks the lookup's own point (node 3, the TwoCpt.m cell), Goldwyn's three
example configurations and the 44 axon nodes' constants, and locates the
45-compartment model: the point whose two-compartment window and threshold
lie closest to the 45-compartment values (relative distance).
"""

import argparse
import json
import math
import pathlib
from typing import Any

import msoaxon
from msoaxon import _bisect, coincidence, constants

MARGIN = 0.03  # inputs this far above the coincident threshold (Myoga: ~3%)
# the 45-compartment model, 10 mV criterion, model EPSG (PYTHON.md, "Comparison
# with the literature": 201.1 us and 48.01 nS, the "about 200 us, 48 nS")
MULTI_WINDOW_US = 201.1
MULTI_THRESHOLD_NS = 48.0
MAX_RATIO = 10.0  # Goldwyn 2019 Methods: couple12 / couple21 does not exceed it
AXIS_LO, AXIS_HI = 0.1, 0.9
DEFAULT_STEP, QUICK_STEP = 0.1, 0.4
DEFAULT_CEILING, QUICK_CEILING = 2000.0, 300.0  # [nS]
NODE = 3  # the lookup node whose constants are the default cell
# Goldwyn 2019 Fig 3's stars; the lookup's point is added in landmarks()
STARS = {
    "weak (Goldwyn)": (0.3, 0.2),
    "forward (Goldwyn)": (0.8, 0.2),
    "strong (Goldwyn)": (0.8, 0.7),
}
STATUSES = ("fires", "silent", "spontaneous", "wide")


def axis(step: float, lo: float = AXIS_LO, hi: float = AXIS_HI) -> list[float]:
    """Returns the coupling values lo, lo + step, ... up to hi.

    Args:
        step: The spacing, > 0.
        lo: The first value, in (0, 1).
        hi: The largest value allowed, in [lo, 1).

    Returns:
        The values, rounded to 6 places so that 0.1 + 0.2 is 0.3.

    Raises:
        ValueError: If step is not > 0 or lo and hi are not 0 < lo <= hi < 1.
    """
    if not step > 0:
        raise ValueError(f"step must be > 0, got {step}")
    if not 0 < lo <= hi < 1:
        raise ValueError(f"need 0 < lo <= hi < 1, got lo={lo}, hi={hi}")
    n = int(math.floor((hi - lo) / step + 1e-9))
    return [round(lo + k * step, 6) for k in range(n + 1)]


def grid(
    step: float,
    lo: float = AXIS_LO,
    hi: float = AXIS_HI,
    max_ratio: float = MAX_RATIO,
) -> list[tuple[float, float]]:
    """Returns the (couple12, couple21) points of Goldwyn's plane.

    Both axes run over axis(step, lo, hi). A point is kept if the forward
    constant is not below the backward one and their ratio is at most
    max_ratio.

    Args:
        step: The spacing of both axes.
        lo: The first value of both axes.
        hi: The largest value of both axes.
        max_ratio: The largest couple12 / couple21.

    Returns:
        The points, forward constant first, then backward, ascending.

    Raises:
        ValueError: If axis() rejects its arguments or max_ratio is < 1.
    """
    if not max_ratio >= 1:
        raise ValueError(f"max_ratio must be >= 1, got {max_ratio}")
    values = axis(step, lo, hi)
    return [
        (c12, c21)
        for c12 in values
        for c21 in values
        if c12 >= c21 and c12 / c21 <= max_ratio * (1 + 1e-9)
    ]


def landmarks() -> dict[str, tuple[float, float]]:
    """Returns the named points: the lookup's own, then Goldwyn's three."""
    default = (
        float(constants.COUPLING1[NODE - 2]),
        float(constants.COUPLING2[NODE - 2]),
    )
    return {f"node {NODE} (TwoCpt.m)": default, **STARS}


def distance(
    window_us: float,
    threshold: float,
    ref_window_us: float = MULTI_WINDOW_US,
    ref_threshold: float = MULTI_THRESHOLD_NS,
) -> float:
    """Returns the relative distance of (window, threshold) to a reference.

    The two quantities have different units, so each is taken as a fraction
    of its reference: sqrt(dw^2 + dth^2).

    Args:
        window_us: The window [us].
        threshold: The threshold [nS].
        ref_window_us: The reference window [us].
        ref_threshold: The reference threshold [nS].

    Returns:
        The distance, or nan if either value is not finite.
    """
    if not (math.isfinite(window_us) and math.isfinite(threshold)):
        return math.nan
    return math.hypot(
        (window_us - ref_window_us) / ref_window_us,
        (threshold - ref_threshold) / ref_threshold,
    )


def locate(
    points: list[dict[str, Any]],
    ref_window_us: float = MULTI_WINDOW_US,
    ref_threshold: float = MULTI_THRESHOLD_NS,
) -> tuple[dict[str, Any], float]:
    """Finds the point closest to a reference (window, threshold).

    Args:
        points: Dicts with at least "window_us" and "threshold" (nan, inf or
            None where there is no value); such points are never closest.
        ref_window_us: The reference window [us].
        ref_threshold: The reference threshold [nS].

    Returns:
        (point, distance): the first of the closest points, and its
        distance (see distance()).

    Raises:
        ValueError: If no point has a finite window and threshold.
    """
    best, best_d = None, math.inf
    for p in points:
        w, th = p["window_us"], p["threshold"]
        d = (
            math.nan
            if w is None or th is None
            else distance(w, th, ref_window_us, ref_threshold)
        )
        if math.isfinite(d) and d < best_d:
            best, best_d = p, d
    if best is None:
        raise ValueError("no point has a finite window and threshold")
    return best, best_d


def measure_point(task: tuple[Any, ...]) -> dict[str, Any]:
    """Measures the threshold and the window of one coupling configuration.

    Args:
        task: (couple12, couple21, mode, ceiling, rel_tol, delay_tol).

    Returns:
        dict(couple12, couple21, status, threshold, window_us). The status
        is "fires" (both found), "silent" (nothing up to the ceiling
        spikes: threshold inf), "spontaneous" (the cell fires with no input)
        or "wide" (the input still spikes at the search's largest delay, so
        no window). Values that were not found are nan (threshold: inf for
        "silent").
    """
    c12, c21, mode, ceiling, rel_tol, delay_tol = task
    kw = {"couple12": c12, "couple21": c21, "mode": mode}
    out = {"couple12": c12, "couple21": c21}
    try:
        th = coincidence.threshold(
            "two", 0.0, model_kw=kw, rel_tol=rel_tol, ceiling=ceiling
        )
    except _bisect.FiresWithoutInputError:
        return {
            **out,
            "status": "spontaneous",
            "threshold": math.nan,
            "window_us": math.nan,
        }
    if not math.isfinite(th):
        return {
            **out,
            "status": "silent",
            "threshold": th,
            "window_us": math.nan,
        }
    width, _ = coincidence.window(
        "two",
        MARGIN,
        model_kw=kw,
        rel_tol=rel_tol,
        delay_tol=delay_tol,
        threshold0=th,
    )
    return {
        **out,
        "status": "fires" if math.isfinite(width) else "wide",
        "threshold": th,
        "window_us": 1000 * width,
    }


def sweep(
    points: list[tuple[float, float]],
    mode: str | None,
    ceiling: float,
    rel_tol: float,
    delay_tol: float,
    pool: Any,
) -> list[dict[str, Any]]:
    """Measures every point in parallel, returning the results in order."""
    tasks = [
        (c12, c21, mode, ceiling, rel_tol, delay_tol) for c12, c21 in points
    ]
    return msoaxon.map_tasks(measure_point, tasks, executor=pool)


def _clean(x: Any) -> Any:
    """Returns x with nan and inf as None, so the JSON is strict."""
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_clean(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def _fmt(x: float | None, spec: str) -> str:
    return "-" if x is None or not math.isfinite(x) else format(x, spec)


def markdown(points: list[dict[str, Any]]) -> str:
    """Returns the points as a Markdown table, by forward then backward."""
    rows = [
        "| couple12 | couple21 | point | status | threshold (nS) "
        "| 3% window (us) | distance to 45-cpt |",
        "|---|---|---|---|---|---|---|",
    ]
    for p in sorted(points, key=lambda p: (p["couple12"], p["couple21"])):
        d = distance(p["window_us"], p["threshold"])
        rows.append(
            f"| {p['couple12']:.4g} | {p['couple21']:.4g} "
            f"| {p.get('label', '')} | {p['status']} "
            f"| {_fmt(p['threshold'], '.2f')} | {_fmt(p['window_us'], '.1f')} "
            f"| {_fmt(d, '.3f')} |"
        )
    return "\n".join(rows)


def extremes(points: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Returns the points of smallest and largest window and threshold.

    Only points with a finite value count.

    Args:
        points: Measured points (see measure_point()).

    Returns:
        {"window_min": p, "window_max": p, "threshold_min": p,
        "threshold_max": p}; a key is missing if no point has the value.
    """
    out = {}
    for key, name in (("window_us", "window"), ("threshold", "threshold")):
        have = [
            p for p in points if p[key] is not None and math.isfinite(p[key])
        ]
        if have:
            out[f"{name}_min"] = min(have, key=lambda p: p[key])
            out[f"{name}_max"] = max(have, key=lambda p: p[key])
    return out


def plot(
    points: list[dict[str, Any]],
    step: float,
    located: dict[str, Any],
    path: pathlib.Path,
) -> None:
    """Draws the window and the threshold over the plane.

    Grid points are coloured squares (grey with a cross where nothing was
    found); the lookup's point, Goldwyn's stars and the located point are
    marked, and the grey dots are the 44 axon nodes' own constants.
    """
    # imported here rather than at the top, like the pool scripts, so importing
    # this module stays cheap
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), sharey=True)
    names = landmarks()
    for ax, key, title in (
        (axes[0], "window_us", "3% coincidence window (us)"),
        (axes[1], "threshold", "Coincident threshold (nS)"),
    ):
        found = [
            p for p in points if p[key] is not None and math.isfinite(p[key])
        ]
        lost = [p for p in points if p not in found]
        vals = [p[key] for p in found]
        sc = ax.scatter(
            [p["couple12"] for p in found],
            [p["couple21"] for p in found],
            c=vals,
            s=900 * step if step < 0.2 else 700,
            marker="s",
            cmap="viridis",
        )
        ax.scatter(
            [p["couple12"] for p in lost],
            [p["couple21"] for p in lost],
            c="0.85",
            s=900 * step if step < 0.2 else 700,
            marker="s",
        )
        ax.scatter(
            [p["couple12"] for p in lost],
            [p["couple21"] for p in lost],
            c="0.4",
            marker="x",
            s=25,
        )
        ax.scatter(
            constants.COUPLING1,
            constants.COUPLING2,
            s=6,
            c="0.55",
            label="the 44 axon nodes' constants",
        )
        for (name, (c12, c21)), colour, marker in zip(
            names.items(),
            ("k", "tab:blue", "tab:green", "tab:red"),
            ("D", "*", "*", "*"),
            strict=True,
        ):
            ax.scatter(
                [c12],
                [c21],
                c=colour,
                marker=marker,
                s=120,
                edgecolors="w",
                label=name,
                zorder=5,
            )
        ax.scatter(
            [located["couple12"]],
            [located["couple21"]],
            facecolors="none",
            edgecolors="tab:orange",
            s=380,
            linewidths=2.2,
            zorder=6,
            label="closest to the 45-compartment model",
        )
        ax.plot([0, 1], [0, 1], "k:", lw=0.8)  # forward = backward
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("forward coupling couple12 (soma to axon)")
        ax.set_title(title)
        fig.colorbar(sc, ax=ax)
    axes[0].set_ylabel("backward coupling couple21 (axon to soma)")
    axes[0].legend(fontsize=7, loc="upper left", frameon=True)
    fig.suptitle(
        "Two-compartment model over Goldwyn et al 2019's coupling plane "
        "(node 3's area ratio, channels and axonal Na); grey: no spike "
        "or no window",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main() -> None:
    """Sweeps the coupling plane and saves the table, figure and JSON."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", default="figures/coupling_plane")
    ap.add_argument(
        "--step",
        type=float,
        help=f"grid spacing of both axes (default {DEFAULT_STEP}; "
        f"--quick {QUICK_STEP})",
    )
    ap.add_argument(
        "--ceiling",
        type=float,
        help="largest EPSG amplitude tried [nS]; a cell with no spike up to "
        f"it is silent (default {DEFAULT_CEILING:.0f}; --quick "
        f"{QUICK_CEILING:.0f})",
    )
    ap.add_argument(
        "--mode",
        choices=("exact", "fast", "native"),
        help="model mode (default: fast); the two-compartment model has "
        "only its numpy path, so the mode changes nothing here",
    )
    ap.add_argument("--workers", type=int, help="processes (default: all)")
    ap.add_argument(
        "--quick",
        action="store_true",
        help="three values per axis and a 300 nS ceiling: "
        "a smoke test, not results",
    )
    a = ap.parse_args()
    step = a.step or (QUICK_STEP if a.quick else DEFAULT_STEP)
    ceiling = a.ceiling or (QUICK_CEILING if a.quick else DEFAULT_CEILING)
    # the window is steep in the amplitude (0.5% above threshold it is a
    # third of the 3% one), so the tolerances stay tight even for --quick
    rel_tol, delay_tol = 1e-4, 1e-4
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    names = landmarks()
    points = grid(step)
    points += [c for c in names.values() if c not in points]
    with msoaxon.process_pool(a.workers) as pool:
        results = sweep(points, a.mode, ceiling, rel_tol, delay_tol, pool)
    by_coords = {v: k for k, v in names.items()}
    for p in results:
        p["label"] = by_coords.get((p["couple12"], p["couple21"]), "")
        p["distance"] = distance(p["window_us"], p["threshold"])
    best, best_d = locate(results)
    default = next(p for p in results if p["label"] == next(iter(names)))

    print(markdown(results))
    print()
    for key, p in extremes(results).items():
        value = p["window_us"] if key.startswith("window") else p["threshold"]
        print(
            f"{key}: {_fmt(value, '.2f')}"
            f" at couple12={p['couple12']:.4g}, couple21={p['couple21']:.4g}"
        )
    print(
        f"default point (node {NODE}): couple12={default['couple12']:.4f}, "
        f"couple21={default['couple21']:.4f}, "
        f"{_fmt(default['window_us'], '.1f')} us, "
        f"{_fmt(default['threshold'], '.2f')} nS, "
        f"distance {_fmt(default['distance'], '.3f')}"
    )
    print(
        f"closest to the 45-compartment model ({MULTI_WINDOW_US} us, "
        f"{MULTI_THRESHOLD_NS} nS): couple12={best['couple12']:.4g}, "
        f"couple21={best['couple21']:.4g}, {best['window_us']:.1f} us, "
        f"{best['threshold']:.2f} nS, distance {best_d:.3f}"
    )
    plot(results, step, best, out / "coupling_plane.png")
    (out / "coupling_plane.json").write_text(
        json.dumps(
            _clean(
                {
                    "mode": a.mode or "default",
                    "quick": a.quick,
                    "margin": MARGIN,
                    "step": step,
                    "ceiling_nS": ceiling,
                    "reference": {
                        "window_us": MULTI_WINDOW_US,
                        "threshold_nS": MULTI_THRESHOLD_NS,
                    },
                    "points": results,
                    "located": {**best, "distance": best_d},
                    "default": default,
                    "extremes": extremes(results),
                }
            ),
            indent=1,
        )
    )
    print(f"saved {out}/coupling_plane.png and .json")


if __name__ == "__main__":
    main()
