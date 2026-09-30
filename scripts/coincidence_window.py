"""Coincidence-detection windows vs Myoga et al 2014, with and without recalibration.

    uv run scripts/coincidence_window.py --out-dir figures/coincidence

Runs three models (45-compartment; two-compartment as in TwoCpt.m; two-compartment
recalibrated to Goldwyn et al 2019's passive targets) with the model's EPSG
(decay 0.18 ms) and Myoga's (0.3 ms), computes each spike-probability half-width,
and checks the threshold shortcut against noisy trials.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from msoaxon._parallel import process_pool
from msoaxon.coincidence import half_width, probability_trials, threshold_curve
from msoaxon.synaptic import EPSG_TAU
from msoaxon.two import GOLDWYN_2019

MYOGA_US = 221  # AP-probability half-width without inhibition (Myoga et al 2014)
MARGINS = (0.005, 0.03)  # inputs this far above coincident threshold (Myoga: "200 pS (~3%)")
KINETICS = {"model EPSG (decay 0.18 ms)": EPSG_TAU, "Myoga EPSG (decay 0.3 ms)": (0.1, 0.3)}
_goldwyn = {k: v for k, v in GOLDWYN_2019.items() if k != "v0"}
CONFIGS = {
    "45-compartment": ("multi", -68.0, {}),
    "2-cpt, TwoCpt.m (10 MOhm, 0.71 ms, -68 mV)": ("two", -68.0, {}),
    "2-cpt, Goldwyn 2019 (8.5 MOhm, 0.34 ms, -58 mV)": ("two", GOLDWYN_2019["v0"], _goldwyn),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="figures/coincidence")
    ap.add_argument("--node", type=int, default=3)
    ap.add_argument("--trials", type=int, default=200, help="noisy trials per delay for the check")
    ap.add_argument("--quick", action="store_true",
                    help="coarse grids and 1%% threshold tolerance: a smoke test, not results")
    a = ap.parse_args()
    with process_pool() as pool:  # one pool for every parallel call below
        run(a, pool)


def run(a, pool):
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    tol = 1e-2 if a.quick else 1e-4  # threshold tolerance
    delays = np.round(np.arange(0, 1.0001, 0.1 if a.quick else 0.02), 4)
    curves, rows = {}, []
    for cname, (model, v0, kw) in CONFIGS.items():
        for kname, tau in KINETICS.items():
            th = threshold_curve(model, delays, node=a.node, v0=v0, epsg_tau=tau, model_kw=kw,
                                 rel_tol=tol, executor=pool)
            curves[(cname, kname)] = th
            widths = [half_width(delays, th, m) * 1e3 for m in MARGINS]
            rows.append((cname, kname, th[0], widths))
            print(f"{cname} | {kname}: threshold(0) {th[0]:.2f}, plateau/zero-delay {th[-1] / th[0]:.2f}, "
                  + ", ".join(f"half-width @{m:.1%} {w:.0f} us" for m, w in zip(MARGINS, widths)), flush=True)

    # check the shortcut with real noisy trials (2-cpt TwoCpt.m, model EPSG, 3% margin)
    cname, kname = list(CONFIGS)[1], list(KINETICS)[0]
    model, v0, kw = CONFIGS[cname]
    th = curves[(cname, kname)]
    amp = (1 + 0.03) * th[0]
    mc_delays = np.round(np.arange(0, 0.3001, 0.1 if a.quick else 0.02), 4)
    # small noise so probability peaks near 100%, as in Myoga's protocol; 1% amplitude
    # jitter plus 5 us onset jitter per EPSG (a noise source the shortcut ignores)
    prob = probability_trials(model, mc_delays, amp, n_trials=a.trials, amp_cv=0.01,
                              jitter=0.005, node=a.node, v0=v0, epsg_tau=KINETICS[kname],
                              model_kw=kw, executor=pool)
    p_half = prob.max() / 2
    below = np.nonzero(prob < p_half)[0]
    if len(below) == 0 or below[0] == 0:  # never falls to half within mc_delays, or never rises
        mc_width = np.nan
    else:
        j = below[0]
        mc_width = 2 * np.interp(p_half, [prob[j], prob[j - 1]], [mc_delays[j], mc_delays[j - 1]]) * 1e3
    shortcut = half_width(delays, th, 0.03) * 1e3
    print(f"noisy-trial check ({cname}, {kname}, 3%): peak probability {prob.max():.2f}, "
          f"half-width {mc_width:.0f} us vs threshold shortcut {shortcut:.0f} us", flush=True)

    # figure: normalised threshold curves, and probability curves
    # imported here rather than at the top: pool workers re-import this script,
    # and only the main process draws
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6.2))
    colours = ["#1f5fa8", "#d2691e", "#2e8b57"]
    for c, cname in zip(colours, CONFIGS):
        for kname, ls in zip(KINETICS, ("-", "--")):
            th = curves[(cname, kname)]
            a1.plot(delays * 1e3, th / th[0], color=c, ls=ls, lw=2,
                    label=f"{cname.split(' (')[0]}, {kname.split(' (')[0]}")
    for m in MARGINS:
        a1.axhline(1 + m, color="0.5", lw=0.8, ls=":")
    a1.set(xlim=(0, 600), ylim=(0.98, 2.05), xlabel="Delay between EPSGs (us)",
           ylabel="Threshold / coincident threshold",
           title="Threshold rises with delay (dotted: 0.5% and 3% input margins)")
    a1.legend(fontsize=7.5, frameon=False, loc="lower right")

    sym = np.concatenate([-mc_delays[:0:-1], mc_delays]) * 1e3
    a2.plot(sym, np.concatenate([prob[:0:-1], prob]), "o-", color=colours[1], lw=1.5,
            label=f"2-cpt TwoCpt.m, noisy trials (FWHM {mc_width:.0f} us)")
    for c, cname, y in zip(colours, CONFIGS, (0.56, 0.5, 0.44)):  # offset so all show
        th = curves[(cname, list(KINETICS)[0])]
        w = half_width(delays, th, 0.03) * 1e3
        a2.plot([-w / 2, w / 2], [y, y], color=c, lw=4, alpha=0.8,
                label=f"{cname.split(' (')[0]}: {w:.0f} us (shortcut, 3%)")
    a2.axvspan(-MYOGA_US / 2, MYOGA_US / 2, color="0.85", zorder=0,
               label=f"Myoga 2014 in vitro: {MYOGA_US} us")
    a2.set(xlim=(-320, 320), ylim=(-0.03, 1.05), xlabel="Delay between EPSGs (us)",
           ylabel="Spike probability",
           title="Coincidence window at half maximum (model EPSG; bars offset around 0.5)")
    a2.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    fig.tight_layout()
    fig.savefig(out / "coincidence_window.png", dpi=140)

    summary = {"delays_ms": delays.tolist(), "margins": MARGINS, "myoga_half_width_us": MYOGA_US,
               "curves": {f"{c} | {k}": v.tolist() for (c, k), v in curves.items()},
               "half_widths_us": [{"model": c, "epsg": k, "threshold0": t0,
                                   **{f"margin_{m}": w for m, w in zip(MARGINS, ws)}}
                                  for c, k, t0, ws in rows],
               "noisy_check": {"half_width_us": mc_width, "shortcut_us": shortcut,
                               "peak_probability": float(prob.max()),
                               "delays_ms": mc_delays.tolist(), "probability": prob.tolist()}}
    (out / "coincidence_window.json").write_text(json.dumps(summary, indent=1))
    print(f"saved {out}/coincidence_window.png and .json")


if __name__ == "__main__":
    main()
