"""The dendritic 45-compartment variant vs the lumped model and the literature.

    uv run scripts/dendrites.py --out-dir figures/dendrites
    uv run scripts/dendrites.py --klt-from-soma --out-dir figures/dendrites-klt-from-soma
    uv run scripts/dendrites.py --dendrite-ra 200 --out-dir figures/dendrites-ra200

Uses with_dendrites() (Lehnert et al 2014, Fig. 8 variant: two 200 x 5 um dendrites,
5 compartments each, soma reduced to 2467 um^2). Checks, in order:
  A. soma input resistance and time constant (lumped vs dendritic)
  B. somatic EPSP from a unitary EPSG at the soma, mid-dendrite and distal dendrite,
     with dendritic KLT on and off (Mathews et al 2010: dendritic Kv1 sharpens EPSPs)
  C. zero-delay threshold: one EPSG on each dendrite vs both on one dendrite
     (Scott et al 2010 predicts bilateral < unilateral)
  D. subthreshold summation: bilateral EPSP vs the sum of the two unilateral ones
     (linear in vivo: van der Heijden et al 2013; Mackenbach & Borst 2023)
  E. Myoga et al 2014-matched coincidence window with bilateral dendritic inputs
  F. EPSG-pair threshold vs delay, lumped vs dendritic (Lehnert: "almost identical")
  G. somatic spike, Scott et al 2005 protocol
A, B and D use rest-rebalanced membranes so sub-mV signals aren't mixed with the
slow resting drift; C, E and F use the default (unbalanced) membranes, like the
lumped results they are compared with.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from msoaxon import mso_axon
from msoaxon._parallel import process_pool
from msoaxon.coincidence import half_width, threshold
from msoaxon.measure import passive_step, soma_on_grid
from msoaxon.multi import LUMPED, membrane, with_dendrites
from msoaxon.somatic import rheobase, spike_amplitude
from msoaxon.synaptic import EPSG_TAU

V0 = -68.0
UNITARY = 26.7  # "unitary EPSG" of the MATLAB code
MYOGA_US = 221


def epsp(morph, site, mem, amp=UNITARY, site2=None, t_end=12.0):
    """Somatic voltage change for an EPSG (or a bilateral pair at 0 delay)."""
    if site2 is None:
        t, y = mso_axon("EPSG", 5, 10, amp, 3, "active-full", t_end, V0, site, morph=morph,
                        mem=mem, max_step=0.01)
    else:
        t, y = mso_axon("EPSGbilateral", 5, 5, amp, 3, "active-full", t_end, V0, site,
                        morph=morph, mem=mem, input_node2=site2, max_step=0.01)
    return soma_on_grid(t, y, V0, 4.9, t_end)


def shape(g, v):
    pk = v.max()
    r10, r90 = g[np.argmax(v >= 0.1 * pk)], g[np.argmax(v >= 0.9 * pk)]
    above = g[v >= pk / 2]
    return dict(amp=float(pk), rise_us=float((r90 - r10) * 1e3),
                half_width_us=float((above[-1] - above[0]) * 1e3))


def rin_tau(morph, mem):
    r = passive_step("multi", V0, morph=morph, mem=mem)
    return dict(rin_steady=float(r["rin_steady"]), rin_peak=float(r["rin_peak"]),
                t63_us=float(r["tau"] * 1e3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="figures/dendrites")
    ap.add_argument("--klt-from-soma", action="store_true",
                    help="start the KLT/h gradient at the soma's density instead of conserving "
                         "the lumped totals (halves total KLT)")
    ap.add_argument("--quick", action="store_true",
                    help="coarse grids and 1%% threshold tolerance: a smoke test, not results")
    ap.add_argument("--dendrite-ra", type=float, default=100.0,
                    help="axial resistivity of the dendrites [Ohm cm] (Mathews et al 2010: 200)")
    a = ap.parse_args()
    with process_pool() as pool:  # one pool for every parallel call below
        run(a, pool)


def run(a, pool):
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    D = with_dendrites(conserve_totals=not a.klt_from_soma, dendrite_ra=a.dendrite_ra)
    mid_l, mid_m, dist_l = int(D.lateral[2]), int(D.medial[2]), int(D.lateral[-1])
    tol = 1e-2 if a.quick else 1e-4  # threshold / rheobase tolerance
    res = {}

    # Every threshold search (C, E, F) and the rheobase (G) is independent of the rest,
    # so all of them go to the pool now, in one batch, and run while A, B and D are
    # computed here. Each search is deterministic, so the results don't depend on this.
    def th_job(delay, tau=EPSG_TAU, **kw):
        return pool.submit(threshold, "multi", float(delay), epsg_tau=tuple(tau), rel_tol=tol, **kw)

    def collect(jobs):
        return np.array([j.result() for j in jobs])

    dend = dict(model_kw=dict(morph=D))
    bil = dict(stim="EPSGbilateral", input_node=mid_l, input_node2=mid_m, **dend)
    mem_d = membrane(V0, morph=D)
    rb_job = pool.submit(rheobase, "multi", mem=mem_d, rel_tol=tol)  # the longest single job
    C_jobs = {"bilateral (one EPSG per dendrite)": th_job(0.0, **bil),
              "unilateral (both on the lateral dendrite)": th_job(
                  0.0, stim="EPSGbilateral", input_node=mid_l, input_node2=mid_l, **dend),
              "both at the soma, dendritic model": th_job(0.0, **dend),
              "both at the soma, lumped model": th_job(0.0)}
    delays = np.round(np.arange(0, 0.6001, 0.1 if a.quick else 0.02), 4)
    kinetics = (("model EPSG (0.18 ms)", EPSG_TAU), ("Myoga EPSG (0.3 ms)", (0.1, 0.3)))
    E_jobs = {kname: [th_job(d, tau, **bil) for d in delays] for kname, tau in kinetics}
    fd = np.round(np.arange(0, 1.0001, 0.2 if a.quick else 0.04), 4)
    F_jobs = {"lumped": [th_job(d) for d in fd], "dendritic": [th_job(d, **dend) for d in fd]}

    # A. passive
    reb = {"lumped": membrane(V0, rebalance_rest=True),
           "dendritic": membrane(V0, rebalance_rest=True, morph=D),
           "dendritic, no dendritic KLT": membrane(V0, rebalance_rest=True, morph=D, dendrite_klt_scale=0)}
    morphs = {"lumped": LUMPED, "dendritic": D, "dendritic, no dendritic KLT": D}
    res["A_rin_tau"] = {k: rin_tau(morphs[k], reb[k]) for k in ("lumped", "dendritic")}
    print("A", res["A_rin_tau"], flush=True)

    # B. EPSP shapes
    B, traces = {}, {}
    for label, morph, site, mem in [
            ("lumped, EPSG at soma", LUMPED, 1, reb["lumped"]),
            ("dendritic, EPSG at soma", D, 1, reb["dendritic"]),
            ("dendritic, EPSG mid-dendrite (100 um)", D, mid_l, reb["dendritic"]),
            ("dendritic, EPSG distal (180 um)", D, dist_l, reb["dendritic"]),
            ("mid-dendrite, dendritic KLT removed", D, mid_l, reb["dendritic, no dendritic KLT"])]:
        g, v = epsp(morph, site, mem)
        B[label] = shape(g, v)
        traces[label] = (g, v)
        print("B", label, B[label], flush=True)
    res["B_epsp"] = B

    # C. bilateral vs unilateral threshold at zero delay (default membranes)
    C = {k: j.result() for k, j in C_jobs.items()}
    res["C_threshold"] = C
    print("C", C, flush=True)

    # D. linear summation, 60% of the bilateral threshold per input
    amp = 0.6 * C["bilateral (one EPSG per dendrite)"]
    g, vl = epsp(D, mid_l, reb["dendritic"], amp=amp)
    _, vm = epsp(D, mid_m, reb["dendritic"], amp=amp)
    _, vb = epsp(D, mid_l, reb["dendritic"], amp=amp, site2=mid_m)
    res["D_summation"] = dict(amplitude=float(amp), unilateral_peaks=[float(vl.max()), float(vm.max())],
                              bilateral_peak=float(vb.max()), sum_of_waveforms_peak=float((vl + vm).max()),
                              ratio=float(vb.max() / (vl + vm).max()))
    print("D", res["D_summation"], flush=True)

    # E. coincidence window with bilateral dendritic inputs
    E, curves = {}, {}
    for kname, _ in kinetics:
        th = curves[kname] = collect(E_jobs[kname])
        E[kname] = {f"margin_{m}": float(half_width(delays, th, m) * 1e3) for m in (0.005, 0.03)}
        print("E", kname, E[kname], flush=True)
    res["E_window_bilateral_dendritic_us"] = E

    # F. EPSG-pair threshold vs delay at the soma, lumped vs dendritic
    F = {k: collect(jobs) for k, jobs in F_jobs.items()}
    res["F_soma_pair_curve"] = {"delays_ms": fd.tolist(), **{k: v.tolist() for k, v in F.items()},
                                "max_rel_diff": float(np.max(np.abs(F["dendritic"] / F["lumped"] - 1)))}
    print("F max relative difference", res["F_soma_pair_curve"]["max_rel_diff"], flush=True)

    # G. somatic spike (Scott protocol), dendritic model
    rb = rb_job.result()
    mults = (1.5, 3.0) if a.quick else (1.5, 2.0, 3.0)
    amp_jobs = {k: pool.submit(spike_amplitude, "multi", rb * k, mem=mem_d) for k in mults}
    res["G_somatic_spike"] = dict(rheobase_pA=float(rb), amplitudes_mV={
        str(k): float(j.result()["amplitude"]) for k, j in amp_jobs.items()})
    print("G", res["G_somatic_spike"], flush=True)

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    for label, (g, v) in traces.items():
        ax[0].plot(g - 5, v / v.max(), lw=1.8, label=f"{label} ({B[label]['amp']:.1f} mV)")
    ax[0].set(xlim=(-0.1, 2.5), xlabel="Time from EPSG onset (ms)", ylabel="Somatic EPSP (normalised)",
              title="B. Unitary EPSP at the soma")
    ax[0].legend(fontsize=7, frameon=False)
    for kname, th in curves.items():
        ax[1].plot(delays * 1e3, th / th[0], lw=2, label=f"bilateral dendritic, {kname}: "
                   f"{E[kname]['margin_0.03']:.0f} us")
    ax[1].axhline(1.03, color="0.5", ls=":", lw=0.8)
    ax[1].set(xlabel="Delay between EPSGs (us)", ylabel="Threshold / coincident threshold",
              title=f"E. Window at 3% margin (Myoga 2014: {MYOGA_US} us)", ylim=(0.98, 2.05))
    ax[1].legend(fontsize=7.5, frameon=False, loc="lower right")
    for k, th in F.items():
        ax[2].plot(fd * 1e3, th, lw=2, ls="-" if k == "lumped" else "--", label=k)
    ax[2].set(xlabel="Delay between EPSGs (us)", ylabel="Threshold (EPSG units)",
              title="F. EPSG pair at the soma: lumped vs dendritic")
    ax[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "dendrites.png", dpi=130)
    (out / "dendrites.json").write_text(json.dumps(res, indent=1))
    print(f"saved {out}/dendrites.png and .json")


if __name__ == "__main__":
    main()
