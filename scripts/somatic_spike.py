"""Somatic spike size vs mature MSO cells (Scott et al 2005 protocol).

    uv run scripts/somatic_spike.py --out-dir figures/somatic

100 ms somatic steps at multiples of rheobase; amplitude from the inflection
point. Includes the dendrotoxin control: KLT removed at soma and AIS, with rest
rebalanced so the cell doesn't fire spontaneously.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from msoaxon.multi import membrane
from msoaxon.somatic import rheobase, spike_amplitude

MULTS = (1.5, 2.0, 3.0)  # below ~1.5x the somatic trace has no distinct inflection
CASES = {
    "45-compartment (msoAxon.m)": ("multi", {}),
    "45-compartment, KLT blocked at soma + AIS (dendrotoxin)":
        ("multi", dict(soma_klt_scale=0, ais_klt_scale=0, rebalance_rest=True)),
    "2-compartment (TwoCpt.m)": ("two", {}),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="figures/somatic")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(CASES), figsize=(15, 4.2), sharey=True)
    results = {}
    for ax, (name, (model, knobs)) in zip(axes, CASES.items()):
        mem = membrane(-68.0, **knobs) if knobs else None
        rb = rheobase(model, mem=mem)
        rows = []
        for k in MULTS:
            r = spike_amplitude(model, rb * k, mem=mem)
            t, x = r.pop("trace")  # the run the amplitude was measured on
            rows.append({"multiple": k, **{key: float(v) for key, v in r.items()}})
            keep = t >= 4.8
            ax.plot(t[keep], x[keep, 0], lw=1.5, label=f"{k}x: {r['amplitude']:.1f} mV")
            ax.plot(r["t_inflection"], r["v_inflection"], "k.", ms=7)
        results[name] = {"rheobase_pA": float(rb), "spikes": rows}
        ax.set_title(f"{name}\nrheobase {rb:.0f} pA", fontsize=9.5)
        ax.set_xlabel("Time (ms)")
        ax.legend(fontsize=8, frameon=False, title="amplitude from inflection (dot)",
                  title_fontsize=7.5)
        print(f"{name}: rheobase {rb:.0f} pA; " + ", ".join(
            f"{r['multiple']}x {r['amplitude']:.1f} mV" for r in rows), flush=True)
    axes[0].set_ylabel("Soma voltage (mV)")
    fig.suptitle("Somatic spike, 100 ms step (Scott et al 2005 protocol). "
                 "Mature MSO: 17 ± 2 mV; dendrotoxin: 15 → 37 mV", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "somatic_spike.png", dpi=130)
    (out / "somatic_spike.json").write_text(json.dumps(results, indent=1))
    print(f"saved {out}/somatic_spike.png and .json")


if __name__ == "__main__":
    main()
