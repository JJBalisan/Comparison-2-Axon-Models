"""Run both models on one stimulus and plot them (port of Combine_all.m).

    uv run scripts/combine_all.py --stim EPSGpair --node 3 --out comparison.png
"""

import argparse

import matplotlib.pyplot as plt

from msoaxon import SynParams, mso_axon, two_cpt
from msoaxon.plotting import graphing
from msoaxon.spiking import count_rising_edges, spiking

# per-stimulus defaults from Combine_all.m: (stop, I)
DEFAULTS = {
    "step": (10.0, 8500.0),
    "ramp": (5.5, 8500.0),  # Combine_all.m has stop == start (5), a zero-length ramp
    "ramp2": (5.5, 8500.0),
    "sine": (5 + 1000 / (2 * 200), 8500.0),
    "EPSG": (10.0, 108.0),
    "EPSGpair": (5.0, 40.0),  # stop is the second EPSG onset
    "Synaptic": (15.0, 0.0),
    "SynapticPair": (15.0, 0.0),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stim", default="EPSGpair", choices=DEFAULTS)
    ap.add_argument("--type", default="active-full",
                    help="passive, active-KLT, active-H, active-KLT+H, active-sodium, active-KHT, active-full")
    ap.add_argument("--node", type=int, default=3, help="axon compartment to compare (1-indexed)")
    ap.add_argument("--input-node", type=int, default=1)
    ap.add_argument("--t-end", type=float, default=20.0)
    ap.add_argument("--v0", type=float, default=-68.0)
    ap.add_argument("--start", type=float, default=5.0)
    ap.add_argument("--stop", type=float, help="override the per-stimulus default")
    ap.add_argument("--I", type=float, help="override the per-stimulus default")
    ap.add_argument("--factor", type=float, default=30.0, help="mV above soma counted as a spike")
    ap.add_argument("--graph", default="ModelComparison",
                    help="ModelComparison, NodeComparison, MultiCompartment, Contour, Input, Multigraph, ModelComparison2")
    ap.add_argument("--out", help="save the figure here instead of showing it")
    a = ap.parse_args()

    stop, I = DEFAULTS[a.stim]
    stop = a.stop if a.stop is not None else stop
    I = a.I if a.I is not None else I
    syn = SynParams(t_end=a.t_end, random_in=1396, diff=1.0, f=200.0)

    args = (a.stim, a.start, stop, I)
    t1, y = mso_axon(*args, a.node, a.type, a.t_end, a.v0, a.input_node, syn)
    t2, x = two_cpt(*args, a.node, a.type, a.t_end, a.v0, a.input_node, syn)
    # Combine_all.m ran TwoCpt twice with identical arguments for the node-3/node-5
    # panels; here the second run uses node 5 so those panels mean what they say.
    t3, z = two_cpt(*args, 5, a.type, a.t_end, a.v0, a.input_node, syn)

    # Combine_all.m counted column 1 of the 2-CPT markers (soma minus soma, always 0);
    # column 2 is the axon compartment it meant.
    n_two = count_rising_edges(spiking(x, a.factor, "Two")[:, 1])
    n_multi = count_rising_edges(spiking(y, a.factor, "Mult")[:, a.node - 1])
    print(f"spikes at node {a.node}: multi-compartment {n_multi}, two-compartment {n_two}")

    graph = dict(Type=a.graph, node=a.node, tEnd=a.t_end, stimType=a.stim,
                 inputNode=a.input_node, start=a.start, stop=stop, I=I, Syn=syn)
    fig = graphing(graph, t1, y, t2, x, t3, z)
    if a.out:
        fig.savefig(a.out, dpi=150)
        print(f"saved {a.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
