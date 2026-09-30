"""Figures (port of Graphing.m).

`graph` is a dict with the same keys as the MATLAB struct: Type, node, tEnd,
stimType, inputNode, start, stop, I, Syn.
"""

import matplotlib.pyplot as plt
import numpy as np

from ._common import stimulus
from .two import applied_current


def _panel(ax, t, V, title, labels, t_end, ylabel="Voltage (mV)"):
    ax.plot(t, V, linewidth=2)
    ax.set_title(title)
    ax.legend(labels, fontsize=8, frameon=False)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, t_end)


def input_trace(graph, t2, x):
    """Applied current along the two-compartment solution.

    Graphing.m recomputes this with its own copy of the stimulus code (and indexes
    the 2-CPT voltage with the multi-compartment time index); here it reuses the
    model's own stimulus function evaluated on the 2-CPT run instead.
    """
    s = stimulus(graph["stimType"], graph["start"], graph["stop"], graph["I"],
                 graph["tEnd"], graph["Syn"])
    return np.array([applied_current(t, V1, graph["stimType"], s) for t, V1 in zip(t2, x[:, 0])])


def graphing(graph, t1, y, t2=None, x=None, t3=None, z=None):
    node, t_end, kind = graph["node"], graph["tEnd"], graph["Type"]
    fig = plt.figure(figsize=(9, 7))

    if kind == "Contour":
        ax = fig.add_subplot(projection="3d")
        X, Y = np.meshgrid(np.arange(1, 46), t1)
        ax.plot_surface(X, Y, y[:, :45], cmap="viridis", linewidth=0)
        ax.set(title="Voltage Propagation", xlabel="Compartment Number (1-45)",
               ylabel="Time", zlabel="Voltage (mV)")

    elif kind == "ModelComparison":
        a1, a2 = fig.subplots(2, 1)
        _panel(a1, t1, y[:, [0, node - 1]], "multi-compartment", ["Soma", f"node {node}"], t_end)
        _panel(a2, t2, x[:, :2], "two-compartment", ["Soma", "node"], t_end)
        fig.suptitle("Comparison of 2 Compartment and Multi-Compartment Models")

    elif kind == "ModelComparison2":
        a1, a2 = fig.subplots(2, 1)
        _panel(a1, t2, x[:, :2], "Two-Compartment (Node 3)", ["Soma", "node"], t_end)
        _panel(a2, t3, z[:, :2], "Two-Compartment (Node 5)", ["Soma", "node"], t_end)
        fig.suptitle("Comparison of 2 Compartment Models")

    elif kind == "NodeComparison":
        a1, a2 = fig.subplots(2, 1)
        labels = ["Multi-Compartment", "2 Compartment"]
        a1.plot(t1, y[:, 0], linewidth=2)
        _panel(a1, t2, x[:, 0], "Soma", labels, t_end)
        a2.plot(t1, y[:, node - 1], linewidth=2)
        _panel(a2, t2, x[:, 1], f"node {node}", labels, t_end)
        fig.suptitle(f"Comparison of Soma and node {node}")

    elif kind == "MultiCompartment":
        _panel(fig.add_subplot(), t1, y[:, [0, node - 1]],
               f"Multi-compartment model, node: {node}", ["Soma", f"node {node}"], t_end)

    elif kind == "TwoCompartment":
        _panel(fig.add_subplot(), t1, y[:, :2], "two-compartment", ["Soma", "node"], t_end)

    elif kind == "Input":
        _panel(fig.add_subplot(), t2, input_trace(graph, t2, x), "Input", [], t_end,
               ylabel="Current (pA)")

    elif kind == "Multigraph":
        gs = fig.add_gridspec(2, 3)
        cols = [1, 3, 5, 9, 13, 21, 29, 37, 45]
        _panel(fig.add_subplot(gs[0, :]), t1, y[:, [c - 1 for c in cols]], "multi-compartment",
               ["Soma"] + [f"node {c}" for c in cols[1:]], t_end)
        _panel(fig.add_subplot(gs[1, 0]), t2, x[:, :2], "Two-Compartment (Node 3)", ["Soma", "node"], t_end)
        _panel(fig.add_subplot(gs[1, 1]), t3, z[:, :2], "Two-Compartment (Node 5)", ["Soma", "node"], t_end)
        _panel(fig.add_subplot(gs[1, 2]), t2, input_trace(graph, t2, x), "Input", [], t_end,
               ylabel="Current (pA)")
        fig.suptitle("Comparison of 2 Compartment and Multi-Compartment Models with Input")

    else:
        raise ValueError(f"unknown graph Type {kind!r}")

    fig.tight_layout()
    return fig


def threshold_figure(x_values, multi, two, title, xlabel):
    """One panel of Making_Threshold_graphs.m."""
    fig, ax = plt.subplots()
    ax.plot(x_values, multi, linewidth=2, label="Multi-Compartment")
    ax.plot(x_values, two, linewidth=2, label="Two-Compartment")
    ax.set(title=title, xlabel=xlabel, ylabel="Threshold")
    ax.legend(fontsize=8, frameon=False)
    return fig
