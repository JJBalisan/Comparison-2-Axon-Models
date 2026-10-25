"""Where each spike starts in the multi-compartment model (issue #80).

Lehnert et al 2014 (J Neurosci 34:5370) found that under sustained,
high-rate input 5-60% of spikes start at node 1 or 2 rather than in the
axon initial segment (AIS). Nothing here changes the model: this reads a
finished run (t, x) and says which compartment fired first.

The criterion. A spike is located by the time of steepest rise (peak dV/dt)
of each spike-generating compartment, and its site is the earliest of them.
Two simpler readouts fail:
- A fixed fraction of each compartment's own peak (or a dV/dt threshold)
  reports the soma. The input enters at the soma, so the EPSP half-rises
  before the AIS has fired (pair at 60 nS: soma 5.15 ms, AIS 5.33 ms), and a
  non-regenerative AIS compartment (peak 35 mV above rest) half-rises on the
  EPSP too.
- A fixed voltage level is meaningless across compartments: nodal spikes are
  ~100 mV, somatic ones ~10-40 mV, and under sustained input the baseline
  moves from spike to spike.
The steepest rise of a compartment that fires is its Na upstroke, which the
EPSP does not match, so the criterion is read on the compartments whose
membrane regenerates (labels in SITE_LABELS: "AIS" and "node") and never on
the soma, internodes or dendrites. Their times are still returned, as the
propagation profile. One more guard: a compartment next to the soma, shunted
by it, can rise faster on the EPSP than on a barely-threshold spike, so a
compartment whose steepest rise is no steeper than the soma's is following
the soma and is not eligible for that spike.

Spikes are found on a reference node, by default the most distal one: every
spike propagates there, so a spike is the reference rising `rise` mV above
its resting value (the first sample). Each spike's window runs from the
midpoint to the previous spike to the midpoint to the next, at most REACH
either side of the reference's rise. The reference's steepest rise in it is
the spike's; each compartment's steepest rise is then sought from LEAD before
that to LAG after it, so neither the EPSP ahead of the spike nor the next
input is mistaken for the upstroke. Times are the solver's samples, with the
vertex of the parabola through the three samples around the maximum.

For the 45-compartment chain the site is the 1-indexed compartment: 2 or 3 is
the AIS, 5 node 1, 7 node 2 (Morphology.labels names them).
"""

import dataclasses
import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from msoaxon import _common, multi, synaptic

SITE_LABELS = ("AIS", "node")  # compartments that can be the initiation site
RISE = 60.0  # a reference spike: this far [mV] above its resting level
# A spike's wide window reaches this far [ms] either side of the reference's
# rise (and no further than the midpoints to the neighbouring spikes). Inside
# it the compartments' steepest rises are sought from LEAD before the
# reference's steepest rise to LAG after it. The AIS leads the last node by
# 0.35-0.4 ms in the 45-compartment chain; much beyond that an EPSP's own peak
# in the AIS starts to pass for the upstroke (the share of spikes that start
# at node 1 or 2 is flat for LEAD 0.4-0.6, then falls: see PYTHON.md).
REACH = 1.5
LEAD = 0.6
LAG = 0.1
PEAK_FRAC = 0.5  # see _steepest
# the groups group_fractions() reports, in its order
GROUPS = ("AIS", "node 1", "node 2", "other")


@dataclasses.dataclass(frozen=True)
class Spike:
    """One spike and where it started.

    Attributes:
        compartment: The 1-indexed initiation site: the eligible
            compartment with the earliest steepest rise.
        time: That compartment's time of steepest rise [ms].
        times: The time of steepest rise of every compartment [ms]
            (n entries, 0-indexed); only the eligible ones compete.
    """

    compartment: int
    time: float
    times: np.ndarray


@dataclasses.dataclass(frozen=True)
class Fractions:
    """The share of spikes starting in each group of compartments.

    Attributes:
        n_spikes: The number of spikes counted.
        ais: Fraction starting in the AIS (labels "AIS").
        node1: Fraction starting at the first node.
        node2: Fraction starting at the second node.
        other: Fraction starting anywhere else.
    """

    n_spikes: int
    ais: float
    node1: float
    node2: float
    other: float

    @property
    def distal(self) -> float:
        """The fraction starting at node 1 or 2: Lehnert 2014's 5-60%."""
        return self.node1 + self.node2


def _sites(morph: multi.Morphology, labels: Sequence[str]) -> np.ndarray:
    """Returns the 0-indexed compartments whose label is in `labels`."""
    return np.array(
        [i for i, name in enumerate(morph.labels) if name in labels], dtype=int
    )


def _crossings(t: np.ndarray, v: np.ndarray, level: float) -> np.ndarray:
    """Returns the times v rises through `level`, interpolated linearly."""
    i = np.flatnonzero((v[:-1] < level) & (v[1:] >= level))
    frac = (level - v[i]) / (v[i + 1] - v[i])
    return t[i] + frac * (t[i + 1] - t[i])


def _steepest(
    t: np.ndarray, v: np.ndarray, last: bool = False
) -> tuple[float, float]:
    """Returns the time [ms] and slope [mV/ms] of a steep rise of v(t).

    The slope on each interval sits at the interval's midpoint; the time is
    the vertex of the parabola through the slope's peak and its two
    neighbours (the peak itself at either end of the window).

    Args:
        t: Strictly increasing times.
        v: The voltage at t.
        last: False takes the steepest rise anywhere. True takes the last
            interior peak of the slope that reaches PEAK_FRAC of the
            highest: the upstroke after an EPSP's own, and never the
            rise that is still climbing at the window's edge. (nan, -inf)
            if there is no interior peak.
    """
    slope = np.diff(v) / np.diff(t)
    tm = (t[:-1] + t[1:]) / 2
    k = int(np.argmax(slope))
    if last:
        inner = np.arange(1, len(slope) - 1)
        up = (slope[inner] > slope[inner - 1]) & (
            slope[inner] >= slope[inner + 1]
        )
        peaks = inner[up & (slope[inner] > 0)]
        if len(peaks) == 0:
            return math.nan, -math.inf
        k = int(peaks[slope[peaks] >= PEAK_FRAC * slope[peaks].max()][-1])
    if k == 0 or k == len(slope) - 1:
        return float(tm[k]), float(slope[k])
    # the vertex of the parabola through three unevenly spaced points; den is
    # positive because y1 > y0 and y1 >= y2
    (x0, x1, x2), (y0, y1, y2) = tm[k - 1 : k + 2], slope[k - 1 : k + 2]
    num = (x1 - x0) ** 2 * (y1 - y2) - (x1 - x2) ** 2 * (y1 - y0)
    den = (x1 - x0) * (y1 - y2) - (x1 - x2) * (y1 - y0)
    return float(np.clip(x1 - num / (2 * den), x0, x2)), float(y1)


def find_spikes(
    t: np.ndarray,
    x: np.ndarray,
    morph: multi.Morphology = multi.LUMPED,
    reference: int | None = None,
    rise: float = RISE,
    labels: Sequence[str] = SITE_LABELS,
    lead: float = LEAD,
) -> list[Spike]:
    """Finds the spikes of a run and the compartment each starts in.

    See the module docstring for the criterion.

    Args:
        t: The run's times [ms].
        x: The run's state, compartment voltages in the first columns.
        morph: The morphology the run used.
        reference: The 1-indexed node spikes are detected on; None is the
            most distal compartment whose label is in `labels`.
        rise: A spike is the reference rising this far above its resting
            level, its first sample [mV].
        labels: The labels of the compartments eligible as the site.
        lead: How long before the reference's steepest rise [ms] a
            compartment's steepest rise is sought.

    Returns:
        One Spike per spike, in time order; empty if there is none.

    Raises:
        ValueError: If no compartment has a label in `labels`, or
            `reference` is out of range.
    """
    n = morph.n
    eligible = _sites(morph, labels)
    if len(eligible) == 0:
        raise ValueError(f"no compartment is labelled any of {tuple(labels)}")
    ref = int(eligible[-1]) if reference is None else reference - 1
    if not 0 <= ref < n:
        raise ValueError(f"reference must be 1..{n}, got {reference}")
    # the solver repeats a time at a restart; a slope needs distinct times
    keep = np.concatenate([[True], np.diff(t) > 0])
    t, v = np.asarray(t)[keep], np.asarray(x)[keep, :n]
    onsets = _crossings(t, v[:, ref], v[0, ref] + rise)
    if len(onsets) == 0:
        return []
    edges = np.concatenate([[t[0]], (onsets[:-1] + onsets[1:]) / 2, [t[-1]]])
    soma = [i for i, name in enumerate(morph.labels) if name == "soma"]
    spikes = []
    for onset, left, right in zip(onsets, edges[:-1], edges[1:], strict=True):
        lo, hi = max(left, onset - REACH), min(right, onset + REACH)
        wide = (t >= lo) & (t <= hi)
        if wide.sum() < 3:
            continue
        # the reference sits far from the input, so its steepest rise is the
        # spike's; the others' are sought just before it
        t_ref = _steepest(t[wide], v[wide, ref])[0]
        w = (t >= max(lo, t_ref - lead)) & (t <= min(hi, t_ref + LAG))
        if w.sum() < 3:
            continue
        found = [_steepest(t[w], v[w, c], last=True) for c in range(n)]
        times, slopes = (np.array(a) for a in zip(*found, strict=True))
        # a compartment that rises no faster than the soma did, anywhere in
        # the wide window, is following the EPSP and not firing
        floor = _steepest(t[wide], v[wide, soma[0]])[1] if soma else -np.inf
        firing = eligible[slopes[eligible] > floor]
        if len(firing) == 0:  # only if `reference` is not itself eligible
            continue
        first = firing[np.argmin(times[firing])]
        spikes.append(Spike(int(first) + 1, float(times[first]), times))
    return spikes


def group_fractions(
    spikes: Sequence[Spike], morph: multi.Morphology = multi.LUMPED
) -> Fractions:
    """Sorts spikes into AIS / node 1 / node 2 / other, by label.

    Node 1 is the first compartment labelled "node" along the axon, node 2
    the second.

    Args:
        spikes: find_spikes' result.
        morph: The morphology the run used.

    Returns:
        The fractions; all nan if there are no spikes.
    """
    ais = set(_sites(morph, ("AIS",)) + 1)
    nodes = _sites(morph, ("node",)) + 1
    where = {int(nodes[0]): "node 1", int(nodes[1]): "node 2"}
    where.update(dict.fromkeys((int(c) for c in ais), "AIS"))
    counts = dict.fromkeys(GROUPS, 0)
    for s in spikes:
        counts[where.get(s.compartment, "other")] += 1
    total = len(spikes)
    share = {g: c / total if total else float("nan") for g, c in counts.items()}
    return Fractions(
        n_spikes=total,
        ais=share["AIS"],
        node1=share["node 1"],
        node2=share["node 2"],
        other=share["other"],
    )


def sustained_input(
    freq: float = 350.0,
    g_e: float = 1.3e7,
    t_end: float = 100.0,
    random_in: int = 1396,
    morph: multi.Morphology | None = None,
    **model_kw: Any,
) -> Fractions:
    """Runs the model on sustained Synaptic input and reads where it fires.

    The input is the gammatone-filtered noise of synaptic.synaptic at
    `freq`, on from 5 ms (_common.START) to the end of the run, into the
    soma.

    Args:
        freq: The input's centre frequency [Hz].
        g_e: The EPSG conductance (SynParams.gE).
        t_end: The run's length [ms]; 100 ms gives ~20 spikes at the
            defaults.
        random_in: The noise seed (SynParams.random_in).
        morph: The morphology; default multi.LUMPED.
        **model_kw: Further mso_axon keywords (mem, mode, ...).

    Returns:
        The fractions of spikes by site.
    """
    morph = morph or multi.LUMPED
    syn = synaptic.SynParams(
        t_end=t_end, freq=freq, gE=g_e, random_in=random_in
    )
    t, x = multi.mso_axon(
        "Synaptic",
        _common.START,
        t_end,
        0.0,
        3,
        "active-full",
        t_end,
        -68.0,
        1,
        syn=syn,
        morph=morph,
        **model_kw,
    )
    return group_fractions(find_spikes(t, x, morph), morph)
