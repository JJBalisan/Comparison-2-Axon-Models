"""Multi-compartment MSO soma + axon model (port of msoAxon.m, 45 compartments).

Equations as in Lehnert et al 2014 unless noted. All quantities are densities
(pA/um^2, nS/um^2).

State layout matches MATLAB's column-major reshape(x, 45, 7): the output y has
315 columns, [V(45), m(45), h(45), p(45), w(45), z(45), a(45)]. So y[:, k-1] is
the voltage of MATLAB compartment k, and y[:, 225] is z1.

Optional dendrites (with_dendrites(), Lehnert et al 2014's Fig. 8 variant) are
appended after the axon, so compartments 1-45 keep their meaning; the state then
has 7 blocks of n = 45 + dendrite compartments.
"""

import dataclasses
from typing import Any

import numpy as np
from scipy import sparse

from msoaxon import _common, _fast, _solve, constants as C, synaptic

N = C.N_CPT

# specific capacitance, uF/cm^2 converted to ms*nS/um^2
CAP = np.concatenate([[0.8, 0.8, 0.8], np.tile([0.01, 0.8], 21)]) * 1e6 * 1e-8
V_NA, V_K, V_H = 69.0, -90.0, -35.0
R_AXIAL = 100.0  # specific axial resistivity [Ohm cm]

# conductance between neighbouring compartments i and i+1 (the (2/R)/(...) term)
_G_AX = (2 / R_AXIAL) / (C.L_CM[:-1] / C.XA_CM[:-1] + C.L_CM[1:] / C.XA_CM[1:])

MSO_STIM_TYPES = _common.STIM_TYPES + ("EPSGbilateral",)


# which gating blocks (m, h, p, w, z, a) evolve for each model type
ACTIVE_GATES = {
    "active-KLT": (0, 0, 0, 1, 1, 0),
    "active-H": (0, 0, 0, 0, 0, 1),
    "active-KLT+H": (0, 0, 0, 1, 1, 1),
    "passive": (0, 0, 0, 0, 0, 0),
    "active-sodium": (1, 1, 0, 0, 0, 0),
    "active-full": (1, 1, 1, 1, 1, 1),
    "active-KHT": (0, 0, 1, 0, 0, 0),
}

# BDF solver settings of mso_axon (msoAxon.m's ode15s options). max_step is
# MAX_STEP_FRACTION * t_end unless given: ode15s's default MaxStep.
RTOL = 1e-8
ATOL = 1e-8
MAX_STEP_FRACTION = 0.1

# Initial gating values (m, h, p, w, z, a), as in msoAxon.m: the steady state
# at -68 mV rounded to two decimals (exact: 0.1192, 0.6699, 0.0006, 0.2868,
# 0.6646, 0.2244; w and z are 0.01 off). The gates therefore start slightly
# off steady state, which is the source of the ~0.9 mV drift of the rest
# potential in membranes built without rebalance_rest. They must stay these
# literals: the golden traces and the Rust port both start from them.
_GATE0 = (0.12, 0.67, 0.0, 0.28, 0.67, 0.22)


def _chain_jac_sparsity() -> sparse.csc_matrix:
    """Returns the Jacobian pattern of msoAxon.m's unbranched chain.

    dV couples to neighbouring V and its own gates; each gate to its own V.
    """
    I = sparse.eye(N)
    tri = sparse.diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(N, N))
    blocks = [[tri] + [I] * 6] + [
        [I] + [I if j == i else None for j in range(6)] for i in range(6)
    ]
    return sparse.bmat(blocks).tocsc()


def _tree_jac_sparsity(
    n: int, parent: np.ndarray, child: np.ndarray
) -> sparse.csc_matrix:
    """Returns _chain_jac_sparsity's pattern for any tree.

    dV couples to its own V, to each V it shares an edge
    (parent[k], child[k]) with, and to its own gates; each gate to its
    own V.
    """
    I = sparse.eye(n, format="csr")
    adj = sparse.coo_matrix(
        (
            np.ones(2 * len(parent)),
            (np.r_[parent, child], np.r_[child, parent]),
        ),
        shape=(n, n),
    )
    couple = (I + adj).tocsr()
    couple.data[:] = 1.0
    blocks = [[couple] + [I] * 6] + [
        [I] + [I if j == i else None for j in range(6)] for i in range(6)
    ]
    return sparse.bmat(blocks).tocsc()


@dataclasses.dataclass(frozen=True, eq=False)
class Morphology:
    """A multi-compartment morphology (LUMPED and with_dendrites()).

    eq=False keeps identity equality and hashing, so the ndarray fields
    never reach a generated __eq__.

    Attributes:
        key: Hashable summary of every setting, part of the pre-stimulus
            cache key.
        n: Number of compartments (45 for LUMPED).
        sa: Membrane area per compartment [um^2].
        cap: Capacitance per compartment [ms nS/um^2].
        g_na: Na density per compartment [nS/um^2]; membrane() may
            override it without touching this.
        g_kht: KHT density per compartment [nS/um^2].
        g_klt: KLT density per compartment [nS/um^2]; membrane() may
            override it.
        g_h: h density per compartment [nS/um^2].
        g_lk: Leak density per compartment [nS/um^2]; membrane() may
            override the leak reversal.
        parent: The axial edges, 0-based: current flows along
            parent[k] -> child[k].
        child: See parent.
        g_ax: Conductance of each axial edge [S].
        chain: True if the edges are exactly 0-1-2-..., so axial_current
            can use msoAxon.m's own float order.
        jac: Sparsity pattern of the 7n x 7n Jacobian for the BDF
            solver.
        labels: A name per compartment ("soma", "AIS", "node", "lateral
            dendrite", ...).
        lateral: 1-indexed compartment numbers of the lateral dendrite,
            proximal to distal; with_dendrites() only.
        medial: As lateral, for the medial dendrite.
    """

    key: str | tuple[Any, ...]
    n: int
    sa: np.ndarray
    cap: np.ndarray
    g_na: np.ndarray
    g_kht: np.ndarray
    g_klt: np.ndarray
    g_h: np.ndarray
    g_lk: np.ndarray
    parent: np.ndarray
    child: np.ndarray
    g_ax: np.ndarray
    chain: bool
    jac: sparse.csc_matrix
    labels: list[str]
    lateral: np.ndarray | None = None
    medial: np.ndarray | None = None


@dataclasses.dataclass(frozen=True, eq=False)
class Membrane:
    """Per-compartment channel overrides for one morphology (membrane()).

    eq=False keeps identity equality and hashing.

    Attributes:
        key: Every setting membrane() was built from, part of the
            pre-stimulus cache key.
        morph_key: The morphology's key.
        morph: The morphology these arrays are for.
        g_na: Na density per compartment [nS/um^2].
        g_klt: KLT density per compartment [nS/um^2].
        na_vhalf: Negated Na inactivation midpoint per compartment
            [mV], as C.hinf's vhalf expects.
        vlk: Leak reversal per compartment [mV].
        y0_gates: Initial (m, h, p, w, z, a) steady-state blocks, each
            of length n; None unless rebalance_rest.
    """

    key: tuple[float | bool, ...]
    morph_key: str | tuple[Any, ...]
    morph: Morphology
    g_na: np.ndarray
    g_klt: np.ndarray
    na_vhalf: np.ndarray
    vlk: np.ndarray
    y0_gates: tuple[np.ndarray, ...] | None


# The unbranched 45-compartment chain of msoAxon.m (soma, 2 AIS, 21
# internode/node pairs)
_PARENT, _CHILD = np.arange(N - 1), np.arange(1, N)
for _a in (
    CAP,
    _G_AX,
    _PARENT,
    _CHILD,
):  # read-only, as in constants.py
    _a.flags.writeable = False
LUMPED = Morphology(
    key="lumped",
    n=N,
    sa=C.SA,
    cap=CAP,
    g_na=C.G_NA,
    g_kht=C.G_KHT,
    g_klt=C.G_KLT,
    g_h=C.G_H,
    g_lk=C.G_LK,
    parent=_PARENT,
    child=_CHILD,
    g_ax=_G_AX,
    chain=True,
    jac=_chain_jac_sparsity(),
    labels=["soma", "AIS", "AIS"] + ["internode", "node"] * 21,
)


def with_dendrites(
    length: float = 200.0,
    diameter: float = 5.0,
    n_seg: int = 5,
    klt_lambda: float = 74.0,
    conserve_totals: bool = True,
    dendrite_ra: float = R_AXIAL,
) -> Morphology:
    """Builds Lehnert et al 2014's dendritic variant (their Fig. 8).

    Two identical unbranched dendrites (lateral = ipsilateral input,
    medial = contralateral) of `length` um and `diameter` um, `n_seg`
    compartments each, attached to the soma and appended after the
    axon. The soma shrinks so total membrane stays 8750 um^2 ("the
    somatic surface was reduced to 2467 um^2"); its Na density is scaled
    up so total Na is unchanged. Dendrites have no Na or KHT; KLT and h
    decay exponentially with distance from the soma (length constant
    74 um, Mathews et al 2010), starting from the soma's densities; leak
    and capacitance as at the soma. The paper doesn't give axial
    resistivity, so the model's 100 Ohm cm is the default.

    The paper doesn't say whether total KLT and h were kept when they
    were spread along the dendrites; it does say so for Na.

    Compartments (1-indexed, as for input_node): 46..45+n_seg lateral
    dendrite, proximal to distal; the next n_seg the medial dendrite.

    Args:
        length: Dendrite length [um].
        diameter: Dendrite diameter [um].
        n_seg: Compartments per dendrite.
        klt_lambda: Length constant of the KLT and h decay [um].
        conserve_totals: True (default) scales soma + dendrite KLT and h
            so their totals equal the lumped model's; it reproduces the
            paper's finding that this variant tunes "almost identically"
            to the lumped model (EPSG-pair thresholds within 12%). False
            starts the gradient at the soma's Table 2 density, which
            halves total KLT (thresholds up to 50% lower, somatic spike
            21-30 mV).
        dendrite_ra: Axial resistivity of the dendrite compartments only
            [Ohm cm] (Mathews et al 2010 used 200 Ohm cm for soma and
            dendrites). The soma and axon keep 100.

    Returns:
        A Morphology, with `lateral` and
        `medial` holding each dendrite's 1-indexed compartment numbers.

    Raises:
        ValueError: If the dendrites need more membrane than the lumped
            soma has.
    """
    n_d = 2 * n_seg
    n = N + n_d
    seg = length / n_seg
    x_mid = (
        np.arange(n_seg) + 0.5
    ) * seg  # distance of each compartment's centre [um]
    sa_d = np.pi * diameter * seg
    soma_sa = C.SA[0] - 2 * n_seg * sa_d
    if soma_sa <= 0:
        raise ValueError(
            "dendrites larger than the lumped soma+dendrite membrane"
        )

    sa = np.concatenate([[soma_sa], C.SA[1:], np.full(n_d, sa_d)])
    L = np.concatenate(
        [[np.sqrt(soma_sa / np.pi)], C.L[1:], np.full(n_d, seg)]
    )  # soma: as Constants.m
    r_cyl = sa / (2 * np.pi * L)
    xa_cm, l_cm = np.pi * r_cyl**2 * 1e-8, L * 1e-4
    decay = np.tile(np.exp(-x_mid / klt_lambda), 2)

    g_na = np.concatenate(
        [[C.G_NA[0] * C.SA[0] / soma_sa], C.G_NA[1:], np.zeros(n_d)]
    )
    g_kht = np.concatenate([C.G_KHT, np.zeros(n_d)])
    g_klt = np.concatenate([C.G_KLT, C.G_KLT[0] * decay])
    g_h = np.concatenate([C.G_H, C.G_H[0] * decay])
    if conserve_totals:
        region = np.r_[0, np.arange(N, n)]  # soma + dendrites
        for g, lumped in ((g_klt, C.G_KLT[0]), (g_h, C.G_H[0])):
            g[region] *= lumped * C.SA[0] / np.sum(g[region] * sa[region])
    g_lk = np.concatenate([C.G_LK, np.full(n_d, C.G_LK[0])])
    cap = np.concatenate([CAP, np.full(n_d, CAP[0])])

    lat = N + np.arange(n_seg)
    med = N + n_seg + np.arange(n_seg)
    parent = np.concatenate([np.arange(N - 1), [0], lat[:-1], [0], med[:-1]])
    child = np.concatenate([np.arange(1, N), lat, med])
    g_ax = (2 / R_AXIAL) / (
        l_cm[parent] / xa_cm[parent] + l_cm[child] / xa_cm[child]
    )
    if dendrite_ra != R_AXIAL:
        # only edges into a dendrite compartment change, so the axon's
        # conductances stay bit-identical; each half-compartment uses its own
        # resistivity
        ra = np.concatenate(
            [np.full(N, R_AXIAL), np.full(n_d, float(dendrite_ra))]
        )
        into = child >= N
        p, c = parent[into], child[into]
        g_ax[into] = 2 / (
            ra[p] * l_cm[p] / xa_cm[p] + ra[c] * l_cm[c] / xa_cm[c]
        )
    labels = (
        LUMPED.labels
        + ["lateral dendrite"] * n_seg
        + ["medial dendrite"] * n_seg
    )
    return Morphology(
        key=(
            "dendrites",
            float(length),
            float(diameter),
            int(n_seg),
            float(klt_lambda),
            bool(conserve_totals),
            float(dendrite_ra),
        ),
        n=n,
        sa=sa,
        cap=cap,
        g_na=g_na,
        g_kht=g_kht,
        g_klt=g_klt,
        g_h=g_h,
        g_lk=g_lk,
        parent=parent,
        child=child,
        g_ax=g_ax,
        chain=False,
        jac=_tree_jac_sparsity(n, parent, child),
        labels=labels,
        lateral=lat + 1,
        medial=med + 1,
    )  # 1-indexed compartment numbers


def axial_current(V: np.ndarray, morph: Morphology = LUMPED) -> np.ndarray:
    """Returns the axial current density, sign convention as in msoAxon.m.

    Args:
        V: Voltage per compartment [mV].
        morph: The morphology.

    Returns:
        The axial current density per compartment [pA/um^2].
    """
    if morph.chain:  # msoAxon.m's own float order
        flow = morph.g_ax * (V[:-1] - V[1:])  # from i to i+1 [mA]
        I = np.zeros_like(V)
        I[:-1] -= flow
        I[1:] += flow
        return -I * 1e9 / morph.sa
    flow = morph.g_ax * (V[morph.parent] - V[morph.child])
    I = np.zeros_like(V)
    np.subtract.at(I, morph.parent, flow)  # the soma is a parent three times
    np.add.at(I, morph.child, flow)
    return -I * 1e9 / morph.sa


# external_current and two.applied_current look like duplicates but each copies
# its own MATLAB file (msoAxon.m / TwoCptODE.m), and they differ on purpose.
# Merging them would change results. The differences:
#
#                 multi (this file)              two (applied_current)
#   units, sign   density [pA/um^2], negative    current [pA], positive
#                 depolarises (enters dV with    depolarises
#                 a minus)
#   where         compartment input_node (but    soma if input_node == 1,
#                 step always into the soma)     else axon
#   area          the soma's for every stimulus  -
#                 but EPSG/EPSGpair, even at
#                 another node
#   window        start < t <= stop (step,       start <= t < stop (step);
#                 sine, Synaptic*, EPSG)         start <= t <= stop (sine,
#                                                Synaptic*)
#   ramp          start <= t <= stop, and only   start <= t <= stop
#                 while 5 < t <= stop + 5
#                 (hardcoded 5); /1000 then
#                 *1e3 round trip
#   ramp2         rises from start               rises from t = 5 (hardcoded)
#   EPSG          cut off at stop                never cut off
#   EPSG, pair    driving force at the input     driving force at the soma,
#                 node                           even with axonal input
#   Synaptic*     driving force at the soma      driving force at the soma
#   SynapticPair  -                              second input zeroed once
#                                                t + diff >= t_end


# One return per stimulus, in msoAxon.m's if/elseif order, so each branch
# can be checked against the table above.
def external_current(  # noqa: PLR0911, PLR0912
    t: float,
    V: np.ndarray,
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    sa: np.ndarray = C.SA,
) -> tuple[int, float]:
    """Returns the input current density as (0-based compartment, value).

    Only one compartment ever receives input. Negative is depolarising
    (it enters dV with a minus). "none" means no stimulus. See the table
    above for how this differs from two.applied_current.

    Args:
        t: Time [ms].
        V: Voltage per compartment [mV].
        stim_type: The stimulus name, or "none".
        s: The stimulus settings (_common.stimulus).
        input_node: Where the stimulus enters, 1-indexed (step always
            enters the soma).
        sa: Membrane area per compartment [um^2].

    Returns:
        (k, current density [pA/um^2]); (0, 0.0) when the stimulus is
        off.
    """
    k = input_node - 1
    if stim_type == "ramp":
        # slope, t_top and t_end_local are derived once, in _common.stimulus
        # (msoAxon.m rebinds tEnd to t_end_local here)
        I0 = (
            float(s.start <= t <= s.t_top)
            * s.I
            * (t - s.start)
            * s.slope
            / 1000
        )
        if _solve.RAMP_RISE < t <= s.t_end_local:  # hardcoded 5 originally
            return k, -I0 * 1e3 / sa[0]
    elif stim_type == "ramp2":
        I0 = min(
            float(t >= s.start) * s.I * (t - s.start) * s.slope / 1000,
            s.I / 1000,
        )
        return k, -I0 * 1e3 / sa[0]
    elif stim_type == "step":
        if s.start < t <= s.stop:
            return 0, -s.I / sa[0]  # always the soma, regardless of input_node
    elif stim_type == "sine":
        wave = np.sin(s.two_pi_f * (t - s.start) / 1000)
        if s.start < t <= s.stop:
            return k, -(s.I * wave * (wave > 0)) / sa[0]
    elif stim_type == "Synaptic":
        if s.start < t <= s.stop:
            g = synaptic.interp_g(s.t_syn, s.g_syn, t)
            return k, g * (V[0] - s.VsynE) / sa[0]
    elif stim_type == "SynapticPair":
        if s.start < t <= s.stop:
            g = synaptic.interp_g(s.t_syn, s.g_syn, t) + synaptic.interp_g(
                s.t_syn, s.g_syn, t + s.diff
            )
            return k, g * (V[0] - s.VsynE) / sa[0]
    elif stim_type == "EPSG":
        te = t - s.start
        if s.start < t <= s.stop:
            return k, s.I * (0 - V[k]) * float(te >= 0) * synaptic.epsg_unitary(
                te, s.epsg_tau
            ) / -sa[k]
    elif stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        wave = float(te >= 0) * synaptic.epsg_unitary(
            te, s.epsg_tau, s.inv_peak
        ) + float(te >= td) * synaptic.epsg_unitary(
            te - td, s.epsg_tau, s.inv_peak
        )
        return k, s.I * (0 - V[k]) * wave / -sa[k]
    return 0, 0.0


def bilateral_current(
    t: float,
    V: np.ndarray,
    s: _common.Stimulus,
    a: int,
    b: int,
    sa: np.ndarray,
) -> tuple[tuple[int, float], tuple[int, float]]:
    """Returns the EPSGbilateral input currents.

    The first EPSG is at compartment a (0-based) at s.start, the second
    at b at s.stop.

    Args:
        t: Time [ms].
        V: Voltage per compartment [mV].
        s: The stimulus settings (_common.stimulus).
        a: The first input's compartment, 0-based.
        b: The second input's compartment, 0-based.
        sa: Membrane area per compartment [um^2].

    Returns:
        ((a, density), (b, density)), densities in pA/um^2.
    """
    te = t - s.start
    td = s.stop - s.start
    ia = (
        s.I
        * (0 - V[a])
        * float(te >= 0)
        * synaptic.epsg_unitary(te, s.epsg_tau, s.inv_peak)
        / -sa[a]
    )
    ib = (
        s.I
        * (0 - V[b])
        * float(te >= td)
        * synaptic.epsg_unitary(te - td, s.epsg_tau, s.inv_peak)
        / -sa[b]
    )
    return (a, ia), (b, ib)


def membrane(
    v0: float,
    soma_klt_scale: float = 1.0,
    ais_klt_scale: float = 1.0,
    soma_na_vhalf: float = -62.5,
    rebalance_rest: bool = False,
    morph: Morphology | None = None,
    dendrite_klt_scale: float = 1.0,
) -> Membrane:
    """Builds per-compartment channel overrides for the multi model.

    Defaults reproduce msoAxon.m exactly. The knobs are the ones the
    mature-MSO literature points to.

    Args:
        v0: The resting potential the leak reversals are built from
            [mV].
        soma_klt_scale: Multiplies gKLT at the soma. Kv1 channels rise
            ~4x from P14 to P23 and set the small mature somatic spike
            (Scott et al 2005, J Neurosci 25:7887).
        ais_klt_scale: Multiplies gKLT at the two AIS compartments.
        soma_na_vhalf: Somatic Na steady-state inactivation midpoint
            [mV]; Scott et al 2010 (J Neurosci 30:2039) measured -77 at
            the soma.
        rebalance_rest: Set each compartment's leak reversal so the net
            current at v0 is zero with gating at steady state, and start
            from that steady state. msoAxon.m does neither, so stronger
            KLT would otherwise shift rest.
        morph: The morphology these arrays are for (default: msoAxon.m's
            chain).
        dendrite_klt_scale: Multiplies gKLT in the dendrite compartments
            (with_dendrites only). Mathews et al 2010 found dendritic Kv1
            sharpens EPSPs.

    Returns:
        A Membrane with key, morph_key, morph, g_na, g_klt, na_vhalf,
        vlk and y0_gates (None unless rebalance_rest).

    Raises:
        ValueError: If dendrite_klt_scale is set for a morphology
            without dendrites.
    """
    morph = morph or LUMPED
    if dendrite_klt_scale != 1.0 and morph.n == N:
        raise ValueError(
            "dendrite_klt_scale needs a morphology with dendrites "
            "(morph=with_dendrites())"
        )
    g_klt = morph.g_klt.copy()
    g_klt[0] *= soma_klt_scale
    g_klt[1:3] *= ais_klt_scale
    g_klt[N:] *= dendrite_klt_scale
    # stored negated, as C.hinf's vhalf (exp((V + vhalf) / 7.77)) expects
    na_vhalf = np.full(morph.n, 62.5)
    na_vhalf[0] = -soma_na_vhalf
    vlk = np.full(morph.n, float(v0))
    y0_gates = None
    if rebalance_rest:
        V = np.full(morph.n, float(v0))
        gates = (
            C.minf(V),
            C.hinf(V, na_vhalf),
            C.pinf(V),
            C.winf(V),
            C.zinf(V),
            C.ainf(V),
        )
        mi, hi, pi, wi, zi, ai = gates
        I_ion = (
            C.na_open(mi, hi, morph.g_na) * (V - V_NA)
            + morph.g_kht * pi * (V - V_K)
            + g_klt * wi**4 * zi * (V - V_K)
            + morph.g_h * ai * (V - V_H)
        )
        vlk = V + I_ion / morph.g_lk
        y0_gates = gates
    return Membrane(
        # v0 is in the key because the leak reversals (vlk) are built from
        # it: without it, membrane(-60) and membrane(-68) would share a
        # pre-stimulus cache entry
        key=(
            float(v0),
            float(soma_klt_scale),
            float(ais_klt_scale),
            float(soma_na_vhalf),
            bool(rebalance_rest),
            float(dendrite_klt_scale),
        ),
        morph_key=morph.key,
        morph=morph,
        g_na=morph.g_na,
        g_klt=g_klt,
        na_vhalf=na_vhalf,
        vlk=vlk,
        y0_gates=y0_gates,
    )


def _dv(
    t: float,
    x: np.ndarray,
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    mem: Membrane,
    morph: Morphology = LUMPED,
    input_node2: int | None = None,
) -> np.ndarray:
    """Computes the voltage derivatives alone: _rhs's first n entries.

    For the settled stop, which needs only the sign of every dV/dt: it
    skips the gates' derivatives, most of _rhs's exp calls. A copy of
    _rhs's voltage block in the same float order rather than a helper
    _rhs calls, which would cost every right-hand side call about a
    microsecond; tests/test_settled.py checks the two agree bit for bit,
    so keep them in step. Which gates evolve doesn't matter here: the
    currents use every gate's value.

    Returns:
        dV/dt per compartment [mV/ms].
    """
    V, m, h, p, w, z, a = x.reshape(7, morph.n)
    INa = C.na_open(m, h, mem.g_na) * (V - V_NA)
    IKHT = morph.g_kht * p * (V - V_K)
    IKLT = mem.g_klt * w**4 * z * (V - V_K)
    Ih = morph.g_h * a * (V - V_H)
    Ilk = morph.g_lk * (V - mem.vlk)
    total = INa + IKHT + IKLT + Ih + Ilk
    if stim_type == "EPSGbilateral":
        for k, iext in bilateral_current(
            t, V, s, input_node - 1, input_node2 - 1, morph.sa
        ):
            total[k] += iext
    else:
        k, iext = external_current(t, V, stim_type, s, input_node, morph.sa)
        total[k] += iext
    return -(total + axial_current(V, morph)) / morph.cap


def _rhs(
    t: float,
    x: np.ndarray,
    stim_type: str,
    s: _common.Stimulus,
    input_node: int,
    active: tuple[bool, ...],
    mem: Membrane,
    morph: Morphology = LUMPED,
    input_node2: int | None = None,
) -> np.ndarray:
    """Computes the right-hand side.

    Arithmetic keeps msoAxon.m's expressions and their order, so the
    golden traces stay stable; the speed comes from writing into one
    output array and skipping gates that are switched off.
    """
    n = morph.n
    V, m, h, p, w, z, a = x.reshape(7, n)
    out = np.empty_like(x)
    d = out.reshape(7, n)

    INa = C.na_open(m, h, mem.g_na) * (V - V_NA)
    IKHT = morph.g_kht * p * (V - V_K)
    IKLT = mem.g_klt * w**4 * z * (V - V_K)
    # linear in activation, as in RM03 (no power given in Lehnert or Baumann)
    Ih = morph.g_h * a * (V - V_H)
    Ilk = morph.g_lk * (
        V - mem.vlk
    )  # leak reversal = resting potential by default

    total = INa + IKHT + IKLT + Ih + Ilk
    if stim_type == "EPSGbilateral":
        for k, iext in bilateral_current(
            t, V, s, input_node - 1, input_node2 - 1, morph.sa
        ):
            total[k] += iext
    else:
        k, iext = external_current(t, V, stim_type, s, input_node, morph.sa)
        total[k] += iext  # the other compartments would add an exact 0.0
    d[0] = -(total + axial_current(V, morph)) / morph.cap

    act_m, act_h, act_p, act_w, act_z, act_a = active
    # sources for each time constant are in constants.py
    d[1] = (C.minf(V) - m) / C.taum(V) if act_m else 0.0
    d[2] = (
        (C.hinf(V, mem.na_vhalf) - h) / C.tauh(V) if act_h else 0.0
    )  # per-compartment midpoint
    d[3] = (C.pinf(V) - p) / C.taup(V) if act_p else 0.0
    d[4] = (C.winf(V) - w) / C.tauw(V) if act_w else 0.0
    d[5] = (C.zinf(V) - z) / C.tauz(V) if act_z else 0.0
    d[6] = (C.ainf(V) - a) / C.taua(V) if act_a else 0.0
    return out


def mso_axon(
    stim_type: str,
    start: float,
    stop: float,
    I: float,
    node: int,
    model_type: str,
    t_end: float,
    v0: float,
    input_node: int,
    syn: synaptic.SynParams | None = None,
    max_step: float | None = None,
    stop_on_spike: float | None = None,
    mem: Membrane | None = None,
    morph: Morphology | None = None,
    input_node2: int | None = None,
    stop_when_settled: bool = False,
    mode: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Runs the multi-compartment model.

    Compartments are 1-indexed and times in ms, as in msoAxon.m.

    Args:
        stim_type: One of MSO_STIM_TYPES. What `start`, `stop` and `I`
            mean depends on it:
            step, sine: current [pA] on from start to stop (step always
                into the soma);
            ramp: rises from 0 at start to I [pA] at stop, then off;
            ramp2: rises with slope I/(stop - start) from start, then
                holds at I [pA];
            EPSG: peak conductance [nS] of the unitary EPSG waveform,
                onset at start (in this model it is also cut off at
                stop);
            EPSGpair: two EPSGs at input_node, onsets at start and stop;
            EPSGbilateral: first EPSG at input_node at start, second at
                input_node2 at stop;
            Synaptic, SynapticPair: I is ignored; the conductance comes
                from `syn`, applied from start to stop.
        start: See stim_type [ms].
        stop: See stim_type [ms].
        I: See stim_type.
        node: The compartment stop_on_spike compares with the soma
            (msoAxon.m itself ignores it). node=1 is the soma, so
            stop_on_spike would never fire.
        model_type: Which channels are active (ACTIVE_GATES).
        t_end: End time [ms].
        v0: Initial and (unless mem says otherwise) leak reversal
            potential [mV].
        input_node: Where the stimulus enters.
        syn: Synaptic waveform settings (SynParams), also the EPSG time
            constants.
        max_step: Largest solver step; defaults to 0.1 * t_end, ode15s's
            default MaxStep.
        stop_on_spike: If given (mV), stop as soon as compartment `node`
            rises that far above the soma; t then ends before t_end.
        mem: Channel overrides from membrane() (default: msoAxon.m's
            own).
        morph: LUMPED or with_dendrites(); defaults to mem's morphology,
            else LUMPED.
        input_node2: Second input site for "EPSGbilateral" (first EPSG
            at input_node at `start`, second at input_node2 at `stop`),
            e.g. the middle compartment of each dendrite.
        stop_when_settled: With stop_on_spike, also stop once the input
            is over and no spike can follow (_solve.settled_event), for
            callers that only need yes/no.
        mode: "exact" runs the numpy right-hand side (_rhs, the
            reference); "fast" runs its Rust port (the `fast` extra),
            which issues a DriftWarning on first use. None (default):
            the MSOAXON_MODE environment variable if set, else "fast"
            when the extension is installed, else "exact". See _fast.

    Returns:
        (t, y) with y shaped (n_times, 7 * n); see the module docstring
        for the state layout.

    Raises:
        ValueError: For unknown names, out-of-range compartments, a
            ramp with stop <= start, a missing input_node2, a mem built
            for another morphology, or an unknown mode.
    """
    morph = morph or (mem.morph if mem is not None else LUMPED)
    _common.check_args(
        stim_type,
        model_type,
        node,
        input_node,
        min_node=1,
        start=start,
        stop=stop,
        stim_types=MSO_STIM_TYPES,
        n_max=morph.n,
    )
    if stim_type == "EPSGbilateral" and not (
        input_node2 and 1 <= input_node2 <= morph.n
    ):
        raise ValueError(
            f"EPSGbilateral needs input_node2 in 1..{morph.n}, "
            f"got {input_node2}"
        )
    if model_type not in ACTIVE_GATES:
        raise ValueError(f"unknown model type {model_type!r}")
    syn = syn or synaptic.SynParams(t_end=t_end)

    s = _common.stimulus(stim_type, start, stop, I, t_end, syn)
    active = tuple(bool(g) for g in ACTIVE_GATES[model_type])
    mem = mem or membrane(v0, morph=morph)
    if mem.morph_key != morph.key:
        raise ValueError(
            "mem was built for a different morphology; "
            "pass morph= to membrane()"
        )
    n = morph.n
    if mem.y0_gates is None:
        y0 = np.concatenate(
            [np.full(n, float(v0))] + [np.full(n, g) for g in _GATE0]
        )
    else:
        y0 = np.concatenate([np.full(n, float(v0)), *mem.y0_gates])
    cuts = _solve.breakpoints(stim_type, start, stop, t_end)
    mode = _fast.resolve(mode)
    if mode == "fast":
        rhs = _fast.build(
            stim_type, s, input_node, active, mem, morph, input_node2
        )
        quiet_rhs = _fast.build("none", s, input_node, active, mem, morph)
        # the settled check slices the Rust RHS: at about 3.5 us a call it
        # costs under 2% of a ramp search, so a dV-only port would not pay
        dv = lambda t, x: rhs(t, x)[:n]  # noqa: E731
    else:
        rhs = lambda t, x: _rhs(  # noqa: E731
            t, x, stim_type, s, input_node, active, mem, morph, input_node2
        )
        quiet_rhs = lambda t, x: _rhs(  # noqa: E731
            t, x, "none", s, input_node, active, mem, morph
        )
        dv = lambda t, x: _dv(  # noqa: E731
            t, x, stim_type, s, input_node, mem, morph, input_node2
        )
    quiet = None
    # unlike two_cpt, step shares the prefix: msoAxon.m switches stimuli on for
    # t > start (the EPSG waveforms are zero at onset), so the equations at
    # t == start are still stimulus-free. The mode is part of the key: the two
    # right-hand sides need not agree in the last bits.
    if _solve.pre_stimulus_is_quiet(cuts, start, stop):
        quiet = (
            ("mso", model_type, float(v0), mem.key, morph.key, mode),
            quiet_rhs,
        )
    return _solve.integrate(
        rhs,
        y0,
        t_end,
        cuts,
        rtol=RTOL,
        atol=ATOL,
        max_step=max_step or MAX_STEP_FRACTION * t_end,
        jac_sparsity=morph.jac,
        quiet=quiet,
        stop_event=_common._stops(
            dv, node - 1, stop_on_spike, stop_when_settled, stim_type, s
        ),
    )
