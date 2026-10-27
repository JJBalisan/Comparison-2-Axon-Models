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
import types
import warnings
from typing import Any

import numpy as np
from scipy import optimize, sparse

from msoaxon import _common, _fast, _solve, constants as C, synaptic

N = C.N_CPT

# specific capacitance, uF/cm^2 converted to ms*nS/um^2
CAP = np.concatenate([[0.8, 0.8, 0.8], np.tile([0.01, 0.8], 21)]) * 1e6 * 1e-8
V_NA, V_K, V_H = 69.0, -90.0, -35.0
R_AXIAL = 100.0  # specific axial resistivity [Ohm cm]

# conductance between neighbouring compartments i and i+1 (the (2/R)/(...) term)
_G_AX = (2 / R_AXIAL) / (C.L_CM[:-1] / C.XA_CM[:-1] + C.L_CM[1:] / C.XA_CM[1:])

MSO_STIM_TYPES = _common.STIM_TYPES + ("EPSGbilateral",)

# msoAxon.m's resting potential [mV]: where membrane(rebalance_rest=True)
# holds the axon (internodes and nodes) unless axon_rest says otherwise
AXON_REST = -68.0
_AXON_LABELS = ("internode", "node")
# largest |dV/dt| [mV/ms] and |dgate/dt| [1/ms] membrane() accepts at the
# steady state it root-finds for a rebalanced rest
_REST_RESIDUAL = 1e-9


class UnstableRestWarning(UserWarning):
    """A rebalanced resting state is unstable: the cell fires on its own.

    membrane(rebalance_rest=True) linearises the full right-hand side at
    the rest it builds and issues this when the leading eigenvalue has a
    positive real part (#128).
    """


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
    never reach a generated __eq__. The arrays are read-only once built
    (_freeze): a Morphology is shared by every run that uses it and the
    quiet-prefix cache and spike heights are keyed by `key`, so a
    mutated bundle would serve stale results.

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
            override the leak reversal, and with klt_to_leak builds a
            copy of the morphology with a larger leak.
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

    def __post_init__(self) -> None:
        """Freezes the arrays: a bundle is shared and keyed by its settings."""
        _freeze(self)

    lateral: np.ndarray | None = None
    medial: np.ndarray | None = None


@dataclasses.dataclass(frozen=True, eq=False)
class Membrane:
    """Per-compartment channel overrides for one morphology (membrane()).

    eq=False keeps identity equality and hashing.

    Attributes:
        key: Every setting membrane() was built from, part of the
            pre-stimulus cache key. Not klt_to_leak: that replaces the
            morphology (morph.key, which morph_key also carries and every
            cache that holds `key` holds too), so the keys stay distinct.
        morph_key: The morphology's key.
        morph: The morphology these arrays are for.
        g_na: Na density per compartment [nS/um^2].
        g_klt: KLT density per compartment [nS/um^2].
        na_vhalf: Negated Na inactivation midpoint per compartment
            [mV], as C.hinf's vhalf expects.
        vlk: Leak reversal per compartment [mV].
        y0_gates: Initial (m, h, p, w, z, a) steady-state blocks, each
            of length n; None unless rebalance_rest.
        v_rest: The root-found resting potential per compartment [mV],
            where initial_state starts V; None when every compartment
            rests at v0 (no rebalance, or the soma and the axon
            rebalanced to the same potential), which starts V at v0.
    """

    key: tuple[float | bool, ...]
    morph_key: str | tuple[Any, ...]
    morph: Morphology
    g_na: np.ndarray
    g_klt: np.ndarray
    na_vhalf: np.ndarray
    vlk: np.ndarray
    y0_gates: tuple[np.ndarray, ...] | None
    v_rest: np.ndarray | None = None

    def __post_init__(self) -> None:
        """Freezes the arrays (see Morphology)."""
        _freeze(self)


def _freeze(bundle: object) -> None:
    """Sets every ndarray field of a dataclass bundle read-only.

    Tuples of arrays (Membrane.y0_gates) are frozen element by element.
    """
    for field in dataclasses.fields(bundle):
        value = getattr(bundle, field.name)
        arrays = value if isinstance(value, tuple) else (value,)
        for a in arrays:
            if isinstance(a, np.ndarray):
                a.flags.writeable = False


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
#   Synaptic*     driving force at the input     driving force at the soma
#                 node (fixed in PR #103; the
#                 area is still the soma's)
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
            return k, g * (V[k] - s.VsynE) / sa[0]
    elif stim_type == "SynapticPair":
        if s.start < t <= s.stop:
            g = synaptic.interp_g(s.t_syn, s.g_syn, t) + synaptic.interp_g(
                s.t_syn, s.g_syn, t + s.diff
            )
            return k, g * (V[k] - s.VsynE) / sa[0]
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


def events_current(
    t: float, V: np.ndarray, s: _common.Stimulus, sa: np.ndarray
) -> np.ndarray:
    """Returns the "Events" input current density per compartment.

    Each event's current (synaptic.event_currents, with the driving
    force at its own compartment, as EPSG and EPSGpair have it) divided
    by its compartment's area, summed per compartment; negative
    depolarises, as external_current's values do.

    Args:
        t: Time [ms].
        V: Voltage per compartment [mV].
        s: The stimulus settings (_common.stimulus), with s.events.
        sa: Membrane area per compartment [um^2].

    Returns:
        The densities [pA/um^2], one per compartment (0.0 where no event
        enters).
    """
    ev = s.events
    k = ev["k"]
    i = synaptic.event_currents(t, V[k], ev, s.I)
    return np.bincount(k, weights=-i / sa[k], minlength=len(sa))


def _klt_to_leak(
    kind: str,
    g_klt: np.ndarray,
    morph: Morphology,
    vlk: np.ndarray,
    v0: float,
    key: tuple[float, ...],
    v_rest: np.ndarray | None = None,
) -> tuple[np.ndarray, Morphology, np.ndarray]:
    """Moves KLT into the leak (membrane's klt_to_leak).

    Args:
        kind: "chord" or "slope".
        g_klt: gKLT per compartment, scale factors applied [nS/um^2].
        morph: The morphology whose leak grows.
        vlk: The leak reversal before the substitution [mV].
        v0: The resting potential [mV].
        key: The scale factors g_klt was built with, for the new
            morphology's key (the leak depends on all of them).
        v_rest: Each compartment's own rest (Membrane.v_rest), where the
            substitution is made instead of v0; None: v0 everywhere.

    Returns:
        (zero gKLT, a copy of morph with the larger leak and its own key,
        the new leak reversal).

    Raises:
        ValueError: For "slope", if the slope conductance of KLT at v0 is
            not positive (v0 well below E_K).
    """
    V = np.full(morph.n, v0) if v_rest is None else v_rest

    def current(v: np.ndarray) -> np.ndarray:
        """Steady-state KLT current density at v [pA/um^2]."""
        return g_klt * C.winf(v) ** 4 * C.zinf(v) * (v - V_K)

    if kind == "chord":
        g_eff, e = g_klt * C.winf(V) ** 4 * C.zinf(V), V_K
    else:
        h = 1e-3  # central difference [mV]
        g_eff = (current(V + h) - current(V - h)) / (2 * h)
        if ((g_klt > 0) & (g_eff <= 0)).any():
            raise ValueError(
                f"klt_to_leak='slope' needs a positive KLT slope conductance "
                f"at v0={v0} mV, which it lacks below about E_K ({V_K} mV)"
            )
        # compartments without KLT have nothing to move (g_eff == 0)
        e = V - current(V) / np.where(g_eff > 0, g_eff, 1.0)
    g_lk = morph.g_lk + g_eff
    # one linear leak standing for the old leak plus g_eff * (V - e)
    vlk = (morph.g_lk * vlk + g_eff * e) / g_lk
    morph = dataclasses.replace(
        morph, key=("klt_to_leak", kind, morph.key, v0, *key), g_lk=g_lk
    )
    return np.zeros(morph.n), morph, vlk


def membrane(
    v0: float,
    soma_klt_scale: float = 1.0,
    ais_klt_scale: float = 1.0,
    soma_na_vhalf: float = -62.5,
    rebalance_rest: bool = False,
    morph: Morphology | None = None,
    dendrite_klt_scale: float = 1.0,
    klt_to_leak: str | None = None,
    axon_rest: float | None = None,
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
        rebalance_rest: Set the leak reversals so the soma rests at v0
            with gating at steady state, and start from that steady
            state. msoAxon.m does neither, so stronger KLT would
            otherwise shift rest. Only the soma, the AIS and any
            dendrites are rebalanced to v0; the axon (internodes and
            nodes) is rebalanced to axon_rest. When the two differ, the
            axial current flows at rest, so the steady state is
            root-found: the soma's leak reversal is solved to hold it at
            v0 exactly, and V starts from the solved profile
            (Membrane.v_rest). The full linearisation at the rest is
            checked, and an UnstableRestWarning is issued when its
            leading eigenvalue has a positive real part: the nodes of
            Ranvier fire on their own once they rest above about -60.5
            mV (#128), which is why the axon is not moved to v0 by
            default.
        morph: The morphology these arrays are for (default: msoAxon.m's
            chain).
        dendrite_klt_scale: Multiplies gKLT in the dendrite compartments
            (with_dendrites only). Mathews et al 2010 found dendritic Kv1
            sharpens EPSPs.
        klt_to_leak: Replace KLT in every compartment by a leak, as the
            dynamic-clamp control of Mikiel-Hunter et al 2016: gKLT is
            set to zero and the compartment's leak conductance grows by
            the conductance g_eff of a linear current g_eff * (V - e)
            with the KLT current at v0 (after the scale factors above),
            so the rest potential is unchanged (exactly, with
            rebalance_rest). None (default) keeps KLT. "chord": the
            KLT conductance at v0 with e = E_K, KLT frozen at rest; its
            steep growth with depolarisation is lost, so the steady
            input resistance rises (2.4 to 9.1 MOhm). "slope": the
            slope of the steady-state KLT current at v0 (the tangent),
            which keeps the steady input resistance. The leak lives in
            the morphology, so the returned Membrane carries a copy
            of it with its own key (pass mem alone, or
            morph=mem.morph). With a rebalanced rest that is not uniform,
            each compartment's KLT is replaced at its own rest.
        axon_rest: With rebalance_rest, the potential the axon's
            internodes and nodes are rebalanced to [mV]. None (default)
            is msoAxon.m's AXON_REST, -68 mV, so at v0 = -68 the whole
            cable is rebalanced to v0 as before. axon_rest=v0 rebalances
            the whole cable to v0, which is unstable above about -60.5
            mV.

    Returns:
        A Membrane with key, morph_key, morph, g_na, g_klt, na_vhalf,
        vlk, y0_gates (None unless rebalance_rest) and v_rest (None
        unless the rebalanced rest is not uniform).

    Raises:
        ValueError: If dendrite_klt_scale is set for a morphology
            without dendrites, klt_to_leak is not None, "chord" or
            "slope", "slope" at a v0 where KLT has no positive slope
            conductance, or axon_rest is set without rebalance_rest.
        RuntimeError: If the steady state of a rebalanced rest cannot
            be found.

    Warns:
        UnstableRestWarning: If a rebalanced rest is unstable.
    """
    morph = morph or LUMPED
    if klt_to_leak not in (None, "chord", "slope"):
        raise ValueError(
            f'klt_to_leak must be None, "chord" or "slope", got {klt_to_leak!r}'
        )
    if axon_rest is not None and not rebalance_rest:
        raise ValueError("axon_rest needs rebalance_rest=True")
    axon_rest = AXON_REST if axon_rest is None else float(axon_rest)
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
    y0_gates = v_rest = None
    if rebalance_rest:
        V = np.full(morph.n, float(v0))
        V[np.isin(morph.labels, _AXON_LABELS)] = axon_rest
        # each compartment balanced on its own, as if V were uniform: exact
        # when it is (no axial current flows), else the starting guess
        gates = _gates_at(V, na_vhalf)
        vlk = V + _ionic(V, gates, morph.g_na, g_klt, morph) / morph.g_lk
        if axon_rest != v0:
            V, vlk = _solve_rest(V, vlk, morph, g_klt, na_vhalf)
            gates, v_rest = _gates_at(V, na_vhalf), V
        y0_gates = gates
    if klt_to_leak:
        g_klt, morph, vlk = _klt_to_leak(
            klt_to_leak,
            g_klt,
            morph,
            vlk,
            float(v0),
            key=(
                float(soma_klt_scale),
                float(ais_klt_scale),
                float(dendrite_klt_scale),
            ),
            v_rest=v_rest,
        )
    mem = Membrane(
        # v0 is in the key because the leak reversals (vlk) are built from
        # it: without it, membrane(-60) and membrane(-68) would share a
        # pre-stimulus cache entry; axon_rest likewise
        key=(
            float(v0),
            float(soma_klt_scale),
            float(ais_klt_scale),
            float(soma_na_vhalf),
            bool(rebalance_rest),
            float(dendrite_klt_scale),
            axon_rest,
        ),
        morph_key=morph.key,
        morph=morph,
        g_na=morph.g_na,
        g_klt=g_klt,
        na_vhalf=na_vhalf,
        vlk=vlk,
        y0_gates=y0_gates,
        v_rest=v_rest,
    )
    if rebalance_rest:
        _warn_if_unstable(mem, float(v0), axon_rest)
    return mem


def _gates_at(V: np.ndarray, na_vhalf: np.ndarray) -> tuple[np.ndarray, ...]:
    """Returns the steady-state gates (m, h, p, w, z, a) at V."""
    return (
        C.minf(V),
        C.hinf(V, na_vhalf),
        C.pinf(V),
        C.winf(V),
        C.zinf(V),
        C.ainf(V),
    )


def _ionic(
    V: np.ndarray,
    gates: tuple[np.ndarray, ...],
    g_na: np.ndarray,
    g_klt: np.ndarray,
    morph: Morphology,
) -> np.ndarray:
    """Returns the ionic current density but the leak [pA/um^2].

    In membrane()'s original float order, which the rebalanced golden
    traces depend on.
    """
    mi, hi, pi, wi, zi, ai = gates
    return (
        C.na_open(mi, hi, g_na) * (V - V_NA)
        + morph.g_kht * pi * (V - V_K)
        + g_klt * wi**4 * zi * (V - V_K)
        + morph.g_h * ai * (V - V_H)
    )


def _solve_rest(
    V: np.ndarray,
    vlk: np.ndarray,
    morph: Morphology,
    g_klt: np.ndarray,
    na_vhalf: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Root-finds a rebalanced rest that is not uniform.

    The unknowns are the soma's leak reversal and every other
    compartment's V; the soma's V stays at V[0] (v0) and the gates
    follow V at steady state, so the zero of _dv is a steady state of
    _rhs with the soma at v0 exactly.

    Args:
        V: The rest each compartment was balanced to [mV]; V[0] is v0.
        vlk: Their leak reversals, balanced compartment by compartment
            (vlk[0] is the starting guess) [mV].
        morph: The morphology.
        g_klt: gKLT per compartment [nS/um^2].
        na_vhalf: Negated Na inactivation midpoint per compartment [mV].

    Returns:
        (V, vlk) at the steady state.

    Raises:
        RuntimeError: If the solver does not reach the steady state.
    """
    v0 = V[0]

    def unpack(u: np.ndarray) -> tuple[np.ndarray, types.SimpleNamespace]:
        """Returns the state and a stand-in Membrane for unknowns u."""
        v = np.concatenate([[v0], u[1:]])
        lk = np.concatenate([u[:1], vlk[1:]])
        x = np.concatenate([v, *_gates_at(v, na_vhalf)])
        mem = types.SimpleNamespace(
            g_na=morph.g_na, g_klt=g_klt, na_vhalf=na_vhalf, vlk=lk
        )
        return x, mem

    def residual(u: np.ndarray) -> np.ndarray:
        """Returns dV/dt at the state the unknowns u describe."""
        x, mem = unpack(u)
        return _dv(0.0, x, "none", None, 1, mem, morph)

    sol = optimize.root(
        residual,
        np.concatenate([vlk[:1], V[1:]]),
        method="hybr",
        options={"xtol": 1e-14},
    )
    x, mem = unpack(sol.x)
    worst = np.abs(_rhs(0.0, x, "none", None, 1, (True,) * 6, mem, morph)).max()
    if not worst <= _REST_RESIDUAL:
        raise RuntimeError(
            f"no steady state for a rebalanced rest at v0={v0} mV "
            f"({sol.message}; largest derivative {worst:.3g})"
        )
    return x[: morph.n], mem.vlk


def rest_jacobian(mem: Membrane) -> tuple[np.ndarray, np.ndarray]:
    """Returns the Jacobian of the quiet right-hand side at mem's rest.

    Every gate evolves (active-full), and the rest is mem's starting
    state (initial_state), which is a steady state only when mem was
    built with rebalance_rest. Central differences, steps of 1e-6
    relative (at least 1e-6 absolute).

    Args:
        mem: The membrane.

    Returns:
        (J, x): the 7n x 7n Jacobian [1/ms] and the state it is taken
        at.
    """
    morph = mem.morph
    x = initial_state(morph, mem, mem.key[0])
    active = (True,) * 6
    J = np.empty((x.size, x.size))
    for j in range(x.size):
        h = 1e-6 * max(1.0, abs(x[j]))
        up, down = x.copy(), x.copy()
        up[j] += h
        down[j] -= h
        J[:, j] = (
            _rhs(0.0, up, "none", None, 1, active, mem, morph)
            - _rhs(0.0, down, "none", None, 1, active, mem, morph)
        ) / (2 * h)
    return J, x


def _warn_if_unstable(mem: Membrane, v0: float, axon_rest: float) -> None:
    """Issues UnstableRestWarning if mem's rebalanced rest is unstable.

    Args:
        mem: A membrane built with rebalance_rest.
        v0: The soma's rest [mV].
        axon_rest: The axon's [mV].
    """
    J, x = rest_jacobian(mem)
    lam, vec = np.linalg.eig(J)
    i = int(np.argmax(lam.real))
    if lam[i].real <= 0:
        return
    n = mem.morph.n
    k = int(np.argmax(np.abs(vec[:n, i])))  # where the growing mode lives
    family = mem.morph.labels[k]
    rests = x[:n][np.array(mem.morph.labels) == family]
    warnings.warn(
        f"the rebalanced rest (soma {v0} mV, axon {axon_rest} mV) is "
        f"unstable: the leading eigenvalue is {lam[i].real:+.3f} "
        f"+/- {abs(lam[i].imag):.3f}j /ms, and its mode is largest in "
        f"the {family} compartments (most at compartment {k + 1}), which "
        f"rest at {rests.min():.1f} to {rests.max():.1f} mV. The cell will "
        f"fire with no input; nodes of Ranvier rebalanced above about "
        f"-60.5 mV do (#128).",
        UnstableRestWarning,
        stacklevel=3,
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
    elif stim_type == "Events":
        total += events_current(t, V, s, morph.sa)
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
    elif stim_type == "Events":
        total += events_current(t, V, s, morph.sa)
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


def initial_state(morph: Morphology, mem: Membrane, v0: float) -> np.ndarray:
    """Returns mso_axon's initial state: V at v0, then the six gates.

    Args:
        morph: The morphology (its n sets the block length).
        mem: The channel overrides; its y0_gates (rebalance_rest) are the
            starting gates, else msoAxon.m's rounded _GATE0, and its
            v_rest, if set, the starting V instead of v0 (a rebalanced
            rest that is not uniform: the root-found steady state).
        v0: The initial potential [mV].

    Returns:
        The state vector, 7 * morph.n long, in the module docstring's
        layout.
    """
    n = morph.n
    if mem.y0_gates is None:
        return np.concatenate(
            [np.full(n, float(v0))] + [np.full(n, g) for g in _GATE0]
        )
    if mem.v_rest is not None:
        return np.concatenate([mem.v_rest, *mem.y0_gates])
    return np.concatenate([np.full(n, float(v0)), *mem.y0_gates])


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
                from `syn`, applied from start to stop;
            Events: the event list in `syn.events` (synaptic.Event: any
                number of EPSG-shaped conductances, each with its own
                onset, compartment, peak, kinetics and reversal, never
                cut off); I scales every peak (1.0 uses them as given);
                start, stop and input_node are ignored.
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
            reference); "fast" runs its Rust port (the `fast` extra)
            under scipy's BDF; "native" also integrates in Rust
            (diffsol's BDF), agreeing with the reference within the
            solver tolerance. Each issues a DriftWarning on first use.
            None (default): the MSOAXON_MODE environment variable if
            set, else "native" when the extension is installed, else
            "exact". See _fast.

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
    if stim_type == "Events":
        synaptic.check_events(tuple(syn.events), morph.n, t_end)

    s = _common.stimulus(stim_type, start, stop, I, t_end, syn)
    if stim_type == "Events":
        start, stop = s.start, s.stop  # the first and last onset
    active = tuple(bool(g) for g in ACTIVE_GATES[model_type])
    mem = mem or membrane(v0, morph=morph)
    if mem.morph_key != morph.key:
        raise ValueError(
            "mem was built for a different morphology; "
            "pass morph= to membrane()"
        )
    n = morph.n
    y0 = initial_state(morph, mem, v0)
    cuts = _solve.breakpoints(
        stim_type,
        start,
        stop,
        t_end,
        onsets=() if s.events is None else s.events["t"],
    )
    mode = _fast.resolve(mode)
    if mode in ("fast", "native"):
        # native: the same Rust right-hand side, its segments solved in Rust
        build = _fast.build if mode == "fast" else _fast.build_native
        rhs = build(stim_type, s, input_node, active, mem, morph, input_node2)
        quiet_rhs = build("none", s, input_node, active, mem, morph)
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
