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

from types import SimpleNamespace

import numpy as np
from scipy.sparse import bmat, coo_matrix, diags, eye

from . import constants as C
from ._common import STIM_TYPES, check_args, stimulus
from ._solve import breakpoints, integrate, pre_stimulus_is_quiet, spike_event
from .synaptic import SynParams, epsg_unitary, interp_g

N = C.N_CPT

# specific capacitance, uF/cm^2 converted to ms*nS/um^2
CAP = np.concatenate([[0.8, 0.8, 0.8], np.tile([0.01, 0.8], 21)]) * 1e6 * 1e-8
V_NA, V_K, V_H = 69.0, -90.0, -35.0
R_AXIAL = 100.0  # specific axial resistivity [Ohm cm]

# conductance between neighbouring compartments i and i+1 (the (2/R)/(...) term)
_G_AX = (2 / R_AXIAL) / (C.L_CM[:-1] / C.XA_CM[:-1] + C.L_CM[1:] / C.XA_CM[1:])

MSO_STIM_TYPES = STIM_TYPES + ("EPSGbilateral",)


# which gating blocks (m, h, p, w, z, a) evolve for each model type
ACTIVE_GATES = {
    "active-KLT":    (0, 0, 0, 1, 1, 0),
    "active-H":      (0, 0, 0, 0, 0, 1),
    "active-KLT+H":  (0, 0, 0, 1, 1, 1),
    "passive":       (0, 0, 0, 0, 0, 0),
    "active-sodium": (1, 1, 0, 0, 0, 0),
    "active-full":   (1, 1, 1, 1, 1, 1),
    "active-KHT":    (0, 0, 1, 0, 0, 0),
}

# initial gating values (steady state at -68 mV, as in msoAxon.m)
_GATE0 = (0.12, 0.67, 0.0, 0.28, 0.67, 0.22)


def _chain_jac_sparsity():
    """dV couples to neighbouring V and its own gates; each gate to its own V."""
    I = eye(N)
    tri = diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(N, N))
    blocks = [[tri] + [I] * 6] + [[I] + [I if j == i else None for j in range(6)]
                                  for i in range(6)]
    return bmat(blocks).tocsc()


def _tree_jac_sparsity(n, parent, child):
    I = eye(n, format="csr")
    adj = coo_matrix((np.ones(2 * len(parent)), (np.r_[parent, child], np.r_[child, parent])),
                     shape=(n, n))
    couple = (I + adj).tocsr()
    couple.data[:] = 1.0
    blocks = [[couple] + [I] * 6] + [[I] + [I if j == i else None for j in range(6)]
                                     for i in range(6)]
    return bmat(blocks).tocsc()


# The unbranched 45-compartment chain of msoAxon.m (soma, 2 AIS, 21 internode/node pairs)
LUMPED = SimpleNamespace(
    key="lumped", n=N, sa=C.SA, cap=CAP, g_na=C.G_NA, g_kht=C.G_KHT, g_klt=C.G_KLT,
    g_h=C.G_H, g_lk=C.G_LK, parent=np.arange(N - 1), child=np.arange(1, N), g_ax=_G_AX,
    chain=True, jac=_chain_jac_sparsity(), labels=["soma", "AIS", "AIS"] + ["internode", "node"] * 21)
for _a in (CAP, _G_AX, LUMPED.parent, LUMPED.child):  # read-only, as in constants.py
    _a.flags.writeable = False


def with_dendrites(length=200.0, diameter=5.0, n_seg=5, klt_lambda=74.0, conserve_totals=True,
                   dendrite_ra=R_AXIAL):
    """Lehnert et al 2014's dendritic variant (their Fig. 8), appended to the axon.

    Two identical unbranched dendrites (lateral = ipsilateral input, medial =
    contralateral) of `length` um and `diameter` um, `n_seg` compartments each,
    attached to the soma. The soma shrinks so total membrane stays 8750 um^2
    ("the somatic surface was reduced to 2467 um^2"); its Na density is scaled up
    so total Na is unchanged. Dendrites have no Na or KHT; KLT and h decay
    exponentially with distance from the soma (length constant 74 um, Mathews et
    al 2010), starting from the soma's densities; leak and capacitance as at the
    soma. The paper doesn't give axial resistivity, so the model's 100 Ohm cm is the
    default; dendrite_ra sets it for the dendrite compartments only (Mathews et al
    2010 used 200 Ohm cm for soma and dendrites). The soma and axon keep 100.

    The paper doesn't say whether total KLT and h were kept when they were spread
    along the dendrites; it does say so for Na. conserve_totals=True (default)
    scales soma + dendrite KLT and h so their totals equal the lumped model's; it
    reproduces the paper's finding that this variant tunes "almost identically" to
    the lumped model (EPSG-pair thresholds within 12%). conserve_totals=False starts
    the gradient at the soma's Table 2 density, which halves total KLT (thresholds
    up to 50% lower, somatic spike 21-30 mV).

    Compartments (1-indexed, as for input_node): 46..45+n_seg lateral dendrite,
    proximal to distal; the next n_seg the medial dendrite.
    """
    n_d = 2 * n_seg
    n = N + n_d
    seg = length / n_seg
    x_mid = (np.arange(n_seg) + 0.5) * seg  # distance of each compartment's centre [um]
    sa_d = np.pi * diameter * seg
    soma_sa = C.SA[0] - 2 * n_seg * sa_d
    if soma_sa <= 0:
        raise ValueError("dendrites larger than the lumped soma+dendrite membrane")

    sa = np.concatenate([[soma_sa], C.SA[1:], np.full(n_d, sa_d)])
    L = np.concatenate([[np.sqrt(soma_sa / np.pi)], C.L[1:], np.full(n_d, seg)])  # soma: as Constants.m
    r_cyl = sa / (2 * np.pi * L)
    xa_cm, l_cm = np.pi * r_cyl ** 2 * 1e-8, L * 1e-4
    decay = np.tile(np.exp(-x_mid / klt_lambda), 2)

    g_na = np.concatenate([[C.G_NA[0] * C.SA[0] / soma_sa], C.G_NA[1:], np.zeros(n_d)])
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
    g_ax = (2 / R_AXIAL) / (l_cm[parent] / xa_cm[parent] + l_cm[child] / xa_cm[child])
    if dendrite_ra != R_AXIAL:
        # only edges into a dendrite compartment change, so the axon's conductances
        # stay bit-identical; each half-compartment uses its own resistivity
        ra = np.concatenate([np.full(N, R_AXIAL), np.full(n_d, float(dendrite_ra))])
        into = child >= N
        p, c = parent[into], child[into]
        g_ax[into] = 2 / (ra[p] * l_cm[p] / xa_cm[p] + ra[c] * l_cm[c] / xa_cm[c])
    labels = LUMPED.labels + ["lateral dendrite"] * n_seg + ["medial dendrite"] * n_seg
    return SimpleNamespace(
        key=("dendrites", float(length), float(diameter), int(n_seg), float(klt_lambda),
             bool(conserve_totals), float(dendrite_ra)),
        n=n, sa=sa, cap=cap, g_na=g_na, g_kht=g_kht, g_klt=g_klt, g_h=g_h, g_lk=g_lk,
        parent=parent, child=child, g_ax=g_ax, chain=False,
        jac=_tree_jac_sparsity(n, parent, child), labels=labels,
        lateral=lat + 1, medial=med + 1)  # 1-indexed compartment numbers


def axial_current(V, morph=LUMPED):
    """Axial current density [pA/um^2], sign convention as in msoAxon.m."""
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


def external_current(t, V, stim_type, s, input_node, sa=C.SA):
    """Input current density [pA/um^2] as (0-based compartment, value).

    Only one compartment ever receives input. Negative is depolarising (it
    enters dV with a minus). "none" means no stimulus.
    """
    k = input_node - 1
    if stim_type == "ramp":
        slope = 1 / (s.stop - s.start)
        t_end_local = s.stop + 5  # msoAxon.m rebinds tEnd here
        t_top = s.start + 1 / slope
        I0 = float(s.start <= t <= t_top) * s.I * (t - s.start) * slope / 1000
        if 5 < t <= t_end_local:  # 5 is hardcoded in the original
            return k, -I0 * 1e3 / sa[0]
    elif stim_type == "ramp2":
        slope = 1 / (s.stop - s.start)
        I0 = min(float(t >= s.start) * s.I * (t - s.start) * slope / 1000, s.I / 1000)
        return k, -I0 * 1e3 / sa[0]
    elif stim_type == "step":
        if s.start < t <= s.stop:
            return 0, -s.I / sa[0]  # always the soma, regardless of input_node
    elif stim_type == "sine":
        wave = np.sin(2 * np.pi * s.f * (t - s.start) / 1000)
        if s.start < t <= s.stop:
            return k, -(s.I * wave * (wave > 0)) / sa[0]
    elif stim_type == "Synaptic":
        if s.start < t <= s.stop:
            g = interp_g(s.t_syn, s.g_syn, t)
            return k, g * (V[0] - s.VsynE) / sa[0]
    elif stim_type == "SynapticPair":
        if s.start < t <= s.stop:
            g = interp_g(s.t_syn, s.g_syn, t) + interp_g(s.t_syn, s.g_syn, t + s.diff)
            return k, g * (V[0] - s.VsynE) / sa[0]
    elif stim_type == "EPSG":
        te = t - s.start
        if s.start < t <= s.stop:
            return k, s.I * (0 - V[k]) * float(te >= 0) * epsg_unitary(te, s.epsg_tau) / -sa[k]
    elif stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        wave = (float(te >= 0) * epsg_unitary(te, s.epsg_tau)
                + float(te >= td) * epsg_unitary(te - td, s.epsg_tau))
        return k, s.I * (0 - V[k]) * wave / -sa[k]
    return 0, 0.0


def bilateral_current(t, V, s, a, b, sa):
    """EPSGbilateral: first EPSG at compartment a (0-based) at s.start, second at b at s.stop."""
    te = t - s.start
    td = s.stop - s.start
    ia = s.I * (0 - V[a]) * float(te >= 0) * epsg_unitary(te, s.epsg_tau) / -sa[a]
    ib = s.I * (0 - V[b]) * float(te >= td) * epsg_unitary(te - td, s.epsg_tau) / -sa[b]
    return (a, ia), (b, ib)


def membrane(v0, soma_klt_scale=1.0, ais_klt_scale=1.0, soma_na_vhalf=-62.5,
             rebalance_rest=False, morph=None, dendrite_klt_scale=1.0):
    """Per-compartment channel overrides for the multi-compartment model.

    Defaults reproduce msoAxon.m exactly. The knobs are the ones the mature-MSO
    literature points to:
    - soma_klt_scale / ais_klt_scale: multiply gKLT at the soma / the two AIS
      compartments. Kv1 channels rise ~4x from P14 to P23 and set the small mature
      somatic spike (Scott et al 2005, J Neurosci 25:7887).
    - soma_na_vhalf: somatic Na steady-state inactivation midpoint [mV]; Scott et
      al 2010 (J Neurosci 30:2039) measured -77 at the soma.
    - rebalance_rest: set each compartment's leak reversal so the net current at
      v0 is zero with gating at steady state, and start from that steady state.
      msoAxon.m does neither, so stronger KLT would otherwise shift rest.
    - morph: the morphology these arrays are for (default: msoAxon.m's chain).
    - dendrite_klt_scale: multiply gKLT in the dendrite compartments (with_dendrites
      only). Mathews et al 2010 found dendritic Kv1 sharpens EPSPs.
    """
    morph = morph or LUMPED
    if dendrite_klt_scale != 1.0 and morph.n == N:
        raise ValueError("dendrite_klt_scale needs a morphology with dendrites (morph=with_dendrites())")
    # v0 is in the key because the leak reversals (vlk) are built from it: without
    # it, membrane(-60) and membrane(-68) would share a pre-stimulus cache entry
    m = SimpleNamespace(key=(float(v0), float(soma_klt_scale), float(ais_klt_scale),
                             float(soma_na_vhalf), bool(rebalance_rest), float(dendrite_klt_scale)),
                        morph_key=morph.key, morph=morph)
    m.g_na = morph.g_na
    m.g_klt = morph.g_klt.copy()
    m.g_klt[0] *= soma_klt_scale
    m.g_klt[1:3] *= ais_klt_scale
    m.g_klt[N:] *= dendrite_klt_scale
    # stored negated, as C.hinf's vhalf (exp((V + vhalf) / 7.77)) expects
    m.na_vhalf = np.full(morph.n, 62.5)
    m.na_vhalf[0] = -soma_na_vhalf
    m.vlk = np.full(morph.n, float(v0))
    m.y0_gates = None
    if rebalance_rest:
        V = np.full(morph.n, float(v0))
        gates = (C.minf(V), C.hinf(V, m.na_vhalf), C.pinf(V),
                 C.winf(V), C.zinf(V), C.ainf(V))
        mi, hi, pi, wi, zi, ai = gates
        I_ion = (m.g_na * mi ** 4 * (0.993 * hi + 0.007) * (V - V_NA)
                 + morph.g_kht * pi * (V - V_K) + m.g_klt * wi ** 4 * zi * (V - V_K)
                 + morph.g_h * ai * (V - V_H))
        m.vlk = V + I_ion / morph.g_lk
        m.y0_gates = gates
    return m


def _rhs(t, x, stim_type, s, input_node, active, mem, morph=LUMPED, input_node2=None):
    """Right-hand side. Arithmetic is kept expression-for-expression identical to
    msoAxon.m's order so results are bit-for-bit stable; the speed comes from
    writing into one output array and skipping gates that are switched off."""
    n = morph.n
    V, m, h, p, w, z, a = x.reshape(7, n)
    out = np.empty_like(x)
    d = out.reshape(7, n)

    INa = mem.g_na * m ** 4 * (0.993 * h + 0.007) * (V - V_NA)
    IKHT = morph.g_kht * p * (V - V_K)
    IKLT = mem.g_klt * w ** 4 * z * (V - V_K)
    # linear in activation, as in RM03 (no power given in Lehnert or Baumann)
    Ih = morph.g_h * a * (V - V_H)
    Ilk = morph.g_lk * (V - mem.vlk)  # leak reversal = resting potential by default

    total = INa + IKHT + IKLT + Ih + Ilk
    if stim_type == "EPSGbilateral":
        for k, iext in bilateral_current(t, V, s, input_node - 1, input_node2 - 1, morph.sa):
            total[k] += iext
    else:
        k, iext = external_current(t, V, stim_type, s, input_node, morph.sa)
        total[k] += iext  # the other compartments would add an exact 0.0
    d[0] = -(total + axial_current(V, morph)) / morph.cap

    act_m, act_h, act_p, act_w, act_z, act_a = active
    # sources for each time constant are in constants.py
    d[1] = (C.minf(V) - m) / C.taum(V) if act_m else 0.0
    d[2] = (C.hinf(V, mem.na_vhalf) - h) / C.tauh(V) if act_h else 0.0  # per-compartment midpoint
    d[3] = (C.pinf(V) - p) / C.taup(V) if act_p else 0.0
    d[4] = (C.winf(V) - w) / C.tauw(V) if act_w else 0.0
    d[5] = (C.zinf(V) - z) / C.tauz(V) if act_z else 0.0
    d[6] = (C.ainf(V) - a) / C.taua(V) if act_a else 0.0
    return out


def mso_axon(stim_type, start, stop, I, node, model_type, t_end, v0, input_node,
             syn: SynParams | None = None, max_step=None, stop_on_spike=None, mem=None,
             morph=None, input_node2=None):
    """Run the multi-compartment model. Returns (t, y) with y shaped (n_times, 7*n).

    Arguments shared with two_cpt (compartments are 1-indexed, times in ms):
    - stim_type: one of MSO_STIM_TYPES. What `start`, `stop` and `I` mean depends
      on it:
        step, sine: current [pA] on from start to stop (step always into the soma)
        ramp: rises from 0 at start to I [pA] at stop, then off
        ramp2: rises with slope I/(stop - start) from start, then holds at I [pA]
        EPSG: peak conductance [nS] of the unitary EPSG waveform, onset at start
          (in this model it is also cut off at stop)
        EPSGpair: two EPSGs at input_node, onsets at start and stop
        EPSGbilateral: first EPSG at input_node at start, second at input_node2 at stop
        Synaptic, SynapticPair: I is ignored; the conductance comes from `syn`,
          applied from start to stop
    - node: the compartment stop_on_spike compares with the soma (msoAxon.m
      itself ignores it). node=1 is the soma, so stop_on_spike would never fire.
    - model_type: which channels are active (ACTIVE_GATES).
    - t_end: end time [ms]; v0: initial and (unless mem says otherwise) leak
      reversal potential [mV]; input_node: where the stimulus enters.
    - syn: synaptic waveform settings (SynParams), also the EPSG time constants.
    max_step defaults to 0.1*t_end, ode15s's default MaxStep.
    stop_on_spike: if given (mV), stop as soon as compartment `node` rises that far
    above the soma; t then ends before t_end.
    mem: channel overrides from membrane() (default: msoAxon.m's own).
    morph: LUMPED or with_dendrites(); defaults to mem's morphology, else LUMPED.
    input_node2: second input site for "EPSGbilateral" (first EPSG at input_node
    at `start`, second at input_node2 at `stop`), e.g. the middle compartment of
    each dendrite.
    """
    morph = morph or (mem.morph if mem is not None else LUMPED)
    check_args(stim_type, model_type, node, input_node, min_node=1,
               stim_types=MSO_STIM_TYPES, n_max=morph.n)
    if stim_type == "EPSGbilateral" and not (input_node2 and 1 <= input_node2 <= morph.n):
        raise ValueError(f"EPSGbilateral needs input_node2 in 1..{morph.n}, got {input_node2}")
    if model_type not in ACTIVE_GATES:
        raise ValueError(f"unknown model type {model_type!r}")
    syn = syn or SynParams(t_end=t_end)

    s = stimulus(stim_type, start, stop, I, t_end, syn)
    active = tuple(bool(g) for g in ACTIVE_GATES[model_type])
    mem = mem or membrane(v0, morph=morph)
    if mem.morph_key != morph.key:
        raise ValueError("mem was built for a different morphology; pass morph= to membrane()")
    n = morph.n
    if mem.y0_gates is None:
        y0 = np.concatenate([np.full(n, float(v0))] + [np.full(n, g) for g in _GATE0])
    else:
        y0 = np.concatenate([np.full(n, float(v0)), *mem.y0_gates])
    cuts = breakpoints(stim_type, start, stop, t_end)
    quiet = None
    # unlike two_cpt, step shares the prefix: msoAxon.m switches stimuli on for
    # t > start (the EPSG waveforms are zero at onset), so the equations at
    # t == start are still stimulus-free
    if pre_stimulus_is_quiet(cuts, start, stop):
        quiet = (("mso", model_type, float(v0), mem.key, morph.key),
                 lambda t, x: _rhs(t, x, "none", s, input_node, active, mem, morph))

    spike_stop = None if stop_on_spike is None else spike_event(node - 1, stop_on_spike)
    return integrate(lambda t, x: _rhs(t, x, stim_type, s, input_node, active, mem,
                                       morph, input_node2),
                     y0, t_end, cuts, rtol=1e-8, atol=1e-8,
                     max_step=max_step or 0.1 * t_end, jac_sparsity=morph.jac,
                     quiet=quiet, stop_event=spike_stop)
