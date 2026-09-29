"""Two-compartment soma + axon model (port of TwoCpt.m and TwoCptODE.m).

State vector (11): V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 — same order as the
MATLAB code, so x[:, 0] is the soma and x[:, 1] the axon compartment.
Currents are in pA, conductances in nS, capacitances in pF, time in ms.
"""

from types import SimpleNamespace

import numpy as np

from . import constants as C
from ._solve import breakpoints, epsg_unitary, integrate
from .synaptic import SynParams, interp_g, synaptic

STIM_TYPES = ("step", "ramp", "ramp2", "sine", "Synaptic", "SynapticPair", "EPSG", "EPSGpair")
MODEL_TYPES = ("passive", "active-KLT", "active-H", "active-KLT+H", "active-sodium",
               "Active-sodium", "active-KHT", "active-full")


def check_args(stim_type, model_type, node, input_node, min_node):
    """Reject inputs that numpy's negative indexing would otherwise accept silently."""
    if stim_type not in STIM_TYPES:
        raise ValueError(f"unknown stimType {stim_type!r}")
    if model_type not in MODEL_TYPES:
        raise ValueError(f"unknown model type {model_type!r}; expected one of {MODEL_TYPES}")
    if not min_node <= node <= C.N_CPT:
        raise ValueError(f"node must be {min_node}..{C.N_CPT} (1-indexed), got {node}")
    if not 1 <= input_node <= C.N_CPT:
        raise ValueError(f"input_node must be 1..{C.N_CPT} (1-indexed), got {input_node}")

# Gating kinetics (from getParam in TwoCpt.m)
_A_TEMP = 3 ** ((32 - 35) / 10)
_P_TEMP = 3 ** ((22 - 35) / 10)


def tauw(V):
    return 21.5 / (6 * np.exp((V + 60) / 7) + 24 * np.exp(-(V + 60) / 50.6)) + 0.35


def taua(V):
    return _A_TEMP * (79 + 417 * np.exp(-(V + 61.5) ** 2 / 800))


def taup(V):
    return _P_TEMP * (100 / (4 * np.exp((V + 60) / 32) + 5 * np.exp(-(V + 60) / 22)) + 5)


def taum(V):  # Rothman-Manis with 35 C adjustment
    return (0.141 + (-0.0826 / (1 + np.exp((-20.5 - V) / 10.8)))) / 3


def tauh(V):
    return (4 + (-3.74 / (1 + np.exp((-40.6 - V) / 5.05)))) / 3


def get_params(v0, node, input_node, model_type):
    """Port of getParam. `node` and `input_node` are 1-indexed like the MATLAB code."""
    P = SimpleNamespace()
    P.couple12 = C.COUPLING1[node - 2]  # forward coupling
    P.couple21 = C.COUPLING2[node - 2]  # backward coupling

    if model_type in ("active-KLT", "active-full", "active-KLT+H"):
        klt1, klt2 = C.KLT_FRAC[0], C.KLT_FRAC[node - 1]
    else:
        klt1 = klt2 = 0.0

    if model_type in ("active-H", "active-full", "active-KLT+H"):
        h1, h2 = C.H_FRAC[0], C.H_FRAC[node - 1]
    else:
        h1 = h2 = 0.0

    # 'Active-sodium' (capital A) is how TwoCpt.m spells it, so 'active-sodium'
    # gets no Na in this model. Kept as-is; see README.
    if model_type in ("Active-sodium", "active-full"):
        na1 = C.NA_FRAC[0]
        P.gNa2 = {3: 119.0, 5: 25.5}.get(node, 140.0)  # 140 is a placeholder
    else:
        na1 = 0.0
        P.gNa2 = 0.0

    P.gKHT = 0.1 / C.SA[0] * 1000 if model_type in ("active-KHT", "active-full") else 0.0

    area_ratio = C.AREA_RATIO[node - 1]
    R1 = 10 * 1e-3  # input resistance to CPT1 [GOhm]
    tau_est = 0.71  # [ms]
    P.Vrest = P.Elk = v0
    P.VK = P.EK = -90.0

    # passive parameters
    P.gC = (1 / R1) * P.couple21 / (1 - P.couple12 * P.couple21)
    P.gTot1 = P.gC * (1 / P.couple21 - 1)
    P.gTot2 = P.gC * (1 / P.couple12 - 1)
    # separation-of-time-scales assumption
    P.tau1 = tau_est * (1 - P.couple12 * P.couple21)
    P.tau2 = P.tau1 * area_ratio * (P.couple12 / P.couple21)
    P.cap1 = P.tau1 * (P.gTot1 + P.gC)
    P.cap2 = P.tau2 * (P.gTot2 + P.gC)

    Vr = P.Vrest
    P.gKLT1 = klt1 * P.gTot1 / (C.winf(Vr) ** 4 * C.zinf(Vr))
    P.gKLT2 = klt2 * P.gTot2 / (C.winf(Vr) ** 4 * C.zinf(Vr))
    P.Vh = -35.0
    P.gh1 = h1 * P.gTot1 / C.ainf(Vr)
    P.gh2 = h2 * P.gTot2 / C.ainf(Vr)
    P.glk1 = (1 - klt1 - h1 - na1) * P.gTot1
    P.glk2 = (1 - klt2 - h2) * P.gTot2
    P.gNa1 = na1 * P.gTot1 / (C.minf(Vr) ** 4 * (0.993 * C.hinf(Vr) + 0.007))
    P.ENa = 69.0
    P.IappLoc = 1 if input_node == 1 else 2

    # resting offsets subtracted so each current is zero at Vrest
    na_rest = C.minf(Vr) ** 4 * (0.993 * C.hinf(Vr) + 0.007) * (Vr - P.ENa)
    klt_rest = C.winf(Vr) ** 4 * C.zinf(Vr) * (Vr - P.EK)
    h_rest = C.ainf(Vr) * (Vr - P.Vh)
    P.INa0_1, P.INa0_2 = P.gNa1 * na_rest, P.gNa2 * na_rest
    P.IKLT0_1, P.IKLT0_2 = P.gKLT1 * klt_rest, P.gKLT2 * klt_rest
    P.Ih0_1, P.Ih0_2 = P.gh1 * h_rest, P.gh2 * h_rest
    P.z_rest = C.zinf(Vr)  # z is frozen at rest in this model
    return P


def applied_current(t, V1, stim_type, s):
    """Input current [pA] at time t (the Iapp branch of TwoCptODE)."""
    if stim_type == "step":
        return s.I if s.start <= t < s.stop else 0.0
    if stim_type == "ramp":
        slope = 1 / (s.stop - s.start)
        t_top = s.start + 1 / slope
        return float(s.start <= t <= t_top) * s.I * (t - s.start) * slope
    if stim_type == "ramp2":
        slope = 1 / (s.stop - s.start)
        return min(s.I, float(t >= 5) * s.I * (t - 5) * slope)  # start hardcoded to 5
    if stim_type == "sine":
        wave = np.sin(2 * np.pi * s.f * (t - s.start) / 1000)
        return float(s.start <= t <= s.stop) * s.I * wave * (wave > 0)
    if stim_type == "Synaptic":
        g = interp_g(s.t_syn, s.g_syn, t)
        return -float(s.start <= t <= s.stop) * g * (V1 - s.VsynE)
    if stim_type == "SynapticPair":
        g1 = interp_g(s.t_syn, s.g_syn, t)
        g2 = interp_g(s.t_syn, s.g_syn, t + s.diff) if t + s.diff < s.t_end else 0.0
        return -float(s.start <= t <= s.stop) * (g1 + g2) * (V1 - s.VsynE)
    if stim_type == "EPSG":
        te = t - s.start
        return s.I * (0 - V1) * float(te >= 0) * epsg_unitary(te)
    if stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        return s.I * (0 - V1) * (float(te >= 0) * epsg_unitary(te)
                                 + float(te >= td) * epsg_unitary(te - td))
    raise ValueError(f"unknown stimType {stim_type!r}")


def _rhs(t, x, P, stim_type, s):
    V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 = x
    z = P.z_rest

    IKHT = P.gKHT * p * (V1 - P.VK)
    INa1 = P.gNa1 * m1 ** 4 * (0.993 * h1 + 0.007) * (V1 - P.ENa) - P.INa0_1
    Ilk1 = P.glk1 * (V1 - P.Elk)
    IKLT1 = P.gKLT1 * w1 ** 4 * z * (V1 - P.EK) - P.IKLT0_1
    Ih1 = P.gh1 * a1 * (V1 - P.Vh) - P.Ih0_1

    Ilk2 = P.glk2 * (V2 - P.Elk)
    INa2 = P.gNa2 * m2 ** 4 * (0.993 * h2 + 0.007) * (V2 - P.ENa) - P.INa0_2
    IKLT2 = P.gKLT2 * w2 ** 4 * z * (V2 - P.EK) - P.IKLT0_2
    Ih2 = P.gh2 * a2 * (V2 - P.Vh) - P.Ih0_2

    IC = P.gC * (V1 - V2)
    Iapp = applied_current(t, V1, stim_type, s)
    Iapp1, Iapp2 = (Iapp, 0.0) if P.IappLoc == 1 else (0.0, Iapp)

    return [
        (-Ilk1 - IKLT1 - IC + Iapp1 - INa1 - Ih1 - IKHT) / P.cap1,
        (-Ilk2 - IKLT2 + IC + Iapp2 - INa2 - Ih2) / P.cap2,
        (C.winf(V1) - w1) / tauw(V1),
        (C.hinf(V1) - h1) / tauh(V1),
        (C.winf(V2) - w2) / tauw(V2),
        (C.minf(V1) - m1) / taum(V1),
        (C.minf(V2) - m2) / taum(V2),
        (C.hinf(V2) - h2) / tauh(V2),
        (C.pinf(V1) - p) / taup(V1),
        (C.ainf(V1) - a1) / taua(V1),
        (C.ainf(V2) - a2) / taua(V2),
    ]


def stimulus(stim_type, start, stop, I, t_end, syn):
    """Bundle stimulus settings (what TwoCpt.m stored on P)."""
    s = SimpleNamespace(start=start, stop=stop, I=I, t_end=t_end)
    if stim_type == "sine":
        s.f = syn.f
    if stim_type in ("Synaptic", "SynapticPair"):
        s.t_syn, s.g_syn = synaptic(syn)
        s.VsynE, s.diff = syn.VsynE, syn.diff
    return s


def two_cpt(stim_type, start, stop, I, node, model_type, t_end, v0, input_node,
            syn: SynParams | None = None):
    """Run the two-compartment model. Returns (t, x) with x shaped (n_times, 11)."""
    check_args(stim_type, model_type, node, input_node, min_node=2)  # node 1 is the soma
    syn = syn or SynParams(t_end=t_end)
    P = get_params(v0, node, input_node, model_type)
    s = stimulus(stim_type, start, stop, I, t_end, syn)
    Vr = P.Vrest
    x0 = [Vr, Vr, C.winf(Vr), C.hinf(Vr), C.winf(Vr), C.minf(Vr), C.minf(Vr),
          C.hinf(Vr), C.pinf(Vr), C.ainf(Vr), C.ainf(Vr)]
    return integrate(lambda t, x: _rhs(t, x, P, stim_type, s), x0, t_end,
                     breakpoints(stim_type, start, stop, t_end),
                     rtol=1e-6, atol=1e-6, max_step=0.1)
