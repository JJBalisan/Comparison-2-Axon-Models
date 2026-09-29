"""45-compartment MSO soma + axon model (port of msoAxon.m).

Equations as in Lehnert et al 2014 unless noted. All quantities are densities
(pA/um^2, nS/um^2).

State layout matches MATLAB's column-major reshape(x, 45, 7): the output y has
315 columns, [V(45), m(45), h(45), p(45), w(45), z(45), a(45)]. So y[:, k-1] is
the voltage of MATLAB compartment k, and y[:, 225] is z1.
"""

import numpy as np
from scipy.sparse import diags, eye, bmat

from . import constants as C
from ._solve import breakpoints, epsg_unitary, integrate, pre_stimulus_is_quiet, spike_event
from .synaptic import SynParams, interp_g
from .two_cpt import check_args, stimulus

N = C.N_CPT

# specific capacitance, uF/cm^2 converted to ms*nS/um^2
CAP = np.concatenate([[0.8, 0.8, 0.8], np.tile([0.01, 0.8], 21)]) * 1e6 * 1e-8
V_NA, V_K, V_H = 69.0, -90.0, -35.0
R_AXIAL = 100.0  # specific axial resistivity [Ohm cm]

# conductance between neighbouring compartments i and i+1 (the (2/R)/(...) term)
_G_AX = (2 / R_AXIAL) / (C.L_CM[:-1] / C.XA_CM[:-1] + C.L_CM[1:] / C.XA_CM[1:])

_P_TEMP = 3 ** ((22 - 35) / 10)
_A_TEMP = 3 ** ((32 - 35) / 10)

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


def axial_current(V):
    """Axial current density [pA/um^2], sign convention as in msoAxon.m."""
    flow = _G_AX * (V[:-1] - V[1:])  # from i to i+1 [mA]
    I = np.zeros_like(V)
    I[:-1] -= flow
    I[1:] += flow
    return -I * 1e9 / C.SA


def external_current(t, V, stim_type, s, input_node):
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
            return k, -I0 * 1e3 / C.SA[0]
    elif stim_type == "ramp2":
        slope = 1 / (s.stop - s.start)
        I0 = min(float(t >= s.start) * s.I * (t - s.start) * slope / 1000, s.I / 1000)
        return k, -I0 * 1e3 / C.SA[0]
    elif stim_type == "step":
        if s.start < t <= s.stop:
            return 0, -s.I / C.SA[0]  # always the soma, regardless of input_node
    elif stim_type == "sine":
        wave = np.sin(2 * np.pi * s.f * (t - s.start) / 1000)
        if s.start < t <= s.stop:
            return k, -(s.I * wave * (wave > 0)) / C.SA[0]
    elif stim_type == "Synaptic":
        if s.start < t <= s.stop:
            g = interp_g(s.t_syn, s.g_syn, t)
            return k, g * (V[0] - s.VsynE) / C.SA[0]
    elif stim_type == "SynapticPair":
        if s.start < t <= s.stop:
            g = interp_g(s.t_syn, s.g_syn, t) + interp_g(s.t_syn, s.g_syn, t + s.diff)
            return k, g * (V[0] - s.VsynE) / C.SA[0]
    elif stim_type == "EPSG":
        te = t - s.start
        if s.start < t <= s.stop:
            return k, s.I * (0 - V[k]) * float(te >= 0) * epsg_unitary(te, s.epsg_tau) / -C.SA[k]
    elif stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        wave = (float(te >= 0) * epsg_unitary(te, s.epsg_tau)
                + float(te >= td) * epsg_unitary(te - td, s.epsg_tau))
        return k, s.I * (0 - V[k]) * wave / -C.SA[k]
    return 0, 0.0


def _rhs(t, x, v0, stim_type, s, input_node, active):
    """Right-hand side. Arithmetic is kept expression-for-expression identical to
    msoAxon.m's order so results are bit-for-bit stable; the speed comes from
    writing into one output array and skipping gates that are switched off."""
    V, m, h, p, w, z, a = x.reshape(7, N)
    out = np.empty_like(x)
    d = out.reshape(7, N)

    INa = C.G_NA * m ** 4 * (0.993 * h + 0.007) * (V - V_NA)
    IKHT = C.G_KHT * p * (V - V_K)
    IKLT = C.G_KLT * w ** 4 * z * (V - V_K)
    # linear in activation, as in RM03 (no power given in Lehnert or Baumann)
    Ih = C.G_H * a * (V - V_H)
    Ilk = C.G_LK * (V - v0)  # leak reversal = resting potential

    total = INa + IKHT + IKLT + Ih + Ilk
    k, iext = external_current(t, V, stim_type, s, input_node)
    total[k] += iext  # the other compartments would add an exact 0.0
    d[0] = -(total + axial_current(V)) / CAP

    act_m, act_h, act_p, act_w, act_z, act_a = active
    # Na: Scott et al 2010, 35 C
    d[1] = (C.minf(V) - m) / ((0.141 + (-0.0826 / (1 + np.exp((-20.5 - V) / 10.8)))) / 3) if act_m else 0.0
    d[2] = (C.hinf(V) - h) / ((4 + (-3.74 / (1 + np.exp((-40.6 - V) / 5.05)))) / 3) if act_h else 0.0
    # KHT: Rothman Manis 2003, 22 C adjusted to 35 C with Q10 of 3
    d[3] = (C.pinf(V) - p) / (_P_TEMP * (100 / (4 * np.exp((V + 60) / 32)
                                                 + 5 * np.exp(-(V + 60) / 22)) + 5)) if act_p else 0.0
    # KLT: Mathews et al 2010, 35 C
    d[4] = (C.winf(V) - w) / (21.5 / (6 * np.exp((V + 60) / 7)
                                      + 24 * np.exp(-(V + 60) / 50.6)) + 0.35) if act_w else 0.0
    d[5] = (C.zinf(V) - z) / (170 / (5 * np.exp((V + 60) / 10)
                                     + np.exp(-(V + 70) / 8)) + 10.7) if act_z else 0.0
    # h: Baumann et al 2013, 32 C adjusted to 35 C
    d[6] = (C.ainf(V) - a) / (_A_TEMP * (79 + 417 * np.exp(-(V + 61.5) ** 2 / 800))) if act_a else 0.0
    return out


def _jac_sparsity():
    """dV couples to neighbouring V and its own gates; each gate to its own V."""
    I = eye(N)
    tri = diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(N, N))
    blocks = [[tri] + [I] * 6] + [[I] + [I if j == i else None for j in range(6)]
                                  for i in range(6)]
    return bmat(blocks).tocsc()


_JAC_SPARSITY = _jac_sparsity()


def mso_axon(stim_type, start, stop, I, node, model_type, t_end, v0, input_node,
             syn: SynParams | None = None, max_step=None, stop_on_spike=None):
    """Run the 45-compartment model. Returns (t, y) with y shaped (n_times, 315).

    `node` is accepted for call parity with two_cpt; msoAxon.m ignores it too.
    max_step defaults to 0.1*t_end, ode15s's default MaxStep.
    stop_on_spike: if given (mV), stop as soon as compartment `node` rises that far
    above the soma; t then ends before t_end.
    """
    check_args(stim_type, model_type, node, input_node, min_node=1)
    if model_type not in ACTIVE_GATES:
        raise ValueError(f"unknown model type {model_type!r}")
    syn = syn or SynParams(t_end=t_end)

    s = stimulus(stim_type, start, stop, I, t_end, syn)
    active = tuple(bool(g) for g in ACTIVE_GATES[model_type])
    y0 = np.concatenate([np.full(N, float(v0))] + [np.full(N, g) for g in _GATE0])
    cuts = breakpoints(stim_type, start, stop, t_end)
    quiet = None
    if pre_stimulus_is_quiet(cuts, start, stop):
        quiet = (("mso", model_type, float(v0)),
                 lambda t, x: _rhs(t, x, v0, "none", s, input_node, active))

    spike_stop = None if stop_on_spike is None else spike_event(node - 1, stop_on_spike)
    return integrate(lambda t, x: _rhs(t, x, v0, stim_type, s, input_node, active),
                     y0, t_end, cuts, rtol=1e-8, atol=1e-8,
                     max_step=max_step or 0.1 * t_end, jac_sparsity=_JAC_SPARSITY,
                     quiet=quiet, stop_event=spike_stop)
