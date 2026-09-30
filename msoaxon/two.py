"""Two-compartment soma + axon model (port of TwoCpt.m and TwoCptODE.m).

State vector (11): V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 — same order as the
MATLAB code, so x[:, 0] is the soma and x[:, 1] the axon compartment.
Currents are in pA, conductances in nS, capacitances in pF, time in ms.
"""

from types import SimpleNamespace

import numpy as np

from . import constants as C
from ._common import check_args, stimulus
from ._solve import breakpoints, integrate, pre_stimulus_is_quiet, spike_event
from .synaptic import SynParams, epsg_unitary, interp_g

# Passive targets of Goldwyn, Remme & Rinzel 2019 (PLoS Comput Biol 15:e1006476), the
# paper this model's coupling-constant framework comes from. TwoCpt.m still carries
# them as comments (%8.5, %-58) beside the values it replaced them with.
GOLDWYN_2019 = dict(r1=8.5, tau_est=0.34, v0=-58.0)


def get_params(v0, node, input_node, model_type, r1=10.0, tau_est=0.71):
    """Port of getParam. `node` and `input_node` are 1-indexed like the MATLAB code.

    r1: soma input resistance of the passive model [MOhm]; tau_est: its soma
    voltage decay time constant [ms]. Defaults are TwoCpt.m's values.
    """
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
    # gets no Na in this model. Kept as-is; see PYTHON.md, "Quirks kept on purpose".
    if model_type in ("Active-sodium", "active-full"):
        na1 = C.NA_FRAC[0]
        # only nodes 3 and 5 were calibrated; TwoCpt.m uses 140 for every other node
        P.gNa2 = {3: 119.0, 5: 25.5}.get(node, 140.0)
    else:
        na1 = 0.0
        P.gNa2 = 0.0

    P.gKHT = 0.1 / C.SA[0] * 1000 if model_type in ("active-KHT", "active-full") else 0.0

    area_ratio = C.AREA_RATIO[node - 1]
    R1 = r1 * 1e-3  # input resistance to CPT1 [GOhm]
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
    """Input current [pA] at time t (the Iapp branch of TwoCptODE). "none" means no stimulus.

    Deliberately not shared with multi.external_current: see the table above that
    function for how the two MATLAB files differ.
    """
    if stim_type == "none":
        return 0.0
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
        return s.I * (0 - V1) * float(te >= 0) * epsg_unitary(te, s.epsg_tau)
    if stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        return s.I * (0 - V1) * (float(te >= 0) * epsg_unitary(te, s.epsg_tau)
                                 + float(te >= td) * epsg_unitary(te - td, s.epsg_tau))
    raise ValueError(f"unknown stimType {stim_type!r}")


def _rhs(t, x, P, stim_type, s):
    """TwoCptODE.m. Compartment 1 is the soma, 2 the axon node; currents in pA."""
    V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 = x
    z = P.z_rest  # KLT inactivation is frozen at rest in this model

    # soma: KHT exists only here. Na, KLT and h have their value at rest subtracted
    # (the *0 terms from get_params); KHT and leak don't
    IKHT = P.gKHT * p * (V1 - P.VK)
    INa1 = P.gNa1 * m1 ** 4 * (0.993 * h1 + 0.007) * (V1 - P.ENa) - P.INa0_1
    Ilk1 = P.glk1 * (V1 - P.Elk)
    IKLT1 = P.gKLT1 * w1 ** 4 * z * (V1 - P.EK) - P.IKLT0_1
    Ih1 = P.gh1 * a1 * (V1 - P.Vh) - P.Ih0_1

    # axon compartment: no KHT
    Ilk2 = P.glk2 * (V2 - P.Elk)
    INa2 = P.gNa2 * m2 ** 4 * (0.993 * h2 + 0.007) * (V2 - P.ENa) - P.INa0_2
    IKLT2 = P.gKLT2 * w2 ** 4 * z * (V2 - P.EK) - P.IKLT0_2
    Ih2 = P.gh2 * a2 * (V2 - P.Vh) - P.Ih0_2

    IC = P.gC * (V1 - V2)  # coupling current, soma to axon
    Iapp = applied_current(t, V1, stim_type, s)
    Iapp1, Iapp2 = (Iapp, 0.0) if P.IappLoc == 1 else (0.0, Iapp)

    # same order as the state vector; the gates relax to steady state with the
    # shared time constants in constants.py
    return [
        (-Ilk1 - IKLT1 - IC + Iapp1 - INa1 - Ih1 - IKHT) / P.cap1,
        (-Ilk2 - IKLT2 + IC + Iapp2 - INa2 - Ih2) / P.cap2,
        (C.winf(V1) - w1) / C.tauw(V1),
        (C.hinf(V1) - h1) / C.tauh(V1),
        (C.winf(V2) - w2) / C.tauw(V2),
        (C.minf(V1) - m1) / C.taum(V1),
        (C.minf(V2) - m2) / C.taum(V2),
        (C.hinf(V2) - h2) / C.tauh(V2),
        (C.pinf(V1) - p) / C.taup(V1),
        (C.ainf(V1) - a1) / C.taua(V1),
        (C.ainf(V2) - a2) / C.taua(V2),
    ]


def two_cpt(stim_type, start, stop, I, node, model_type, t_end, v0, input_node,
            syn: SynParams | None = None, stop_on_spike=None, r1=10.0, tau_est=0.71,
            max_step=0.1):
    """Run the two-compartment model. Returns (t, x) with x shaped (n_times, 11).

    Arguments are as in mso_axon (see its docstring for what start, stop and I
    mean for each stim_type), with these differences, all from TwoCptODE.m:
    - currents enter compartment 1 if input_node == 1, else compartment 2 (the
      axon); step, ramp and sine are on for start <= t < stop (or <= stop)
      rather than start < t <= stop;
    - ramp2 always starts at t = 5, whatever `start` is;
    - EPSG is not cut off at stop.
    - node (2..45) picks the axon compartment whose coupling constants and
      channel fractions set the second compartment; only nodes 3 and 5 have a
      calibrated axonal Na density (see get_params).
    - max_step is a fixed 0.1 ms, not a fraction of t_end.

    r1 [MOhm] and tau_est [ms] set the passive calibration (TwoCpt.m: 10 and 0.71;
    Goldwyn et al 2019: 8.5 and 0.34, with v0 = -58, see GOLDWYN_2019).

    stop_on_spike: if given (mV), stop as soon as the axon compartment rises that
    far above the soma; t then ends before t_end.
    """
    check_args(stim_type, model_type, node, input_node, min_node=2)  # node 1 is the soma
    syn = syn or SynParams(t_end=t_end)
    P = get_params(v0, node, input_node, model_type, r1=r1, tau_est=tau_est)
    s = stimulus(stim_type, start, stop, I, t_end, syn)
    Vr = P.Vrest
    x0 = [Vr, Vr, C.winf(Vr), C.hinf(Vr), C.winf(Vr), C.minf(Vr), C.minf(Vr),
          C.hinf(Vr), C.pinf(Vr), C.ainf(Vr), C.ainf(Vr)]
    cuts = breakpoints(stim_type, start, stop, t_end)
    quiet = None
    # step and the synaptic inputs are already on at t == start (TwoCptODE uses
    # t >= start), and the first segment's last step evaluates there, so sharing a
    # stimulus-free prefix would change their results slightly. They skip the cache.
    if (stim_type not in ("step", "Synaptic", "SynapticPair")
            and pre_stimulus_is_quiet(cuts, start, stop)):
        quiet = (("two", node, model_type, float(v0), input_node, r1, tau_est),
                 lambda t, x: _rhs(t, x, P, "none", s))
    spike_stop = None if stop_on_spike is None else spike_event(1, stop_on_spike)
    return integrate(lambda t, x: _rhs(t, x, P, stim_type, s), x0, t_end, cuts,
                     rtol=1e-6, atol=1e-6, max_step=max_step, quiet=quiet,
                     stop_event=spike_stop)
