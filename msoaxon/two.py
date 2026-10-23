"""Two-compartment soma + axon model (port of TwoCpt.m and TwoCptODE.m).

State vector (11): V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 — same order as the
MATLAB code, so x[:, 0] is the soma and x[:, 1] the axon compartment.
Currents are in pA, conductances in nS, capacitances in pF, time in ms.
"""

import dataclasses

import numpy as np

from msoaxon import _common, _solve, constants as C, synaptic

# Passive targets of Goldwyn, Remme & Rinzel 2019 (PLoS Comput Biol
# 15:e1006476), the paper this model's coupling-constant framework comes from.
# TwoCpt.m still carries them as comments (%8.5, %-58) beside the values it
# replaced them with.
GOLDWYN_2019 = {"r1": 8.5, "tau_est": 0.34, "v0": -58.0}

# BDF solver settings of two_cpt (TwoCpt.m's ode15s options). Unlike
# multi.mso_axon, whose max_step is a fraction of t_end, this model's is a
# fixed 0.1 ms.
RTOL = 1e-6
ATOL = 1e-6
MAX_STEP = 0.1  # [ms]


@dataclasses.dataclass(frozen=True, eq=False)
class TwoCptParams:
    """TwoCpt.m's parameter struct P, from get_params.

    eq=False keeps identity equality and hashing. Compartment 1 is the
    soma, 2 the axon node; conductances in nS, capacitances in pF.

    Attributes:
        couple12: Forward coupling constant, soma to axon.
        couple21: Backward coupling constant, axon to soma.
        gNa2: Axonal Na conductance [nS].
        gKHT: Somatic KHT conductance [nS].
        Vrest: Resting potential [mV].
        Elk: Leak reversal potential [mV].
        VK: K reversal potential of the KHT current [mV].
        EK: K reversal potential of the KLT current [mV].
        gC: Coupling conductance between the compartments [nS].
        gTot1: Total resting conductance of the soma [nS].
        gTot2: Total resting conductance of the axon [nS].
        tau1: Soma time constant [ms].
        tau2: Axon time constant [ms].
        cap1: Soma capacitance [pF].
        cap2: Axon capacitance [pF].
        gKLT1: Somatic KLT conductance [nS].
        gKLT2: Axonal KLT conductance [nS].
        Vh: h reversal potential [mV].
        gh1: Somatic h conductance [nS].
        gh2: Axonal h conductance [nS].
        glk1: Somatic leak conductance [nS].
        glk2: Axonal leak conductance [nS].
        gNa1: Somatic Na conductance [nS].
        ENa: Na reversal potential [mV].
        IappLoc: The compartment (1 or 2) the stimulus enters.
        INa0_1: Somatic Na current at rest, subtracted [pA].
        INa0_2: Axonal Na current at rest, subtracted [pA].
        IKLT0_1: Somatic KLT current at rest, subtracted [pA].
        IKLT0_2: Axonal KLT current at rest, subtracted [pA].
        Ih0_1: Somatic h current at rest, subtracted [pA].
        Ih0_2: Axonal h current at rest, subtracted [pA].
        z_rest: KLT inactivation at rest, frozen in this model.
    """

    couple12: float
    couple21: float
    gNa2: float
    gKHT: float
    Vrest: float
    Elk: float
    VK: float
    EK: float
    gC: float
    gTot1: float
    gTot2: float
    tau1: float
    tau2: float
    cap1: float
    cap2: float
    gKLT1: float
    gKLT2: float
    Vh: float
    gh1: float
    gh2: float
    glk1: float
    glk2: float
    gNa1: float
    ENa: float
    IappLoc: int
    INa0_1: float
    INa0_2: float
    IKLT0_1: float
    IKLT0_2: float
    Ih0_1: float
    Ih0_2: float
    z_rest: float


def get_params(
    v0: float,
    node: int,
    input_node: int,
    model_type: str,
    r1: float = 10.0,
    tau_est: float = 0.71,
) -> TwoCptParams:
    """Ports getParam: the two-compartment parameters.

    Args:
        v0: Resting (and leak reversal) potential [mV].
        node: The axon compartment, 1-indexed like the MATLAB code,
            whose coupling constants and channel fractions set the
            second compartment.
        input_node: Where the stimulus enters, 1-indexed.
        model_type: Which channels are active.
        r1: Soma input resistance of the passive model [MOhm].
        tau_est: The passive model's soma voltage decay time constant
            [ms]. Defaults are TwoCpt.m's values.

    Returns:
        A TwoCptParams with TwoCpt.m's field names.
    """
    couple12 = C.COUPLING1[node - 2]  # forward coupling
    couple21 = C.COUPLING2[node - 2]  # backward coupling

    if model_type in ("active-KLT", "active-full", "active-KLT+H"):
        klt1, klt2 = C.KLT_FRAC[0], C.KLT_FRAC[node - 1]
    else:
        klt1 = klt2 = 0.0

    if model_type in ("active-H", "active-full", "active-KLT+H"):
        h1, h2 = C.H_FRAC[0], C.H_FRAC[node - 1]
    else:
        h1 = h2 = 0.0

    # 'Active-sodium' (capital A) is how TwoCpt.m spells it, so 'active-sodium'
    # gets no Na in this model. Kept as-is; see PYTHON.md, "Quirks kept on
    # purpose".
    if model_type in ("Active-sodium", "active-full"):
        na1 = C.NA_FRAC[0]
        # only nodes 3 and 5 were calibrated; TwoCpt.m uses 140 for every other
        # node
        gNa2 = {3: 119.0, 5: 25.5}.get(node, 140.0)
    else:
        na1 = 0.0
        gNa2 = 0.0

    gKHT = (
        0.1 / C.SA[0] * 1000
        if model_type in ("active-KHT", "active-full")
        else 0.0
    )

    area_ratio = C.AREA_RATIO[node - 1]
    R1 = r1 * 1e-3  # input resistance to CPT1 [GOhm]
    Vrest = Elk = v0
    VK = EK = -90.0

    # passive parameters
    gC = (1 / R1) * couple21 / (1 - couple12 * couple21)
    gTot1 = gC * (1 / couple21 - 1)
    gTot2 = gC * (1 / couple12 - 1)
    # separation-of-time-scales assumption
    tau1 = tau_est * (1 - couple12 * couple21)
    tau2 = tau1 * area_ratio * (couple12 / couple21)
    cap1 = tau1 * (gTot1 + gC)
    cap2 = tau2 * (gTot2 + gC)

    Vr = Vrest
    gKLT1 = klt1 * gTot1 / (C.winf(Vr) ** 4 * C.zinf(Vr))
    gKLT2 = klt2 * gTot2 / (C.winf(Vr) ** 4 * C.zinf(Vr))
    Vh = -35.0
    gh1 = h1 * gTot1 / C.ainf(Vr)
    gh2 = h2 * gTot2 / C.ainf(Vr)
    glk1 = (1 - klt1 - h1 - na1) * gTot1
    glk2 = (1 - klt2 - h2) * gTot2
    gNa1 = na1 * gTot1 / C.na_open(C.minf(Vr), C.hinf(Vr))
    ENa = 69.0

    # resting offsets subtracted so each current is zero at Vrest
    na_rest = C.na_open(C.minf(Vr), C.hinf(Vr)) * (Vr - ENa)
    klt_rest = C.winf(Vr) ** 4 * C.zinf(Vr) * (Vr - EK)
    h_rest = C.ainf(Vr) * (Vr - Vh)
    return TwoCptParams(
        couple12=couple12,
        couple21=couple21,
        gNa2=gNa2,
        gKHT=gKHT,
        Vrest=Vrest,
        Elk=Elk,
        VK=VK,
        EK=EK,
        gC=gC,
        gTot1=gTot1,
        gTot2=gTot2,
        tau1=tau1,
        tau2=tau2,
        cap1=cap1,
        cap2=cap2,
        gKLT1=gKLT1,
        gKLT2=gKLT2,
        Vh=Vh,
        gh1=gh1,
        gh2=gh2,
        glk1=glk1,
        glk2=glk2,
        gNa1=gNa1,
        ENa=ENa,
        IappLoc=1 if input_node == 1 else 2,
        INa0_1=gNa1 * na_rest,
        INa0_2=gNa2 * na_rest,
        IKLT0_1=gKLT1 * klt_rest,
        IKLT0_2=gKLT2 * klt_rest,
        Ih0_1=gh1 * h_rest,
        Ih0_2=gh2 * h_rest,
        z_rest=C.zinf(Vr),  # z is frozen at rest in this model
    )


# One return per stimulus, in TwoCptODE.m's if/elseif order, so each branch
# can be checked against the table above multi.external_current.
def applied_current(  # noqa: PLR0911
    t: float, V1: float, stim_type: str, s: _common.Stimulus
) -> float:
    """Returns the input current at time t (the Iapp branch of TwoCptODE).

    Deliberately not shared with multi.external_current: see the table
    above that function for how the two MATLAB files differ.

    Args:
        t: Time [ms].
        V1: Soma voltage [mV].
        stim_type: The stimulus name; "none" means no stimulus.
        s: The stimulus settings (_common.stimulus).

    Returns:
        The current [pA]; positive depolarises.

    Raises:
        ValueError: For an unknown stim_type.
    """
    if stim_type == "none":
        return 0.0
    if stim_type == "step":
        return s.I if s.start <= t < s.stop else 0.0
    if stim_type == "ramp":
        return float(s.start <= t <= s.t_top) * s.I * (t - s.start) * s.slope
    if stim_type == "ramp2":
        return min(
            s.I,
            float(t >= _solve.RAMP_RISE)
            * s.I
            * (t - _solve.RAMP_RISE)
            * s.slope,
        )  # start hardcoded to 5
    if stim_type == "sine":
        wave = np.sin(s.two_pi_f * (t - s.start) / 1000)
        return float(s.start <= t <= s.stop) * s.I * wave * (wave > 0)
    if stim_type == "Synaptic":
        g = synaptic.interp_g(s.t_syn, s.g_syn, t)
        return -float(s.start <= t <= s.stop) * g * (V1 - s.VsynE)
    if stim_type == "SynapticPair":
        g1 = synaptic.interp_g(s.t_syn, s.g_syn, t)
        g2 = (
            synaptic.interp_g(s.t_syn, s.g_syn, t + s.diff)
            if t + s.diff < s.t_end
            else 0.0
        )
        return -float(s.start <= t <= s.stop) * (g1 + g2) * (V1 - s.VsynE)
    if stim_type == "EPSG":
        te = t - s.start
        return (
            s.I
            * (0 - V1)
            * float(te >= 0)
            * synaptic.epsg_unitary(te, s.epsg_tau, s.inv_peak)
        )
    if stim_type == "EPSGpair":
        te = t - s.start
        td = s.stop - s.start
        return (
            s.I
            * (0 - V1)
            * (
                float(te >= 0)
                * synaptic.epsg_unitary(te, s.epsg_tau, s.inv_peak)
                + float(te >= td)
                * synaptic.epsg_unitary(te - td, s.epsg_tau, s.inv_peak)
            )
        )
    raise ValueError(f"unknown stimType {stim_type!r}")


def _dv(
    t: float,
    x: np.ndarray,
    P: TwoCptParams,
    stim_type: str,
    s: _common.Stimulus,
) -> list[float]:
    """Computes the voltage derivatives alone: _rhs's first two entries.

    For the settled stop, which needs only their signs. A copy of _rhs's
    current terms in the same float order (see multi._dv for why not a
    shared helper); tests/test_settled.py checks the two agree bit for
    bit, so keep them in step.

    Returns:
        [dV1/dt, dV2/dt] [mV/ms].
    """
    V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 = x
    z = P.z_rest
    IKHT = P.gKHT * p * (V1 - P.VK)
    INa1 = C.na_open(m1, h1, P.gNa1) * (V1 - P.ENa) - P.INa0_1
    Ilk1 = P.glk1 * (V1 - P.Elk)
    IKLT1 = P.gKLT1 * w1**4 * z * (V1 - P.EK) - P.IKLT0_1
    Ih1 = P.gh1 * a1 * (V1 - P.Vh) - P.Ih0_1
    Ilk2 = P.glk2 * (V2 - P.Elk)
    INa2 = C.na_open(m2, h2, P.gNa2) * (V2 - P.ENa) - P.INa0_2
    IKLT2 = P.gKLT2 * w2**4 * z * (V2 - P.EK) - P.IKLT0_2
    Ih2 = P.gh2 * a2 * (V2 - P.Vh) - P.Ih0_2
    IC = P.gC * (V1 - V2)
    Iapp = applied_current(t, V1, stim_type, s)
    Iapp1, Iapp2 = (Iapp, 0.0) if P.IappLoc == 1 else (0.0, Iapp)
    return [
        (-Ilk1 - IKLT1 - IC + Iapp1 - INa1 - Ih1 - IKHT) / P.cap1,
        (-Ilk2 - IKLT2 + IC + Iapp2 - INa2 - Ih2) / P.cap2,
    ]


def _rhs(
    t: float,
    x: np.ndarray,
    P: TwoCptParams,
    stim_type: str,
    s: _common.Stimulus,
) -> list[float]:
    """Computes TwoCptODE.m's right-hand side.

    Compartment 1 is the soma, 2 the axon node; currents in pA.
    """
    V1, V2, w1, h1, w2, m1, m2, h2, p, a1, a2 = x
    z = P.z_rest  # KLT inactivation is frozen at rest in this model

    # soma: KHT exists only here. Na, KLT and h have their value at rest
    # subtracted (the *0 terms from get_params); KHT and leak don't
    IKHT = P.gKHT * p * (V1 - P.VK)
    INa1 = C.na_open(m1, h1, P.gNa1) * (V1 - P.ENa) - P.INa0_1
    Ilk1 = P.glk1 * (V1 - P.Elk)
    IKLT1 = P.gKLT1 * w1**4 * z * (V1 - P.EK) - P.IKLT0_1
    Ih1 = P.gh1 * a1 * (V1 - P.Vh) - P.Ih0_1

    # axon compartment: no KHT
    Ilk2 = P.glk2 * (V2 - P.Elk)
    INa2 = C.na_open(m2, h2, P.gNa2) * (V2 - P.ENa) - P.INa0_2
    IKLT2 = P.gKLT2 * w2**4 * z * (V2 - P.EK) - P.IKLT0_2
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


def two_cpt(
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
    stop_on_spike: float | None = None,
    r1: float = 10.0,
    tau_est: float = 0.71,
    max_step: float = MAX_STEP,
    stop_when_settled: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Runs the two-compartment model.

    Arguments are as in mso_axon (see its docstring for what start, stop
    and I mean for each stim_type), with these differences, all from
    TwoCptODE.m:

    - currents enter compartment 1 if input_node == 1, else compartment
      2 (the axon); step, ramp and sine are on for start <= t < stop (or
      <= stop) rather than start < t <= stop;
    - ramp2 always starts at t = 5, whatever `start` is;
    - EPSG is not cut off at stop.

    Args:
        stim_type: One of _common.STIM_TYPES.
        start: Stimulus start [ms].
        stop: Stimulus stop [ms].
        I: Stimulus amplitude.
        node: The axon compartment (2..45) whose coupling constants and
            channel fractions set the second compartment; only nodes 3
            and 5 have a calibrated axonal Na density (see get_params).
        model_type: Which channels are active.
        t_end: End time [ms].
        v0: Initial and leak reversal potential [mV].
        input_node: 1 for the soma; anything else means the axon.
        syn: Synaptic waveform settings (SynParams).
        stop_on_spike: If given (mV), stop as soon as the axon
            compartment rises that far above the soma; t then ends before
            t_end.
        r1: Passive calibration: soma input resistance [MOhm] (TwoCpt.m:
            10; Goldwyn et al 2019: 8.5, with v0 = -58, see GOLDWYN_2019).
        tau_est: Passive calibration: soma time constant [ms] (TwoCpt.m:
            0.71; Goldwyn et al 2019: 0.34).
        max_step: Largest solver step: a fixed 0.1 ms, not a fraction of
            t_end.
        stop_when_settled: With stop_on_spike, also stop once no spike
            can follow (see mso_axon).

    Returns:
        (t, x) with x shaped (n_times, 11); see the module docstring for
        the state layout.

    Raises:
        ValueError: For unknown names, out-of-range compartments, or a
            ramp with stop <= start.
    """
    _common.check_args(
        stim_type,
        model_type,
        node,
        input_node,
        min_node=2,  # node 1 is the soma
        start=start,
        stop=stop,
    )
    syn = syn or synaptic.SynParams(t_end=t_end)
    P = get_params(v0, node, input_node, model_type, r1=r1, tau_est=tau_est)
    s = _common.stimulus(stim_type, start, stop, I, t_end, syn)
    Vr = P.Vrest
    # unlike multi's rounded _GATE0, every gate starts at its exact steady state
    # at Vrest (z is frozen there), and get_params subtracts each current's
    # value at rest, so the model has no rest drift
    x0 = [
        Vr,
        Vr,
        C.winf(Vr),
        C.hinf(Vr),
        C.winf(Vr),
        C.minf(Vr),
        C.minf(Vr),
        C.hinf(Vr),
        C.pinf(Vr),
        C.ainf(Vr),
        C.ainf(Vr),
    ]
    cuts = _solve.breakpoints(stim_type, start, stop, t_end)
    quiet = None
    # step and the synaptic inputs are already on at t == start (TwoCptODE uses
    # t >= start), and the first segment's last step evaluates there, so sharing
    # a stimulus-free prefix would change their results slightly. They skip the
    # cache.
    if stim_type not in (
        "step",
        "Synaptic",
        "SynapticPair",
    ) and _solve.pre_stimulus_is_quiet(cuts, start, stop):
        quiet = (
            ("two", node, model_type, float(v0), input_node, r1, tau_est),
            lambda t, x: _rhs(t, x, P, "none", s),
        )
    rhs = lambda t, x: _rhs(t, x, P, stim_type, s)  # noqa: E731
    dv = lambda t, x: _dv(t, x, P, stim_type, s)  # noqa: E731
    return _solve.integrate(
        rhs,
        x0,
        t_end,
        cuts,
        rtol=RTOL,
        atol=ATOL,
        max_step=max_step,
        quiet=quiet,
        stop_event=_common._stops(
            dv, 1, stop_on_spike, stop_when_settled, stim_type, s
        ),
    )
