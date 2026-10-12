"""Synaptic input conductance (port of Synaptic.m)."""

import dataclasses

import numpy as np

EPSG_TAU = (
    0.1,
    0.18,
)  # (rise, decay) [ms], the Franken et al 2015 waveform used by the MATLAB


def epsg_unitary(t, tau=EPSG_TAU):
    """Returns the unitary EPSG waveform (difference of exponentials).

    The peak is 1. The default kinetics keep the MATLAB expression,
    peak-normalised by 0.21317, so results are unchanged; other
    (rise, decay) pairs are normalised exactly.

    Args:
        t: Time since onset [ms]; zero for t <= 0.
        tau: The (rise, decay) time constants [ms].

    Returns:
        The waveform at t.
    """
    if tau == EPSG_TAU:
        return (1 / 0.21317) * (t > 0) * (np.exp(-t / 0.18) - np.exp(-t / 0.1))
    rise, decay = tau
    t_peak = np.log(decay / rise) * rise * decay / (decay - rise)
    peak = np.exp(-t_peak / decay) - np.exp(-t_peak / rise)
    return (1 / peak) * (t > 0) * (np.exp(-t / decay) - np.exp(-t / rise))


@dataclasses.dataclass
class SynParams:
    """Stand-in for the MATLAB `Syn` struct.

    Only the fields the chosen stimulus uses need to be meaningful.

    Attributes:
        t_end: End of the conductance table [ms].
        freq: Input frequency [Hz].
        gE: EPSG conductance (non-meaningful units).
        VsynE: Excitatory reversal potential [mV].
        random_in: RNG seed.
        diff: Time offset of the second input for SynapticPair [ms].
        f: Sine frequency [Hz].
        epsg_tau: EPSG/EPSGpair (rise, decay) [ms]; Myoga et al 2014
            used decay 0.3.
    """

    t_end: float = 20.0
    freq: float = 350.0  # input frequency [Hz]
    gE: float = 1.3e7  # epsg conductance (non-meaningful units)
    VsynE: float = 0.0  # excitatory reversal potential [mV]
    random_in: int = 1396  # RNG seed
    diff: float = 1.0  # time offset of the second input for SynapticPair [ms]
    f: float = 200.0  # sine frequency [Hz]
    epsg_tau: tuple = EPSG_TAU  # EPSG/EPSGpair (rise, decay) [ms]; Myoga et al 2014 used decay 0.3


def synaptic(syn: SynParams):
    """Makes the synaptic conductance: filtered noise convolved with an EPSG.

    The noise is gammatone-filtered and rectified. numpy's RNG differs
    from MATLAB's rng/randn, so the same seed gives a different
    (statistically equivalent) noise realisation.

    Args:
        syn: The settings (t_end, freq, gE, random_in).

    Returns:
        (tSyn, gSyn): the time grid [ms] and the conductance on it.
    """
    dt = 0.01
    t = np.arange(0, syn.t_end + dt / 2, dt)
    nt = len(t)
    n = np.random.default_rng(syn.random_in).standard_normal(nt)

    wc = syn.freq * 2 * np.pi / 1000  # center frequency
    gw = 24.7 * (4.37 * wc / (2 * np.pi) + 1)
    gam_filter = (
        (t / 1000) ** 4 * np.exp(-2 * np.pi * t * gw / 1000) * np.cos(t * wc)
    )
    rect_gam = np.maximum(np.convolve(gam_filter, n)[:nt], 0)

    # 17.57 normalises so max is G, 37 is the unitary epsg amplitude (Lehnert et al)
    epsg = (
        syn.gE
        * 37
        * 17.57
        * ((1 - np.exp(-t / 1.0)) ** 1.3 * np.exp(-t / 0.27))
    )
    g = np.convolve(rect_gam, epsg)[:nt]
    return t, g


def interp_g(t_syn, g_syn, t):
    """Interpolates the conductance table like interp1q, but 0 outside it.

    MATLAB's interp1q gives NaN outside the table; TwoCptODE's guard
    against that is applied everywhere.

    Args:
        t_syn: The table's times [ms].
        g_syn: The table's conductances.
        t: The time [ms].

    Returns:
        The conductance at t as a float.
    """
    return float(np.interp(t, t_syn, g_syn, left=0.0, right=0.0))
