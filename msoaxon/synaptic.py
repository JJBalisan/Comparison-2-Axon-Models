"""Synaptic input conductance (port of Synaptic.m)."""

from dataclasses import dataclass

import numpy as np


@dataclass
class SynParams:
    """Stand-in for the MATLAB `Syn` struct.

    Only the fields the chosen stimulus uses need to be meaningful.
    """

    t_end: float = 20.0
    freq: float = 350.0  # input frequency [Hz]
    gE: float = 1.3e7  # epsg conductance (non-meaningful units)
    VsynE: float = 0.0  # excitatory reversal potential [mV]
    random_in: int = 1396  # RNG seed
    diff: float = 1.0  # time offset of the second input for SynapticPair [ms]
    f: float = 200.0  # sine frequency [Hz]


def synaptic(syn: SynParams):
    """Return (tSyn, gSyn): gammatone-filtered, rectified noise convolved with an EPSG.

    numpy's RNG differs from MATLAB's rng/randn, so the same seed gives a
    different (statistically equivalent) noise realisation.
    """
    dt = 0.01
    t = np.arange(0, syn.t_end + dt / 2, dt)
    nt = len(t)
    n = np.random.default_rng(syn.random_in).standard_normal(nt)

    wc = syn.freq * 2 * np.pi / 1000  # center frequency
    gw = 24.7 * (4.37 * wc / (2 * np.pi) + 1)
    gam_filter = (t / 1000) ** 4 * np.exp(-2 * np.pi * t * gw / 1000) * np.cos(t * wc)
    rect_gam = np.maximum(np.convolve(gam_filter, n)[:nt], 0)

    # 17.57 normalises so max is G, 37 is the unitary epsg amplitude (Lehnert et al)
    epsg = syn.gE * 37 * 17.57 * ((1 - np.exp(-t / 1.0)) ** 1.3 * np.exp(-t / 0.27))
    g = np.convolve(rect_gam, epsg)[:nt]
    return t, g


def interp_g(t_syn, g_syn, t):
    """interp1q, but 0 instead of NaN outside the table (TwoCptODE's guard, applied everywhere)."""
    return float(np.interp(t, t_syn, g_syn, left=0.0, right=0.0))
