"""Geometry, channel densities and conductance fractions (port of Constants.m).

Constants.m wrote Area.mat and Fractions.mat; here they are computed on import
instead. Coupling constants have no generator in the original, so they are read
from data/coupling.json (extracted from Coupling.mat).

All arrays are 0-indexed numpy arrays of length N_CPT. Compartment k in the
MATLAB code (1-indexed) is index k-1 here.
"""

import json
from pathlib import Path

import numpy as np

N_CPT = 45  # 1 soma, 2 AIS, 21 internodes, 21 nodes
V_REST_REF = -68.0  # resting potential used to compute the fractions [mV]


def _alternating(internode, node):
    """The 42 axon compartments after soma + 2 AIS: internode, node, internode, ..."""
    return np.tile([internode, node], 21)


# --- Channel densities, Lehnert et al 2014 Table 2 [nS/um^2] ---------------------
G_NA = np.concatenate([[0.2, 4, 4], _alternating(0, 4)]).astype(float)
G_KHT = np.zeros(N_CPT)
G_KHT[0] = 0.1
G_KLT = np.concatenate([[1.55, 1.55, 1.55], _alternating(0, 1.55)])
G_H = np.concatenate([[0.02, 0.02, 0.02], np.zeros(42)])
G_LK = np.concatenate([[0.0005, 0.0005, 0.0005], _alternating(0.0002, 0.05)])


# --- Steady-state gating functions shared by both models -------------------------
# V in mV; each returns the open fraction at steady state. Sources as msoAxon.m gives
# them (Lehnert et al 2014 used the same channel models).


def minf(V):
    """Na activation: Scott et al 2010."""
    return 1.0 / (1.0 + np.exp((V + 46.0) / -11.0))


def hinf(V, vhalf=62.5):
    """Na inactivation: Scott et al 2010. vhalf is the negated midpoint (default -62.5 mV)."""
    return 1.0 / (1.0 + np.exp((V + vhalf) / 7.77))


def pinf(V):
    """KHT activation: Rothman & Manis 2003, without their slow activation gate."""
    return 1.0 / (1.0 + np.exp(-(V + 23.0) / 6.0))


def winf(V):
    """KLT activation: Mathews et al 2010."""
    return 1.0 / (1.0 + np.exp((V + 57.34) / -11.7))


def zinf(V):
    """KLT inactivation: Mathews et al 2010; never below 0.27 (non-inactivating part)."""
    return (1 - 0.27) / (1.0 + np.exp((V + 67.0) / 6.16)) + 0.27


def ainf(V):
    """h activation: Baumann et al 2013."""
    return 1.0 / (1.0 + np.exp(0.1 * (V + 80.4)))


# --- Gating time constants [ms] at 35 C, shared by both models --------------------
# The expressions are msoAxon.m's and TwoCpt.m's, term for term, so results stay
# bit-identical. The 2-CPT model holds z at rest and never uses tauz.
_P_TEMP = 3 ** ((22 - 35) / 10)  # Q10 of 3 from 22 C
_A_TEMP = 3 ** ((32 - 35) / 10)  # Q10 of 3 from 32 C


def taum(V):
    """Na activation: Scott et al 2010 in msoAxon.m (TwoCpt.m credits Rothman-Manis)."""
    return (0.141 + (-0.0826 / (1 + np.exp((-20.5 - V) / 10.8)))) / 3


def tauh(V):
    """Na inactivation: Scott et al 2010."""
    return (4 + (-3.74 / (1 + np.exp((-40.6 - V) / 5.05)))) / 3


def taup(V):
    """KHT activation: Rothman & Manis 2003, 22 C adjusted to 35 C."""
    return _P_TEMP * (100 / (4 * np.exp((V + 60) / 32) + 5 * np.exp(-(V + 60) / 22)) + 5)


def tauw(V):
    """KLT activation: Mathews et al 2010."""
    return 21.5 / (6 * np.exp((V + 60) / 7) + 24 * np.exp(-(V + 60) / 50.6)) + 0.35


def tauz(V):
    """KLT inactivation: Mathews et al 2010."""
    return 170 / (5 * np.exp((V + 60) / 10) + np.exp(-(V + 70) / 8)) + 10.7


def taua(V):
    """h activation: Baumann et al 2013, 32 C adjusted to 35 C."""
    return _A_TEMP * (79 + 417 * np.exp(-(V + 61.5) ** 2 / 800))


# --- Conductance fractions at rest -----------------------------------------------
def _fractions(v_rest=V_REST_REF):
    h = G_H * ainf(v_rest)
    na = G_NA * minf(v_rest) ** 4 * (0.993 * hinf(v_rest) + 0.007)
    kht = G_KHT * pinf(v_rest)
    klt = G_KLT * winf(v_rest) ** 4 * zinf(v_rest)
    total = klt + kht + na + h + G_LK
    safe = np.where(total != 0, total, 1.0)
    frac = lambda g: np.where(total != 0, g / safe, 0.0)
    return frac(klt), frac(h), frac(na)


KLT_FRAC, H_FRAC, NA_FRAC = _fractions()


# --- Geometry ----------------------------------------------------------------------
# membrane surface area [um^2]; conical frustum for the tapered AIS, inner area for
# myelinated internodes. The 1.66/3 term is reproduced as written in Constants.m
# (likely meant 0.66/2) because the stored Area.mat was built from it.
SA = np.concatenate([
    [8750.0,
     np.pi * (1.64 / 2 + 0.66 / 2) * np.sqrt((1.64 / 2 - 1.66 / 3) ** 2 + 10 ** 2),
     np.pi * 10 * 0.66],
    _alternating(np.pi * 100 * 0.66, np.pi * 1 * 0.66),
])
SA_CM = SA * 1e-8  # [cm^2]
AREA_RATIO = SA_CM / SA_CM[0]

# compartment length [um]; soma treated as a cylinder as long as the sphere's diameter
L = np.concatenate([[np.sqrt(8750 / np.pi), 10.0, 10.0], _alternating(100.0, 1.0)])
L_CM = L * 1e-4  # [cm]

# cross-sectional area [um^2], every compartment treated as a cylinder with the same SA
R_CYL = SA / (2 * np.pi * L)
XA = np.pi * R_CYL ** 2
XA_CM = XA * 1e-8  # [cm^2]


# --- Coupling constants (Coupling.mat) ----------------------------------------------
_coupling = json.loads((Path(__file__).parent / "data" / "coupling.json").read_text())
COUPLING1 = np.array(_coupling["coupling1"])  # forward, length 44, index node-2
COUPLING2 = np.array(_coupling["coupling2"])  # backward


# The arrays above are shared by reference (e.g. membrane() hands them to the
# model unchanged), so an in-place edit would change every later run while the
# pre-stimulus cache keys stay the same. Make them read-only.
for _a in [v for v in list(globals().values()) if isinstance(v, np.ndarray)]:
    _a.flags.writeable = False
