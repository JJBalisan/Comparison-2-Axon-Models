"""The stimulus current of multi._rhs: which compartment it drives."""

import numpy as np
import pytest

from msoaxon import _common, multi, synaptic


@pytest.mark.parametrize("stim", ["Synaptic", "SynapticPair"])
def test_synaptic_input_drives_its_own_compartment(stim):
    # the conductance into node 3 (compartment 2) sees compartment 2's
    # voltage, not the soma's: dV[2] does not change with V[0], with which
    # compartment 2 shares no edge in the chain; its neighbour, compartment
    # 1, does (through the axial current). Until 1.4.0 the driving force
    # was the soma's, V[0] - VsynE, whatever input_node said.
    syn = synaptic.SynParams(t_end=20.0, diff=1.0)
    s = _common.stimulus(stim, 5.0, 15.0, 0.0, 20.0, syn)
    mem = multi.membrane(-68.0)
    active = tuple(bool(g) for g in multi.ACTIVE_GATES["active-full"])
    n = multi.LUMPED.n
    t = 8.0
    assert synaptic.interp_g(s.t_syn, s.g_syn, t) > 0  # the input is on
    x = np.concatenate([np.full(n, -68.0), np.full(6 * n, 0.3)])
    x_soma = x.copy()
    x_soma[0] = -40.0
    d = multi._rhs(t, x, stim, s, 3, active, mem)
    d_soma = multi._rhs(t, x_soma, stim, s, 3, active, mem)
    assert d[2] == d_soma[2]
    assert d[1] != d_soma[1]
    # and the input does drive compartment 2: its own voltage matters
    x_node = x.copy()
    x_node[2] = -40.0
    d_node = multi._rhs(t, x_node, stim, s, 3, active, mem)
    assert d_node[2] != d[2]
