"""Python port of the 2-compartment vs 45-compartment MSO axon comparison.

The models live in msoaxon.multi and msoaxon.two (named after the
model="multi" / "two" switch the analysis functions use); the package
re-exports their entry points, mso_axon and two_cpt, and DriftWarning, the
warning mso_axon's Rust right-hand side (mode="fast") issues on first use.
It also re-exports process_pool and map_tasks, the worker-pool helpers the
analysis scripts share, and the initiation-site readout (msoaxon.initiation):
find_spikes, group_fractions, sustained_input and their result types.
"""

import importlib.metadata

# the package's public API: re-exported names, so `msoaxon.mso_axon` works
from msoaxon._fast import DriftWarning
from msoaxon._parallel import map_tasks, process_pool
from msoaxon.initiation import (
    Fractions,
    Spike,
    find_spikes,
    group_fractions,
    sustained_input,
)
from msoaxon.multi import mso_axon
from msoaxon.synaptic import Event, SynParams
from msoaxon.two import two_cpt

__version__ = importlib.metadata.version("msoaxon")

__all__ = [
    "DriftWarning",
    "Event",
    "Fractions",
    "Spike",
    "SynParams",
    "__version__",
    "find_spikes",
    "group_fractions",
    "map_tasks",
    "mso_axon",
    "process_pool",
    "sustained_input",
    "two_cpt",
]
