"""Python port of the 2-compartment vs 45-compartment MSO axon comparison.

The models live in msoaxon.multi and msoaxon.two (named after the
model="multi" / "two" switch the analysis functions use); the package
re-exports their entry points, mso_axon and two_cpt.
"""

import importlib.metadata

# the package's public API: re-exported names, so `msoaxon.mso_axon` works
from msoaxon.multi import mso_axon
from msoaxon.synaptic import SynParams
from msoaxon.two import two_cpt

__version__ = importlib.metadata.version("msoaxon")

__all__ = ["SynParams", "__version__", "mso_axon", "two_cpt"]
