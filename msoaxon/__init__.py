"""Python port of the two-compartment vs 45-compartment MSO axon model comparison.

The models live in msoaxon.multi and msoaxon.two (named after the model="multi" /
"two" switch the analysis functions use); the package re-exports their entry
points, mso_axon and two_cpt.
"""

from .multi import mso_axon
from .synaptic import SynParams
from .two import two_cpt

__all__ = ["mso_axon", "two_cpt", "SynParams"]
