# Comparison-of-2-compartment-and-muti-compartment-MSO-Axon-Models

Combine_all Inputs

type
node
Input_node
stimType
tEnd
v0
GraphType
M


Combine_all.m This is a "starter" program that lets the user set model parameters and run the two-compartment model. It is called in combine all.

Graphing.m graphs results from Combine_all.m

TwoCpt.m sets various parameters (for 2-CPT model), equations used in TwoCptODE.m

TwoCptODE.m runs the 2-CPT mode... System of ordinary differential equations that define the dynamics of the two-compartment model in response to current input.

msoAxon.m sets runs the 45-CPT model

Making_Threshold_graphs.m creates spiking current thresholds... runs BinarySearch.m

BinarySearch.m runs a binary search to create the spiking current thresholds... passes values back to Making_Threshold_graphs.m

Coupling.mat contains our coupling constants

Synaptic.m creates the Synpatic input type

Spiking.m spike testing code

fitExp.m is used to calculate the inverse of our tau constant.



Programs marked with asterisk require the Zilany,Bruce, Carney "updated" Auditory Nerve model (2014).
Code for this model is available at https://www.urmc.rochester.edu/labs/carney/publications/auditory-models.aspx Using the AN model requires that the Zilany et al code by mex-compiled, and placed in a path to which MATLAB has access.

make sure klt and 

Na and KHT are for spiking

try to build