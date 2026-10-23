# Comparison of 2-compartment and multi-compartment MSO axon models

This repo compares two models of the axon of a medial superior olive (MSO)
neuron: the 45-compartment model of Lehnert et al. 2014 (`msoAxon.m`) and the
two-compartment model of Goldwyn et al. 2019 (`TwoCpt.m`). Both are run on the
same stimuli (steps, ramps, sines, EPSGs and EPSG pairs) to compare their
spiking thresholds. The comparison began in 2019 as MATLAB code. In 2026 it was
ported to Python, the `msoaxon` package, which is now the reference: its
results are pinned by the golden traces in `tests/`. The `.m` files stay as the
history. [PYTHON.md](PYTHON.md) has the details: validation, fast mode,
calibration presets, the somatic spike, the coincidence window and benchmarks.

## Quick start

```bash
uv sync                  # numpy/scipy only
uv sync --extra fast     # also builds the Rust right-hand side (needs cargo)
uv run pytest            # the tests, in parallel
```

Fast mode makes the 45-compartment runs about 2.7 times faster. Without cargo,
install a prebuilt wheel from the assets of a GitHub release (steps in
PYTHON.md, "Fast mode").

| Script | Produces |
|---|---|
| `scripts/combine_all.py` | One stimulus on both models, plotted together |
| `scripts/making_threshold_graphs.py` | Spiking-current thresholds vs stimulus |
| `scripts/coincidence_window.py` | Coincidence-detection windows vs Myoga et al. 2014 |
| `scripts/dendrites.py` | The dendritic 45-compartment variant vs the lumped model |
| `scripts/somatic_spike.py` | Somatic spike size vs Scott et al. 2005 |
| `scripts/benchmark_search.py` | Simulation counts and timings of the threshold searches |

## The MATLAB files

They are unchanged from 2019 and no longer have to match the Python.

| File | What it did | Python |
|---|---|---|
| `Combine_all.m` | Sets the parameters and runs both models on one stimulus | `scripts/combine_all.py` |
| `Graphing.m` | Plots the results of `Combine_all.m` | `msoaxon/plotting.py` |
| `TwoCpt.m`, `TwoCptODE.m` | Two-compartment model: parameters, and its ODEs | `msoaxon/two.py` |
| `msoAxon.m` | The 45-compartment model | `msoaxon/multi.py` |
| `Constants.m`, `Area.mat`, `Fractions.mat`, `Coupling.mat` | Areas, fractions and coupling constants | `msoaxon/constants.py`, `msoaxon/data/` |
| `Synaptic.m` | The synaptic input type | `msoaxon/synaptic.py` |
| `Spiking.m` | Spike counting | `msoaxon/spiking.py` |
| `BinarySearch.m` | Binary search for the spiking current threshold | `msoaxon/threshold.py` |
| `Making_Threshold_graphs.m` | Threshold graphs, via `BinarySearch.m` | `scripts/making_threshold_graphs.py` |
| `fitExp.m` | Fits the inverse of the time constant | Not ported (see PYTHON.md) |

`msoAxon.m~` is an earlier draft of `msoAxon.m`, kept as part of that record.
`EPSGpair_Thresholds*.jpg` are the 2019 threshold graphs that the Python is
compared with.

## Licence

GNU General Public License v3.0; see [LICENSE](LICENSE).
