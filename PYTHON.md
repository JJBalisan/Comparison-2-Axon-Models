# Python port

A numpy/scipy port of the MATLAB model comparison. The `.m` files are untouched and
remain the reference.

```bash
uv sync
uv run pytest                                              # constants + sanity checks
uv run scripts/combine_all.py --stim EPSGpair --node 3     # Combine_all.m
uv run scripts/making_threshold_graphs.py --node 3         # Making_Threshold_graphs.m
uv run scripts/making_threshold_graphs.py --only EPSGpair --jpg-grid   # the repo's jpgs
```

| MATLAB | Python |
|---|---|
| `Constants.m` → `Area.mat`, `Fractions.mat` | `msoaxon/constants.py` (computed on import) |
| `Coupling.mat` | `msoaxon/data/coupling.json` |
| `TwoCpt.m`, `TwoCptODE.m` | `msoaxon/two_cpt.py` → `two_cpt(...)` |
| `msoAxon.m` | `msoaxon/mso_axon.py` → `mso_axon(...)` |
| `Synaptic.m` | `msoaxon/synaptic.py` (`SynParams` replaces the `Syn` struct) |
| `Spiking.m` | `msoaxon/spiking.py` |
| `BinarySearch.m` | `msoaxon/threshold.py` → `binary_search(...)` |
| `Graphing.m` | `msoaxon/plotting.py` → `graphing(graph, ...)` |
| `Combine_all.m`, `Making_Threshold_graphs.m` | `scripts/` |

`two_cpt` and `mso_axon` take the MATLAB argument order
`(stimType, start, stop, I, node, type, tEnd, v0, inputNode, syn)`. `node` and
`inputNode` stay **1-indexed** like the MATLAB code. Output arrays are 0-indexed:
`y[:, node-1]` is compartment `node` of the 45-CPT model, and its 315 state columns
keep MATLAB's column-major layout `[V(45), m, h, p, w, z, a]`, so `z1` is `y[:, 225]`.

## Validation

There is no MATLAB here, so nothing was compared run-for-run. What was checked:

- **Constants:** `Area.mat`, `Fractions.mat` and `Coupling.mat` match the Python values to 1e-12 (`tests/`).
- **Rest:** with no input the 2-CPT model stays at `v0` to 1e-6 mV. The 45-CPT model settles about 0.9 mV away from it. The same happens in MATLAB, since `msoAxon.m` starts from rounded gating values and doesn't subtract resting currents.
- **EPSGpair thresholds vs `EPSGpair_Thresholds*.jpg`:** every plotted value is a multiple of q = 150/128, i.e. raw output of the halving search. Re-running that grid (`--jpg-grid`, delays 0–0.6 ms, nodes 3 and 5, both models) matches 20 of 28 values exactly. The other 8 are one possible output away (2q ≈ 2.3 units, 2.5–5%). The curve shapes, the steep rise between 0.2 and 0.5 ms, and the plateau positions all agree. The misses aren't random: node-5 two-compartment sits one step low at four consecutive delays, and the node-3 multi-compartment plateau is 94.9 vs 97.3.
  - Tightening the tolerances 1000× (1e-9, `max_step` 0.02) doesn't move any Python value, so the Python side is converged. The gap is either ode15s error in the original runs or a code change after the jpgs were made (they came from an older 11-point `BinarySearch.m`). Telling those apart needs MATLAB.
- **Not validated numerically:** `step`, `ramp`, `ramp2`, `sine`, `EPSG`, and the synaptic inputs. Every stimulus, graph type and threshold sweep runs, and the traces look physiologically sensible, but there's no reference output to compare against.

## Deliberate differences

- **Solver:** `ode15s` → scipy `solve_ivp(method="BDF")` with the same tolerances (2-CPT 1e-6, `max_step` 0.1; 45-CPT 1e-8, `max_step` 0.1·tEnd, ode15s's default). Integration restarts at stimulus on/off times so the solver can't step over a narrow EPSG. The 45-CPT model gets a sparse Jacobian pattern.
- **Threshold search stops each run at its first spike** (`stop_on_spike`), since it only needs yes or no. The stop is found on the continuous solution, not the saved solver steps, so in principle it could count a crossing that falls between two steps. Across every production sweep (184 thresholds) the results are identical to running each simulation to `tEnd`. Sweep points also run in parallel (`workers`, default all CPUs).
- **Synaptic noise:** numpy's RNG can't reproduce MATLAB's `rng(seed); randn`, so `Synaptic`/`SynapticPair` runs are statistically equivalent, not identical.
- **`interp1q` out of range** returns 0, not NaN. `TwoCptODE.m` already guarded this for `SynapticPair`. `msoAxon.m` didn't, which is likely the "Problem including Synaptic pair" comment in `Combine_all.m`.
- **`BinarySearch.m`:** the two-compartment loop tested `location1` in its `while` condition. It now tests its own `location2`.
- **`Combine_all.m`:** the 2-CPT spike count read column 1 (soma minus soma, always 0); it now reads the axon column. The second `TwoCpt` run used identical arguments; it now uses node 5 so the "Node 5" panels mean what they say. The default `ramp` stop was equal to start, a zero-length ramp; it's now 5.5.
- **`Graphing.m`:** the `Input` panel reuses the model's own stimulus function instead of a hand-copied version that indexed 2-CPT voltages with multi-compartment time steps. The `z1`/`a` line fits (broken: `graph.graph.tEnd`, length mismatch) were not ported.
- **EPSGpair threshold x-axis** shows the delay between the two EPSGs, like the jpgs do, rather than the second EPSG's absolute time (`5 + (i-1)/25`) that `Making_Threshold_graphs.m` plots.
- **`fitExp.m`** was not ported. It contains an unresolved merge conflict and fits 301 samples against a 49-point x-axis.

## Quirks kept on purpose

These are in the MATLAB and are reproduced unchanged, because changing them would change results:

- `TwoCpt.m` enables sodium for `'Active-sodium'` (capital A), so `active-sodium` means *no* Na in the 2-CPT model but Na-only in the 45-CPT model.
- The 2-CPT model holds KLT inactivation `z` fixed at its resting value.
- In `msoAxon.m`, `step` always injects into the soma whatever `inputNode` is; `ramp` only applies for `t > 5` (hardcoded); `EPSG`/`EPSGpair` scale by the input compartment's area while the other stimuli use the soma's.
- `Spiking.m` doesn't clear `reset` between columns.
- `SA(2)` uses `1.66/3` where `0.66/2` was probably meant. The stored `Area.mat` was built with it.
- `BinarySearch.m` reports the ceiling (`max`) when nothing fires, so a returned 15000 means "no threshold found", not "threshold ≈ 15000". It also assumes spiking gets easier as input grows. That fails for `EPSG` with factor 30: a very large conductance holds the soma near 0 mV, so the axon never gets 30 mV above it. The `ramp` sweep's first point (stop = 5.1) returns the ceiling in both models for the same reason.
