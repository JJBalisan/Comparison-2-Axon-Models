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

## Comparison with the literature

`uv run scripts/coincidence_window.py` reproduces the numbers below; `figures/coincidence/` gets the plot and a JSON of every curve.

### Lineage

The 45-compartment model follows Lehnert et al. 2014 ([doi](https://doi.org/10.1523/JNEUROSCI.4038-13.2014)). The two-compartment model is the coupling-constant framework of Goldwyn, Remme & Rinzel 2019 ([doi](https://doi.org/10.1371/journal.pcbi.1006476)).

`TwoCpt.m` changes Goldwyn's passive targets:

| | Goldwyn 2019 | `TwoCpt.m` |
|---|---|---|
| Input resistance | 8.5 MΩ | 10 MΩ |
| Soma time constant | 0.34 ms | 0.71 ms |
| Resting potential | −58 mV | −68 mV (from Lehnert) |

The old values survive as comments beside the new ones (`%8.5`, `%-58`). `GOLDWYN_2019` in `two_cpt.py` restores them: `two_cpt(..., v0=-58, r1=8.5, tau_est=0.34)`.

### Coincidence window

Measured the way Myoga et al. 2014 did in adult gerbil MSO at 35 °C ([doi](https://doi.org/10.1038/ncomms4790)):
- two identical EPSGs, their relative timing stepped in 20 µs;
- inputs slightly above the coincident threshold, so spike probability peaks near 100%;
- result: the full width of the spike-probability curve at half maximum. They measured **221 µs**.

| Model | Model EPSG (decay 0.18 ms) | Myoga EPSG (decay 0.3 ms) |
|---|---|---|
| 45-compartment | 201 µs | 218 µs |
| 2-compartment, `TwoCpt.m` | 199 µs | 216 µs |
| 2-compartment, Goldwyn 2019 calibration | 166 µs | 181 µs |

These widths use inputs 3% above threshold. At 0.5% they shrink to 65–89 µs. The margin matters a lot, and the paper's "200 pS (~3%)" is ambiguous: 200 pS is closer to 0.5% of their 43 nS EPSGs.

How the widths are computed:
- They come from thresholds resolved to 1e-4 (`msoaxon/coincidence.py`). Probability is 50% where the threshold rises to the input level.
- Noisy trials check this: 1% amplitude jitter plus 5 µs onset jitter give 200 µs against the estimate's 199 µs. If the noise is large enough to keep peak probability near 85%, the half-maximum width comes out about 10% wider.

Findings:
- With matched methods, the original calibration reproduces the measured window closely.
- Goldwyn's faster membrane narrows it by about 17%.
- An earlier estimate here (half-way point of the threshold curve, ~340 µs) measured a different quantity.

### Other benchmarks

Values are for the Python port.

| | Model | Literature |
|---|---|---|
| EPSP half-width at soma | 523–549 µs | 0.52–0.6 ms |
| Spike initiation site | AIS | AIS (Lehnert 2014) |
| Input resistance, active, at rest | 2.5 MΩ steady, ~4.8 MΩ peak | 5 MΩ (Lehnert's model); ~7 MΩ measured |
| Somatic spike, from inflection (100 ms step) | 7.6 / 13.4 / 19.4 mV at 1.5× / 2× / 3× rheobase | 17 ± 2 mV mature, 5–15 mV near threshold, graded |
| Best frequency for half-wave sine input | 400–500 Hz | subthreshold resonance 242–300 Hz |

An earlier version of this table listed a ~40 mV somatic spike as the clearest mismatch. That figure was measured from rest during an EPSG pair at 2× threshold, so it counted the synaptic depolarisation as spike. Measured like the experiments (next section), the spike matches mature cells.

### Missing relative to current models

- dendrites
- glycinergic inhibition
- binaural or history-dependent input
- variation in cell properties along the frequency map

### Somatic spike

`uv run scripts/somatic_spike.py` reproduces this; `figures/somatic/` gets the traces and a JSON.

Scott et al. 2005 ([doi](https://doi.org/10.1523/JNEUROSCI.1016-05.2005)) evoked spikes with 100 ms somatic current steps and measured amplitude from the inflection point. Mature cells (≥P21) give **17 ± 2 mV**: 5–15 mV near threshold, growing with the stimulus. Blocking Kv1 (KLT) with dendrotoxin in P20–21 cells raised it **from 15 to 37 mV**.

`msoaxon/somatic.py` measures the model the same way. The inflection is the peak of d²V/dt² after the first 0.15 ms of charging: MSO neurons fire at step onset, while the membrane is still charging at 50–80 mV/ms, so fixed dV/dt criteria pick up the charging rather than the spike.

| Somatic spike amplitude | 1.5× rheobase | 2× | 3× |
|---|---|---|---|
| 45-compartment (`msoAxon.m`) | 7.6 mV | 13.4 mV | 19.4 mV |
| 45-compartment, KLT removed at soma + AIS, rest rebalanced (dendrotoxin) | 25.1 mV | 27.9 mV | 31.9 mV |
| 2-compartment (`TwoCpt.m`) | 6.2 mV | 15.6 mV | 17.6 mV |

- **The unmodified model already matches mature cells.** Lehnert et al. 2014 report "∼10 mV" for their model too.
- **The KLT block reproduces the developmental mechanism:** the spike roughly doubles, and rheobase falls from 3784 to 200 pA. No tuning was needed.
- **Below ~1.5× rheobase** the somatic response is a smooth hump with no distinct inflection, so no value is given.
- **Spikes driven by EPSG pairs** measure 9–11 mV at 2× threshold, against 8.5 ± 1.3 mV in vivo (van der Heijden et al. 2013). Near threshold the EPSP and spike merge and the inflection is ambiguous.

`mso_axon(..., mem=membrane(v0, ...))` exposes the knobs used for the block:
- `soma_klt_scale` and `ais_klt_scale`
- `soma_na_vhalf`: Scott et al. 2010 measured −77 mV for somatic Na inactivation; the model uses −62.5.
- `rebalance_rest`: sets each leak reversal so rest stays at `v0`. Without it, removing KLT makes the cell fire with no input.

The defaults reproduce `msoAxon.m` exactly.
