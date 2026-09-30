# Python port

A numpy/scipy port of the MATLAB model comparison. The `.m` files are untouched and
remain the reference.

```bash
uv sync
uv run pytest                                              # the tests, in parallel (-n 0: serially)
uv run pytest --cov                                        # plus line + branch coverage (as CI)
uv run scripts/combine_all.py --stim EPSGpair --node 3     # Combine_all.m
uv run scripts/making_threshold_graphs.py --node 3         # Making_Threshold_graphs.m
uv run scripts/making_threshold_graphs.py --only EPSGpair --jpg-grid   # the repo's jpgs
```

`dendrites.py`, `coincidence_window.py` and `somatic_spike.py` take `--quick` (coarse
delay grids, 1% threshold tolerance, fewer multiples), and `making_threshold_graphs.py`
takes `--points N`. The same code runs in seconds instead of minutes; the numbers are
for smoke tests, not results. `tests/test_scripts.py` runs every script this way; those
tests are marked `slow` (`uv run pytest -m "not slow"` skips them), and CI runs them on
Python 3.12 only.

| MATLAB | Python |
|---|---|
| `Constants.m` → `Area.mat`, `Fractions.mat` | `msoaxon/constants.py` (computed on import) |
| `Coupling.mat` | `msoaxon/data/coupling.json` |
| `TwoCpt.m`, `TwoCptODE.m` | `msoaxon/two.py` → `two_cpt(...)` |
| `msoAxon.m` | `msoaxon/multi.py` → `mso_axon(...)` |
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

Code that isn't a port of a MATLAB file:

| Module | What it holds |
|---|---|
| `msoaxon/_common.py` | stimulus bundle and argument checks shared by both models |
| `msoaxon/_solve.py` | the ode15s stand-in: segmented BDF, pre-stimulus cache, spike stop |
| `msoaxon/_dispatch.py` | `run_model` / `spikes`: one entry point for `model="multi"` or `"two"` |
| `msoaxon/_parallel.py` | `map_tasks`; pass one `executor=` to reuse a process pool across calls |
| `msoaxon/_bisect.py` | the bracketing bisection behind `coincidence.threshold` and `somatic.rheobase` |
| `msoaxon/coincidence.py` | Myoga-style coincidence windows; `stim=`, `input_node=`, `input_node2=` pick the input sites, `model_kw` passes model options such as `morph` |
| `msoaxon/somatic.py` | rheobase and the Scott et al. 2005 somatic spike amplitude |
| `msoaxon/measure.py` | input resistance and time constant from a current step; soma trace on a grid |

## Validation

There is no MATLAB here, so nothing was compared run-for-run. What was checked:

- **Constants:** `Area.mat`, `Fractions.mat` and `Coupling.mat` match the Python values to 1e-12 (`tests/`).
- **Rest:** with no input the 2-CPT model stays at `v0` to 1e-6 mV. The 45-CPT model settles about 0.9 mV away from it. The same happens in MATLAB, since `msoAxon.m` starts from rounded gating values and doesn't subtract resting currents.
- **EPSGpair thresholds vs `EPSGpair_Thresholds*.jpg`:** every plotted value is a multiple of q = 150/128, i.e. raw output of the halving search. Re-running that grid (`--jpg-grid`, delays 0–0.6 ms, nodes 3 and 5, both models) matches 20 of 28 values exactly. The other 8 are one possible output away (2q ≈ 2.3 units, 2.5–5%). The curve shapes, the steep rise between 0.2 and 0.5 ms, and the plateau positions all agree. The misses aren't random: node-5 two-compartment sits one step low at four consecutive delays, and the node-3 multi-compartment plateau is 94.9 vs 97.3.
  - Tightening the tolerances 1000× (1e-9, `max_step` 0.02) doesn't move any Python value, so the Python side is converged. The gap is either ode15s error in the original runs or a code change after the jpgs were made (they came from an older 11-point `BinarySearch.m`). Telling those apart needs MATLAB.
- **Not validated numerically:** `step`, `ramp`, `ramp2`, `sine`, `EPSG`, and the synaptic inputs. Every stimulus, graph type and threshold sweep runs, and the traces look physiologically sensible, but there's no reference output to compare against.

## Deliberate differences

- **Solver:** `ode15s` → scipy `solve_ivp(method="BDF")` with the same tolerances (2-CPT 1e-6, `max_step` 0.1; 45-CPT 1e-8, `max_step` 0.1·tEnd, ode15s's default). Integration restarts at stimulus on/off times so the solver can't step over a narrow EPSG. The 45-CPT model gets a sparse Jacobian pattern.
- **Threshold search stops each run at its first spike** (`stop_on_spike`), since it only needs yes or no. The stop is found on the continuous solution, not the saved solver steps, so in principle it could count a crossing that falls between two steps. Across every production sweep (184 thresholds) the results are identical to running each simulation to `tEnd`. The search used to also count spikes in the saved steps, as `BinarySearch.m` does; that never decided a result, since a sampled crossing is one the stop already caught, so it was dropped. Sweep points also run in parallel (`workers`, default all CPUs).
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
- `TwoCptODE.m` starts `ramp2` at t = 5 whatever `start` is; `msoAxon.m` starts it at `start`.
- `TwoCpt.m` has calibrated axonal Na (`gNa2`) only for nodes 3 (119) and 5 (25.5). Every other node gets the placeholder 140, so 2-CPT runs at those nodes use an uncalibrated value.
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

The old values survive as comments beside the new ones (`%8.5`, `%-58`). `GOLDWYN_2019` in `msoaxon/two.py` restores them: `two_cpt(..., v0=-58, r1=8.5, tau_est=0.34)`.

### Coincidence window

Measured the way Myoga et al. 2014 did in adult gerbil MSO at 35 °C ([doi](https://doi.org/10.1038/ncomms4790)):
- two identical EPSGs, their relative timing stepped in 20 µs;
- inputs slightly above the coincident threshold, so spike probability peaks near 100%;
- result: the full width of the spike-probability curve at half maximum. They measured **221 µs**.

| Model | Model EPSG (decay 0.18 ms) | Myoga EPSG (decay 0.3 ms) |
|---|---|---|
| 45-compartment | 201 µs | 219 µs |
| 2-compartment, `TwoCpt.m` | 199 µs | 217 µs |
| 2-compartment, Goldwyn 2019 calibration | 167 µs | 182 µs |

These widths use inputs 3% above threshold. At 0.5% they shrink to 70–91 µs. The margin matters a lot, and the paper's "200 pS (~3%)" is ambiguous: 200 pS is closer to 0.5% of their 43 nS EPSGs.

How the widths are computed:
- `coincidence.window()` does what the experiment does: it fixes the input at (1 + margin) × the coincident threshold (resolved to 1e-4) and searches over delay, to 0.1 µs, for where that input stops spiking. Probability is 50% there. About 20 simulations per window, within 0.5 µs of the same search at 1e-7.
- Both searches predict rather than halve. At the 10 mV criterion the peak of axon − soma rises smoothly through 10 mV (7.7, 9.7, 10.0, 10.3, 17.8 mV at −1%, −0.1%, 0, +0.1%, +1% of threshold), so each run that doesn't spike says how close it came. That takes 8–12 simulations per threshold instead of 15.
- Earlier versions read the crossing off a threshold curve on a 20 µs grid. Linear interpolation across a curve that bends upward read 0.1–4 µs low, most at the 0.5% margin.
- Noisy trials check this: 1% amplitude jitter plus 5 µs onset jitter give 200 µs against `window()`'s 199 µs. If the noise is large enough to keep peak probability near 85%, the half-maximum width comes out about 10% wider.

The spike criterion (axon − soma, 10 mV here, as in the MATLAB EPSG-pair sweeps) sits below where the spike takes off, 1–3% above that threshold. It barely affects the windows. Set as a fraction of each model's own spike height instead (about 38.5 mV of axon − soma in the 45-compartment and `TwoCpt.m` models, 20 mV with Goldwyn's calibration), anywhere from 25% to 75% of it, the 3% windows change by at most 1.7% (201 → 198 µs for the 45-compartment model) and 3.5% for Goldwyn. The coincident thresholds move more: up to 3%, and 13% for Goldwyn at 25%, where the criterion is below its takeoff. A fixed 20 or 30 mV criterion is out of reach for Goldwyn's calibration altogether.

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

## Dendrites

`mso_axon(..., morph=with_dendrites())` adds two dendrites, following the dendritic variant in Lehnert et al. 2014 (their Fig. 8). The lumped model's 8750 µm² "soma" stands for soma plus dendrites: in the paper's words it "combines the somatic and dendritic membrane surface", sized for 70 pF. So the dendrites take their membrane from it rather than being added on top:
- two unbranched dendrites, each 200 µm × 5 µm, five compartments each;
- the soma reduced to 2467 µm², so total membrane stays 8750 µm²;
- somatic Na density scaled up so total Na is unchanged;
- no Na or KHT in the dendrites;
- KLT and h decaying along the dendrites with a 74 µm length constant (Mathews et al. 2010, [doi](https://doi.org/10.1038/nn.2530)).

How to use it:
- Dendrite compartments are appended after the axon, so compartments 1–45 keep their meaning. Numbers 46–50 are the lateral (ipsilateral) dendrite, 51–55 the medial (contralateral), proximal to distal.
- New stimulus `EPSGbilateral`: one EPSG at `input_node` at `start`, the other at `input_node2` at `stop`.
- `membrane(..., morph=D, dendrite_klt_scale=...)` scales dendritic KLT.
- `uv run scripts/dendrites.py` runs every check below.

Two things the paper doesn't state:
- **Axial resistivity.** The model's 100 Ω·cm is kept by default; Mathews used 200. `with_dendrites(dendrite_ra=200)` (script: `--dendrite-ra 200`) changes it for the dendrites only, and the results are in the section after this one.
- **Whether total KLT and h were conserved.** `conserve_totals=True` (default) keeps soma + dendrite totals equal to the lumped model's. `False` starts the gradient at the soma's Table 2 density, halving total KLT.

| | Lumped | Dendrites, totals conserved (default) | Dendrites, KLT from soma density | Literature |
|---|---|---|---|---|
| Soma input resistance, steady / peak | 2.4 / 4.7 MΩ | 2.6 / 4.8 MΩ | 4.9 / 7.7 MΩ | 5 MΩ (Lehnert's model); ~7 MΩ measured |
| Somatic EPSP half-width, unitary EPSG at soma | 547 µs | 527 µs | 691 µs | 0.52–0.6 ms |
| Same EPSG at mid-dendrite (100 µm): size, rise, half-width | — | 4.1 mV, 221 µs, 536 µs | 5.3 mV, 265 µs, 697 µs | |
| Mid-dendrite half-width with dendritic KLT removed | — | 689 µs | 910 µs | wider without it (Mathews 2010) |
| Threshold at 0 delay: bilateral vs unilateral | — | 58.5 vs 82.2 | 37.1 vs 44.0 | bilateral lower (Scott 2010) |
| Bilateral / summed unilateral EPSP | — | 0.98 | 1.00 | linear in vivo |
| Coincidence window at 3%, model EPSG / Myoga EPSG | 201 / 219 µs (soma input) | 165 / 181 µs | 223 / 249 µs | 221 µs (Myoga 2014) |
| EPSG-pair threshold curve at the soma vs lumped (max difference) | — | 12% | 50% | "almost identical" (Lehnert) |
| Somatic spike at 1.5 / 2 / 3× rheobase | 7.6 / 13.4 / 19.4 mV | 16.3 / 19.5 / 26.6 mV | 21.4 / 26.1 / 29.6 mV | 17 ± 2 mV mature |

Conserving the totals is the default because it reproduces what Lehnert et al. reported for their own dendritic variant: tuning almost identical to the lumped model. It also fixes the EPSP width and keeps the somatic spike near mature size.

Neither variant matches everything:
- **Conserving totals** narrows the coincidence window below Myoga's value.
- **The soma-density gradient** matches the input resistance and the window, but leaves the spike too large.

Both reproduce the qualitative dendritic results: attenuation along the cable, EPSP sharpening by dendritic KLT, and the bilateral advantage.

### Dendritic axial resistivity 200 Ω·cm

`with_dendrites(dendrite_ra=200)` uses Mathews et al.'s 200 Ω·cm in the dendrites. The soma and axon keep 100, so their conductances stay bit-identical. Values below are for Ra = 100 → 200.

| | Totals conserved, 100 → 200 | KLT from soma density, 100 → 200 |
|---|---|---|
| Soma input resistance, steady / peak | 2.6 / 4.8 → 2.8 / 5.0 MΩ | 4.9 / 7.7 → 5.0 / 7.8 MΩ |
| Unitary EPSG at soma: size, half-width | 5.1 mV, 527 µs → 5.8 mV, 476 µs | 6.2 mV, 691 µs → 7.0 mV, 624 µs |
| Same EPSG at mid-dendrite: size, half-width | 4.1 mV, 536 µs → 3.5 mV, 529 µs | 5.3 mV, 697 µs → 4.7 mV, 687 µs |
| Mid / soma EPSP size | 0.79 → 0.60 | 0.84 → 0.67 |
| Threshold at 0 delay: bilateral vs unilateral | 58.5 vs 82.2 → 76.6 vs 206.7 | 37.1 vs 44.0 → 41.5 vs 63.5 |
| Bilateral / summed unilateral EPSP | 0.98 → 1.00 | 1.00 → 1.02 |
| Coincidence window at 3%, model / Myoga EPSG | 165 / 181 → 139 / 149 µs | 223 / 249 → 209 / 230 µs |
| Pair threshold curve at the soma vs lumped (max difference) | 12% → 24% | 50% → 54% |
| Somatic spike at 1.5 / 2 / 3× rheobase | 16.3 / 19.5 / 26.6 → 19.3 / 24.5 / 29.9 mV | 21.4 / 26.1 / 29.6 → 25.5 / 29.7 / 35.8 mV |

What the higher resistivity does:
- **It isolates the soma from the dendrites.** Somatic EPSPs get larger and narrower, dendritic EPSPs are attenuated more, and the soma sees less of the dendritic membrane.
- **It makes the bilateral advantage much stronger.** With conserved totals, two EPSGs on one dendrite need 2.7× the bilateral threshold, up from 1.4×. The local depolarisation saturates the synaptic driving force, and dendritic KLT shunts it (Scott et al. 2010's argument).
- **Summation stays linear**, within 2%.
- **It moves each variant away from what it already matched.** With conserved totals the window narrows further from Myoga's 221 µs, and the curve drifts from Lehnert's "almost identical". With the soma-density gradient the window is still close (230 µs with the Myoga EPSG), but the spike grows further past 17 mV.

So 200 Ω·cm doesn't reconcile the two variants. 100 stays the default.
