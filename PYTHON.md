# Python port

A numpy/scipy port of the MATLAB model comparison. The Python package is the
reference: its results, pinned by `tests/test_golden.py`, are the ones the
analyses below use, on any hardware. The `.m` files are the historical original
it was ported from. They stay in the repo, untouched, but the Python no longer
has to match them: where this page says the port reproduces the MATLAB, that
describes how it was built, not a promise for future changes.

```bash
uv sync                                                    # or: uv sync --extra fast (Rust, below)
uv run pytest                                              # the tests, in parallel (-n 0: serially)
uv run pytest --cov                                        # plus line + branch coverage (as CI)
uv run scripts/combine_all.py --stim EPSGpair --node 3     # Combine_all.m
uv run scripts/making_threshold_graphs.py --node 3         # Making_Threshold_graphs.m
uv run scripts/making_threshold_graphs.py --only EPSGpair --jpg-grid   # the repo's jpgs
```

`dendrites.py`, `coincidence_window.py` and `somatic_spike.py` take `--quick`
(coarse delay grids, 1% threshold tolerance, fewer multiples), and
`making_threshold_graphs.py` takes `--points N`. The same code runs in seconds
instead of minutes; the numbers are for smoke tests, not results.
`tests/test_scripts.py` runs every script this way; those tests are marked
`slow` (`uv run pytest -m "not slow"` skips them), and CI runs them on Python
3.12 only.

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
`inputNode` stay **1-indexed** like the MATLAB code. Output arrays are
0-indexed: `y[:, node-1]` is compartment `node` of the 45-CPT model, and its 315
state columns keep MATLAB's column-major layout `[V(45), m, h, p, w, z, a]`, so
`z1` is `y[:, 225]`.

Code that isn't a port of a MATLAB file:

| Module | What it holds |
|---|---|
| `msoaxon/_common.py` | stimulus bundle and argument checks shared by both models |
| `msoaxon/_solve.py` | the ode15s stand-in: segmented BDF, pre-stimulus cache, spike stop |
| `msoaxon/_dispatch.py` | `run_model` / `spikes`: one entry point for `model="multi"` or `"two"` |
| `msoaxon/_parallel.py` | `map_tasks`; pass one `executor=` to reuse a process pool across calls |
| `msoaxon/_bisect.py` | the predictive search behind `coincidence.threshold`, `coincidence.window` and `somatic.rheobase` |
| `msoaxon/coincidence.py` | Myoga-style coincidence windows; `stim=`, `input_node=`, `input_node2=` pick the input sites, `model_kw` passes model options such as `morph` |
| `msoaxon/somatic.py` | rheobase and the Scott et al. 2005 somatic spike amplitude |
| `msoaxon/measure.py` | input resistance and time constant from a current step; soma trace on a grid |
| `msoaxon/_fast.py` | the optional Rust right-hand side: mode selection, self-check, `DriftWarning` |
| `rust/` | the `msoaxon-fast` extension (PyO3, built by maturin): `multi._rhs` in Rust |

## Fast mode: the 45-compartment right-hand side in Rust

`multi._rhs` spends nearly all of its time in numpy's per-call overhead on
45-element arrays, and the 45-compartment runs are where the scripts spend
their time. `rust/` ports that one function to Rust. scipy's BDF solver, the
breakpoints, the events and the pre-stimulus cache are unchanged; the solver
just calls a compiled right-hand side.

```bash
uv sync --extra fast        # builds rust/ with maturin; needs cargo (rustup.rs)
uv sync                     # numpy only: no Rust toolchain needed
```

**Without cargo.** Each release carries prebuilt `msoaxon-fast` wheels (Linux
x86_64 and aarch64, macOS arm64 and x86_64, Windows x86_64; one wheel per
platform covers Python 3.11-3.13). Download the one for your machine, or all of
them, from the release's assets and install it next to `uv sync`:

```bash
gh release download vX.Y.Z -p 'msoaxon_fast-*' -D wheels   # or save them from the release page
uv sync
uv pip install --find-links wheels msoaxon-fast
```

`uv run` leaves the wheel in place; a later plain `uv sync` removes it again
(use `uv sync --inexact`). The path source stays for development: it is what
`uv sync --extra fast` builds when you change `rust/`.

**On PyPI (prepared, not yet enabled).** The `publish` job in
`.github/workflows/wheels.yml` uploads the five wheels to PyPI by trusted
publishing after a `v*` tag builds them, but it is skipped until the
repository variable `PUBLISH_PYPI` is `true`, and no `msoaxon-fast` project
exists on PyPI yet. Once one does, `uv sync --extra fast` resolves the wheel
from the index and needs neither cargo nor a download; `[tool.uv.sources]`
stays for development only (it makes uv build `rust/` instead of using the
index). To enable it:

1. On pypi.org, create the project `msoaxon-fast` (a pending publisher does
   this on the first upload) with a trusted publisher: owner `JJBalisan`,
   repository `Comparison-2-Axon-Models`, workflow `wheels.yml`, environment
   left blank.
2. Set the repository variable: `gh variable set PUBLISH_PYPI --body true`
   (or Settings > Secrets and variables > Actions > Variables).
3. The next `v*` tag publishes; the crate's version is the package's, so
   each release has one `msoaxon-fast` version on PyPI.

The first upload of a version is final: PyPI never accepts the same file
name twice, so a re-run of a published tag fails at the upload.

`mso_axon(..., mode=...)` picks the right-hand side; `_dispatch.run_model`
passes `mode` on, so it also works in `model_kw`
(`coincidence.threshold("multi", 0.1, model_kw={"mode": "exact"})`) and
`measure.passive_step(..., mode=...)`. The two-compartment model has only its
numpy path and ignores it.

| `mode` | Right-hand side |
|---|---|
| `"exact"` | `multi._rhs`, the numpy reference |
| `"fast"` | the Rust port; warns and runs `"exact"` if the extension is missing |
| `"native"` | the Rust port under a Rust BDF integrator (opt-in, see below); warns and runs `"exact"` if the extension is missing |
| `None` (default) | `MSOAXON_MODE` if set (`exact`, `fast` or `native`), else `"fast"` when the extension is installed, else `"exact"` (silently) |

`MSOAXON_MODE=exact uv run scripts/...` runs a whole script on the reference,
worker processes included (pools hand their workers the parent's setting).

**Drift.** The port keeps numpy's floating-point evaluation order expression for
expression (`m**4` is `powf(4.0)`, not `powi`; `np.interp`'s fused multiply-add
on arm64 is reproduced). On darwin-arm64, where numpy and Rust call the same
libm, it is bit-identical to `multi._rhs`, and so are whole runs, thresholds and
windows. On other hardware numpy may use its own SIMD `exp`/`pow` while Rust
calls the system libm, so results can differ in the last bits; through a spike
upstroke that grows to the size of the cross-platform drift the golden traces
already tolerate. The first fast run in each process therefore checks the port
against `multi._rhs` on 24 fixed states (about 5 ms) and issues a
`msoaxon.DriftWarning` saying whether it was bit-identical on this machine.
Silence it with

```python
import warnings

import msoaxon

warnings.filterwarnings("ignore", category=msoaxon.DriftWarning)
```

or pass `mode="exact"` to use the reference. Everything `multi._rhs` handles is
ported: every stimulus (Synaptic, SynapticPair and Events included), the
lumped chain and `with_dendrites()` trees, `membrane()` overrides and every model type, so
`mode="fast"` never falls back to `"exact"` for a configuration.

**Speed** (Apple M1 Max, Python 3.12, best of 3; the machine was busy, so read
the ratios rather than the seconds):

| | `"exact"` | `"fast"` | Speed-up |
|---|---|---|---|
| one right-hand side call (315 states) | 63.5 µs | 3.5 µs | 18× |
| one run (EPSG pair, 20 ms) | 0.42 s | 0.15 s | 2.8× |
| `coincidence.threshold("multi", 0.1)` | 3.16 s | 1.19 s | 2.7× |
| `coincidence.window("multi", 0.03, delay_tol=1e-3)` | 3.78 s | 1.41 s | 2.7× |
| `scripts/benchmark_search.py --check` (all 20 searches, both models) | 274 s | 124 s | 2.2× |

Every value above came out identical in both modes, and `benchmark_search.py
--check` reported the same values, run counts and 0 mismatches in each. The rest
of a run is scipy's own Python (step control, the sparse finite-difference
Jacobian, LU), which caps the gain near 3×; the two-compartment model is
unchanged.

`tests/test_fast.py` compares the two right-hand sides on sampled states across
every covered configuration (to 1e-12 relative, and bit for bit on
darwin-arm64), and `tests/test_golden.py` runs every multi-compartment trace in
both modes. CI runs the suite without the extension (the fallback) on Linux and
macOS, and with it on both (`pytest-fast`), which also prints the self-check.

### Native mode: the BDF integrator in Rust

`mode="native"` (or `MSOAXON_MODE=native`) takes the rest of the time loop out
of Python: the same Rust right-hand side, integrated by a Rust BDF solver
([diffsol](https://github.com/martinjrobins/diffsol) 0.18.0, pure Rust, no
SUNDIALS or SuiteSparse) instead of scipy's. It is opt-in; the default stays
`"fast"` until native has been judged on more than the checks below.

What runs where: `mso_axon` still splits the run at the stimulus breakpoints,
keeps the pre-stimulus cache (keyed with the mode, so a native prefix never
serves a fast run) and builds the two event closures; `_solve._segment` hands
each segment to `_fast.Native.segment`, which calls
`msoaxon_fast.Integrator.run` (`rust/src/integrate.rs`). Rust then does every
step, Newton iteration, Jacobian, LU and event test of the segment and returns
the whole trajectory at once: the segment start, every accepted step, and the
end (`t_end` exactly) or the event. The integrator is built once per right-hand
side and reused for every segment of the run.

**Agreement.** Native's steps are diffsol's, not scipy's, so it is not bit for
bit any other mode; at the same `RTOL = ATOL = 1e-8` it agrees with the
reference within the solver tolerance. It is judged by the golden traces
(`tests/test_golden.py` runs every multi case in native too, at rtol 1e-4 /
atol 5e-3 mV and `t_last` exact, dendrite trees included) and by
`MSOAXON_MODE=native uv run scripts/benchmark_search.py --check` (0 flipped
spike decisions). The first native run in a process issues its own
`msoaxon.DriftWarning` saying so.

**Jacobian.** diffsol is given the pattern explicitly, not left to detect it:
`core::tree_jac_sparsity` of the run's own morphology (the Rust twin of
`multi._tree_jac_sparsity`; every stimulus couples a compartment only to
itself). The Jacobian is filled by forward-difference
products `J v` on two right-hand-side calls (step `1e-7 (1 + max |x|)`), one
per colour of the pattern.

**Newton solve.** The Newton matrix I - cJ has a fixed shape (each gate
couples only to itself and its own compartment's voltage; the voltage block is
a tree), so it is solved directly on that structure: every gate is eliminated
into its compartment's diagonal, then the tree is solved leaves to root, in
O(n) with no ordering or symbolic step. A linearisation whose pattern is not
exactly that, or whose tree pivot falls under 1e-3 of its diagonal, goes to
faer's sparse LU instead (`MSOAXON_NATIVE_LOG_FALLBACK=1` logs each one; none
occurs on the golden cases or the benchmark).

**max_step.** diffsol has no step cap, so `max_step` is emulated: before every
step the controller's growth bounds are narrowed so that no step is proposed
past the cap, and should one still be longer the integrator's stop time is
moved to one `max_step` ahead. Only the step length is bounded; the tolerances
are never changed.

**Events.** `spike_event` and `settled_event` carry a `spec` the native segment
reads; Rust evaluates the same functions as root functions, stops on a rising
crossing only (solve_ivp's `direction = 1`) and returns the state interpolated
at the root, which lies inside the last accepted step (that step's end is not
output, as in solve_ivp). The settled stop keeps the strict "every dV < 0".
With no input it can therefore fire right at `t_off` in a run that scipy
integrated to `t_end`: the same answer, an earlier stop.

`"Events"` runs under native like every other stimulus: its event table goes
to Rust whole and each event is evaluated at its own compartment, so the
Jacobian pattern above is unchanged.

## Port contract

[Issue #28](https://github.com/JJBalisan/Comparison-2-Axon-Models/issues/28)
plans a Rust port with its own integrator; Python stays the reference. What the
tests hold a new integrator to:

**Survive any integrator** (Python against Python): `tests/test_settled.py`
(`_dv` against `_rhs`, bit for bit), the step-time prefix checks in
`tests/test_analysis.py` and `tests/test_coverage.py` (events must never alter
the steps, as scipy's do not), the default-membrane check in
`tests/test_somatic.py`, and `tests/test_solve.py`.

**Retire or move:** the bitwise tests in `tests/test_fast.py` compare the Rust
RHS with the numpy one through scipy, so they become an RHS-only comparison.
`tests/test_golden.py`'s exact `t_last ==` needs the port to land on `t_end`.

**Can flip with another step sequence** (JJ decides, nobody loosens): the 28
jpg-grid halvings in `tests/test_coverage.py`, the one-flip budget at rel 2e-4
in `tests/test_analysis.py` and `tests/test_scripts.py`, and the somatic table
in `tests/test_somatic.py` (per-point tolerances from the `max_step` spread, see
"How stable the numbers are").

**Python-only behaviour to reproduce, or keep in Python:**
- `spike_amplitude`: a not-a-knot `CubicSpline` on a float `np.arange` grid at
  1 µs, after the `np.diff(t) > 0` filter.
- Events: brentq on BDF's dense output. `settled_event` jumps from -1 to +1;
  `spikes()` tells the stops apart by the final state, at `0.75 * factor`.
- `breakpoints()`: a set, strictly inside `(0, t_end)`, plus 5.0 and
  `stop + 5.0` for ramps. `input_end()`: one rule per stimulus.
  `pre_stimulus_is_quiet()`: `cuts[0] <= min(start, stop)`, and `two_cpt` keeps
  step, Synaptic and SynapticPair out of the shared prefix (`mso_axon` does
  not). All pinned in `tests/test_solve.py`.
- `interp_g` is `np.interp` with left = right = 0. `epsg_unitary` gives `-0.0`
  for -70 ms < t < 0, which the bitwise tests rely on.
- `matlab_round`, `half_width`, `passive_step`, `soma_on_grid`, the golden
  sampler; no `first_step`, `jac_sparsity` with column grouping, every accepted
  step as output.

**Scope:** the Synaptic tables come from numpy's PCG64 and `np.convolve`, whose
streams are not stable across releases, so the port takes the table as an input.
`two._rhs` has no Rust counterpart: whether `two_cpt` is in scope is open.
`uv.lock` pins two dependency sets (numpy 2.4.6 / scipy 1.17.1 below Python
3.12, 2.5.3 / 1.18.1 above), so "the Python reference" is itself two.

## Validation

There is no MATLAB here, so nothing was compared run-for-run. What was checked:

- **Constants:** `Area.mat`, `Fractions.mat` and `Coupling.mat` match the Python
  values to 1e-12 (`tests/`).
- **Rest:** with no input the 2-CPT model stays at `v0` to 1e-6 mV. The 45-CPT
  model settles about 0.9 mV away from it. The same happens in MATLAB, since
  `msoAxon.m` starts from rounded gating values and doesn't subtract resting
  currents.
- **EPSGpair thresholds vs `EPSGpair_Thresholds*.jpg`:** every plotted value is
  a multiple of q = 150/128, i.e. raw output of the halving search. Re-running
  that grid (`--jpg-grid`, delays 0–0.6 ms, nodes 3 and 5, both models) matches
  20 of 28 values exactly. The other 8 are one possible output away (2q ≈ 2.3
  units, 2.5–5%). The curve shapes, the steep rise between 0.2 and 0.5 ms, and
  the plateau positions all agree. The misses aren't random: node-5
  two-compartment sits one step low at four consecutive delays, and the node-3
  multi-compartment plateau is 94.9 vs 97.3. `tests/test_coverage.py` pins all
  28 Python values, with the jpg values alongside.
  - Tightening the tolerances 1000× (1e-9, `max_step` 0.02) doesn't move any
    Python value, so the Python side is converged. The gap is either ode15s
    error in the original runs or a code change after the jpgs were made (they
    came from an older 11-point `BinarySearch.m`). Telling those apart needs
    MATLAB.
- **Not validated numerically:** `step`, `ramp`, `ramp2`, `sine`, `EPSG`, and
  the synaptic inputs. Every stimulus, graph type and threshold sweep runs, and
  the traces look physiologically sensible, but there's no reference output to
  compare against.

## Deliberate differences

Where the port departs from the MATLAB, and why. Now that the Python is the
reference these are history; they are not reverted.

- **Solver:** `ode15s` → scipy `solve_ivp(method="BDF")` with the same
  tolerances (2-CPT 1e-6, `max_step` 0.1; 45-CPT 1e-8, `max_step` 0.1·tEnd,
  ode15s's default). Integration restarts at stimulus on/off times so the solver
  can't step over a narrow EPSG. The 45-CPT model gets a sparse Jacobian
  pattern.
- **Threshold search stops each run at its first spike** (`stop_on_spike`),
  since it only needs yes or no. The stop is found on the continuous solution,
  not the saved solver steps, so in principle it could count a crossing that
  falls between two steps. Across every production sweep (184 thresholds) the
  results are identical to running each simulation to `tEnd`. The search used to
  also count spikes in the saved steps, as `BinarySearch.m` does; that never
  decided a result, since a sampled crossing is one the stop already caught, so
  it was dropped. Sweep points also run in parallel (`workers`, default all
  CPUs).
- **Synaptic noise:** numpy's RNG can't reproduce MATLAB's `rng(seed); randn`,
  so `Synaptic`/`SynapticPair` runs are statistically equivalent, not identical.
- **`interp1q` out of range** returns 0, not NaN. `TwoCptODE.m` already guarded
  this for `SynapticPair`. `msoAxon.m` didn't, which is likely the "Problem
  including Synaptic pair" comment in `Combine_all.m`.
- **`BinarySearch.m`:** the two-compartment loop tested `location1` in its
  `while` condition. It now tests its own `location2`.
- **`Combine_all.m`:** the 2-CPT spike count read column 1 (soma minus soma,
  always 0); it now reads the axon column. The second `TwoCpt` run used
  identical arguments; it now uses node 5 so the "Node 5" panels mean what they
  say. The default `ramp` stop was equal to start, a zero-length ramp; it's now
  5.5.
- **`Graphing.m`:** the `Input` panel reuses the model's own stimulus function
  instead of a hand-copied version that indexed 2-CPT voltages with
  multi-compartment time steps. The `z1`/`a` line fits (broken:
  `graph.graph.tEnd`, length mismatch) were not ported.
- **EPSGpair threshold x-axis** shows the delay between the two EPSGs, like the
  jpgs do, rather than the second EPSG's absolute time (`5 + (i-1)/25`) that
  `Making_Threshold_graphs.m` plots.
- **Thresholds and coincidence windows** come from the predictive search in
  `msoaxon/_bisect.py`, not a MATLAB-style bisection. Full reruns against plain
  bisection move windows by under 0.4 µs and thresholds by under 1e-4
  (relative); one rounded value changed (Goldwyn, Myoga EPSG: 183 → 182 µs).
  These values are the reference. Single runs never depended on the search.
  Rheobase has used it since issue #19: the multi and two step responses rise
  smoothly through 10 mV like the EPSG pairs, so it takes 10–11 runs instead of
  15 and moves by under 5e-5 (relative); with KLT blocked the response is all or
  none and the search makes bisection's runs, to the same value.
- **Settled stop:** searches end a non-spiking run once the input is over and no
  spike can follow (`STOP_WHEN_SETTLED` in `msoaxon/_dispatch.py`; the
  single-run models leave it off). It stays: in the search benchmark it changed
  0 of 365 spike decisions, for about 6% less solver time.
- **`fitExp.m`** was not ported and won't be. It contains an unresolved merge
  conflict and fits 301 samples against a 49-point x-axis, and with MATLAB
  compatibility no longer a goal nothing here needs the fit. The `.m` file stays
  as part of the historical record.

## Quirks kept on purpose

These are in the MATLAB and are reproduced unchanged. They were kept so that the
port matched the MATLAB; they are now kept for continuity with the results on
this page. Any of them could be changed by a deliberate decision, which would
move the golden traces (`tests/golden_traces.json`) and the numbers here:

- `TwoCpt.m` enables sodium for `'Active-sodium'` (capital A), so
  `active-sodium` means *no* Na in the 2-CPT model but Na-only in the 45-CPT
  model.
- The 2-CPT model holds KLT inactivation `z` fixed at its resting value.
- In `msoAxon.m`, `step` always injects into the soma whatever `inputNode` is;
  `ramp` only applies for `t > 5` (hardcoded); `EPSG`/`EPSGpair` scale by the
  input compartment's area while the other stimuli use the soma's.
- `TwoCptODE.m` starts `ramp2` at t = 5 whatever `start` is; `msoAxon.m` starts
  it at `start`.
- `TwoCpt.m` has calibrated axonal Na (`gNa2`) only for nodes 3 (119) and 5
  (25.5). Every other node gets the placeholder 140, so 2-CPT runs at those
  nodes use an uncalibrated value. `two_cpt(..., gNa2=...)`
  replaces it.
- `Spiking.m` doesn't clear `reset` between columns.
- `SA(2)` uses `1.66/3` where `0.66/2` was probably meant. The stored `Area.mat`
  was built with it.
- `BinarySearch.m` reports the ceiling (`max`) when nothing fires, so a returned
  15000 means "no threshold found", not "threshold ≈ 15000". It also assumes
  spiking gets easier as input grows. That fails for `EPSG` with factor 30: a
  very large conductance holds the soma near 0 mV, so the axon never gets 30 mV
  above it. The `ramp` sweep's first point (stop = 5.1) returns the ceiling in
  both models for the same reason.

## Comparison with the literature

`uv run scripts/coincidence_window.py` reproduces the numbers below;
`figures/coincidence/` gets the plot and a JSON of every curve.

### Lineage

The 45-compartment model follows Lehnert et al. 2014
([doi](https://doi.org/10.1523/JNEUROSCI.4038-13.2014)). The two-compartment
model is the coupling-constant framework of Goldwyn, Remme & Rinzel 2019
([doi](https://doi.org/10.1371/journal.pcbi.1006476)).

`TwoCpt.m` changes Goldwyn's passive targets:

| | Goldwyn 2019 | `TwoCpt.m` |
|---|---|---|
| Input resistance | 8.5 MΩ | 10 MΩ |
| Soma time constant | 0.34 ms | 0.71 ms |
| Resting potential | −58 mV | −68 mV (from Lehnert) |

The old values survive as comments beside the new ones (`%8.5`, `%-58`).
`GOLDWYN_2019` in `msoaxon/two.py` restores them:
`two_cpt(..., v0=-58, r1=8.5, tau_est=0.34)`.

The coupling `node` looks up (forward and backward constants, the area ratio,
the KLT and h fractions and the axonal Na conductance) is a default, not a
limit: keyword-only arguments of `two_cpt` and `two.get_params` replace any of
it, each `None` keeping the lookup's value, so every call without them gives the
same bits as before.

| Argument | Meaning | Range |
|---|---|---|
| `gC` | coupling conductance [nS]; replaces what `r1` and the constants give, so `r1` stays at its default | > 0 |
| `couple12`, `couple21` | forward (soma → axon) and backward coupling constants | (0, 1) |
| `area_ratio` | axon to soma area ratio (sets the axon's time constant and capacitance) | > 0 |
| `klt1`, `klt2`, `h1`, `h2`, `na1` | fraction of the compartment's resting conductance that is KLT, h and (soma) Na; a compartment's fractions sum to at most 1 | [0, 1] |
| `gNa2` | axonal Na conductance [nS] | ≥ 0 |
| `soma_na_vhalf`, `axon_na_vhalf` | Na inactivation midpoint [mV] (default −62.5); see "Somatic spike" | finite |

They are validated before any run. The coupling plane of Goldwyn et al. 2019 is
a sweep of `couple12` and `couple21`; searches take them in `model_kw`
(`coincidence.threshold("two", 0.1, model_kw={"couple12": 0.5})`). Jing et al.
2025 give their soma-axon coupling as gax = 56.4 nS; read as `gC`, on node 3's
constants, that is `two_cpt(..., 3, "active-full", ..., gC=56.4)`. That figure comes from
issue #84; the paper was not read here, so whether their gax is this `gC` is
unconfirmed.

### Coincidence window

Measured the way Myoga et al. 2014 did in adult gerbil MSO at 35 °C
([doi](https://doi.org/10.1038/ncomms4790)):
- two identical EPSGs, their relative timing stepped in 20 µs;
- inputs slightly above the coincident threshold, so spike probability peaks
  near 100%;
- result: the full width of the spike-probability curve at half maximum. They
  measured **221 µs**.

| Model | Model EPSG (decay 0.18 ms) | Myoga EPSG (decay 0.3 ms) |
|---|---|---|
| 45-compartment | 201 µs | 219 µs |
| 2-compartment, `TwoCpt.m` | 199 µs | 217 µs |
| 2-compartment, Goldwyn 2019 calibration | 167 µs | 182 µs |

These widths use inputs 3% above threshold. At 0.5% they shrink to 70–91 µs. The
margin matters a lot, and the paper's "200 pS (~3%)" is ambiguous: 200 pS is
closer to 0.5% of their 43 nS EPSGs.

How the widths are computed:
- `coincidence.window()` does what the experiment does: it fixes the input at
  (1 + margin) × the coincident threshold (resolved to 1e-4) and searches over
  delay, to 0.1 µs, for where that input stops spiking. Probability is 50%
  there. About 20 simulations per window, within 0.5 µs of the same search at
  1e-7.
- Both searches predict rather than halve. At the 10 mV criterion the peak of
  axon − soma rises smoothly through 10 mV (7.7, 9.7, 10.0, 10.3, 17.8 mV at
  −1%, −0.1%, 0, +0.1%, +1% of threshold), so each run that doesn't spike says
  how close it came. That takes about 10 simulations per threshold on average
  (125 for the twelve benchmark thresholds, at most 13) instead of 15.
- Earlier versions read the crossing off a threshold curve on a 20 µs grid.
  Linear interpolation across a curve that bends upward read 0.1–4 µs low, most
  at the 0.5% margin.
- Noisy trials check this: 1% amplitude jitter plus 5 µs onset jitter give 200
  µs against `window()`'s 199 µs. If the noise is large enough to keep peak
  probability near 85%, the half-maximum width comes out about 10% wider.

Findings:
- With matched methods, the original calibration reproduces the measured window
  closely.
- Goldwyn's faster membrane narrows it by about 17%.
- An earlier estimate here (half-way point of the threshold curve, ~340 µs)
  measured a different quantity.
- Inhibition is not modelled (#77). The targets once it is: Myoga et al. 2014
  found an AP half-width of 208 µs with inhibition 0.2 ms before the pair,
  and the optimum moved by −73 µs with inhibition leading by 0.1 ms and by
  +41 µs with it lagging by 0.6 ms. Roberts et al. 2013
  ([doi](https://doi.org/10.1016/j.neuron.2013.04.022), dynamic clamp),
  Franken et al. 2015 ([doi](https://doi.org/10.1038/nn.3948), in vivo
  strychnine) and Casarez et al. 2026
  ([doi](https://doi.org/10.7554/eLife.111364), reviewed preprint) instead
  find that inhibition narrows the ITD function without shifting it. Myoga's
  model had no Kv1 history effect; both models here have KLT, so they can test
  which outcome follows.

#### Spike criterion

A run counts as a spike once axon − soma exceeds 10 mV, as in the MATLAB
EPSG-pair sweeps (`coincidence.FACTOR`, `somatic.FACTOR`). That fixed value
means different things in different models:
- **45-compartment and `TwoCpt.m`:** 10 mV is about 26% of their spike (38 mV of
  axon − soma). It sits on the foot of the takeoff, which comes 1–3% above that
  threshold.
- **Goldwyn's calibration:** 10 mV is about 50% of its spike (20 mV). Its
  response is graded rather than all-or-none: the peak grows smoothly from 5.8
  to 19.6 mV between 0.9× and 1.5× threshold. So for Goldwyn the criterion
  *defines* the threshold rather than locating it. A 5 or 15 mV criterion moves
  its coincident threshold by −13% or +8%. Its windows are robust, though: they
  move by at most 4%.
- A fixed 20 or 30 mV criterion is out of reach for Goldwyn's calibration
  altogether.

The criterion is an option (issue #21). `coincidence.threshold`,
`threshold_curve`, `window` and `probability_trials`, and `somatic.fires`,
`rheobase` and `spike_amplitude` take `criterion=`:
- a float is absolute, in mV; the default is 10.0, so nothing changes unless you
  ask;
- `coincidence.Relative(0.5)` is 50% of that configuration's own spike height.

`coincidence.spike_height()` measures the height as the largest axon − soma of a
coincident EPSG pair at 2× the 10 mV threshold, run with no spike stop. That
costs one threshold search plus one run. The height is flat (within 2.5%) from
1.5× to 3× threshold, so the 10 mV seed doesn't bias it. It is cached per
configuration: model, node, `v0`, EPSG kinetics, stimulus sites and every model
keyword. Membranes and morphologies are keyed by their settings, so the
dendritic and dendrotoxin variants, whose spikes differ, each get their own
height. `somatic` measures the height with an EPSG pair on the same membrane.
The MATLAB ports (`threshold.binary_search`, `spiking`,
`combine_all.py --factor`) keep their absolute 10 and 30 mV: those values are
part of what they reproduce.

At 50% of each model's height (model EPSG; `Relative(0.5)` against the 10 mV
default):

| Model | Spike height | 50% | Coincident threshold | 3% window | Rheobase (100 ms step) |
|---|---|---|---|---|---|
| 45-compartment | 38.4 mV | 19.2 mV | 48.01 → 48.81 nS (+1.7%) | 201.1 → 199.7 µs (−0.7%) | 3785 → 3852 pA (+1.8%) |
| 2-compartment, `TwoCpt.m` | 38.0 mV | 19.0 mV | 47.72 → 48.22 nS (+1.1%) | 199.0 → 198.0 µs (−0.5%) | 3848 → 3877 pA (+0.8%) |
| 2-compartment, Goldwyn 2019 calibration | 20.2 mV | 10.1 mV | 26.35 → 26.39 nS (+0.1%) | 167.5 → 166.9 µs (−0.4%) | 2070 → 2072 pA (+0.1%) |

- With the Myoga EPSG the shifts are the same size: the 3% windows go 219 →
  217.5, 216.8 → 215.5 and 182.4 → 182.5 µs.
- The 0.5% windows narrow by 0.5–0.7% for the 45-compartment and `TwoCpt.m`
  models.
- At 50% Goldwyn's criterion is today's, so its results don't move. The other
  two models' criterion moves onto the steep upstroke, which raises their
  thresholds 1.0–1.8% and narrows their windows by at most 0.75%: about 2 µs on
  200, inside the 0.5% vs 3% margin ambiguity above. Goldwyn still narrows the
  window by about 17% (−16.7% at 10 mV, −16.4% at 50%).
- At 25% or 75% of the height, the 3% windows change by at most 1.8% (4% for
  Goldwyn). The coincident thresholds move more: up to 3%, and Goldwyn's by −13%
  at 25%, where the criterion is below its takeoff.
- Goldwyn's rheobase is the same 100 ms step run with its calibration
  (`v0=-58, r1=8.5, tau_est=0.34`). `somatic.fires`, `response`, `rheobase` and
  `spike_amplitude` take a `model_kw` dict like `coincidence.threshold`, so the
  2070 pA is
  `somatic.rheobase("two", v0=-58, model_kw={"r1": 8.5, "tau_est": 0.34})`.
  `model_kw` also carries `mode="exact"`, `morph`, `model_type` and `max_step`;
  `mem=` stays as a shorthand for `model_kw={"mem": ...}`.
- The searches take about 7% more simulations at 50%: the interpolation still
  converges on the upstroke, a little less efficiently.

### Other benchmarks

Values are for the Python port.

| | Model | Literature |
|---|---|---|
| EPSP half-width at soma | 523–549 µs | 0.52–0.6 ms; PSP half-width 467 µs with the 0.3 ms EPSG (Myoga et al. 2014) |
| Spike initiation site | AIS (compartment 3) for an EPSG pair at 1.25–3× threshold; under sustained input 0–36% of spikes start at node 1 or 2 where 50 or more spikes fired, 8% at 350 Hz and the default `gE` (`initiation.py`, #80) | AIS; 5–60% of spikes start at node 1 or 2 under sustained input (Lehnert 2014) |
| Input resistance, active, at rest | 2.4 MΩ steady, ~4.7 MΩ peak with the leak rebalanced at −68 mV (`rebalance_rest=True`, as in the dendrites table); 2.5 / 4.8 MΩ for the default membrane, which settles at −68.13 mV where the conductances differ a little (`measure.passive_step`, #68) | 5 MΩ (Lehnert's model); ~7 MΩ measured |
| Somatic spike, from inflection (100 ms step) | 7.4 / 13.3 / 18.2 mV at 1.5× / 2× / 3× rheobase | 17 ± 2 mV mature, 5–15 mV near threshold, graded |
| Best frequency for half-wave sine input | 400–500 Hz | spike resonance 285 Hz, EPSC-train threshold minimum 317 Hz; subthreshold resonance 242 Hz (Mikiel-Hunter et al. 2016, [doi](https://doi.org/10.1371/journal.pcbi.1005166)); a suprathreshold measure belongs against the first two (#83) |

The initiation site is read from a finished run by `msoaxon.initiation`. It
cannot use a fraction of each compartment's own peak: the input enters at the
soma, so the EPSP half-rises before the AIS has fired (60 nS pair: soma 5.15 ms,
AIS 5.33 ms). Instead each spike gets the time of the steepest rise (peak dV/dt)
of every AIS and node compartment, taken just before the steepest rise of the
last node, and the site is the earliest. A compartment that rises no faster than
the soma did is following the EPSP and does not count. The AIS leads node 1 by
about 20 µs when it initiates.

Fraction of spikes starting at node 1 or 2 (spikes pooled over 3 noise seeds of
200 ms; `uv run scripts/initiation_site.py`):

| `gE` | 200 Hz | 350 Hz | 500 Hz | 800 Hz |
|---|---|---|---|---|
| 0.65e7 | 10% (103) | 36% (80) | no spikes | no spikes |
| 1.3e7 (default) | 1% (111) | 8% (165) | 50% (10) | no spikes |
| 2.6e7 | 1% (114) | 3% (191) | 11% (111) | no spikes |
| 5.2e7 | 0% (116) | 1% (197) | 4% (223) | no spikes |
| 1.04e8 | 0% (117) | 1% (201) | 1% (267) | 1 spike |

Spike counts are in brackets. Distal initiation grows as the drive weakens, and
only the weakly driven cells (350 Hz at 0.65e7, 500 Hz at 1.3e7) reach into
Lehnert's 5–60%; 800 Hz hardly fires at all. Of 39 spikes that started beyond
the AIS in the 350 Hz / 0.65e7 (3 seeds), 500 Hz / 1.3e7 (3 seeds) and 200 Hz /
0.65e7 (1 seed) runs, 14 had no AIS upstroke at all (the AIS was only dragged
along), and 12 of the other 25 led the AIS by under 10 µs. Node 1 and node 2 are
not told apart that finely. The pair protocol is AIS from 1.25× threshold up; at
1.02–1.05× node 1 leads the AIS by 0–11 µs, a tie at the criterion's resolution.

The 36% (350 Hz, 0.65e7) depends on how far ahead of the last node's upstroke
the others are searched: 38% for 0.4 ms, 36% for 0.6 ms (used), 23% for 0.8 ms
and 5% for 1.0 ms. The AIS leads the last node by 0.35–0.4 ms, and past that an
EPSP's own peak in the AIS starts to pass for its upstroke. The share does not
move with the other settings tried (the 60 mV rise that marks a spike, 40–80 mV;
how high a dV/dt peak must be to count, 30–90% of the highest).

An earlier version of this table listed a ~40 mV somatic spike as the clearest
mismatch. That figure was measured from rest during an EPSG pair at 2×
threshold, so it counted the synaptic depolarisation as spike. Measured like the
experiments (next section), the spike matches mature cells.

### Missing relative to current models

Each item has an issue (label `literature-2026-10`):

- glycinergic inhibition as a protocol: an inhibitory event is one `Event`
  with a reversal below rest (below), but no `dt_inh` axis on the window yet
  (#77)
- protocols on the event-list input (below): fibre counts, Poisson arrivals,
  pre-pulse trains, and reaching it from the searches, which build their own
  `SynParams` (#78)
- axon geometry as parameters: internode length, diameter, node length, AIS
  length and start and myelin capacitance are folded into the constants
  (#81)
- variation in cell properties along the frequency map: Ih density, V½ and
  time constant, Na oscillation parameters (#82)
- temperature as a parameter rather than two Q10 literals (#82)
- a Kv7 / M-type current (#85)

Dendrites, listed here until October 2026, are in (see below), and so is
the initiation-site readout (#80, `msoaxon.initiation`, table above), and so is
the free two-compartment coupling (#84, see "Lineage").

### Somatic spike

`uv run scripts/somatic_spike.py` reproduces this; `figures/somatic/` gets the
traces and a JSON.

Scott et al. 2005 ([doi](https://doi.org/10.1523/JNEUROSCI.1016-05.2005)) evoked
spikes with 100 ms somatic current steps and measured amplitude from the
inflection point. Mature cells (≥P21) give **17 ± 2 mV**: 5–15 mV near
threshold, growing with the stimulus. Blocking Kv1 (KLT) with dendrotoxin in
P20–21 cells raised it **from 15 to 37 mV**.

`msoaxon/somatic.py` measures the model the same way. The inflection is the peak
of d²V/dt² (of a cubic spline through the solver's samples) after the first 0.15
ms of charging: MSO neurons fire at step onset, while the membrane is still
charging at 50–80 mV/ms, so fixed dV/dt criteria pick up the charging rather
than the spike.

| Somatic spike amplitude | 1.5× rheobase | 2× | 3× |
|---|---|---|---|
| 45-compartment (`msoAxon.m`) | 7.4 mV | 13.3 mV | 18.2 mV |
| 45-compartment, KLT removed at soma + AIS, rest rebalanced (dendrotoxin) | 25.4 mV | 28.2 mV | 31.7 mV |
| 2-compartment (`TwoCpt.m`) | 6.4 mV | 13.0 mV | 17.8 mV |

- **The unmodified model already matches mature cells.** Lehnert et al. 2014
  report "∼10 mV" for their model too.
- **The KLT block reproduces the developmental mechanism:** the spike roughly
  doubles, and rheobase falls from 3784 to 200 pA. No tuning was needed.
- **Valid range: about 1.35× to 4.5× rheobase** (multi and two models, 10 mV
  criterion). `spike_amplitude` returns `distinct=True` there. Outside it the
  step still fires but `distinct` is `False`, and `amplitude`, `t_inflection`
  and `v_inflection` are `None` (`peak_above_rest`, `t_spike` and `trace` stay
  valid). The inflection is distinct when the d²V/dt² peak in the search window
  is the largest before the spike peak (so the takeoff was not left in the
  skipped charging window), is not at either end of the window, and is at least
  `somatic.D2_MIN` = 100 mV/ms². Below the range the soma is a smooth hump with
  no upstroke curvature (1.05× and 1.2× used to return 9–12 mV from a flat
  d²V/dt² ≈ 0 plateau). Above it the takeoff falls in the skipped window, and
  the old picker returned ≈0 mV from the upstroke next to the peak at 10×
  rheobase. The limits move with the membrane: with KLT removed every step from
  1.05× to 10× is distinct.
- **How stable the numbers are:** the pick is the maximum of the spline's
  d²V/dt², which lies on a solver sample, so it depends on where the solver put
  its steps. Rerunning with `max_step` 3, 4, 6 and 7 µs instead of 5 moves the
  dendrotoxin row by 0.00 mV and the multi rows by 0.00 mV; the 2-compartment
  rows move by up to 0.06 mV (1.5×), 0.00 (2×) and 0.20 mV (3×) (remeasured in
  exact mode: 6.372–6.435 and 17.576–17.775 mV; the other rows move by under
  0.0001 mV). `tests/test_somatic.py` gives each table point a tolerance from
  its spread: 0.06 mV, but 0.1 and 0.25 mV for those two. Before
  issue #56 the pick was the argmax of d²V/dt² of the linearly interpolated
  trace, a comb with a tooth at every solver step. The dendrotoxin row then
  moved 0.6–1.9 mV under the same changes (25.08 mV on macOS against 25.97 mV on
  Linux CI at 1.5×), and so did the controls (2-compartment 2× spread 12.8–17.6
  mV). The old control values (7.6 / 13.4 / 19.4 and 6.2 / 15.6 / 17.6 mV) were
  single teeth of that comb: the 2-compartment 2× value is 13.0 mV, not 15.6.
- **Spikes driven by EPSG pairs** measure 9–11 mV at 2× threshold, against 8.5 ±
  1.3 mV in vivo (van der Heijden et al. 2013). Near threshold the EPSP and
  spike merge and the inflection is ambiguous.

`mso_axon(..., mem=membrane(v0, ...))` exposes the knobs used for the block:
- `soma_klt_scale` and `ais_klt_scale`
- `soma_na_vhalf`: Scott et al. 2010
  ([doi](https://doi.org/10.1523/JNEUROSCI.2385-09.2010)) measured −77 mV
  for somatic Na inactivation, with about 92% of channels inactivated at
  rest (−58 mV); the model uses −62.5. With the model's own slope (7.77 mV),
  `constants.hinf` gives:

  | Rest | V½ −62.5 (default) | V½ −77 (Scott 2010) |
  |---|---|---|
  | −68 mV (this model) | 33% inactivated | 76% |
  | −58 mV (Goldwyn 2019, Scott 2010) | 64% | 92% |

  So the two numbers reproduce the measurement only together, and the default
  rests with about 8× the available somatic Na that was measured. Whether that
  is a deliberate choice in Lehnert et al. 2014 is open. `two_cpt` has the same
  knob as `soma_na_vhalf` (and `axon_na_vhalf`), without a rebalance: it
  subtracts the resting currents, so rest holds at `v0`. Its somatic Na
  conductance stays calibrated at −62.5 mV, so a shift changes how many
  channels are available at rest, not how many there are.

  **Proposed table, for JJ's decision (#79).** Nothing here is a default or a
  published number: `uv run scripts/na_availability.py` prints it and saves
  `figures/na_availability/`. Rest {−68, −58} mV × somatic V½ {−62.5, −77} mV,
  both models, V½ at the soma only, exact mode (no `native` shortcut). The
  45-compartment cells use `membrane(rest, soma_na_vhalf=..., rebalance_rest=
  True)` throughout, so its default cell differs a little from the published
  `msoAxon.m` row (rheobase 3793 against 3784 pA, spike 7.2 / 13.1 / 17.8
  against 7.4 / 13.3 / 18.2 mV). The EPSP is one somatic EPSG of a quarter of
  that cell's coincident threshold; the window is at 3% above that threshold.

  | 45-compartment: rest, V½ | inactivated | rheobase | Rin | EPSP half-width | 3% window | somatic spike 1.5× / 2× / 3× |
  |---|---|---|---|---|---|---|
  | −68, −62.5 | 33% | 3793 pA | 2.43 MΩ | 572 µs | 200 µs | 7.2 / 13.1 / 17.8 mV |
  | −68, −77 | 76% | 4405 pA | 2.39 MΩ | 561 µs | 193 µs | no inflection (peak 23.8 / 30.1 / 38.9) |
  | −58, −62.5 | 64% | 5289 pA | no stable rest | 329 µs | 137 µs | no inflection (peak 14.9 / 19.9 / 28.0) |
  | −58, −77 | 92% | 6438 pA | no stable rest | 318 µs | 133 µs | no inflection (peak 14.2 / 18.2 / 25.0) |

  | 2-compartment: rest, V½ | inactivated | rheobase | Rin | EPSP half-width | 3% window | somatic spike 1.5× / 2× / 3× |
  |---|---|---|---|---|---|---|
  | −68, −62.5 | 33% | 3846 pA | 2.43 MΩ | 588 µs | 199 µs | 6.4 / 13.0 / 17.8 mV |
  | −68, −77 | 76% | 4426 pA | 2.39 MΩ | 577 µs | 191 µs | no inflection (peak 23.5 / 29.1 / 37.8) |
  | −58, −62.5 | 64% | 1961 pA | 2.36 MΩ | 595 µs | 218 µs | no inflection (peak 13.2 / 16.2 / 21.5) |
  | −58, −77 | 92% | 2033 pA | 2.33 MΩ | 588 µs | 215 µs | no inflection (peak 13.0 / 15.9 / 21.0) |

  - **Scott's V½ moves the rheobase and little else:** at −68 mV, −77 raises
    rheobase by 16% (45-compartment) and 15% (2-compartment), and moves Rin,
    the EPSP half-width and the 3% window by 2–4%.
  - **The somatic spike has no distinct inflection** (`spike_amplitude`'s
    `distinct` is False) in any cell but the default one, at every multiple
    scanned (1.1× to 10×, 2-compartment): with less somatic Na, or a depolarized
    rest, the soma is a smooth hump that carries the axonal spike. The amplitude
    protocol of Scott 2005 does not apply there; the peak above rest is listed
    instead.
  - **The 45-compartment model has no stable rest at −58 mV:** with the leak
    rebalanced to hold −58, a 2 pA kick sets off a 0.18 mV oscillation with a
    1.4 ms period that does not decay, whatever the kick, solver mode or step
    size (the same at −60 mV, steady at −63 mV). The two-compartment model is
    steady there. Rin is only reported where the rest is steady.
  - **The models disagree at −58 mV:** rheobase rises 39% in the
    45-compartment model and falls 49% in the 2-compartment one, and the window
    narrows (137 against 218 µs). The two hold rest at `v0` differently (subtracted
    resting currents against a rebalanced leak reversal), which may be part
    of it; this was not isolated.
  - **A preset is proposed, not added:** a named Scott-consistent membrane
    (like `GOLDWYN_2019`: rest −58 mV, `soma_na_vhalf=-77`) would give the
    rows above as one option. It needs the 45-compartment instability at −58
    resolved first, and the full texts of Scott 2010 and Lehnert 2014's methods
    to say which V½ goes where.
- `rebalance_rest`: sets each leak reversal so rest stays at `v0`. Without it,
  removing KLT makes the cell fire with no input.

The defaults are `msoAxon.m`'s membrane, unchanged.

## Event-list input

`stim_type="Events"` drives either model with any number of synaptic
conductance events, each an `msoaxon.Event(t, node, g, tau, e_rev)`: an
EPSG-shaped conductance (`epsg_unitary`, peak `g` nS) from onset `t`,
entering compartment `node` with reversal `e_rev` and its own `(rise, decay)`.
Events add, never end, and may share a compartment or an onset. The run's `I`
scales every peak, so `I=1.0` uses them as given; `start`, `stop` and
`input_node` are ignored. Pass them in `SynParams(events=(...))`:

```python
syn = msoaxon.SynParams(
    t_end=20.0,
    events=(
        msoaxon.Event(4.7, 1, 30.0, (0.28, 1.5), -80.0),  # a leading IPSG
        msoaxon.Event(5.0, 1, 50.0),  # the EPSG pair
        msoaxon.Event(5.2, 1, 50.0),
    ),
)
t, x = msoaxon.two_cpt("Events", 0, 0, 1.0, 3, "active-full", 20, -68, 1, syn)
```

What it shares with the other stimuli: every onset is a breakpoint, the
stimulus-free prefix before the first onset is cached, the settled stop uses
the slowest event's decay, and `two_cpt` maps node 1 to the soma and any other
node to the axon. Two things it does differently, on purpose:

- **It is not bit-identical to `EPSGpair`.** The pair multiplies amplitude and
  driving force after summing its two waveforms; a list sums per event. The
  right-hand sides agree to 1e-12, and a single event is bit-identical to the
  matching `EPSG` (`tests/test_events.py`). The old names keep their own code
  paths, which the golden traces pin.
- **Each event drives its own compartment**, in `two_cpt` too, where `EPSG`
  keeps `TwoCptODE.m`'s soma driving force even for axonal input.

**Rust.** `mode="fast"` and `mode="native"` run it in the Rust right-hand side
(`Stim::Events`, #107): the event arrays of `synaptic.event_arrays` are handed
over once, and each call adds the events' currents per compartment in list
order, as `np.bincount` does, so fast stays bit-identical to `"exact"` on
darwin-arm64 (`tests/test_fast.py`). Native agrees with the reference to the
golden tolerance (`tests/test_events.py`); there is no Events golden trace.

Not yet: the searches (`coincidence`, `somatic` build their own `SynParams`),
and a breakpoint policy for dense lists, where a Casarez-style hundred-synapse
drive would mean hundreds of solver restarts (#78).

## Dendrites

`mso_axon(..., morph=with_dendrites())` adds two dendrites, following the
dendritic variant in Lehnert et al. 2014 (their Fig. 8). The lumped model's 8750
µm² "soma" stands for soma plus dendrites: in the paper's words it "combines the
somatic and dendritic membrane surface", sized for 70 pF. So the dendrites take
their membrane from it rather than being added on top:
- two unbranched dendrites, each 200 µm × 5 µm, five compartments each;
- the soma reduced to 2467 µm², so total membrane stays 8750 µm²;
- somatic Na density scaled up so total Na is unchanged;
- no Na or KHT in the dendrites;
- KLT and h decaying along the dendrites with a 74 µm length constant (Mathews
  et al. 2010, [doi](https://doi.org/10.1038/nn.2530)).

How to use it:
- Dendrite compartments are appended after the axon, so compartments 1–45 keep
  their meaning. Numbers 46–50 are the lateral (ipsilateral) dendrite, 51–55 the
  medial (contralateral), proximal to distal.
- New stimulus `EPSGbilateral`: one EPSG at `input_node` at `start`, the other
  at `input_node2` at `stop`.
- `membrane(..., morph=D, dendrite_klt_scale=...)` scales dendritic KLT.
- `uv run scripts/dendrites.py` runs every check below.

Two things the paper doesn't state:
- **Axial resistivity.** The model's 100 Ω·cm is kept by default; Mathews
  used 200. `with_dendrites(dendrite_ra=200)` (script: `--dendrite-ra 200`)
  changes it for the dendrites only, and the results are in the section after
  this one.
- **Whether total KLT and h were conserved.** `conserve_totals=True` (default)
  keeps soma + dendrite totals equal to the lumped model's. `False` starts the
  gradient at the soma's Table 2 density, halving total KLT.

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
| Somatic spike at 1.5 / 2 / 3× rheobase | 7.4 / 13.3 / 18.2 mV | 16.7 / 20.9 / 25.9 mV | 21.1 / 26.1 / 30.8 mV | 17 ± 2 mV mature |

Conserving the totals is the default because it reproduces what Lehnert et al.
reported for their own dendritic variant: tuning almost identical to the lumped
model. It also fixes the EPSP width and keeps the somatic spike near mature
size.

Neither variant matches everything:
- **Conserving totals** narrows the coincidence window below Myoga's value.
- **The soma-density gradient** matches the input resistance and the window, but
  leaves the spike too large.

Both reproduce the qualitative dendritic results: attenuation along the cable,
EPSP sharpening by dendritic KLT, and the bilateral advantage.

Not yet compared: Winters et al. 2017
([doi](https://doi.org/10.1523/JNEUROSCI.3110-16.2017)) found dendritic
spontaneous EPSPs growing with distance from the soma (length constant 53.6
µm) yet arriving there at a uniform ~0.2 mV, and Callan et al. 2020
([doi](https://doi.org/10.1523/JNEUROSCI.1055-20.2020)) found single fibres
contacting a dendrite at a median of four sites. Both are calibration targets
for synaptic strength against distance once inputs can be placed as fibres
(#78).

### Dendritic axial resistivity 200 Ω·cm

`with_dendrites(dendrite_ra=200)` uses Mathews et al.'s 200 Ω·cm in the
dendrites. The soma and axon keep 100, so their conductances stay bit-identical.
Values below are for Ra = 100 → 200.

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
| Somatic spike at 1.5 / 2 / 3× rheobase | 16.7 / 20.9 / 25.9 → 19.8 / 24.5 / 29.5 mV | 21.1 / 26.1 / 30.8 → 25.8 / 30.9 / 35.4 mV |

What the higher resistivity does:
- **It isolates the soma from the dendrites.** Somatic EPSPs get larger and
  narrower, dendritic EPSPs are attenuated more, and the soma sees less of the
  dendritic membrane.
- **It makes the bilateral advantage much stronger.** With conserved totals, two
  EPSGs on one dendrite need 2.7× the bilateral threshold, up from 1.4×. The
  local depolarisation saturates the synaptic driving force, and dendritic KLT
  shunts it (Scott et al. 2010's argument).
- **Summation stays linear**, within 2%.
- **It moves each variant away from what it already matched.** With conserved
  totals the window narrows further from Myoga's 221 µs, and the curve drifts
  from Lehnert's "almost identical". With the soma-density gradient the window
  is still close (230 µs with the Myoga EPSG), but the spike grows further past
  17 mV.

So 200 Ω·cm doesn't reconcile the two variants. 100 stays the default.
