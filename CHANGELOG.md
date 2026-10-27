# Changelog

Versions follow [semantic versioning](https://semver.org/). The Python
package is the reference: a release pins the numbers the analysis produces,
and a change that moves them beyond the golden-trace tolerance is a new
minor version at least. The `msoaxon-fast` crate in `rust/` carries the
same version as the package, bumped in the same release PR.

## 1.8.0 (27 Oct 2026)

Phase 1's third round, three lanes reviewed and merged in parallel: native
mode becomes the default, the two-compartment model can trade its KLT for a
leak that keeps rest, and the coupling plane of Goldwyn et al. 2019 gets its
sweep. The golden traces, the 28 pinned thresholds and the benchmark's 300
runs hold in every mode and no tolerance moved. Under the new default the
results move within the solver tolerance (at most 0.31 of the golden
tolerance over the 15 golden cases, with every last-time equal); `mode="fast"`
or `MSOAXON_MODE=fast` gives the previous numbers bit for bit.

### Added

- `klt_to_leak` (`"chord"` | `"slope"` | `None`) on `two_cpt` and
  `two.get_params`, replacing the soma's and the axon's KLT by a leak that
  keeps rest; no existing result changes (#132).
- `scripts/coupling_plane.py`: the two-compartment window and threshold over
  the (forward, backward) coupling plane of Goldwyn et al. 2019, with the
  45-compartment model located on it (#84). A cell counts as spiking only
  with a spike of at least 30 mV, since graded responses also clear the
  10 mV criterion. `--quick`, `--step`, `--ceiling`, `--mode`, `--workers`;
  saves `figures/coupling_plane/coupling_plane.{png,json}`. With
  `tests/test_coupling_plane.py` and the plane-sweep protocol in PYTHON.md
  under the free two-compartment coupling (#130).

### Changed

- `mode="native"` (diffsol's BDF in Rust) is now the default whenever the
  `msoaxon-fast` extension is installed (#106). Results move within the
  solver tolerance; `mode="fast"` and `mode="exact"`, or
  `MSOAXON_MODE=fast|exact`, give the previous results. Without the
  extension the default is still `"exact"` (#131).
- A diffsol failure in native mode raises
  `RuntimeError("integration failed on [a, b]: diffsol: ...")`, like
  scipy's in the other modes, instead of `ValueError`. Both DriftWarnings
  name the new default (#131).
- `scripts/inhibition_window.py` builds its equal-leak column from
  `membrane(klt_to_leak="chord")` instead of its own helper (printed tables
  unchanged; JSON digits differ at 1e-7) and adds the two-compartment
  equal-leak column; `scripts/resonance.py` adds two-compartment
  substitution arms (#132).

## 1.7.0 (26 Oct 2026)

Phase 1's second round and the three Python lanes from the 25 October
literature review, reviewed and merged in parallel: native mode fills its
Jacobian analytically, the searches take event lists, the two-compartment
model's coupling and Na availability are free parameters, inhibition is an
event with a timing axis, and the resonance protocol replaces the old
half-wave-sine row. Model results are unchanged: fast stays bit-identical
to the reference on darwin-arm64, the golden traces, the 28 pinned
thresholds and the benchmark's 300 runs hold in every mode, no tolerance
moved, and the default mode is still `"fast"`.

### Added

- `events=` on `coincidence.threshold`, `threshold_curve`, `window`,
  `probability_trials` and `somatic.fires`, `response`, `rheobase`: an
  event tuple or an event source (a callable giving the events for a delay
  and an amplitude); calls without it are unchanged. On top of it:
  `coincidence.window_edges`, `best_delay`, `trial_spikes`,
  `bilateral_probability`, `fibre_window_trials`, `itd_sweep`,
  `prepulse_thresholds`; `synaptic.PairSource`, `Fixed`, `fibre_events`,
  `bilateral_events`, `current_event`, `epsc_train`;
  `measure.curve_summary`, `train_rin`; and `scripts/event_protocols.py`
  for the fibre-drive, pre-pulse, EPSC-train and ITD protocols (#78).
- Inhibition as an event: `synaptic.ipsg`, `synaptic.Inhibition`,
  `PairSource(inh=...)`, and `dt_inh` / `inh` on `coincidence.window` and
  `probability_trials` (`inh_jitter` on the latter), with
  `scripts/inhibition_window.py` (#77).
- Keyword-only `gC`, `couple12`, `couple21`, `area_ratio`, `klt1`, `klt2`,
  `h1`, `h2`, `na1` and `gNa2` on `two_cpt` and `two.get_params`, each
  defaulting to the `node` lookup's value, validated before any run (#84);
  `soma_na_vhalf` and `axon_na_vhalf` on the same (default −62.5 mV),
  `TwoCptParams.na_vhalf1` and `na_vhalf2`, and `scripts/na_availability.py`
  for the rest × Na V½ 2 × 2 of both models (#79).
- `msoaxon.resonance` (the soma's impedance, and the sine-train and
  EPSG-train thresholds over frequency), `scripts/resonance.py`, and
  `membrane(klt_to_leak="chord" | "slope")`, which replaces KLT by a leak
  that keeps the rest potential; the default membrane is unchanged (#83).

### Changed

- Native mode fills the BDF Jacobian analytically (`MultiRhs::jacobian`)
  instead of by forward differences: the same steps on the EPSGpair and
  step harness cases (the Synaptic goldens move by a few steps, within
  tolerance), about 8–10% less CPU per run, and an exact Jacobian where the
  old fill was off by up to 1.5e-4 relative (#105).
- PYTHON.md: the "Best frequency for half-wave sine input" row is replaced
  by subthreshold, spike and EPSG-train resonance rows against
  Mikiel-Hunter et al. 2016, with a Resonance section that says why the
  model's frequencies sit 1.3–1.5× above the paper's; the two-compartment
  keyword table and the proposed Na-availability 2 × 2 are documented.
- CI: the `pytest` job's timeout is 20 minutes; the 3.12 slow-test run,
  which runs every script's `--quick` smoke test without the Rust
  extension, outgrew 10.

## 1.6.0 (26 Oct 2026)

Phase 1 of the Rust port (#28), three parallel lanes: the `"Events"` stimulus
joins the Rust right-hand side, native mode meets the 8x bar with a structured
Newton solve, and CI builds the crate on every wheel platform with a cache.
Model results are unchanged: fast stays bit-identical to the reference on
darwin-arm64, the golden traces, the 28 pinned thresholds and the benchmark's
300 runs hold in every mode, and the default mode is still `"fast"`.

### Added

- `mode="fast"` and `mode="native"` run the `"Events"` stimulus in the Rust
  right-hand side instead of falling back to `"exact"` with a warning: the
  event table goes to Rust whole and each event is evaluated at its own
  compartment in numpy's summation order, so fast is bit-identical to
  `"exact"` on darwin-arm64 and native agrees within the golden tolerance
  (#107).
- Native mode solves the Newton matrix directly on the model's
  tree-plus-diagonal structure, with faer's sparse LU as the fallback (none
  occurs on the golden cases or the benchmark). The integrator is about 1.4x
  faster and a threshold search about 8x faster than fast mode (0.121 s
  against the 0.963 s baseline; the coincidence window 8.8x). PYTHON.md's
  "Native mode" section describes the solve and the `max_step` emulation
  (#105).

### Changed

- Native mode bounds `max_step` through the step controller's growth limits
  rather than a stop-time hop at every capped step, so a capped 110 ms step
  refactorises 1492 times instead of 13008; uncapped runs take the same steps
  as before (#105).
- `rust/Cargo.toml` depends on `diffsol-la =0.3.0` directly (the linear-solver
  trait diffsol does not re-export).
- CI: `Swatinem/rust-cache` in the fast and wheel jobs (second-run job times
  about halved), the `prototype/diffsol-real` crate built and tested on
  ubuntu, `prototype/diffsol-probe` frozen as a snapshot, and the fast jobs
  run on Windows and linux-aarch64 too, so native mode has run on every wheel
  target (#106, #109, #113, #114). On linux-aarch64 Rust fuses the
  multiply-add in `interp` and numpy's build does not, one ulp apart; the
  non-bitwise comparison allows 1e-15 absolute for it.

### Fixed

- Native mode: a non-terminal root crossed in the step before a stop time no
  longer panics inside diffsol (`StopTimeBeforeCurrentTime`); the stop time is
  re-armed after every root (#106).
- `prototype/diffsol-real` builds against the merged `core.rs` again
  (its stimulus spec lacked the new `events` field), which is what its
  new CI step exists to catch (#113).

## 1.5.0 (26 Oct 2026)

The first release with the Rust integrator (#28, Phase 0 and lane D): an
opt-in `mode="native"` that runs the whole 45-compartment integration in Rust,
a fix to the synaptic inputs' driving force that changes results for off-soma
synaptic runs only, and the Phase 0 tooling. The golden traces, the 28 pinned
thresholds, the somatic table and the benchmark's 300 runs are unchanged in
every mode, and the default mode stays `"fast"`.

### Added

- `mode="native"` (or `MSOAXON_MODE=native`): the Rust right-hand side under a
  Rust BDF integrator (diffsol 0.18.0, pure Rust) takes the whole time loop
  out of Python; `mso_axon` still splits the run at the stimulus breakpoints
  and keeps the pre-stimulus cache. It agrees with the reference within the
  solver tolerance on every golden trace (`tests/test_golden.py` runs each
  multi case in native too, at most 31% of the tolerance consumed, `t_last`
  exact), `benchmark_search --check` flips none of its 300 spike decisions,
  and thresholds and windows run 7.2x / 7.9x faster than fast mode. Opt-in:
  the default stays `"fast"` until #106; `"Events"` falls back to `"exact"`
  in native as in fast (#107). PYTHON.md gains "Native mode: the BDF
  integrator in Rust" (#102).
- `scripts/phase0_harness.py`: golden-margin, search-timing and step-count
  measurements per mode, the baseline the integrator is judged by (#101).
- `prototype/diffsol-probe/` and `prototype/diffsol-real/`: the Phase 0 probe
  (a toy system with the model's sparsity and its REPORT.md) and the E1
  experiment crate (diffsol over the real right-hand side). Reference code,
  not shipped behaviour (#100, #104).

### Changed

- The Rust right-hand side lives in a pyo3-free `rust/src/core.rs` with a
  thin binding in `lib.rs`, and the crate also builds as an rlib. No
  behaviour change (#99).
- The `msoaxon-fast` crate depends on diffsol, pinned at `=0.18.0` until
  #110: the extension grows from 0.70 to 1.85 MB (about 1.15 MB more per
  wheel), a clean build takes about ten times longer (about 100 s on an M1
  Max), and `rust/Cargo.lock` holds 162 crates, all pure Rust with
  permissive licences.

### Fixed

- `Synaptic` and `SynapticPair` inputs at `input_node != 1` drove their
  compartment with the soma's driving force, `V[0] - VsynE`, in the
  reference and, bit for bit, in the Rust port. Both now use the input
  compartment's voltage, in the same evaluation order, so fast stays
  bit-identical to exact. Every golden synaptic case has `input_node=1`, so
  the golden traces are unchanged; results change only for off-soma synaptic
  runs. The native Jacobian pattern drops the entry the bug had required
  (#103).

## 1.4.0 (25 Oct 2026)

The last batch before the Rust port (#28): a general synaptic input, a
measurement of where spikes start, two corrected helpers, documentation
brought up to the 25 Oct 2026 literature review, and a CI that can hold long
port runs. Model results are unchanged: the golden traces, the 28 pinned
thresholds, the somatic table and the benchmark's 300 runs all hold. Only
`measure.passive_step` on a membrane that is not rest-balanced returns
different, now correct, numbers.

### Added

- `stim_type="Events"` for both models: a list of synaptic conductance events
  (`msoaxon.Event`, in `SynParams.events`), each with its own onset,
  compartment, peak, kinetics and reversal, so inhibitory and multi-fibre
  inputs can be composed; every onset is a breakpoint and the settled stop and
  quiet-prefix cache apply. It sits beside the EPSG stimuli rather than
  replacing them, whose float order the golden traces pin; a single event is
  bit-identical to `EPSG`. `mode="fast"` falls back to `"exact"` with a
  warning until the Rust side follows (#78).
- `msoaxon.initiation`: reads which compartment of the 45-compartment model
  fires first for each spike (`find_spikes`, `group_fractions`,
  `sustained_input`), with `scripts/initiation_site.py` sweeping input
  frequency and conductance. PYTHON.md's spike initiation row is now a
  measurement: the AIS for an EPSG pair, and 1–8% of spikes at node 1 or 2
  under sustained input at the default drive (#80).

### Changed

- The Rust toolchain is pinned to 1.99.1 (`rust/rust-toolchain.toml`,
  `rust-version`, and the CI and wheels workflows) so compiler releases cannot
  change CI results or the extension's bits unasked.
- CI: every job has a `timeout-minutes`, runners are pinned (`ubuntu-24.04`,
  `macos-26`, `windows-2025-vs2026`), and the suite also runs on Windows and
  Python 3.14. The two tests that need a forkserver context are skipped on
  Windows (#70).
- PYTHON.md follows the 25 Oct 2026 literature review: the "Missing relative
  to current models" list names each remaining gap with its issue; the
  input-resistance benchmark names the membrane behind each value (2.4 / 4.7
  MΩ rebalanced, 2.5 / 4.8 MΩ default); the best-frequency row compares
  against Mikiel-Hunter 2016's 285 / 317 Hz; the coincidence-window findings,
  the `soma_na_vhalf` note (33% of somatic Na inactivated at rest against
  Scott 2010's 92%) and the dendrites section carry their literature targets;
  and a section documents the event-list input (#87).

### Fixed

- `measure.passive_step` gave wrong input resistance and time constant on a
  membrane that is not rest-balanced (61 / 69 MΩ and τ 0 for the default
  membrane, instead of about 2.5 / 4.8 MΩ and 0.12 ms) and for a positive
  step current: it subtracted `v0` rather than the drift-free baseline and
  assumed a negative step. It now subtracts an identical run without the
  step, works for either sign of `I`, and raises `ValueError` for `I == 0`.
  Results on rest-balanced membranes, and so the dendrites table, are
  unchanged (#68).
- `coincidence.threshold`, `threshold_curve`, `window`, `somatic.rheobase`
  and `threshold.binary_search` raise `ValueError` for a `guess`, `rel_tol`,
  `delay_tol` or `zoom` that is not > 0 (and `_bisect.smallest_crossing` for
  a `descent` not > 1), instead of searching for ever (#69).

## 1.3.3 (23 Oct 2026)

Preparation for the Rust port (#28): every number the port must copy has one
source, the solver driver's rules have direct tests, and the Rust right-hand
side takes its constants from Python. Results are bit-identical; one NaN edge
far outside any script's reach is fixed.

### Changed

- The magic numbers shared by the Python models and the Rust port are named
  once: `_solve.RAMP_RISE`, `EPSG_TAIL`, `SETTLED_FRACTION` and
  `SPIKE_TEST_FRACTION`; the solver tolerances and step sizes (`multi.RTOL`,
  `ATOL`, `MAX_STEP_FRACTION`; `two.RTOL`, `ATOL`, `MAX_STEP`); the search
  constants in `_bisect` (`PREDICT_MARGIN`, `STALL_PROBES`, `PEAK_FLOOR`,
  `LOG_SLOPE`, `MIN_STEP_UP`); `synaptic.EPSG_NORM`; and the analysis
  defaults `START`, `T_END`, `FACTOR`, now defined in `_common` and
  re-exported by `coincidence`, `threshold` and `somatic`. Each holds exactly
  the literal it replaced (#28).
- `_common.Stimulus` carries the ramp `slope`, `t_top` and `t_end_local`,
  the sine `two_pi_f` and the EPSG `inv_peak`, computed once per run instead
  of on every right-hand side call, with the operand order the inline
  expressions had. `multi.initial_state` and `two.initial_state` build the
  initial state; `_dispatch.axon_column` replaces three copies of the
  compared-column rule, and `spikes` and `response` share one body (#28).
- The Rust right-hand side (`msoaxon-fast`) no longer holds its own copies of
  the reversal potentials or the gating-kinetics constants: `multi.V_NA`,
  `V_K`, `V_H`, `constants._P_TEMP` / `_A_TEMP` and the literals of the
  gating functions are passed in by `_fast` when a run is built, so a change
  in Python can no longer be mistaken for hardware drift by the self-check.
  `msoaxon_fast.interp` exposes the port's `np.interp` for the tests (#28,
  #71).
- `cargo test` runs in the `pytest-fast` CI job: nine crate tests, among them
  a bit-exact guard on one right-hand-side evaluation and the interpolation
  edge cases (#71).

### Fixed

- `synaptic.epsg_unitary` returned `nan` for onsets more than about 71 ms in
  the past (35 ms for a 0.05 ms rise), which crashed both models for an
  EPSGpair whose second onset was 80 ms after the first. It now returns
  `-0.0`, as it already did for slightly negative times; every reachable
  value is unchanged bit for bit, and `Epsg::unitary` in the Rust port does
  the same (#71).

### Removed

- The unreachable midpoint branch in `_bisect.crossing` (#71).

### Tests

- `tests/test_solve.py` pins the solver driver's rules: `breakpoints`,
  `input_end` for every stimulus, `pre_stimulus_is_quiet`, `spike_event`,
  `settled_event`, `interp_g`, and the two/multi quiet-prefix difference
  (#71).
- `tests/test_analysis.py` pins `epsg_unitary` after, just before and far
  before onset, and an EPSGpair 80 ms apart on both models (#71).
- `tests/test_somatic.py` gives each somatic table point its own tolerance
  from the measured `max_step` spread: 0.1 and 0.25 mV for the 2-compartment
  1.5x and 3x points, 0.06 mV elsewhere (#28).
- PYTHON.md gains a "Port contract" section: which tests survive a new
  integrator, which retire, which can flip, and the Python-only behaviour the
  port must reproduce or keep (#28).

## 1.3.2 (23 Oct 2026)

Housekeeping; results are unchanged. The parameter bundles become typed
dataclasses, the worker pool and `spiking()` get a public, checked
interface, and the Rust crate takes the package's version from here on.

### Changed

- The parameter bundles are frozen dataclasses with typed, documented
  fields: `_common.Stimulus`, `multi.Morphology`, `multi.Membrane` and
  `two.TwoCptParams` replace the `types.SimpleNamespace` objects. Attribute
  names and the cache keys are unchanged, so results are bit-identical
  (golden traces, `tests/test_settled.py`, the fast-mode BITWISE tests and
  `_fast.self_check`); `two._rhs` is about 4% faster. `tests/test_params.py`
  pins frozenness, identity hashing and the key tuples (#42).
- The Na conductance expression `g * m**4 * (0.993 * h + 0.007)` is written
  once, as `constants.na_open(m, h, g)`, with the evaluation order every copy
  had (#45).
- `msoaxon.process_pool` and `msoaxon.map_tasks` are public, re-exported
  from the package; the analysis scripts use them instead of the private
  `msoaxon._parallel` (#45).
- `spiking.spiking` accepts `model="two"` or `"multi"` (and the `"Two"` /
  `"Mult"` aliases of the port) and raises `ValueError` for any other name
  instead of treating it as the 45-compartment model (#45).
- `msoaxon-fast` carries the same version as `msoaxon`, bumped in the same
  release PR, so release wheels are named after the release;
  `tests/test_version.py` keeps `rust/pyproject.toml`, `rust/Cargo.toml` and
  `rust/Cargo.lock` in step with the package. The wheels workflow can publish
  to PyPI by trusted publishing, but stays off until the `PUBLISH_PYPI`
  repository variable is `true` (PYTHON.md "Fast mode" has the steps) (#57).
- Internal: `_fast.effective` and `_fast.resolve` lose their always-true
  `covered` parameter, left over from #37 (#45).

### Fixed

- Threshold searches that climb to their ceiling no longer emit numpy
  `overflow` RuntimeWarnings: `_bisect._predict`'s log-amplitude
  extrapolation along a flat curve overflows `np.exp`, which is now
  suppressed where it happens; the returned `inf` ("does not reach the
  criterion") is unchanged, and so are the searches (#45).

### Notes

- Those overflow warnings (e.g. a 30 mV criterion with Goldwyn's
  calibration) never came from the model: `smallest_crossing` clamps the
  prediction's `inf` to 2x, so every result before the fix was already
  right (#45).

## 1.3.1 (23 Oct 2026)

The published somatic spike table moves: `spike_amplitude`'s inflection
pick is now stable across platforms, and the old values turn out to have
been artefacts of the old pick. Everything else (golden traces, thresholds,
windows, script outputs) is unchanged.

### Fixed

- `somatic.spike_amplitude` picks the inflection from the second derivative
  of a cubic spline through the solver's own samples, evaluated on the 1 us
  grid, instead of the argmax of a twice-differentiated linear interpolant.
  The old pick turned every solver step into a tooth of d2V/dt2 and the
  argmax into a lottery over the teeth, so it followed the solver's step
  placement: the dendrotoxin row at 1.5x rheobase gave 25.08 mV on macOS and
  25.97 mV on the Linux runners, and under small `max_step` changes the
  2-compartment 2x value spread 12.8-17.6 mV. The published table changes
  with the pick (at 1.5x / 2x / 3x rheobase): 45-compartment 7.4 / 13.3 /
  18.2 mV (was 7.6 / 13.4 / 19.4), 2-compartment 6.4 / 13.0 / 17.8 mV (was
  6.2 / 15.6 / 17.6), dendrotoxin 25.4 / 28.2 / 31.7 mV (was 25.1 / 27.9 /
  31.9), and the dendritic variants' rows in PYTHON.md likewise. The new
  values move by at most 0.06 mV under the same step changes (0.2 mV for the
  2-compartment 3x point), and the dendrotoxin test is pinned to 0.1 mV
  instead of bounded (#56).

### Changed

- `scripts/coincidence_window.py` and `scripts/dendrites.py` search each
  delay-0 threshold once and share it between the windows (`threshold0=`)
  and the curves: 18 searches become 6, and 12 become 5, for 130 and 82
  fewer simulations per full run. `coincidence_window.py` and
  `scripts/making_threshold_graphs.py` queue all their curve points and
  sweeps on the shared pool at once instead of draining it between batches.
  Every JSON output is byte-identical; the wall-clock gain is within noise,
  since the windows now wait for the delay-0 job (#43).
- `synaptic.synaptic()` caches its table on `(t_end, freq, gE, random_in)`,
  the fields that fix the waveform, in a 16-entry cache, and returns
  read-only arrays (#43).

### Docs

- README.md is a short landing page: what the repo is, quick start, the
  scripts, and a table of the historical `.m` files with their Python
  replacements. It replaces the 2019 to-do fragments and the unused
  Zilany/Bruce/Carney note (#44).
- Stale descriptions corrected: the `coincidence` module docstring
  (predictive search, not bisection), the "simulations per point" counts in
  `making_threshold_graphs.py` (about 8 for EPSGpair, about 14 for the other
  sweeps) and PYTHON.md (about 10 per threshold on average, at most 13).
  MATLAB-era claims in `constants.py`, `multi.py` and `_bisect.py` now say
  that the golden traces, not bit-identity with MATLAB, pin the results
  (#44).
- PYTHON.md prose is wrapped to about 80 columns, with no wording changes
  (#44).

## 1.3.0 (23 Oct 2026)

Results are unchanged: the golden traces, thresholds, windows and the
somatic table keep their numbers. `spike_amplitude` gains a flag, the
`somatic` and BinarySearch.m searches gain model keywords, three bad inputs
fail early, and fast mode ships as wheels.

### Added

- `model_kw` on `somatic.fires`, `response`, `rheobase` and
  `spike_amplitude`, and on `threshold.binary_search`, passed to the models
  per call the way `coincidence.threshold` does: `mode`, `morph`,
  `model_type`, `max_step` and, for "two", `r1`/`tau_est`. `mem=` stays as
  a shorthand for `model_kw={"mem": ...}` and raises if given both ways.
  Goldwyn's rheobase (2070 pA) is now reproducible through the public API,
  `somatic.rheobase("two", v0=-58, model_kw={"r1": 8.5, "tau_est": 0.34})`
  (#40, #29).
- Prebuilt `msoaxon-fast` wheels (Linux x86_64 and aarch64, macOS arm64 and
  x86_64, Windows x86_64; one abi3 wheel per platform covers Python
  3.11-3.13) are built, smoke-tested and attached to each GitHub release by
  `.github/workflows/wheels.yml`, so fast mode no longer needs cargo;
  PYTHON.md says how to install one. The 0.1.0 wheels are also on the
  1.2.1 release, whose crate is identical (#39).

### Fixed

- `somatic.spike_amplitude` no longer returns a number when the soma has
  no distinct inflection: about 0 mV at 10x rheobase, where the axon crosses
  the criterion before the inflection window opens, and inflated 9-12 mV at
  1.05x and 1.2x, a smooth hump. The result carries a `distinct` flag; when
  it is False, `amplitude`, `t_inflection` and `v_inflection` are None while
  `peak_above_rest`, `t_spike` and `trace` stay valid. The inflection counts
  as distinct when the d2V/dt2 peak in the search window is the largest
  before the spike peak, is not at a window edge, and is at least
  `somatic.D2_MIN` = 100 mV/ms^2. PYTHON.md states the valid range, about
  1.35x to 4.5x rheobase for the multi and two models; the published 1.5x,
  2x and 3x values are unchanged (#38).
- `two_cpt` and `mso_axon` reject a `ramp` or `ramp2` with `stop <= start`
  with `ValueError("ramp needs stop > start")` before any mode dispatch;
  `stop == start` used to raise `ZeroDivisionError` and `stop < start` gave
  a negative slope silently (#37).
- `coincidence.window()` returns `(nan, th0)` when the threshold is not
  finite (nothing fires up to the ceiling, e.g. a fixed 30 mV criterion
  with Goldwyn's calibration) instead of failing inside scipy with
  `ValueError: array must not contain infs or NaNs` (#37).

### Changed

- `threshold.sweep_setting("Synaptic", ...)`, and so
  `binary_search("Synaptic", ...)`, raise "not supported" like
  `SynapticPair`; MATLAB's pinned `I = 0` only ever reported the ceiling
  (#37).
- Fast mode has no ramp exception any more: `_fast.covers()` is gone and
  every configuration `multi._rhs` handles runs in `mode="fast"` (#37).

### Tests

- `tests/test_threshold_search.py` pins `threshold._search`'s probe
  sequence, returned value and run count with a fake predicate, including
  the never-fires case (#46).
- `tests/conftest.py` adds a `fresh_quiet_cache` fixture; the quiet-prefix
  cache tests use it instead of clearing `_solve._QUIET_CACHE` by hand
  (#46).
- pytest runs with `--strict-markers`, so a misspelt marker fails instead
  of silently running in the `-m "not slow"` jobs (#46).
- The dendrotoxin spike-amplitude test bounds the 1.5x and 3x amplitudes
  (20-40 mV) instead of pinning them: Linux gives 25.97 mV at 1.5x where
  macOS gives 25.08, so that membrane's inflection pick is platform
  sensitive at about 1 mV; the PYTHON.md dendrotoxin row (25.1, 31.9 mV)
  is from macOS (#38).

## 1.2.1 (13 Oct 2026)

Tests and housekeeping only; results are unchanged.

### Tests

- A slow regression test pins all 28 EPSG-pair thresholds of the jpg grid
  (`making_threshold_graphs.py --only EPSGpair --jpg-grid`: delays 0-0.6
  ms, nodes 3 and 5, both models) exactly, as multiples of q = 150/128,
  guarding `threshold.binary_search`'s halving output. It replaces the test
  that pinned two of them. The table also records each point's value in
  the repo's MATLAB jpgs: 20 of 28 match and the other 8 are one output
  (2q) away, as history rather than a requirement. The script's EPSG-pair
  grid moved into `epsgpair_grid()` so the test runs the same grid.

### Notes

- 1.2.0's "13 runs instead of 11" for the 45-compartment, model-EPSG 3%
  window counts the delay search's probes; counting `window()`'s
  `max_delay` probe as well gives 14. #41 tracks whether the hedge should
  apply to that search.
- `msoAxon.m~`, a tracked editor backup that differs from `msoAxon.m`, is
  kept as part of the historical record (an earlier draft of the 45-CPT
  model); README.md now says so, and `*~` is ignored so no new backups get
  committed (#45).

## 1.2.0 (13 Oct 2026)

### Changed

- The predictive (B4) search no longer bisects right after a prediction
  overshoots: it probes a quarter of the bracket back from the overshoot
  (`_bisect.HEDGE`), which usually lands just short of the threshold and
  gives the next prediction a close point. The bilateral-dendrite threshold
  at 0.3 ms delay takes 13 runs instead of 17 (bisection: 15); the twelve
  benchmark thresholds take 125 instead of 145, and none takes more. The
  twelve windows of coincidence_window.py take 113 instead of 120, though
  one (45-compartment, model EPSG, 3%) takes 13 instead of 11.
  Thresholds and windows move within their search tolerances (benchmark
  thresholds by at most 3.2e-5 relative, rel_tol 1e-4; the twelve windows
  of coincidence_window.py by at most 0.08 us, delay_tol 0.1 us); the
  golden traces and single runs are unchanged (#18).
- The settled stop computes only the voltage derivatives it checks
  (`multi._dv`, `two._dv`), not the whole right-hand side. In exact mode a
  check costs about 20 µs instead of 81 µs (45-CPT) and 4 µs instead of
  13 µs (2-CPT). The settled checks fall from 11% to 3% of the 45-CPT
  BinarySearch.m ramp point (factor 30), which runs about 7% faster; the
  2-CPT ramp point gains about 3% and an EPSG-pair threshold search about
  1%. Event values, stop times, trajectories, thresholds and run counts are
  identical. `_rhs` itself is unchanged. Fast mode keeps slicing the Rust
  right-hand side: all its settled checks cost under 2% of that ramp point,
  too little for a Rust dV-only entry point to save a measurable amount.
- `somatic.rheobase` uses the predictive (B4) search
  (`_bisect.smallest_crossing`) instead of yes/no bisection. Near rheobase
  the multi and two step responses rise smoothly through the 10 mV
  criterion (9.96 mV at 0.9999x, 10.02 at 1.00005x; the spike takes off
  1-2% above), as the EPSG pairs do: with the overshoot hedge above,
  rheobase takes 10 runs instead of 15 for both multi and two. With KLT blocked the response is all
  or none (0.3 mV just below, a 15 mV spike just above); there the search
  makes exactly bisection's 19 runs, to the same value. Rheobase moves
  within rel_tol (multi 3784.180 -> 3784.101 pA, -2.1e-5; two 3845.947 ->
  3845.990 pA, +1.1e-5; KLT blocked unchanged), and rheobase now depends on
  the search, so this is a minor version (#19).

### Added

- `somatic.response()`: a step run's (spiked, peak axon - soma).
- `descent=` on `_bisect.smallest_crossing`: the first step down when the
  guess spikes (default 1.25, unchanged; rheobase uses 2) (#19).

## 1.1.0 (12 Oct 2026)

### Added

- `mode="fast"`: the 45-compartment right-hand side (`multi._rhs`) in Rust,
  as the optional `fast` extra (`uv sync --extra fast`, needs cargo; the
  `msoaxon-fast` package in `rust/`). scipy's BDF solver is unchanged. Single
  runs, thresholds and windows are about 2.7 times faster on an M1 Max. The
  port covers every stimulus, the chain and dendrite morphologies,
  `membrane()` overrides and every model type (#30).
- `mode=` on `mso_axon` and `_dispatch.run_model` (so also in `model_kw`):
  `"exact"` runs the numpy reference, `"fast"` the Rust port. The default is
  `"fast"` when the extension is installed, else `"exact"`; the
  `MSOAXON_MODE` environment variable overrides it, in pool workers too
  (#30).
- `msoaxon.DriftWarning`, issued once per process at the first fast run.
  The run self-checks the port against numpy on 24 states and says whether
  it was bit-identical on this machine: it is on darwin-arm64, and on other
  hardware it may differ in the last bits (#30).
- CI: a `pytest-fast` job on Linux and macOS lints the crate (`cargo fmt`,
  `clippy -D warnings`), builds it and runs the suite in fast mode (#30).

## 1.0.0 (12 Oct 2026)

The first versioned release, and the first in which the Python package,
not the MATLAB code, is the reference.

### Decisions

- Python is the reference. The `.m` files stay as the historical original,
  and MATLAB compatibility is no longer a requirement (#17, #25).
- The predictive (B4) search's thresholds and windows are the reference
  values. They differ from MATLAB's bisection by less than 0.4 µs and
  1e-4 relative (#17).
- The settled stop stays on (#20). `fitExp.m` will not be ported (#22).
- Results are hardware agnostic: tests compare against relative
  tolerances, and the golden traces use one tier on every platform (#26).

### Added

- `criterion=` on the coincidence and somatic searches. A float is an
  absolute criterion in mV (default 10.0), and `coincidence.Relative(f)` is
  a fraction of each configuration's own spike height, from the new cached
  `coincidence.spike_height()`. PYTHON.md explains why 10 mV stays the
  default, with the Goldwyn caveat (#21, #29).
- `_bisect.FiresWithoutInputError`: a search on a cell that fires with no
  input raises, naming the configuration and suggesting
  `rebalance_rest=True`. Before, the threshold search looped for ever and
  the rheobase search returned about 0 (#27).
- `msoaxon.__version__`.
- CI: `ruff check` is enforced (#25), and a macOS (arm64) job runs next to
  Linux (#26).

### Changed

- Full type annotations across `msoaxon/` and `scripts/` (#25).
- Reassigning `coincidence.FACTOR` or `somatic.FACTOR` at runtime no longer
  changes the searches; pass `criterion=` instead (#29).

## 0.1.0 (Sep 2026, untagged)

The Python port of the MATLAB model comparison (#1–#24): both models,
every stimulus, the threshold sweeps and graphs, the literature comparison,
dendrites, the somatic spike, the B4 search, and the Google-style refactor
with golden-trace tests.
