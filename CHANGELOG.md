# Changelog

Versions follow [semantic versioning](https://semver.org/). The Python
package is the reference: a release pins the numbers the analysis produces,
and a change that moves them beyond the golden-trace tolerance is a new
minor version at least.

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
