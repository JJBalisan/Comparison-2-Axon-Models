# Changelog

Versions follow [semantic versioning](https://semver.org/). The Python
package is the reference: a release pins the numbers the analysis produces,
and a change that moves them beyond the golden-trace tolerance is a new
minor version at least.

## Unreleased

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
  golden traces, single runs and rheobase are unchanged (#18).
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
