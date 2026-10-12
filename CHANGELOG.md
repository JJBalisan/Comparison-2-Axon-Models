# Changelog

Versions follow [semantic versioning](https://semver.org/). The Python
package is the reference: a release pins the numbers the analysis produces,
and a change that moves them beyond the golden-trace tolerance is a new
minor version at least.

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
