# diffsol-real

Experiment E1 of #28 (PR #104): replays a captured `mso_axon` case under
diffsol 0.18 and times the integrator alone, with no pyo3. It includes
`rust/src/core.rs` by `#[path]`, so every change to `core.rs` compiles into
it; `REPORT.md` has the findings and the method.

## In CI

`test.yml`'s `pytest-fast` job, on ubuntu-24.04 only, runs here (so a change
to `core.rs` cannot silently break the crate):

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --release   # uses the committed fixtures/epsgpair.json
```

## By hand (not in CI)

`dump_cases.py` and `compare.py` monkeypatch and import `msoaxon` and
`tests/test_golden.py`; they are run by hand from the repository root and
are not exercised by CI (ruff lints them):

```bash
uv run python prototype/diffsol-real/dump_cases.py OUT_DIR [--fixture F]
uv run python prototype/diffsol-real/compare.py CASE_DIR OUT_DIR
```

Keep this crate's `Cargo.lock` diffsol pin in step with `rust/Cargo.lock`.
