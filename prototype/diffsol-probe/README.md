> **Frozen snapshot.** Written against `cfe0ec8` (PR #100, 25 Oct 2026) with
> diffsol 0.18.0 and rustc 1.99.1. CI does not build it and nothing promises
> it still builds. `rust/src/integrate.rs` is the source of truth;
> `src/pattern.rs` here is a verbatim copy of `core.rs`'s sparsity code at
> that commit and will drift. New measurements belong in `../diffsol-real`.

# diffsol-probe

The Phase 0 toy probe of #28: a 315-state cable with the model's sparsity
pattern, solved by diffsol 0.18 the way `_solve._segment` solves it with
scipy. `REPORT.md` has the per-item verdicts and the explicit-sparsity
(E2) findings; it is the record of why the sparse path, the explicit
pattern, `interpolate` at the root and the `set_stop_time` hop were chosen.

By hand, if it still builds (it has its own `Cargo.lock`):

```bash
cargo run --release   # prints every number in REPORT.md
cargo test            # asserts the behaviour (tests/probe.rs)
```
