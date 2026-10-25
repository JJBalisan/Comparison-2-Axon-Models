# diffsol 0.18 probe for `_solve._segment` (Phase 0, lane A, issue #28)

Question: can diffsol 0.18 replace scipy's BDF in `_solve._segment`?
Pinned: `diffsol = "=0.18.0"` (diffsol-la 0.3.0, diffsol-nl 0.3.0, faer
0.24.4, nalgebra 0.35.0; see the committed `Cargo.lock`), rustc 1.99.1,
Apple silicon arm64, release build. Nothing outside `prototype/diffsol-probe/`
changed.

## Toy system

`src/model.rs`: 315 states = voltage block `V[0..45]` then six gate blocks.
Tridiagonal (sealed cable) V block, diagonal gate blocks, V-gate coupling per
compartment, i.e. the `multi._chain_jac_sparsity` pattern: 943 non-zeros of
99225 (0.95%; the solver's auto-detected pattern, equal to 133 + 3*270).
Axial coupling d = 100 /ms (largest cable eigenvalue about 400 /ms), gate time
constants 0.01 to 20 ms: stiffness ratio about 4e4. Gate 0 (E = +50 mV) gives
an upstroke and gates 1-2 repolarise, so V[0] has a genuine peak (23.64 mV)
inside a 0.3 ms current pulse (t = 2.0-2.3 ms, piecewise constant, set at
segment starts like msoaxon's stimulus). Two terminal roots, both negative at
rest: `V[0] - (-20)` and `V[22] - 0`. The toy's RHS costs 1.0 us per call, so
per-step numbers are mostly solver overhead, which is the point.

Driver: `src/solve.rs`, one `solve` per backend (macro `backend!`), mirroring
`_segment`: one problem, then per segment a fresh `BdfState` and `Bdf`
solver. rtol 1e-6, atol 1e-8. `cargo run --release` prints every number
below; `cargo test` asserts the behaviour (`tests/probe.rs`).

## Results

### 1. Terminal events: PASS

- API: `OdeBuilder::root(|x,p,t,g| ..., nroots)` (src/solve.rs:135);
  `solver.step()` returns `OdeSolverStopReason::RootFound(t_root, idx)`
  (:168); state at root is `solver.interpolate(t_root)` (:185). Tests:
  `item1_*`.
- The step containing the root is accepted, so `solver.state().t` is the END
  of that step (2.051037), not the root (2.049634). The state at the root must
  come from `interpolate`, not `state().y`. The solver stops when the caller
  stops calling `step()`.
- Placement: Sundials-style modified secant on the interpolant
  (nonlinear_solver/root.rs), tolerance about 100 eps. Measured
  |g(interpolate(t_root))|/|threshold| = 8e-13 (target 1e-10: met). Against a
  rtol-1e-12 reference the root time differs by 3.5e-9 relative; that is ODE
  error at rtol 1e-6, not root-finder error, the same behaviour as scipy
  `events` (root of the dense output).
- `idx` says which root function fired (0 here), so `spike_event` vs
  `settled_event` can be told apart. Direction is not filtered and not needed.
- Gotcha for a continuing (non-terminal) run: `RootFound` returns before the
  solver processes `tstop`, so a step that both crosses a root and ends at the
  stop time never reports `TstopReached`, and the next `step()` panics inside
  diffsol (`StopTimeBeforeCurrentTime`, bdf.rs:1675 `.unwrap()`). Fix used:
  call `set_stop_time` again after any `RootFound` (src/solve.rs:205-209).
  Irrelevant for terminal events.

### 2. Every accepted step as output: PASS

- API: the `step()` loop; read `solver.state().t` / `.y` after each call
  (src/solve.rs:177-181). Test: `item2_every_accepted_step_is_collected`.
- 603 accepted steps over t = 0..12 ms, all 315 states each, strictly
  increasing t. Rejected steps are not reported (same as scipy). The initial
  point is not emitted by `step()`; push it by hand. Cloning y per step is
  small next to the step cost.

### 3. Re-initialisation at a new t0 from a given state: PASS

- API (src/solve.rs:147-157): `BdfState::new_without_initialise(&problem)`,
  then `state.as_mut()`: set `*t`, write y with `set_index`, set dy with
  `problem.eqn.rhs().call_inplace`,
  `set_step_size(h0, &problem.atol, problem.rtol, &problem.eqn, 1)`; then
  `problem.bdf_solver::<LS>(state)`. The problem (and its sparsity/colouring)
  is built once and reused for every segment. Forcing that changes at
  breakpoints is a `Cell` read by the RHS closure, set at each segment start
  (the `init` closure only supplies the very first state). `set_stop_time(b)`
  lands the segment end (`TstopReached`; t can differ from b by roundoff, e.g.
  2.049999999999999, so compare with a tolerance). Tests: `item3_*`.
- Correctness: 24 segments on a quiet model reproduce a single-segment solve
  to <1e-6; a stimulus change at a breakpoint changes the dynamics.
- Cost (sparse backend, quiet model, t = 0..12, median of 5): 1 segment
  1516 us / 96 steps; 240 segments 63605 us / 2318 steps, about 260 us per
  re-init (about 9 extra order-1 restart steps at about 14 us, plus new solver
  and sparse-LU symbolic set-up). Negligible for msoaxon's few stimulus edges
  per run. Dense backend: 23 ms per re-init.

### 4. FD Jacobian + supplied sparsity + FaerSparseLU vs default dense: PASS (sparse); dense FAILS on speed

Per-step cost, 315 states, staged run t = 0..12 (3 segments, 603 steps, 84 LU
set-ups), release, median of 9:

| backend | Jv | wall (us) | us/step |
|---|---|---|---|
| FaerSparseMat + FaerSparseLU | analytic | 8282 | 13.7 |
| FaerSparseMat + FaerSparseLU | forward difference | 8364 | 13.9 |
| NalgebraMat + NalgebraLU (dense) | analytic | 228982 | 380 |
| NalgebraMat + NalgebraLU (dense) | forward difference | 229243 | 380 |
| FaerMat + FaerLU (dense) | analytic | 337099 | 559 |
| FaerMat + FaerLU (dense) | forward difference | 340797 | 565 |

- Sparse is about 28x faster per step than dense nalgebra. FD costs nothing
  measurable because the Jacobian is rebuilt only 84 times in 603 steps and
  colouring needs only a handful of Jv calls. Step counts are identical across
  all six; final states agree (`item4_backends_agree`, 1e-3 mV).
- diffsol has no built-in finite-difference Jacobian: the builder wants a
  Jacobian-vector product closure (`OdeBuilder::rhs_implicit(f, jv)`,
  src/solve.rs:124). The FD run supplies `jv_fd` (src/model.rs:149), a forward
  difference of the RHS alone (eps = 1e-7(1+max|x|), no norm of v), which is
  what `MultiRhs::eval` can provide. Two RHS calls per Jv; the closure API
  gives no hook to reuse f(x).
- Sparsity is not passed explicitly: with a sparse matrix type the builder
  detects the pattern at build time by probing Jv with NaN-tainted unit
  vectors (`find_jacobian_non_zeros`) and colours it (`JacobianColoring`); it
  found exactly the expected 943 entries. Consequences: Jv must propagate NaN
  structurally (no `if v == 0 {skip}`), and the pattern is fixed at t0/y0.
  Supplying `jacobian_sparsity` by hand would need a hand-written
  `OdeEquations` impl via `build_from_eqn`; not tried. Dendrite trees need no
  change in principle (pattern is detected, not assumed); not tested.
- Context: scipy BDF with `jac_sparsity` on the same toy (numpy RHS at
  14 us/call, same tolerances, 3 segments, 578 steps) took 79500 us,
  138 us/step, about 105 us/step of Python overhead after subtracting its RHS
  calls (script not committed). diffsol sparse is about 9.6x faster end to
  end, about 8x per step against scipy's overhead alone, with a 1 us RHS. Not
  the final speed-up: the real RHS is heavier and the issue's 10-20x is for
  the real model.

### 5. Max step size: WORKAROUND

- diffsol 0.18 has no max-step option: `OdeBuilder` has `h0` only;
  `OdeSolverOptions` has `min_timestep` and growth/shrink factors but no upper
  bound.
- Workaround A (scipy `max_step` semantics): arm `set_stop_time(t + max_step)`
  whenever the solver's next `state().h` exceeds the cap
  (src/solve.rs:200-209, `Cfg::max_step`). Test
  `item5_stop_time_hops_bound_the_step` asserts every accepted step <= 0.005
  (+1e-9). Toy cost (cap 0.005 ms, t = 0..12): 2620 steps, 41154 us vs
  603 steps, 8387 us uncapped (about 5x; the cap binds on this slow toy).
  Re-arming after every step (`hop_always`) is worse: 5322 steps, 76410 us.
  Peak V[0] 23.641520 vs 23.641797 reference.
- Workaround B (dense output only, no step cap): `interpolate` on a 0.005 grid
  after every accepted step (`Cfg::grid`, `emit_grid!`). 603 steps, 9357 us,
  peak 23.640035, 1.8e-3 mV low (uncapped step values: 3e-4 low). Test
  `item5_dense_output_grid_finds_peak` (tolerance 0.05 mV).
- msoaxon needs the cap in one place (`spike_amplitude`'s second run): use A
  there, nothing elsewhere (0.1*t_end rarely binds).

## Not settled here

- diffsol's BDF is not scipy's (different order/step policy and
  interpolant), so traces differ at solver-tolerance level. Golden tests and
  `benchmark_search.py --check` are lane D's gate; not tested here.
- Timings: one arm64 Mac, normal desktop load; repeat before quoting.
- The toy has no stimulus-edge discontinuity inside a segment, no dendrite
  tree, and a cheap RHS.
- diffsol pulls faer 0.24 and nalgebra 0.35; confirm the pyo3/maturin build
  (abi3-py311, CI targets) and wheel-size impact.

## Verdict

Lane D can proceed with diffsol 0.18. Terminal events (root placed to 8e-13
on the interpolant), every-accepted-step output and re-initialisation at a
breakpoint from a given state all work with documented APIs. The sparse path
(`FaerSparseMat` + `FaerSparseLU`, pattern auto-detected from a Jv closure, FD
Jv built on the existing RHS) costs about 14 us per step on the 315-state toy,
28x less than the default dense path and roughly 8-10x less than scipy's
per-step Python overhead; the dense path must not be used. The only missing
feature is a max-step setting; it is emulated with `set_stop_time` (about 5x
more steps on the toy where it binds, which matters only for
`spike_amplitude`'s second run) or replaced by dense output where 2e-3 mV peak
resolution is acceptable. Two pitfalls to build in: the state at a root must
come from `interpolate`, not `state()`, and a continuing run must re-arm
`set_stop_time` after a `RootFound`.
