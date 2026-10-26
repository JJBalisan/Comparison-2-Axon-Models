> **Frozen snapshot.** Written against `cfe0ec8` (PR #100, 25 Oct 2026) with
> diffsol 0.18.0 and rustc 1.99.1. CI does not build it and nothing promises
> it still builds; `rust/src/integrate.rs` is the source of truth and
> `src/pattern.rs` is a copy that will drift. New measurements belong in
> `../diffsol-real`.

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
  [Addendum (E2), 25 Oct 2026: the hand-written `build_from_eqn` path is now
  implemented (`src/explicit.rs`) and reproduces the auto path bit for bit;
  `if v == 0 {skip}` turns out to be safe (NaN == 0 is false), the real
  hazard is a guard NaN falls out of (`abs() > 0`, `max`/`min`, zero-
  coefficient skips), which made the solve fail; dendrite trees tested. See
  the addendum below.]
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
  [E2: a 45-compartment tree is now tested (addendum, 4c); the explicit
  sparsity path is no longer "not tried" (4a). Still open: the real
  `MultiRhs::eval` Jv's NaN behaviour was only read, not probed (4b), and LU
  fill-in was measured by time, not by counting L+U non-zeros (4c).]
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

## Addendum (E2), 25 Oct 2026: explicit sparsity, NaN-probe hazards, trees

Lane D open question Q4. Same crate, same pins; new code in
`src/pattern.rs`, `src/explicit.rs`, `src/lu.rs`, the tree and Jv variants
in `src/model.rs`, `sparse_explicit` in `src/solve.rs`; `cargo run --release
-- 4a|4b|4c` prints the numbers below, tests `item4a_*`, `item4b_*`,
`item4c_*` assert them. Timings: release, each a median of 9 solves, then the
median of 5 invocations; the machine was heavily loaded (`uptime` load
average 33 to 50 from two parallel experiments), so absolute us/step is
about 20% above item 4's 13.7. Compare within this table only.

### 4a. Explicit pattern via `build_from_eqn`: works, bit-identical

- Why hand-written: `OdeBuilder` has no method that takes a pattern, and
  `Closure`'s `sparsity`/`coloring` fields are private with
  `calculate_sparsity` (the NaN probe) as their only setter
  (diffsol `op/closure.rs:54-66`). So a supplied pattern needs an
  `OdeEquations` of one's own, about 150 lines (`src/explicit.rs`).
- Traits diffsol 0.18.0 needs for `build_from_eqn` + `bdf_solver`
  (`OdeEquationsImplicit` is a blanket impl over these):
  `Op` on the equations and each operator (`impl_op!`, explicit.rs:70-95,
  :130); `OdeEquationsRef<'a>` naming `Rhs`, `Init`, `Root` and
  `ParameterisedOp<UnitCallable<M>>` for the unused `Mass`/`Out`/`Reset`
  (:149); `OdeEquations` (`rhs`, `mass` = None, `root`, `init`, no-op
  `set_params`/`get_params`, :158); `NonLinearOp` on the RHS (:97);
  `NonLinearOpJacobian` on the RHS (:104) with `jac_mul_inplace` = the Jv,
  `jacobian_sparsity` returning the supplied pattern (:112) and
  `jacobian_inplace` filling it through a `JacobianColoring` built from the
  same pairs (:109, :51). Without that override diffsol's default
  `jacobian_inplace` makes N = 315 Jv calls per Jacobian; with it, one call
  per colour. `ConstantOp` on the initial state (:117); `NonLinearOp` on the
  root function, nout 2 (:123). Build: `OdeBuilder::<M>::new().t0().h0()
  .rtol().atol().build_from_eqn(eqn)` (src/solve.rs:263-269), `p` left
  empty (nparams 0).
- Conversion (src/pattern.rs:63): `core.rs::tree_jac_sparsity`'s
  `(rows, cols)` zipped into `(row, col)` pairs, then
  `<M::Sparsity as MatrixSparsity<M>>::try_from_indices(n, n, pairs)`
  (`diffsol::matrix::sparsity::MatrixSparsity`; for `FaerSparseMat` the type
  is faer's `SymbolicSparseColMat<usize>`). The same pairs go to
  `JacobianColoring::new(&sparsity, &pairs, ctx)`. Test
  `item4a_pattern_converts_to_the_diffsol_sparsity`: nnz 943 and
  `indices()` gives the pairs back in the same CSC order.
  `src/pattern.rs` is a verbatim copy of `chain_jac_sparsity` /
  `tree_jac_sparsity` from origin/main `rust/src/core.rs:761-804`.
- Result, staged run t = 0..12, chain, 315 states:

  | path | Jv | nnz | steps | LU set-ups | build (us) | solve (us) | us/step |
  |---|---|---|---|---|---|---|---|
  | auto-detect | analytic | 943 | 603 | 84 | 811 | 10079 | 16.7 |
  | explicit | analytic | 943 | 603 | 84 | 97 | 10190 | 16.9 |
  | auto-detect | FD | 943 | 603 | 84 | 1199 | 10292 | 17.1 |
  | explicit | FD | 943 | 603 | 84 | 93 | 10224 | 17.0 |

  Every accepted step time and state vector is bitwise equal between the two
  paths (test `item4a_explicit_pattern_matches_auto_detection`): the probe
  emits non-zeros in CSC order, the same order as `tree_jac_sparsity`, so the
  greedy colouring and every Jacobian are identical. Per-step cost is the
  same within noise; the only difference is build time: the NaN probe plus
  colouring cost 0.7 to 1.1 ms per problem build, the explicit path 0.1 ms.

### 4b. NaN-probe hazard: real, but not where expected

How detection works (diffsol `jacobian/mod.rs` `probe_columns` /
`probe_triplets`, called from `op/closure.rs:55`): one batched Jv call with
N lanes, lane j's v all zeros except NaN at j; an output entry is in the
pattern iff it comes back NaN. Since 0 * NaN = NaN, any branch-free
arithmetic propagates the probe whatever the coefficients are. Detected nnz
(explicit pattern: 943):

| Jv | at rest y0 | at "edge" y0 (all V = E_K = -80, gates exactly 0 or 1) |
|---|---|---|
| analytic, branch-free | 943 | 943 |
| forward difference of the RHS | 943 | 943 |
| `if v == 0.0 { continue }` on gate inputs (`jv_skip_zero`) | 943 | 943 |
| `if v.abs() > 0.0 { use it }` on gate inputs (`jv_guarded`) | 673 | 673 |
| skip a term when its coefficient is exactly 0 (`jv_coef_skip`) | 943 | 853 |

- The `if v == 0 {skip}` shortcut the original item 4 warned about is
  harmless: NaN == 0.0 is false, so the probe takes the "use it" branch.
  (Item 4's sentence "no `if v == 0 {skip}`" was wrong.) The real rule: the
  probe is lost by any test NaN falls out of: `> 0`, `< 0`, `abs() > 0`,
  `is_finite()`, `f64::max`/`min`/`clamp` (they return the non-NaN operand),
  and value-dependent skips of zero coefficients. `== 0.0 -> skip` and
  `!= 0.0 -> include` are safe. The same holds for the RHS when Jv is a
  forward difference of it: a state clamp like `g.max(0.0)` in the RHS would
  drop that column.
- The guarded Jv is exact for every finite v, yet the probe misses all
  6 x 45 = 270 dV/dgate entries (673 detected). Consequence, staged run from
  rest: the auto path FAILS ("Exceeded maximum number of nonlinear solver
  failures (51) at time = 0.102"); the explicit path with the very same Jv
  runs 603 steps / 84 set-ups, final state 5e-7 mV from an rtol-1e-10
  reference, identical to the analytic run. Test
  `item4b_guarded_jv_undercounts_and_breaks_the_auto_path`.
- y0 dependence: none for branch-free Jv (analytic and FD give 943 at the
  edge state with gates exactly 0 and 1 and V exactly at E_K). With a
  zero-coefficient skip it appears: at V = E_K the dV/dg coefficient
  gbar (E_K - V) is exactly 0 for the two gates reversing at -80 mV, so
  2 x 45 = 90 entries vanish (853), and because the pattern is fixed at t0
  they stay missing after V leaves E_K. Staged run from the edge state: auto
  path FAILS at t = 2.189 (inside the stimulus pulse); explicit path 707
  steps / 91 set-ups, 9.5e-7 mV from the reference. Test
  `item4b_detection_depends_on_y0_only_through_branches`. In this toy gates
  enter dV linearly, so a gate at exactly 0 or 1 cannot remove an entry; in
  the real model dV/dh carries m^4 (`core.rs` `eval`), so `m = 0` at t0 would
  remove dV/dh under such a skip.
- The opposite failure also exists, by reasoning rather than test: an FD Jv
  that scales its step by a norm of v turns the whole output NaN, so the
  probe reports a dense pattern (99225 entries) and the sparse path silently
  becomes the 28x slower dense one. `jv_fd` avoids it (no norm of v).
- The real RHS, read not probed: `core.rs::eval` (origin/main :653-758)
  branches only on `active` flags and on t (`Epsg::unitary`, stimulus
  windows), and uses `powf`, which propagates NaN, so an FD Jv over it should
  detect correctly. Not run here.

Conclusion: lane D should pass the pattern explicitly. It is
`core.rs::tree_jac_sparsity`, the same single source of truth as Python's
`multi._tree_jac_sparsity`/`_chain_jac_sparsity`; it costs one ~150-line
`OdeEquations` impl, gives bit-identical results where detection is right,
removes the dependence on how Jv treats NaN and on y0, and saves ~1 ms per
problem build (one build per `_segment` call; re-initialisations reuse it).
Keep a debug assertion that the detected pattern is a subset of the
explicit one (what `detected_pattern` in src/solve.rs computes) as a cheap
guard against a Jv that touches entries the pattern lacks.

### 4c. Dendrite tree: same pattern size, same cost, no fill-in trouble

`Model::tree()` (src/model.rs:136): 45 compartments, a 25-compartment trunk
`0..25` plus two 10-compartment branches `25..35` and `35..45` both off
compartment 5, built from an explicit parent/child edge list; the axial term
is d * sum over neighbours (x_j - x_i) (`Model::axial`), and the pattern is
`tree_jac_sparsity(45, parent, child)`. The chain keeps its original
arithmetic, so items 1-5 are unchanged (603 steps, 84 set-ups, as before).

- Pattern: 943 entries for the tree too (a tree on n nodes has n - 1 edges,
  so 45 + 2 x 44 + 18 x 45 whatever the branching). Auto-detected = explicit,
  entry for entry, with analytic and FD Jv; it differs from the chain's
  (column 5 holds rows 4, 5, 6, 25, 35). Test
  `item4c_tree_patterns_agree_and_the_solve_completes` (both paths complete
  t = 0..12 with equal steps and final state).
- Cost, staged run t = 0..12 (the tree takes fewer steps: different
  dynamics, the stimulus charges more membrane):

  | model | path | Jv | steps | set-ups | build (us) | solve (us) | us/step |
  |---|---|---|---|---|---|---|---|
  | chain | auto | analytic | 603 | 84 | 810 | 9956 | 16.5 |
  | chain | explicit | analytic | 603 | 84 | 97 | 9984 | 16.6 |
  | chain | auto | FD | 603 | 84 | 1204 | 10117 | 16.8 |
  | chain | explicit | FD | 603 | 84 | 93 | 10106 | 16.8 |
  | tree | auto | analytic | 504 | 75 | 834 | 8137 | 16.1 |
  | tree | explicit | analytic | 504 | 75 | 100 | 8160 | 16.2 |
  | tree | auto | FD | 504 | 75 | 1245 | 8324 | 16.5 |
  | tree | explicit | FD | 504 | 75 | 98 | 8272 | 16.4 |

- `FaerSparseLU` alone on the BDF matrix `J - I/c` at rest (src/lu.rs,
  median of 201): symbolic set-up (`set_problem`: COLAMD ordering,
  elimination tree, allocation) 37.1 us chain / 36.2 us tree; numeric
  factorisation (`set_linearisation`) 24.0 / 23.9 us at c = 1e-3 and
  29.5 / 29.8 us at c = 0.1; one solve 4.8 / 4.7 us. No sign of fill-in
  trouble: faer orders columns with COLAMD (`factorize_symbolic_lu`), and a
  tree's elimination graph has no cycles to fill. L+U non-zeros were not
  counted (faer's numeric LU keeps them private); time is the measure here.

### 4d. Pattern fixed at t0: confirmed

- Detection happens once, in `OdeBuilder::build` (diffsol
  `ode_solver/builder.rs:1828` evaluates y0 = init(t0); `:1861-1862` calls
  `rhs.calculate_sparsity(&y0, self.t0, &p)` only when `M::is_sparse()` or
  `use_coloring`; that is `op/closure.rs:54-66`, probe then
  `JacobianColoring::new`). `build_from_eqn` (`:1953`) never probes.
- It is read, not recomputed: `BdfCallable::new` (`op/bdf.rs:126-150`)
  takes `rhs().jacobian_sparsity()` and unions it with the identity for the
  iteration matrix; `Bdf::new` (`ode_solver/bdf.rs:331-333`) then runs the
  linear solver's `set_problem`, the symbolic LU (37 us above) that every
  segment's re-initialisation repeats; that is not a re-detection. Every
  Jacobian update fills the same pattern through the same colouring.
- For `active` flags: an inactive gate's derivative is the constant 0
  (`core.rs::eval`, `if act_x {...} else {0.0}`), so its row has true
  structural zeros. The explicit pattern stores them as explicit zeros:
  harmless, a few extra stored entries, and the identity union keeps the
  diagonal of I - cJ non-zero so the LU is fine. The auto path would instead
  drop that gate's row (no NaN passes a constant), so the two paths' nnz
  differ on the real model even with a clean Jv: harmless too, but one more
  reason to use one pattern. The hazard is only values that are zero by
  coincidence at t0 behind a branch (4b): structural zeros never hurt, and
  a branch-free Jv cannot produce coincidental zeros under the NaN probe.
