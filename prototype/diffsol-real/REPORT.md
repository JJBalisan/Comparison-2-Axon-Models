# diffsol 0.18 on the real model (lane D, experiment E1, issue #28)

**Setup.**
- Pinned: diffsol =0.18.0 (Cargo.lock committed, faer 0.24.4), rustc 1.99.1,
  Apple M1 Max, release build, numpy 2.5.3, scipy 1.18.1.
- The real model is `rust/src/core.rs` included verbatim
  (`#[path = "../../../rust/src/core.rs"]`), with no pyo3.
- The driver is the probe's `src/solve.rs`, adapted:
  - FaerSparseMat with FaerSparseLU.
  - Forward-difference Jv built on two `MultiRhs::eval` calls.
  - Sparsity pattern auto-detected.
  - A new `BdfState` for each segment, and every accepted step is output.
  - `max_step` emulated with lazy `set_stop_time` hops.
  - rtol = atol = 1e-8, as in `mso_axon`.
  - The quiet prefix [0, cuts[0]] is integrated with the "none" right-hand
    side, as `_solve.integrate` does.

**How the cases were captured.**
- `dump_cases.py` wraps `_fast._ext.MultiRhs`, `_solve.integrate` and
  `_common._stops`, and runs `mso_axon(mode="fast")`. Each case JSON
  therefore holds exactly what reaches the solver, plus scipy's result and
  wall time.
- It sets `_fast._STATE["warned"]` first, so the self-check's three extra
  MultiRhs objects are not captured.
- Scipy's headline counts are reproduced: EPSGpair = 1292 steps, 110 ms
  somatic step = 1576.

**Load and timing method.**
- The consolidated timing pass ran 21:01-21:03 GMT. Load averages fell from
  13.5 to 9.2 (1-minute) over that pass. Earlier passes at load 25-85 gave
  identical step counts and margins, with noisier wall times.
- Python walls are the median of 3 with the quiet-prefix cache cleared
  before each run (cold, as the harness does).
- Rust walls are the median of 5. They include the problem build (sparsity
  detection and colouring, 2.9 ms) and the quiet prefix.

## Per-case table (Q1, Q2, Q8)

| case | scipy steps | diffsol steps | diffsol wall ms | python fast wall ms | ratio | margin | scipy margin | t_last exact (raw) |
|---|---|---|---|---|---|---|---|---|
| multi/EPSG-subthreshold | 365 | 412 | 12.3 | 50.6 | 4.1 | 0.126 | 0 | True (True) |
| multi/EPSG | 1252 | 1316 | 33.0 | 200.6 | 6.1 | 0.0614 | 0 | True (True) |
| multi/EPSGbilateral | 1320 | 1389 | 33.9 | 208.6 | 6.2 | 0.0478 | 0 | True (True) |
| multi/EPSGpair-input-node-3 | 1306 | 1394 | 34.0 | 209.3 | 6.2 | 0.0506 | 0 | True (True) |
| multi/EPSGpair-subthreshold | 447 | 484 | 13.5 | 64.2 | 4.8 | 0.0479 | 0 | True (True) |
| multi/EPSGpair | 1292 | 1363 | 33.6 | 208.5 | 6.2 | 0.166 | 0 | True (True) |
| multi/Synaptic | 670 | 673 | 19.2 | 105.2 | 5.5 | 0.172 | 0 | True (True) |
| multi/SynapticPair | 766 | 781 | 21.7 | 121.1 | 5.6 | 0.0505 | 0 | True (True) |
| multi/ramp | 1304 | 1389 | 34.4 | 215.9 | 6.3 | 0.242 | 0 | True (True) |
| multi/ramp2 | 1286 | 1317 | 33.6 | 215.4 | 6.4 | 0.0845 | 0 | True (True) |
| multi/sine | 1314 | 1330 | 33.6 | 209.6 | 6.2 | 0.313 | 0 | True (True) |
| multi/step-input-node-3 | 1507 | 1597 | 38.8 | 239.7 | 6.2 | 0.146 | 0 | True (True) |
| multi/step-subthreshold | 463 | 514 | 14.4 | 66.3 | 4.6 | 0.0122 | 0 | True (True) |
| multi/step | 1507 | 1597 | 38.7 | 241.7 | 6.2 | 0.146 | 0 | True (True) |
| multi/type/active-full | 1234 | 1298 | 32.3 | 200.2 | 6.2 | 0.0614 | 0 | True (True) |

"Scipy margin" is fast-mode scipy against the golden file. It is 0, so
diffsol's margin against golden equals its margin against scipy.

**FD epsilon.** Two step choices were tried:
- The probe's: 1e-7·(1+max|x|).
- Standard: sqrt(eps)·max(1, max|x_j| over the probed entries).

Both give max margin 0.313 and identical margins to 3 digits, the same step
counts in 13 of 15 cases (Synaptic 673 vs 682, SynapticPair 781 vs 779), and
0 Newton failures. The probe's choice is kept.

**Sparsity.**
- NaN probing detected 943 non-zeros, equal to `core::chain_jac_sparsity(45)`.
  `eval` propagates NaN structurally for active-full.
- Caveat 1: the pattern is detected once, at (t0, y0). Inactive gates write an
  exact 0.0, so other model types detect fewer entries. That is correct, but
  the count will not be 943.
- Caveat 2: a stimulus coupling outside the chain pattern would be missed.
  `Synaptic` with k≠0 puts V[0] into dV[k]; all 15 cases have k=0.
- An explicit pattern (E2) avoids both caveats and also saves the 2.9 ms
  (630 RHS calls) that detection costs on every run.

**Side result (not a recommendation).** At rtol = atol = 1e-7, diffsol still
passes all 15 cases (max margin 0.666, 1.3x faster). At 1e-6 it fails 9 of 15
(max 2.47).

## Q1 Golden tolerance

Yes, all 15 cases pass.
- At most 31% of the margin is consumed: sine 0.313, ramp 0.242, Synaptic
  0.172, EPSGpair 0.166, subthreshold cases 0.01-0.13. That leaves a factor
  of 3 to the edge.
- t_last equals t_end exactly in 15 of 15 cases. The raw `state().t` at the
  final `TstopReached` was already exact; the driver also snaps to the
  breakpoint, but that never had to act.

## Q2 Step count

diffsol takes 0.4-13% more accepted steps than scipy:
- Spiking cases: 1-7% more. EPSGpair is 1363 vs 1292; the 110 ms somatic step
  is 1655 vs 1576.
- Subthreshold cases: 8-13% more. EPSG-subthreshold is 412 vs 365.

The default cap (0.1·t_end) never binds: runs without it give identical
counts. On EPSGpair there were 62 error-test failures, and there were 0
Newton failures anywhere.

## Q3 Real RHS cost

- `MultiRhs::eval` costs 4.09 µs per call: 10k calls on y0, median of 9
  batches, six repeats between 4.08 and 4.09 µs, the same with the stimulus
  on (t=5.1) or off (t=1.0).
- On EPSGpair there are 3057 RHS calls in total:
  - 630 for sparsity detection (315 NaN Jv probes);
  - 2427 over 1362 accepted steps, which is 1.78 calls per step including the
    FD Jv calls. During integration there were 261 Jv calls and 182 LU setups.
- Per step: (33.2 − 2.9 ms build) / 1362 = 22.2 µs. About 7.3 µs of that is
  the RHS and about 15 µs is solver overhead.
- The toy measured 14 µs/step total with a 1 µs RHS, so solver overhead is
  about the same on the real model.
- Fast-mode scipy costs 161 µs per step (208.5 ms / 1292), so diffsol is 7x
  cheaper per step.

## Q5 max_step

110 ms somatic step (5-105 ms, 8500 pA, t_end 110):

| run | steps | wall | LU setups |
|---|---|---|---|
| scipy, max_step 0.005 | 22831 | 2622 ms | n/a |
| diffsol, set_stop_time hops | 22906 | 624 ms (27 µs/step) | 13008 |
| scipy, default max_step (11 ms) | 1576 | 253 ms | n/a |
| diffsol, default (cap never binds) | 1655 | 40 ms | 228 |

- Hop emulation works: the largest step is 0.005 + 3e-15 ms, and the final V
  is within 4e-9 mV of scipy's.
- It is costly, though: 57% of steps trigger an LU setup. Each hop shortens
  the step the controller wanted, which rescales the BDF history and drops
  the factorisation. So the capped run is 4.2x faster than scipy rather than
  6x.
- **It is acceptable for spike_amplitude.** That function's real second run
  goes only from 0 to ts+2 ms (7.24 ms here), not 110 ms.
- `amplitude.py` ran `spike_amplitude("multi", 8500)` with its second run on
  diffsol:
  - Scipy: 2219 steps, 307 ms.
  - diffsol with hops: 2276 steps, 55 ms. Amplitude 14.752766 vs 14.752761 mV
    (5e-6 mV apart), the same t_inflection (5.271 ms), distinct=True.
  - diffsol with no cap at all: 1243 steps, 31 ms, amplitude 14.752767 mV
    with the same inflection.
- So diffsol's own steps at 1e-8 are dense enough for the cubic spline in
  this case. Keeping the cap for parity with scipy costs little.

## Q6 Root state

`spike_event` is g = x[2] − x[0] − 10: node 3 (column 2) and FACTOR = 10 mV,
the threshold search's criterion. Tested on the golden EPSGpair run:
- Scipy stops at t = 5.338068594 ms (269 steps); diffsol at 5.338068543 ms
  (292 steps). That is 5.0e-8 ms apart.
- The state at the event agrees to 1.3e-5 mV in V and 4e-7 in the gates. |g|
  at diffsol's root is 2.2e-11.
- The root falls inside the last accepted step (step end 5.339816 ms). The
  state therefore comes from `interpolate(t_root)`, not `state().y`, and the
  trace ends at the root, as solve_ivp's does.
- The settled event is also implemented: piecewise −1/+1, with the dV check
  done through `eval`. The same root finder handles it.

**Replay of the threshold search.** All 13 runs of
`coincidence.threshold("multi", 0.1)` were replayed in diffsol:
- The spike/settled decision agreed in 13 of 13 runs.
- Stop times agree within 2e-5 ms, except the run at the final amplitude
  (1.5e-4 ms apart, a near-critical trajectory).
- One difference: in the zero-input check (amp 0), scipy ran to t_end while
  diffsol's settled event fired at t_off = 8.7. At rest the sign of every dV
  is round-off. The answer (no spike) is unchanged and the run stops earlier.

**Pitfalls from the probe.** The re-arm-stop-time pitfall does not arise
here: both events are terminal and the driver ends the run on `RootFound`. A
continuing run would still need to re-arm.

## Q8 Overall gain

**Per full run** (cold, everything included):
- 6.1-6.4x faster on spiking cases.
- 4.1-5.6x on subthreshold and synaptic cases. On short runs the fixed 2.9 ms
  build and 3.4 ms quiet prefix weigh more.

**Re-initialisation cost per breakpoint on the real model:**
- 638 µs, measured on EPSGpair with 50 extra breakpoints in 10-20 ms:
  +31.9 ms, +1083 steps, +470 LU setups.
- 515 µs with 200 extra breakpoints.
- About 21 extra restart steps each, against the probe's 260 µs on the toy.
  With msoaxon's two breakpoints per run that is about 1.3 ms of a 33 ms run.

**Threshold search** (13 runs, threshold 49.43):
- Python fast mode: 1.229 s best of 3 at load ~10 (0.966 s in the quiet
  baseline).
- diffsol replay, with every run rebuilding the problem and integrating the
  quiet prefix: 249 ms total, 4.9x faster.
- Caching the quiet prefix as Python already does (12 of 13 runs skip it):
  about 208 ms, 5.9x. This is the like-for-like comparison.
- Also hoisting the problem build (or using an explicit pattern, E2): about
  173 ms, 7.1x.
- Scaled to the 0.966 s baseline: about 0.16-0.20 s, or about 0.14 s with
  both savings. That excludes pyo3 call and output-array overhead, which
  should be small (about 1300×315 floats per run).
- This is below the issue's hoped-for 10-20x. diffsol's per-step solver
  overhead (about 15 µs) is now twice the RHS cost, and it takes about 5% more
  steps.

## Not settled here

- The machine was loaded. The ratios were stable across passes, but the
  absolute times should not be quoted.
- Only lumped active-full was run. Dendrite trees, other model types and
  Events were not tested.
- The search replay compared only the spike/settled decision.
  `_bisect.smallest_crossing` also interpolates on the peak axon−soma of
  non-spiking runs, and those traces were not stored. So the threshold value
  49.43 itself was not reproduced end to end on diffsol: unchecked.
- The search logic outside the integrator was not timed separately.

## Verdict

diffsol 0.18 on the real model:
- Passes the golden gate with a factor-3 margin to spare and lands exactly
  on t_end.
- Finds spike and settled stops that agree with scipy's.
- Needs no change to core.rs.

It is about 6x faster per full run and about 5-7x on the threshold search
(0.97 s → about 0.14-0.20 s). Caching the problem build and the quiet prefix
accounts for about 1.4x of that. Hop emulation of max_step is acceptable for
spike_amplitude: the same amplitude to 5e-6 mV, 5.6x faster than scipy on
that run.

**Judgement calls.** The diffsol Jv is forward difference only; no analytic
Jacobian was added. Initial step h0 = 1e-6 in every segment, as the probe
used. The settled event in Rust is the piecewise ±1 function solved by
diffsol's root finder, not a port of solve_ivp's event handling.
