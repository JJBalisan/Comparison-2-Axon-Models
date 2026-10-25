//! diffsol 0.18 BDF on msoaxon's real 45-compartment right-hand side.
//!
//! Lane D experiment E1 (issue #28). The model is `rust/src/core.rs`,
//! included verbatim (no pyo3). A [`Case`] is one `mso_axon` run as
//! `dump_cases.py` captured it: the constructor arguments of the stimulus
//! and the stimulus-free ("none") right-hand sides, y0, the breakpoints,
//! t_end, the tolerances, max_step and the stop events. [`solve`] mirrors
//! `_solve.integrate`: one BDF solver per segment, re-initialised from the
//! previous segment's last state at every breakpoint, every accepted step
//! recorded. The driver is the probe's (`origin/phase0-diffsol-probe`,
//! `prototype/diffsol-probe/src/solve.rs`), adapted.

#[path = "../../../rust/src/core.rs"]
pub mod core;

use std::cell::Cell;
use std::collections::HashMap;
use std::time::{Duration, Instant};

use diffsol::{
    FaerSparseLU, FaerSparseMat, NonLinearOp, NonLinearOpJacobian, OdeBuilder, OdeEquations,
    OdeSolverMethod, OdeSolverState, OdeSolverStopReason, Vector,
};
use serde::{Deserialize, Serialize};

type M = FaerSparseMat<f64>;
type LS = FaerSparseLU<f64>;

/// The keyword arguments of one `msoaxon_fast.MultiRhs(...)` call.
#[derive(Clone, Debug, Deserialize)]
pub struct RhsArgs {
    pub g_na: Vec<f64>,
    pub g_kht: Vec<f64>,
    pub g_klt: Vec<f64>,
    pub g_h: Vec<f64>,
    pub g_lk: Vec<f64>,
    pub vlk: Vec<f64>,
    pub na_vhalf: Vec<f64>,
    pub cap: Vec<f64>,
    pub sa: Vec<f64>,
    pub g_ax: Vec<f64>,
    pub parent: Vec<i64>,
    pub child: Vec<i64>,
    pub active: [bool; 6],
    pub gating: HashMap<String, f64>,
    pub stim_type: String,
    #[serde(default)]
    pub start: f64,
    #[serde(default)]
    pub stop: f64,
    #[serde(default)]
    pub amp: f64,
    #[serde(default)]
    pub k: i64,
    #[serde(default)]
    pub k2: i64,
    #[serde(default)]
    pub f: f64,
    pub epsg: Option<(f64, f64, f64)>,
    pub t_syn: Option<Vec<f64>>,
    pub g_syn: Option<Vec<f64>>,
    pub vsyn: Option<f64>,
    pub diff: Option<f64>,
}

impl RhsArgs {
    /// Builds the core right-hand side, as `lib.rs`'s binding does.
    pub fn build(&self) -> Result<core::MultiRhs, core::CoreError> {
        let params = core::Params {
            g_na: self.g_na.clone(),
            g_kht: self.g_kht.clone(),
            g_klt: self.g_klt.clone(),
            g_h: self.g_h.clone(),
            g_lk: self.g_lk.clone(),
            vlk: self.vlk.clone(),
            na_vhalf: self.na_vhalf.clone(),
            cap: self.cap.clone(),
            sa: self.sa.clone(),
            g_ax: self.g_ax.clone(),
            parent: self.parent.clone(),
            child: self.child.clone(),
            active: self.active,
            gating: core::Gating::from_map(&self.gating)?,
        };
        let spec = core::StimSpec {
            stim_type: self.stim_type.clone(),
            start: self.start,
            stop: self.stop,
            amp: self.amp,
            k: self.k,
            k2: self.k2,
            f: self.f,
            epsg: self.epsg.unwrap_or((1.0, 1.0, 1.0)),
            t_syn: self.t_syn.clone(),
            g_syn: self.g_syn.clone(),
            vsyn: self.vsyn.unwrap_or(0.0),
            diff: self.diff,
        };
        core::MultiRhs::new(params, &spec)
    }
}

/// The terminal events of a run (`_common._stops`).
#[derive(Clone, Debug, Deserialize)]
pub struct Stops {
    pub axon_col: usize,
    pub factor: f64,
    pub settled: bool,
    pub t_off: f64,
    pub settled_fraction: f64,
}

/// One `mso_axon` run as `dump_cases.py` wrote it.
#[derive(Clone, Debug, Deserialize)]
pub struct Case {
    pub name: String,
    pub rhs: RhsArgs,
    pub quiet_rhs: Option<RhsArgs>,
    pub y0: Vec<f64>,
    pub t_end: f64,
    pub cuts: Vec<f64>,
    pub rtol: f64,
    pub atol: f64,
    pub max_step: f64,
    pub use_quiet: bool,
    pub stops: Option<Stops>,
    pub axon_col: usize,
    pub scipy: serde_json::Value,
}

/// How the forward-difference `Jv` picks its step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fd {
    /// The probe's: `1e-7 * (1 + max |x|)` over the whole state.
    Probe,
    /// `sqrt(eps) * max(1, max |x_j|)` over the probed entries `v_j != 0`.
    Scaled,
}

/// Run configuration.
#[derive(Clone, Debug)]
pub struct Cfg {
    pub fd: Fd,
    /// Step cap emulated with `set_stop_time` hops; None: no cap.
    pub max_step: Option<f64>,
    /// Further breakpoints (to measure the re-initialisation cost).
    pub extra_cuts: Vec<f64>,
    /// Apply the case's stop events.
    pub stops: bool,
    /// Initial step handed to every segment.
    pub h0: f64,
}

impl Cfg {
    /// The configuration that mirrors `mso_axon`: the case's max_step.
    pub fn for_case(case: &Case) -> Self {
        Self {
            fd: Fd::Probe,
            max_step: Some(case.max_step),
            extra_cuts: Vec::new(),
            stops: true,
            h0: 1e-6,
        }
    }
}

/// Summed solver statistics over the segments.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Stats {
    pub accepted_steps: usize,
    pub lu_setups: usize,
    pub error_test_failures: usize,
    pub newton_iterations: usize,
    pub newton_failures: usize,
}

/// A located root.
#[derive(Clone, Debug, Serialize)]
pub struct Root {
    pub t: f64,
    pub index: usize,
    /// State from `interpolate(t)`.
    pub y: Vec<f64>,
    /// End of the accepted step that contained the root.
    pub step_end: f64,
}

/// Result of one run.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Trace {
    pub name: String,
    /// Every output time: the initial point, then every accepted step (a
    /// root ends the trace at the root, like solve_ivp's terminal events).
    pub t: Vec<f64>,
    pub soma: Vec<f64>,
    pub axon: Vec<f64>,
    pub final_y: Vec<f64>,
    /// `len(t)`, comparable with scipy's.
    pub steps: usize,
    pub stats: Stats,
    /// Right-hand side calls: by the solver (`f` and inside `Jv`).
    pub evals: usize,
    pub jv_calls: usize,
    /// Of `evals`, those made while building the problem (sparsity
    /// detection and colouring).
    pub build_evals: usize,
    /// Calls made by the settled event's dV check.
    pub root_evals: usize,
    pub segments: usize,
    /// Wall time of each segment, re-initialisation included [us].
    pub segment_us: Vec<f64>,
    /// Segment ends where `state().t` was not exactly the breakpoint (the
    /// trace is snapped to the breakpoint there).
    pub inexact_segment_ends: usize,
    /// Raw `state().t` at the end of the last segment.
    pub t_last_raw: f64,
    pub root: Option<Root>,
    pub jac_nnz: usize,
    /// Problem build (sparsity detection and colouring) [us].
    pub build_us: f64,
    /// Whole run including the build [us].
    pub wall_us: f64,
}

/// The two right-hand sides of a run and the counters the closures bump.
struct Model {
    rhs: core::MultiRhs,
    quiet: Option<core::MultiRhs>,
    use_quiet: Cell<bool>,
    roots_on: Cell<bool>,
    evals: Cell<usize>,
    jv_calls: Cell<usize>,
    root_evals: Cell<usize>,
}

impl Model {
    fn f(&self, x: &[f64], t: f64, y: &mut [f64]) {
        self.evals.set(self.evals.get() + 1);
        match (&self.quiet, self.use_quiet.get()) {
            (Some(q), true) => q.eval(t, x, y),
            _ => self.rhs.eval(t, x, y),
        }
    }

    fn jv(&self, fd: Fd, x: &[f64], t: f64, v: &[f64], y: &mut [f64]) {
        self.jv_calls.set(self.jv_calls.get() + 1);
        let h = match fd {
            Fd::Probe => 1e-7 * (1.0 + x.iter().fold(0.0_f64, |a, b| a.max(b.abs()))),
            Fd::Scaled => {
                // v_j != 0 holds for the NaN probes of sparsity detection.
                let xm = x
                    .iter()
                    .zip(v)
                    .filter(|(_, b)| **b != 0.0)
                    .fold(1.0_f64, |a, (xj, _)| a.max(xj.abs()));
                f64::EPSILON.sqrt() * xm
            }
        };
        let xp: Vec<f64> = x.iter().zip(v).map(|(a, b)| a + h * b).collect();
        let mut f0 = vec![0.0; x.len()];
        self.f(x, t, &mut f0);
        self.f(&xp, t, y);
        for (a, b) in y.iter_mut().zip(&f0) {
            *a = (*a - b) / h;
        }
    }

    /// `spike_event` (g[0]) and `settled_event` (g[1]); -1 when off.
    fn roots(&self, stops: Option<&Stops>, x: &[f64], t: f64, g: &mut [f64]) {
        g[0] = -1.0;
        g[1] = -1.0;
        let Some(s) = stops else { return };
        if !self.roots_on.get() {
            return;
        }
        let d = x[s.axon_col] - x[0];
        g[0] = d - s.factor;
        if s.settled && t >= s.t_off && d < s.factor * s.settled_fraction {
            self.root_evals.set(self.root_evals.get() + 1);
            let mut dx = vec![0.0; x.len()];
            self.rhs.eval(t, x, &mut dx);
            let n = self.rhs.n();
            if dx[..n].iter().all(|v| *v < 0.0) {
                g[1] = 1.0;
            }
        }
    }
}

/// Runs `case` under `cfg` and returns its trace.
pub fn solve(case: &Case, cfg: &Cfg) -> Result<Trace, Box<dyn std::error::Error>> {
    let start = Instant::now();
    let model = Model {
        rhs: case.rhs.build()?,
        quiet: case.quiet_rhs.as_ref().map(|q| q.build()).transpose()?,
        use_quiet: Cell::new(false),
        roots_on: Cell::new(false),
        evals: Cell::new(0),
        jv_calls: Cell::new(0),
        root_evals: Cell::new(0),
    };
    let n = case.y0.len();
    let mut cuts = case.cuts.clone();
    cuts.extend(cfg.extra_cuts.iter().copied());
    cuts.sort_by(|a, b| a.total_cmp(b));
    cuts.dedup();
    let mut edges = vec![0.0];
    edges.extend(cuts.iter().copied().filter(|c| *c > 0.0 && *c < case.t_end));
    edges.push(case.t_end);
    let quiet_first = case.use_quiet && model.quiet.is_some() && edges.len() > 2;
    let stops = if cfg.stops { case.stops.as_ref() } else { None };
    let fd = cfg.fd;
    let ystart = case.y0.clone();
    let m = &model;
    model.use_quiet.set(quiet_first);
    let problem = OdeBuilder::<M>::new()
        .t0(0.0)
        .h0(cfg.h0)
        .rtol(case.rtol)
        .atol([case.atol])
        .rhs_implicit(
            |x, _p, t, y| m.f(x, t, y),
            move |x, _p, t, v, y| m.jv(fd, x, t, v, y),
        )
        .init(|_p, _t, y| y.copy_from_slice(&ystart), n)
        .root(move |x, _p, t, g| m.roots(stops, x, t, g), 2)
        .build()?;
    let mut tr = Trace {
        name: case.name.clone(),
        jac_nnz: problem
            .eqn
            .rhs()
            .jacobian_sparsity()
            .map_or(0, |s| s.compute_nnz()),
        build_us: start.elapsed().as_secs_f64() * 1e6,
        build_evals: model.evals.get(),
        ..Trace::default()
    };
    let ac = case.axon_col;
    let push = |tr: &mut Trace, t: f64, y: &[f64]| {
        tr.t.push(t);
        tr.soma.push(y[0]);
        tr.axon.push(y[ac]);
    };
    push(&mut tr, 0.0, &case.y0);
    let mut ycur = case.y0.clone();
    'segments: for (i, w) in edges.windows(2).enumerate() {
        let (a, b) = (w[0], w[1]);
        model.use_quiet.set(quiet_first && i == 0);
        model.roots_on.set(!(quiet_first && i == 0));
        let seg_start = Instant::now();
        // Re-initialise at (a, ycur): the `_segment` step.
        let mut state = diffsol::BdfState::new_without_initialise(&problem)?;
        {
            let mut s = state.as_mut();
            *s.t = a;
            for (j, v) in ycur.iter().enumerate() {
                s.y.set_index(j, *v);
            }
            problem.eqn.rhs().call_inplace(s.y, a, s.dy);
            s.set_step_size(cfg.h0, &problem.atol, problem.rtol, &problem.eqn, 1);
        }
        let mut solver = problem.bdf_solver::<LS>(state)?;
        let mut stop = cfg.max_step.map_or(b, |h| (a + h).min(b));
        solver.set_stop_time(stop)?;
        tr.segments += 1;
        loop {
            let reason = solver.step()?;
            let tt = solver.state().t;
            let root = match reason {
                OdeSolverStopReason::RootFound(t, idx) => Some((t, idx)),
                _ => None,
            };
            if let Some((t, index)) = root {
                // The step containing the root is accepted; its end is not
                // output (solve_ivp ends a terminal run at the event).
                let y = solver.interpolate(t)?.clone_as_vec();
                push(&mut tr, t, &y);
                tr.root = Some(Root {
                    t,
                    index,
                    y: y.clone(),
                    step_end: tt,
                });
                ycur = y;
                tr.t_last_raw = t;
                add_stats(&mut tr.stats, solver.get_statistics());
                tr.segment_us.push(seg_start.elapsed().as_secs_f64() * 1e6);
                break 'segments;
            }
            let at_stop = tt >= stop - 1e-9;
            if at_stop && stop >= b - 1e-9 {
                // Snap to the breakpoint, as solve_ivp clips to t_bound.
                if tt != b {
                    tr.inexact_segment_ends += 1;
                }
                tr.t_last_raw = tt;
                ycur = solver.state().y.clone_as_vec();
                push(&mut tr, b, &ycur);
                break;
            }
            push(&mut tr, tt, &solver.state().y.clone_as_vec());
            let want = match cfg.max_step {
                Some(h) if solver.state().h > h => (tt + h).min(b),
                _ => b,
            };
            if at_stop || want != stop {
                stop = want;
                solver.set_stop_time(stop)?;
            }
        }
        add_stats(&mut tr.stats, solver.get_statistics());
        tr.segment_us.push(seg_start.elapsed().as_secs_f64() * 1e6);
    }
    tr.final_y = ycur;
    tr.steps = tr.t.len();
    tr.evals = model.evals.get();
    tr.jv_calls = model.jv_calls.get();
    tr.root_evals = model.root_evals.get();
    tr.wall_us = start.elapsed().as_secs_f64() * 1e6;
    Ok(tr)
}

fn add_stats(s: &mut Stats, st: &diffsol::ode_solver::OdeSolverStatistics) {
    s.accepted_steps += st.number_of_steps;
    s.lu_setups += st.number_of_linear_solver_setups;
    s.error_test_failures += st.number_of_error_test_failures;
    s.newton_iterations += st.number_of_nonlinear_solver_iterations;
    s.newton_failures += st.number_of_nonlinear_solver_fails;
}

/// Median wall time per `MultiRhs::eval` call [us] over `batches` batches
/// of `calls` calls on the fixed state `x` at time `t`.
pub fn bench_eval(case: &Case, t: f64, x: &[f64], calls: usize, batches: usize) -> f64 {
    let rhs = case.rhs.build().expect("valid case");
    let mut out = vec![0.0; x.len()];
    let mut per: Vec<f64> = (0..batches)
        .map(|_| {
            let t0 = Instant::now();
            for _ in 0..calls {
                rhs.eval(t, std::hint::black_box(x), &mut out);
                std::hint::black_box(&out);
            }
            t0.elapsed().as_secs_f64() * 1e6 / calls as f64
        })
        .collect();
    per.sort_by(|a, b| a.total_cmp(b));
    per[per.len() / 2]
}

/// Median of `f()`'s wall times over `n` runs.
pub fn median_wall<T>(n: usize, mut f: impl FnMut() -> T) -> (T, Duration) {
    let mut out = None;
    let mut walls: Vec<Duration> = (0..n)
        .map(|_| {
            let t0 = Instant::now();
            out = Some(f());
            t0.elapsed()
        })
        .collect();
    walls.sort();
    (out.expect("n > 0"), walls[n / 2])
}
