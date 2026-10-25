//! diffsol drivers mirroring `_solve._segment`: one BDF solver per segment,
//! re-initialised from a given state at every breakpoint.
//!
//! The same driver is instantiated for three matrix/linear-solver
//! combinations by `backend!`: sparse (`FaerSparseMat` + `FaerSparseLU`),
//! dense nalgebra (`NalgebraMat` + `NalgebraLU`) and dense faer
//! (`FaerMat` + `FaerLU`).

use std::time::{Duration, Instant};

use diffsol::{
    FaerLU, FaerMat, FaerSparseLU, FaerSparseMat, NalgebraLU, NalgebraMat, NonLinearOp,
    NonLinearOpJacobian, OdeBuilder, OdeEquations, OdeSolverMethod, OdeSolverState,
    OdeSolverStopReason, Vector,
};

use diffsol::matrix::sparsity::MatrixSparsity;
use diffsol::Op;

use crate::explicit::ExplicitEqn;
use crate::model::{Jv, Model, N};

/// Run configuration.
#[derive(Clone)]
pub struct Cfg {
    /// Which `Jv` the solver gets (analytic, forward difference, or one of
    /// the shortcut variants of item 4b).
    pub jv: Jv,
    pub rtol: f64,
    pub atol: f64,
    /// Segment edges: `breaks[0]` is t0, the last is t_end.
    pub breaks: Vec<f64>,
    /// Stop at the first root (terminal event).
    pub terminal: bool,
    /// Emulated max step, enforced with `set_stop_time` (diffsol 0.18 has
    /// no max-step option): the next stop is at most this far ahead.
    pub max_step: Option<f64>,
    /// With `max_step`: re-arm the stop after every step (naive) instead of
    /// only when the solver's next step `h` would exceed the cap (lazy).
    pub hop_always: bool,
    /// Collect `(t, y)` for every accepted step.
    pub record: bool,
    /// Also sample the dense output on a uniform grid of this spacing.
    pub grid: Option<f64>,
    /// Initial step size handed to the first step of every segment.
    pub h0: f64,
}

impl Default for Cfg {
    fn default() -> Self {
        Self {
            jv: Jv::Analytic,
            rtol: 1e-6,
            atol: 1e-8,
            breaks: vec![0.0, 10.0],
            terminal: true,
            max_step: None,
            hop_always: false,
            record: false,
            grid: None,
            h0: 1e-6,
        }
    }
}

/// A located root.
#[derive(Clone, Debug)]
pub struct Root {
    pub t: f64,
    pub index: usize,
    /// State from `interpolate(t)`.
    pub y: Vec<f64>,
    /// Time of the accepted step that contained the root (its end).
    pub step_end: f64,
}

/// Result of a run.
#[derive(Default)]
pub struct Trace {
    pub t: Vec<f64>,
    pub y: Vec<Vec<f64>>,
    /// Grid times and states sampled with `interpolate`.
    pub gt: Vec<f64>,
    pub gy: Vec<Vec<f64>>,
    pub root: Option<Root>,
    pub final_y: Vec<f64>,
    pub final_t: f64,
    pub steps: usize,
    pub wall: Duration,
    pub segments: usize,
    pub lu_setups: usize,
    /// Non-zeros of the Jacobian pattern in use (sparse backends only).
    pub jac_nnz: Option<usize>,
    /// Time to build the problem (includes sparsity detection and
    /// colouring on the auto path); not part of `wall`.
    pub build: Duration,
}

/// Sample `interpolate` at every grid point up to `$upto` (inclusive).
macro_rules! emit_grid {
    ($solver:expr, $tr:expr, $next:expr, $upto:expr, $h:expr) => {
        while $next <= $upto + 1e-12 {
            let tg = $next.min($upto);
            $tr.gt.push(tg);
            $tr.gy.push($solver.interpolate(tg)?.clone_as_vec());
            $next += $h;
        }
    };
}

/// The `_segment` loop shared by every backend: per segment a fresh
/// `BdfState` re-initialised at the breakpoint, then `step()` to its end.
macro_rules! drive {
    ($problem:ident, $LS:ty, $model:ident, $y0:ident, $cfg:ident, $nnz:expr, $build:expr) => {{
        let start = Instant::now();
        let mut tr = Trace::default();
        tr.jac_nnz = $problem.eqn.rhs().jacobian_sparsity().and_then($nnz);
        tr.build = $build;
        let mut ycur = $y0.to_vec();
        let mut next_grid = $cfg.breaks[0];
        for w in $cfg.breaks.windows(2) {
            let (a, b) = (w[0], w[1]);
            $model.stim.set($model.stim_at(a));
            // Re-initialise at (a, ycur): the `_segment` step.
            let mut state = diffsol::BdfState::new_without_initialise(&$problem)?;
            {
                let mut s = state.as_mut();
                *s.t = a;
                for (i, v) in ycur.iter().enumerate() {
                    s.y.set_index(i, *v);
                }
                $problem.eqn.rhs().call_inplace(s.y, a, s.dy);
                s.set_step_size($cfg.h0, &$problem.atol, $problem.rtol, &$problem.eqn, 1);
            }
            let mut solver = $problem.bdf_solver::<$LS>(state)?;
            let mut stop = b;
            if let Some(h) = $cfg.max_step {
                stop = (a + h).min(b);
            }
            solver.set_stop_time(stop)?;
            tr.segments += 1;
            loop {
                let reason = solver.step()?;
                tr.steps += 1;
                let (tt, root) = match reason {
                    OdeSolverStopReason::RootFound(t, i) => (solver.state().t, Some((t, i))),
                    _ => (solver.state().t, None),
                };
                if let Some(h) = $cfg.grid {
                    let upto = root.map_or(tt, |r| r.0);
                    emit_grid!(solver, tr, next_grid, upto, h);
                }
                if $cfg.record {
                    tr.t.push(tt);
                    tr.y.push(solver.state().y.clone_as_vec());
                }
                if let Some((t, index)) = root {
                    tr.root = Some(Root {
                        t,
                        index,
                        y: solver.interpolate(t)?.clone_as_vec(),
                        step_end: tt,
                    });
                    if $cfg.terminal {
                        let st = solver.get_statistics();
                        tr.lu_setups += st.number_of_linear_solver_setups;
                        tr.wall = start.elapsed();
                        return Ok(tr);
                    }
                }
                let at_stop = tt >= stop - 1e-9;
                if at_stop && stop >= b - 1e-9 {
                    ycur = solver.state().y.clone_as_vec();
                    tr.final_t = tt;
                    break;
                }
                let want = match $cfg.max_step {
                    Some(h) if $cfg.hop_always || solver.state().h > h => (tt + h).min(b),
                    _ => b,
                };
                // RootFound returns before the solver looks at
                // tstop, so a continuing run must re-arm it.
                if at_stop || root.is_some() || want != stop {
                    stop = want;
                    solver.set_stop_time(stop)?;
                }
            }
            let st = solver.get_statistics();
            tr.lu_setups += st.number_of_linear_solver_setups;
        }
        tr.final_y = ycur;
        tr.wall = start.elapsed();
        Ok(tr)
    }};
}

macro_rules! backend {
    ($name:ident, $M:ty, $LS:ty, $nnz:expr) => {
        pub mod $name {
            use super::*;

            type M = $M;
            type LS = $LS;

            /// Run `cfg` from `y0` with the pattern auto-detected from `Jv`.
            pub fn solve(
                model: &Model,
                y0: &[f64],
                cfg: &Cfg,
            ) -> Result<Trace, diffsol::DiffsolError> {
                let jv = cfg.jv;
                let ystart = y0.to_vec();
                let t_build = Instant::now();
                let problem = OdeBuilder::<M>::new()
                    .t0(cfg.breaks[0])
                    .h0(cfg.h0)
                    .rtol(cfg.rtol)
                    .atol([cfg.atol])
                    .rhs_implicit(
                        |x, _p, t, y| model.rhs(x, t, y),
                        move |x, _p, t, v, y| model.jv_by(jv, x, t, v, y),
                    )
                    .init(|_p, _t, y| y.copy_from_slice(&ystart), N)
                    .root(|x, _p, t, g| model.root(x, t, g), 2)
                    .build()?;
                let build = t_build.elapsed();
                drive!(problem, LS, model, y0, cfg, $nnz, build)
            }
        }
    };
}

backend!(sparse, FaerSparseMat<f64>, FaerSparseLU<f64>, |s| Some(
    s.compute_nnz()
));
backend!(dense, NalgebraMat<f64>, NalgebraLU<f64>, |_| None);
backend!(faer_dense, FaerMat<f64>, FaerLU<f64>, |_| None);

/// Sparse backend with the pattern supplied by hand (`build_from_eqn` on
/// [`ExplicitEqn`]) instead of detected.
pub mod sparse_explicit {
    use super::*;

    type M = FaerSparseMat<f64>;
    type LS = FaerSparseLU<f64>;

    /// Run `cfg` from `y0` with the model's own pattern
    /// (`core.rs::tree_jac_sparsity` on its edges).
    pub fn solve(model: &Model, y0: &[f64], cfg: &Cfg) -> Result<Trace, diffsol::DiffsolError> {
        let (rows, cols) = model.pattern();
        solve_with_pattern(model, y0, cfg, &rows, &cols)
    }

    /// Run `cfg` from `y0` with the pattern `(rows, cols)`.
    pub fn solve_with_pattern(
        model: &Model,
        y0: &[f64],
        cfg: &Cfg,
        rows: &[usize],
        cols: &[usize],
    ) -> Result<Trace, diffsol::DiffsolError> {
        let t_build = Instant::now();
        let eqn = ExplicitEqn::<M>::new(model, cfg.jv, y0, rows, cols)?;
        let problem = OdeBuilder::<M>::new()
            .t0(cfg.breaks[0])
            .h0(cfg.h0)
            .rtol(cfg.rtol)
            .atol([cfg.atol])
            .build_from_eqn(eqn)?;
        let build = t_build.elapsed();
        drive!(
            problem,
            LS,
            model,
            y0,
            cfg,
            |s| Some(s.compute_nnz()),
            build
        )
    }
}

/// The pattern diffsol auto-detects (NaN probe of `jv` at `y0`, t = 0), as
/// sorted `(row, col)` pairs.
pub fn detected_pattern(model: &Model, y0: &[f64], jv: Jv) -> Vec<(usize, usize)> {
    let ystart = y0.to_vec();
    let problem = OdeBuilder::<FaerSparseMat<f64>>::new()
        .rhs_implicit(
            |x, _p, t, y| model.rhs(x, t, y),
            move |x, _p, t, v, y| model.jv_by(jv, x, t, v, y),
        )
        .init(|_p, _t, y| y.copy_from_slice(&ystart), N)
        .build()
        .unwrap();
    let s = problem.eqn.rhs().jacobian_sparsity().unwrap();
    let mut nz = MatrixSparsity::<FaerSparseMat<f64>>::indices(&s);
    nz.sort_unstable();
    nz
}

/// The Jacobian at `(y0, t = 0)` on the model's explicit pattern, filled
/// through the colouring (what each BDF Jacobian update computes).
pub fn jacobian_at(model: &Model, y0: &[f64], jv: Jv) -> FaerSparseMat<f64> {
    let (rows, cols) = model.pattern();
    let eqn = ExplicitEqn::<FaerSparseMat<f64>>::new(model, jv, y0, &rows, &cols).unwrap();
    let x = <FaerSparseMat<f64> as diffsol::MatrixCommon>::V::from_vec(y0.to_vec(), *eqn.context());
    eqn.rhs().jacobian(&x, 0.0)
}
