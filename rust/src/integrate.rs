//! The BDF integrator of `mode="native"`: diffsol 0.18 over [`MultiRhs`].
//!
//! Plain Rust, no pyo3. [`Integrator`] is built once per right-hand side
//! (the Jacobian pattern, its colouring and the diffsol problem) and
//! [`Integrator::run`] then solves one segment of `_solve.integrate`, the
//! job `_solve._segment` gives scipy's BDF in the other modes: from a fresh
//! state at the segment start, every accepted step recorded, stopped by the
//! two terminal events of `_solve` (`spike_event`, `settled_event`).
//!
//! The pattern is `core::tree_jac_sparsity` of the run's own morphology
//! plus the stimulus's own coupling (`MultiRhs::stim_jac_entries`), handed to diffsol through a hand-written `OdeEquations` (diffsol's
//! builder can only detect one, by NaN probes); the Jacobian is filled by
//! forward-difference `Jv` products, one per colour of that pattern.
//! `max_step` is not a diffsol option: it is emulated with stop-time hops,
//! set only when the next step would be longer.

use diffsol::matrix::sparsity::MatrixSparsity;
use diffsol::{
    BdfState, ConstantOp, DiffsolError, FaerContext, FaerSparseLU, FaerSparseMat, FaerVec,
    JacobianColoring, Matrix, NonLinearOp, NonLinearOpJacobian, OdeBuilder, OdeEquations,
    OdeEquationsRef, OdeSolverMethod, OdeSolverProblem, OdeSolverState, OdeSolverStopReason, Op,
    ParameterisedOp, UnitCallable, Vector,
};

use crate::core::{tree_jac_sparsity, CoreError, MultiRhs};

type M = FaerSparseMat<f64>;
type V = FaerVec<f64>;
type C = FaerContext;
type LS = FaerSparseLU<f64>;

/// The initial step of every segment [ms] (diffsol's own guess would
/// differ from scipy's; a small fixed one lets the error control grow it).
pub const H0: f64 = 1e-6;

/// `_solve.spike_event`: stop when x[axon_col] - x[0] rises to `factor`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spike {
    pub axon_col: usize,
    pub factor: f64,
}

/// `_solve.settled_event`: stop at the first time after `t_off` when
/// x[axon_col] - x[0] is below `threshold` (factor * SETTLED_FRACTION, which
/// Python computes) and every dV/dt is negative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settled {
    pub t_off: f64,
    pub axon_col: usize,
    pub threshold: f64,
}

/// The terminal events of one segment; both stop the run on a rising
/// crossing, as `direction = 1` does in solve_ivp.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Events {
    pub spike: Option<Spike>,
    pub settled: Option<Settled>,
}

/// One segment's output.
#[derive(Clone, Debug, PartialEq)]
pub struct Run {
    /// The segment start, every accepted step, then t1 (or the event time).
    pub t: Vec<f64>,
    /// The states at `t`, row-major (`t.len()` rows of `n` entries).
    pub y: Vec<f64>,
    /// Whether an event ended the segment.
    pub stopped: bool,
}

fn err(e: impl Into<DiffsolError>) -> CoreError {
    CoreError(format!("diffsol: {}", e.into()))
}

/// The equations diffsol solves: the right-hand side, its pattern and
/// colouring, and the events of the current run.
struct Eqn {
    rhs: MultiRhs,
    events: Events,
    nstates: usize,
    y0: V,
    sparsity: <M as Matrix>::Sparsity,
    coloring: JacobianColoring<M>,
    ctx: C,
}

impl Eqn {
    /// `J v` by one forward difference on two `eval` calls, with the step
    /// 1e-7 (1 + max |x|) of the diffsol prototypes (issue #28).
    fn jv(&self, x: &[f64], t: f64, v: &[f64], y: &mut [f64]) {
        let h = 1e-7 * (1.0 + x.iter().fold(0.0_f64, |a, b| a.max(b.abs())));
        let xp: Vec<f64> = x.iter().zip(v).map(|(a, b)| a + h * b).collect();
        let mut f0 = vec![0.0; x.len()];
        self.rhs.eval(t, x, &mut f0);
        self.rhs.eval(t, &xp, y);
        for (a, b) in y.iter_mut().zip(&f0) {
            *a = (*a - b) / h;
        }
    }

    /// The two root functions: g[0] the spike's, g[1] the settled one's,
    /// -1 when that event is off (as `settled_event` returns until it fires).
    fn roots(&self, x: &[f64], t: f64, g: &mut [f64]) {
        g[0] = match self.events.spike {
            Some(s) => x[s.axon_col] - x[0] - s.factor,
            None => -1.0,
        };
        g[1] = -1.0;
        if let Some(s) = self.events.settled {
            if t >= s.t_off && x[s.axon_col] - x[0] < s.threshold {
                let mut dx = vec![0.0; x.len()];
                self.rhs.eval(t, x, &mut dx);
                if dx[..self.rhs.n()].iter().all(|d| *d < 0.0) {
                    g[1] = 1.0;
                }
            }
        }
    }
}

/// The right-hand side operator.
struct Rhs<'a>(&'a Eqn);
/// The initial-state operator (unused: every segment sets its own state).
struct Init<'a>(&'a Eqn);
/// The two root functions.
struct Root<'a>(&'a Eqn);

macro_rules! impl_op {
    ($ty:ident, $nout:expr) => {
        impl Op for $ty<'_> {
            type T = f64;
            type V = V;
            type M = M;
            type C = C;
            fn context(&self) -> &C {
                &self.0.ctx
            }
            fn nstates(&self) -> usize {
                self.0.nstates
            }
            fn nout(&self) -> usize {
                $nout(self.0)
            }
            fn nparams(&self) -> usize {
                0
            }
        }
    };
}

impl_op!(Rhs, |e: &Eqn| e.nstates);
impl_op!(Init, |e: &Eqn| e.nstates);
impl_op!(Root, |_: &Eqn| 2);

impl NonLinearOp for Rhs<'_> {
    fn call_inplace(&self, x: &V, t: f64, y: &mut V) {
        let rhs = &self.0.rhs;
        y.for_each_batch_host([x], |y, [x], _| rhs.eval(t, x, y));
    }
}

impl NonLinearOpJacobian for Rhs<'_> {
    fn jac_mul_inplace(&self, x: &V, t: f64, v: &V, y: &mut V) {
        let e = self.0;
        y.for_each_batch_host([x, v], |y, [x, v], _| e.jv(x, t, v, y));
    }
    fn jacobian_inplace(&self, x: &V, t: f64, y: &mut M) {
        self.0.coloring.jacobian_inplace(self, x, t, y);
    }
    fn jacobian_sparsity(&self) -> Option<<M as Matrix>::Sparsity> {
        Some(self.0.sparsity.clone())
    }
}

impl ConstantOp for Init<'_> {
    fn call_inplace(&self, _t: f64, y: &mut V) {
        y.copy_from(&self.0.y0);
    }
}

impl NonLinearOp for Root<'_> {
    fn call_inplace(&self, x: &V, t: f64, y: &mut V) {
        let e = self.0;
        y.for_each_batch_host([x], |g, [x], _| e.roots(x, t, g));
    }
}

impl Op for Eqn {
    type T = f64;
    type V = V;
    type M = M;
    type C = C;
    fn context(&self) -> &C {
        &self.ctx
    }
    fn nstates(&self) -> usize {
        self.nstates
    }
    fn nout(&self) -> usize {
        self.nstates
    }
    fn nparams(&self) -> usize {
        0
    }
}

impl<'a> OdeEquationsRef<'a> for Eqn {
    type Rhs = Rhs<'a>;
    type Init = Init<'a>;
    type Root = Root<'a>;
    type Mass = ParameterisedOp<'a, UnitCallable<M>>;
    type Out = ParameterisedOp<'a, UnitCallable<M>>;
    type Reset = ParameterisedOp<'a, UnitCallable<M>>;
}

impl OdeEquations for Eqn {
    fn rhs(&self) -> Rhs<'_> {
        Rhs(self)
    }
    fn mass(&self) -> Option<ParameterisedOp<'_, UnitCallable<M>>> {
        None
    }
    fn root(&self) -> Option<Root<'_>> {
        // no root function at all when no event is asked for
        (self.events != Events::default()).then_some(Root(self))
    }
    fn init(&self) -> Init<'_> {
        Init(self)
    }
    fn set_params(&mut self, _p: &V) {}
    fn get_params(&self, _p: &mut V) {}
}

/// diffsol's BDF over one right-hand side, built once and reused for every
/// segment of every run that uses it.
pub struct Integrator {
    problem: OdeSolverProblem<Eqn>,
}

impl Integrator {
    /// Builds the problem: the pattern of `rhs`'s own tree, its colouring
    /// and the diffsol problem around them (tolerances are set per run).
    pub fn new(rhs: MultiRhs) -> Result<Self, CoreError> {
        let nstates = 7 * rhs.n();
        let (parent, child) = rhs.edges();
        let (rows, cols) = tree_jac_sparsity(rhs.n(), parent, child);
        let mut pairs: Vec<(usize, usize)> = rows.into_iter().zip(cols).collect();
        // plus the stimulus's own coupling, kept in CSC order
        pairs.extend(rhs.stim_jac_entries());
        pairs.sort_unstable_by_key(|&(r, c)| (c, r));
        pairs.dedup();
        let ctx = C::default();
        let sparsity = <<M as Matrix>::Sparsity as MatrixSparsity<M>>::try_from_indices(
            nstates,
            nstates,
            pairs.clone(),
        )
        .map_err(err)?;
        let coloring = JacobianColoring::new(&sparsity, &pairs, ctx);
        let eqn = Eqn {
            rhs,
            events: Events::default(),
            nstates,
            y0: V::zeros(nstates, ctx),
            sparsity,
            coloring,
            ctx,
        };
        let problem = OdeBuilder::<M>::new()
            .h0(H0)
            .build_from_eqn(eqn)
            .map_err(err)?;
        Ok(Self { problem })
    }

    /// The number of entries in the Jacobian pattern.
    pub fn nnz(&self) -> usize {
        self.problem.eqn.sparsity.compute_nnz()
    }

    /// The number of states (7n).
    pub fn nstates(&self) -> usize {
        self.problem.eqn.nstates
    }

    /// Solves dx/dt = rhs(t, x) on [t0, t1] from y0.
    ///
    /// The output starts with (t0, y0) (diffsol's `step` does not emit
    /// it), then has every accepted step. Without an event it ends at t1
    /// exactly (the last stop time is snapped to t1, as solve_ivp clips to
    /// t_bound); an event ends it at the root, with the state interpolated
    /// there (the root lies inside the last accepted step, whose end is not
    /// output), and sets `stopped`. A falling crossing does not stop it.
    ///
    /// `max_step` (none if not finite and positive) bounds every step: when
    /// the solver's next step would be longer, the stop time is moved to
    /// one `max_step` ahead; the tolerances are never touched.
    #[allow(clippy::too_many_arguments)]
    pub fn run(
        &mut self,
        y0: &[f64],
        t0: f64,
        t1: f64,
        rtol: f64,
        atol: f64,
        max_step: f64,
        h0: f64,
        events: &Events,
    ) -> Result<Run, CoreError> {
        let n = self.nstates();
        if y0.len() != n {
            return Err(CoreError(format!(
                "y0 has {} entries, expected {n}",
                y0.len()
            )));
        }
        if t1.partial_cmp(&t0) != Some(std::cmp::Ordering::Greater) {
            return Err(CoreError(format!("need t0 < t1, got {t0}, {t1}")));
        }
        let ctx = self.problem.eqn.ctx;
        self.problem.rtol = rtol;
        self.problem.atol = V::from_element(n, atol, ctx);
        self.problem.h0 = h0;
        self.problem.t0 = t0;
        self.problem.eqn.events = *events;
        let problem = &self.problem;
        let cap = (max_step.is_finite() && max_step > 0.0).then_some(max_step);

        let mut out = Run {
            t: vec![t0],
            y: y0.to_vec(),
            stopped: false,
        };
        // a fresh state at (t0, y0): `_segment`'s re-initialisation
        let mut state = BdfState::new_without_initialise(problem).map_err(err)?;
        {
            let mut s = state.as_mut();
            *s.t = t0;
            for (j, v) in y0.iter().enumerate() {
                s.y.set_index(j, *v);
            }
            problem.eqn.rhs().call_inplace(s.y, t0, s.dy);
            s.set_step_size(h0, &problem.atol, problem.rtol, &problem.eqn, 1);
        }
        let mut solver = problem.bdf_solver::<LS>(state).map_err(err)?;
        let mut stop = cap.map_or(t1, |h| (t0 + h).min(t1));
        solver.set_stop_time(stop).map_err(err)?;
        let mut g_prev = [0.0; 2];
        loop {
            let reason = solver.step().map_err(err)?;
            let tt = solver.state().t;
            if let OdeSolverStopReason::RootFound(t_root, idx) = reason {
                // the step started at the last row; a rising crossing is
                // one that started below zero
                let start = &out.y[out.y.len() - n..];
                problem
                    .eqn
                    .roots(start, out.t[out.t.len() - 1], &mut g_prev);
                if g_prev[idx] < 0.0 {
                    let y = solver.interpolate(t_root).map_err(err)?.clone_as_vec();
                    out.t.push(t_root);
                    out.y.extend_from_slice(&y);
                    out.stopped = true;
                    return Ok(out);
                }
            }
            let at_stop = matches!(reason, OdeSolverStopReason::TstopReached) || tt >= stop - 1e-9;
            if at_stop && stop >= t1 - 1e-9 {
                out.t.push(t1);
                out.y.extend_from_slice(&solver.state().y.clone_as_vec());
                return Ok(out);
            }
            out.t.push(tt);
            out.y.extend_from_slice(&solver.state().y.clone_as_vec());
            let want = match cap {
                Some(h) if solver.state().h > h => (tt + h).min(t1),
                _ => t1,
            };
            if at_stop || want != stop {
                stop = want;
                solver.set_stop_time(stop).map_err(err)?;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::tests::{model, state, step_model, synaptic_model};

    const N: usize = 7 * 45;
    const TOL: f64 = 1e-8;

    fn row(r: &Run, i: usize) -> &[f64] {
        &r.y[i * N..(i + 1) * N]
    }

    /// x[col] - x[0] at every output row.
    fn diff(r: &Run, col: usize) -> Vec<f64> {
        (0..r.t.len())
            .map(|i| row(r, i)[col] - row(r, i)[0])
            .collect()
    }

    fn run(it: &mut Integrator, y0: &[f64], t1: f64, max_step: f64, ev: &Events) -> Run {
        it.run(y0, 0.0, t1, TOL, TOL, max_step, H0, ev).unwrap()
    }

    #[test]
    fn integrator_is_send() {
        // not Sync (diffsol keeps RefCell scratch space): the binding holds
        // it in a Mutex
        fn check<T: Send>() {}
        check::<Integrator>();
    }

    #[test]
    fn the_chain_pattern_has_943_entries() {
        let it = Integrator::new(model()).unwrap();
        assert_eq!(it.nnz(), 943);
        assert_eq!(it.nstates(), N);
    }

    #[test]
    fn a_synaptic_input_off_the_soma_adds_its_soma_column_entry() {
        // (1, 0) is a chain entry already; (2, 0) is not
        assert_eq!(Integrator::new(synaptic_model(0)).unwrap().nnz(), 943);
        assert_eq!(Integrator::new(synaptic_model(1)).unwrap().nnz(), 943);
        assert_eq!(Integrator::new(synaptic_model(2)).unwrap().nnz(), 944);
    }

    #[test]
    fn a_run_starts_at_y0_and_ends_on_t1_exactly() {
        let mut it = Integrator::new(model()).unwrap();
        let y0 = state();
        let r = it
            .run(&y0, 0.25, 7.3, TOL, TOL, 1.0, H0, &Events::default())
            .unwrap();
        assert!(!r.stopped);
        assert_eq!(r.y.len(), r.t.len() * N);
        assert_eq!((r.t[0], row(&r, 0)), (0.25, &y0[..]));
        assert_eq!(*r.t.last().unwrap(), 7.3);
        assert!(r.t.windows(2).all(|w| w[0] < w[1]));
        assert!(r.y.iter().all(|v| v.is_finite()));
        // the same build serves a second run unchanged
        let again = it
            .run(&y0, 0.25, 7.3, TOL, TOL, 1.0, H0, &Events::default())
            .unwrap();
        assert_eq!(again, r);
    }

    #[test]
    fn max_step_bounds_every_step() {
        let mut it = Integrator::new(model()).unwrap();
        let r = run(&mut it, &state(), 0.5, 0.005, &Events::default());
        assert!(r.t.len() > 100);
        assert!(r.t.windows(2).all(|w| w[1] - w[0] <= 0.005 + 1e-12));
        assert_eq!(*r.t.last().unwrap(), 0.5);
    }

    #[test]
    fn a_spike_root_stops_on_the_crossing() {
        // the hyperpolarised soma lets x[10] - x[0] rise
        let mut it = Integrator::new(step_model(-5000.0)).unwrap();
        let y0 = state();
        let free = run(&mut it, &y0, 5.0, 0.5, &Events::default());
        let d = diff(&free, 10);
        let max = d.iter().copied().fold(f64::MIN, f64::max);
        assert!(max > d[0] + 1.0, "no rise to stop on: {} -> {max}", d[0]);
        let factor = 0.5 * (d[0] + max);
        let spike = Spike {
            axon_col: 10,
            factor,
        };
        let ev = Events {
            spike: Some(spike),
            settled: None,
        };
        let r = run(&mut it, &y0, 5.0, 0.5, &ev);
        assert!(r.stopped);
        let t_root = *r.t.last().unwrap();
        let last = row(&r, r.t.len() - 1);
        assert!((last[10] - last[0] - factor).abs() < 1e-9);
        // the event-free run crosses between the same two steps
        let k = d.iter().position(|v| *v >= factor).unwrap();
        assert!(free.t[k - 1] <= t_root && t_root <= free.t[k]);
        // up to the stop the two runs are the same steps
        assert_eq!(r.t[..r.t.len() - 1], free.t[..r.t.len() - 1]);
    }

    #[test]
    fn a_settled_root_stops_after_t_off() {
        // every V above its rest: all fall from the start
        let mut y0 = state();
        y0[..45].iter_mut().for_each(|v| *v = -55.0);
        let mut it = Integrator::new(model()).unwrap();
        let settled = Settled {
            t_off: 1.0,
            axon_col: 10,
            threshold: 1e9,
        };
        let ev = Events {
            spike: None,
            settled: Some(settled),
        };
        let r = run(&mut it, &y0, 4.0, 0.5, &ev);
        assert!(r.stopped);
        let t_stop = *r.t.last().unwrap();
        assert!((1.0..4.0).contains(&t_stop), "stopped at {t_stop}");
        // the condition holds where it stopped
        let mut dx = vec![0.0; N];
        model().eval(t_stop, row(&r, r.t.len() - 1), &mut dx);
        assert!(dx[..45].iter().all(|v| *v < 0.0));
    }

    #[test]
    fn a_falling_crossing_does_not_stop_the_run() {
        // x[10] starts 30 mV above its neighbours and relaxes
        let mut it = Integrator::new(model()).unwrap();
        let mut y0 = state();
        y0[10] += 30.0;
        let free = run(&mut it, &y0, 2.0, 0.5, &Events::default());
        let d = diff(&free, 10);
        let min = d.iter().copied().fold(f64::MAX, f64::min);
        assert!(min < d[0] - 0.1, "no fall to cross: {} -> {min}", d[0]);
        let spike = Spike {
            axon_col: 10,
            factor: 0.5 * (d[0] + min),
        };
        let ev = Events {
            spike: Some(spike),
            settled: None,
        };
        let r = run(&mut it, &y0, 2.0, 0.5, &ev);
        assert!(!r.stopped);
        assert_eq!(r.t, free.t);
    }
}
