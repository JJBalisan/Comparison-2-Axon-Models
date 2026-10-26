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
//! (every stimulus couples a compartment only to itself), handed to diffsol
//! through a hand-written `OdeEquations` (diffsol's builder can only detect
//! one, by NaN probes); the Jacobian is filled by
//! forward-difference `Jv` products, one per colour of that pattern, and
//! the Newton matrix I - cJ is solved directly on that pattern's structure
//! ([`TreeLu`]), with faer's sparse LU as the fall-back.
//! `max_step` is not a diffsol option: it is emulated by narrowing the
//! controller's growth bounds before every step, so that no step is
//! proposed past it, with a stop-time hop as the fallback when a step
//! still would be longer.
//!
//! `prototype/diffsol-probe` is the design record (why diffsol, the
//! pitfalls, the costs); this file is the source of truth.

use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};

use diffsol::matrix::sparsity::MatrixSparsity;
use diffsol::{
    BdfState, ConstantOp, Context, DiffsolError, FaerContext, FaerSparseLU, FaerSparseMat, FaerVec,
    JacobianColoring, LaError, Matrix, MatrixCommon, NonLinearOp, NonLinearOpJacobian, OdeBuilder,
    OdeEquations, OdeEquationsRef, OdeSolverMethod, OdeSolverProblem, OdeSolverState,
    OdeSolverStopReason, Op, ParameterisedOp, UnitCallable, Vector,
};
use diffsol_la::{LinearOp, LinearSolver as LaLinearSolver};

use crate::core::{tree_jac_sparsity, CoreError, MultiRhs};

type M = FaerSparseMat<f64>;
type V = FaerVec<f64>;
type C = FaerContext;
type LS = TreeLu;

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
    /// Two state-sized buffers for `jv` and `roots`, so neither allocates.
    scratch: RefCell<[Vec<f64>; 2]>,
}

impl Eqn {
    /// `J v` by one forward difference on two `eval` calls, with the step
    /// 1e-7 (1 + max |x|) of the diffsol prototypes (issue #28).
    fn jv(&self, x: &[f64], t: f64, v: &[f64], y: &mut [f64]) {
        let h = 1e-7 * (1.0 + x.iter().fold(0.0_f64, |a, b| a.max(b.abs())));
        let [xp, f0] = &mut *self.scratch.borrow_mut();
        for ((p, a), b) in xp.iter_mut().zip(x).zip(v) {
            *p = a + h * b;
        }
        self.rhs.eval(t, x, f0);
        self.rhs.eval(t, xp, y);
        for (a, b) in y.iter_mut().zip(f0.iter()) {
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
                let [dx, _] = &mut *self.scratch.borrow_mut();
                self.rhs.eval(t, x, dx);
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

/// No parent: a root of the V tree.
const ROOT: usize = usize::MAX;

/// The smallest tree pivot [`TreeLu`] takes, as a fraction of the V
/// diagonal entry it came from; a smaller one hands that linearisation to
/// faer.
///
/// The tree is eliminated without pivoting, and an elimination with pivot
/// p multiplies the entries it updates by at most |a / p|. Here the V
/// diagonal of I - cJ is 1 + c * (conductances + axial) / C, never below
/// the sum of its off-diagonals, so a pivot down to 1e-3 of it bounds that
/// growth near 1e3: at most 3 of 16 digits lost in a Newton correction
/// that only needs a few. A pivot gets that small only where the Na
/// current's own feedback (eliminated from the gates into the V diagonal)
/// nearly cancels it, on a long step through an upstroke; there faer's
/// partial pivoting is the better solve.
const PIVOT_FLOOR: f64 = 1e-3;

/// Linearisations [`TreeLu`] has handed to faer, over the process.
static FALLBACKS: AtomicUsize = AtomicUsize::new(0);

/// The number of linearisations solved by faer's LU rather than the tree
/// elimination since the process started (0 on every model `core.rs`
/// builds; a log line per fallback goes to stderr when
/// `MSOAXON_NATIVE_LOG_FALLBACK` is set).
pub fn tree_lu_fallbacks() -> usize {
    FALLBACKS.load(Ordering::Relaxed)
}

/// The structure of I - cJ read off its pattern, and one factorisation.
struct Tree {
    n: usize,
    /// The V tree: every node's parent (ROOT for a root) and a
    /// breadth-first order (every parent before its children).
    parent: Vec<usize>,
    order: Vec<usize>,
    /// For every stored entry, in CSC order, where its value goes in `m`.
    slot: Vec<usize>,
    /// The entries: diag (n), d (6n), e (6n), f (6n), up (n), down (n),
    /// with d = A[g, g], e = A[g, V_i] and f = A[V_i, g] for gate g of
    /// compartment i, and up[i] = A[V_p, V_i], down[i] = A[V_i, V_p] for
    /// i's parent p (both 0 at a root).
    m: Vec<f64>,
    /// 1 / d per gate, and the tree pivots.
    rd: Vec<f64>,
    piv: Vec<f64>,
    /// The matrix diffsol fills.
    a: M,
}

impl Tree {
    /// Reads the structure off `a`'s pattern: None unless it is exactly
    /// msoaxon's (`core::tree_jac_sparsity` of a forest, with every
    /// diagonal stored), every stored entry accounted for once.
    fn new(a: M) -> Option<Self> {
        let nstates = a.nrows();
        if nstates == 0 || !nstates.is_multiple_of(7) || a.ncols() != nstates {
            return None;
        }
        let n = nstates / 7;
        let entries: Vec<(usize, usize)> = a.triplet_iter().0.collect();
        let mut nbrs = vec![Vec::new(); n];
        for &(r, c) in &entries {
            if r < n && c < n && r != c {
                nbrs[c].push(r);
            }
        }
        // breadth-first from every unvisited node; a neighbour seen before
        // that is not the parent closes a cycle
        let mut parent = vec![ROOT; n];
        let mut seen = vec![false; n];
        let mut order = Vec::with_capacity(n);
        for root in 0..n {
            if seen[root] {
                continue;
            }
            seen[root] = true;
            let mut head = order.len();
            order.push(root);
            while head < order.len() {
                let i = order[head];
                head += 1;
                for &j in &nbrs[i] {
                    if !seen[j] {
                        seen[j] = true;
                        parent[j] = i;
                        order.push(j);
                    } else if parent[i] != j {
                        return None;
                    }
                }
            }
        }
        let (dg, eg, fg, up, down) = (n, 7 * n, 13 * n, 19 * n, 20 * n);
        let mut slot = Vec::with_capacity(entries.len());
        let mut hit = vec![false; 21 * n];
        for &(r, c) in &entries {
            let s = match (r < n, c < n) {
                (true, true) if r == c => r,
                (true, true) if parent[c] == r => up + c,
                (true, true) if parent[r] == c => down + r,
                (false, true) if r % n == c => eg + r - n,
                (true, false) if c % n == r => fg + c - n,
                (false, false) if r == c => dg + r - n,
                _ => return None,
            };
            if std::mem::replace(&mut hit[s], true) {
                return None;
            }
            slot.push(s);
        }
        // a coupling may be absent (a zero), a diagonal may not
        if !(0..n).all(|i| hit[i]) || !hit[dg..eg].iter().all(|h| *h) {
            return None;
        }
        Some(Self {
            n,
            parent,
            order,
            slot,
            m: vec![0.0; 21 * n],
            rd: vec![0.0; 6 * n],
            piv: vec![0.0; n],
            a,
        })
    }

    /// Factorises the values in `self.a`: false if a gate or tree pivot is
    /// zero, not finite or under the floor.
    fn factor(&mut self) -> bool {
        let n = self.n;
        self.m.fill(0.0);
        for (&s, v) in self.slot.iter().zip(self.a.triplet_iter().1) {
            self.m[s] = v;
        }
        let (diag, rest) = self.m.split_at(n);
        let (d, rest) = rest.split_at(6 * n);
        let (e, rest) = rest.split_at(6 * n);
        let (f, rest) = rest.split_at(6 * n);
        let (up, down) = rest.split_at(n);
        let s = &mut self.piv;
        s.copy_from_slice(diag);
        // each gate's row eliminates on its own, into its V's diagonal
        for g in 0..6 * n {
            if !(d[g] != 0.0 && d[g].is_finite()) {
                return false;
            }
            self.rd[g] = 1.0 / d[g];
            s[g % n] -= f[g] * e[g] * self.rd[g];
        }
        // then the tree, leaves first, each node into its parent
        for &i in self.order.iter().rev() {
            if !(s[i].abs() >= PIVOT_FLOOR * diag[i].abs() && s[i].is_finite()) {
                return false;
            }
            let p = self.parent[i];
            if p != ROOT {
                s[p] -= up[i] / s[i] * down[i];
            }
        }
        true
    }

    /// Solves A x = b, or A^T x = b, in place.
    fn solve(&self, b: &mut [f64], transpose: bool) {
        let n = self.n;
        let (_, rest) = self.m.split_at(7 * n);
        let (e, rest) = rest.split_at(6 * n);
        let (f, rest) = rest.split_at(6 * n);
        let (up, down) = rest.split_at(n);
        // A^T swaps e with f and up with down, with the same pivots
        let (e, f, up, down) = if transpose {
            (f, e, down, up)
        } else {
            (e, f, up, down)
        };
        let (bv, bg) = b.split_at_mut(n);
        for g in 0..6 * n {
            bv[g % n] -= f[g] * self.rd[g] * bg[g];
        }
        for &i in self.order.iter().rev() {
            let p = self.parent[i];
            if p != ROOT {
                bv[p] -= up[i] / self.piv[i] * bv[i];
            }
        }
        for &i in &self.order {
            let p = self.parent[i];
            let coupled = if p == ROOT { 0.0 } else { down[i] * bv[p] };
            bv[i] = (bv[i] - coupled) / self.piv[i];
        }
        for g in 0..6 * n {
            bg[g] = self.rd[g] * (bg[g] - e[g] * bv[g % n]);
        }
    }
}

/// The linear solver of the Newton iteration: I - cJ solved directly on
/// its structure, every gate on its own and then the V tree leaves to root
/// (O(n), no ordering, no symbolic step per segment), where the pattern is
/// msoaxon's and every pivot clears [`PIVOT_FLOOR`]; faer's sparse LU
/// otherwise. On EPSGpair faer's general LU was ~37% of the run.
#[derive(Default)]
struct TreeLu {
    tree: Option<Tree>,
    faer: Option<FaerSparseLU<f64>>,
    /// Whether the current linearisation is faer's.
    use_faer: bool,
}

impl LaLinearSolver<M> for TreeLu {
    fn set_sparsity<O: LinearOp<T = f64, V = V, M = M, C = C>>(&mut self, op: &O) {
        let (nrows, ncols) = (op.nrows(), op.ncols());
        self.tree = match op.sparsity() {
            Some(s) if op.context().nbatch() == 1 => {
                Tree::new(M::new_from_sparsity(nrows, ncols, Some(s), *op.context()))
            }
            _ => None,
        };
        self.faer = None;
        self.use_faer = false;
    }

    fn set_linearisation<O: LinearOp<T = f64, V = V, M = M, C = C>>(&mut self, op: &O) {
        if let Some(t) = self.tree.as_mut() {
            op.matrix_inplace(&mut t.a);
            if t.factor() {
                self.use_faer = false;
                return;
            }
        }
        FALLBACKS.fetch_add(1, Ordering::Relaxed);
        if std::env::var_os("MSOAXON_NATIVE_LOG_FALLBACK").is_some() {
            let why = if self.tree.is_some() {
                "a small pivot"
            } else {
                "the pattern"
            };
            eprintln!("msoaxon native: faer LU for this linearisation ({why})");
        }
        let lu = self.faer.get_or_insert_with(|| {
            let mut lu = FaerSparseLU::default();
            lu.set_sparsity(op);
            lu
        });
        lu.set_linearisation(op);
        self.use_faer = true;
    }

    fn solve_in_place(&self, x: &mut V) -> Result<(), LaError> {
        self.solve_impl(x, false)
    }

    fn solve_transpose_in_place(&self, x: &mut V) -> Result<(), LaError> {
        self.solve_impl(x, true)
    }
}

impl TreeLu {
    fn solve_impl(&self, x: &mut V, transpose: bool) -> Result<(), LaError> {
        match (&self.tree, &self.faer) {
            (Some(t), _) if !self.use_faer => {
                x.for_each_batch_host([], |b, [], _| t.solve(b, transpose));
                Ok(())
            }
            (_, Some(lu)) if transpose => lu.solve_transpose_in_place(x),
            (_, Some(lu)) => lu.solve_in_place(x),
            // diffsol always linearises before it solves
            _ => panic!("TreeLu: solve before set_linearisation"),
        }
    }
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
        let pairs: Vec<(usize, usize)> = rows.into_iter().zip(cols).collect();
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
            scratch: RefCell::new([vec![0.0; nstates], vec![0.0; nstates]]),
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
    /// `max_step` (none if not finite and positive) bounds every step:
    /// before each step the controller's growth bounds are narrowed so the
    /// step after it can be at most `max_step`, and should the next step
    /// still be longer the stop time is moved to one `max_step` ahead; the
    /// tolerances are never touched.
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
        // diffsol's growth bounds, narrowed below while max_step binds. A
        // hop shortens the step the controller chose after the Newton
        // matrix was already refactorised for it, so every hop cost two LU
        // set-ups; a step that is never proposed past the cap costs none.
        let (grow_max, grow_min) = {
            let c = solver.config();
            (c.maximum_timestep_growth, c.minimum_timestep_growth)
        };
        let next_stop = |h: f64, t: f64| match cap {
            Some(c) if h > c => (t + c).min(t1),
            _ => t1,
        };
        let mut stop = next_stop(h0, t0);
        solver.set_stop_time(stop).map_err(err)?;
        // the vector the output rows are copied out of (diffsol's vectors
        // lend their slice only to a closure on a mutable operand)
        let mut row = V::zeros(n, ctx);
        let push = |out: &mut Run, row: &mut V, t: f64, y: &V| {
            out.t.push(t);
            row.for_each_batch_host([y], |_, [y], _| out.y.extend_from_slice(y));
        };
        let mut g_prev = [0.0; 2];
        loop {
            if let Some(c) = cap {
                // the next step is h * factor with factor <= the growth
                // bound; lowering the dead zone's floor with it lets the
                // step grow to the cap itself, as scipy's clip does
                let g = grow_max.min(c / solver.state().h);
                let cfg = solver.config_mut();
                cfg.maximum_timestep_growth = g;
                cfg.minimum_timestep_growth = grow_min.min(g);
            }
            let reason = solver.step().map_err(err)?;
            let tt = solver.state().t;
            let root = matches!(reason, OdeSolverStopReason::RootFound(..));
            if let OdeSolverStopReason::RootFound(t_root, idx) = reason {
                // the step started at the last row; a rising crossing is
                // one that started below zero
                let start = &out.y[out.y.len() - n..];
                problem
                    .eqn
                    .roots(start, out.t[out.t.len() - 1], &mut g_prev);
                if g_prev[idx] < 0.0 {
                    let y = solver.interpolate(t_root).map_err(err)?;
                    push(&mut out, &mut row, t_root, &y);
                    out.stopped = true;
                    return Ok(out);
                }
            }
            let at_stop = matches!(reason, OdeSolverStopReason::TstopReached) || tt >= stop - 1e-9;
            if at_stop && stop >= t1 - 1e-9 {
                push(&mut out, &mut row, t1, solver.state().y);
                return Ok(out);
            }
            push(&mut out, &mut row, tt, solver.state().y);
            // a RootFound step returns before diffsol checks the stop time,
            // so the next step could pass it (StopTimeBeforeCurrentTime, a
            // panic inside diffsol): re-arm it after every root
            let want = next_stop(solver.state().h, tt);
            if at_stop || root || want != stop {
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
    fn a_synaptic_input_off_the_soma_adds_no_entry() {
        // the input drives its own compartment, a diagonal entry already
        for k in 0..3 {
            assert_eq!(Integrator::new(synaptic_model(k)).unwrap().nnz(), 943);
        }
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

    /// The relaxing x[10] of `a_falling_crossing_does_not_stop_the_run`:
    /// its y0 and event-free run.
    fn relaxing(it: &mut Integrator, t1: f64, max_step: f64, h0: f64) -> (Vec<f64>, Run) {
        let mut y0 = state();
        y0[10] += 30.0;
        let free = it
            .run(&y0, 0.0, t1, TOL, TOL, max_step, h0, &Events::default())
            .unwrap();
        (y0, free)
    }

    /// A spike event that `r` crosses, falling, inside its step k.
    fn falling_in_step(r: &Run, k: usize) -> Events {
        let d = diff(r, 10);
        assert!(d[k] < d[k - 1], "step {k} does not fall");
        let spike = Spike {
            axon_col: 10,
            factor: 0.5 * (d[k - 1] + d[k]),
        };
        Events {
            spike: Some(spike),
            settled: None,
        }
    }

    #[test]
    fn a_falling_crossing_just_before_t1_still_ends_on_t1() {
        // diffsol returns RootFound before it shortens the next step to the
        // stop time: without a re-arm the next step passes t1 and diffsol
        // panics (StopTimeBeforeCurrentTime)
        let mut it = Integrator::new(model()).unwrap();
        let (y0, free) = relaxing(&mut it, 2.0, 0.5, H0);
        let k = free.t.len() / 2;
        let ev = falling_in_step(&free, k);
        let t1 = 0.5 * (free.t[k] + free.t[k + 1]);
        let r = it.run(&y0, 0.0, t1, TOL, TOL, 0.5, H0, &ev).unwrap();
        assert!(!r.stopped);
        assert_eq!(*r.t.last().unwrap(), t1);
        assert_eq!(r.t[..=k], free.t[..=k]);
    }

    #[test]
    fn a_falling_crossing_inside_a_hop_does_not_stop_the_run() {
        // h0 above max_step: the first step is cut to a stop time one
        // max_step ahead, and the crossing lies inside it
        let (t1, cap) = (5e-5, 1e-6);
        let mut it = Integrator::new(model()).unwrap();
        let (y0, free) = relaxing(&mut it, t1, cap, 2.0 * cap);
        assert!((free.t[1] - cap).abs() < 1e-18, "no hop: {}", free.t[1]);
        let ev = falling_in_step(&free, 1);
        let r = it.run(&y0, 0.0, t1, TOL, TOL, cap, 2.0 * cap, &ev).unwrap();
        assert!(!r.stopped);
        assert_eq!(r.t, free.t);
        assert!(r.t.windows(2).all(|w| w[1] - w[0] <= cap * (1.0 + 1e-9)));
    }

    /// A fixed matrix as the linear operator a solver is handed.
    struct Fixed(M, <M as Matrix>::Sparsity, C);

    impl LinearOp for Fixed {
        type T = f64;
        type V = V;
        type M = M;
        type C = C;
        fn nrows(&self) -> usize {
            self.0.nrows()
        }
        fn ncols(&self) -> usize {
            self.0.ncols()
        }
        fn context(&self) -> &C {
            &self.2
        }
        fn matrix_inplace(&self, y: &mut M) {
            y.copy_from(&self.0);
        }
        fn sparsity(&self) -> Option<<M as Matrix>::Sparsity> {
            Some(self.1.clone())
        }
    }

    /// I - cJ of `rhs` at (t, x), on the pattern's entries plus `extra`
    /// (each 0.25), as the operator diffsol hands over.
    fn newton_matrix(rhs: MultiRhs, x: &[f64], t: f64, c: f64, extra: &[(usize, usize)]) -> Fixed {
        let it = Integrator::new(rhs).unwrap();
        let eqn = &it.problem.eqn;
        let ctx = eqn.ctx;
        let mut jac = M::new_from_sparsity(N, N, Some(eqn.sparsity.clone()), ctx);
        Rhs(eqn).jacobian_inplace(&V::from_slice(x, ctx), t, &mut jac);
        let (idx, vals) = jac.triplet_iter();
        let mut pairs: Vec<(usize, usize)> = idx.collect();
        let mut values: Vec<f64> = pairs
            .iter()
            .zip(vals)
            .map(|(&(r, k), v)| f64::from(r == k) - c * v)
            .collect();
        for &p in extra {
            pairs.push(p);
            values.push(0.25);
        }
        let a = M::try_from_triplets(N, N, pairs.clone(), values, ctx).unwrap();
        let sparsity =
            <<M as Matrix>::Sparsity as MatrixSparsity<M>>::try_from_indices(N, N, pairs).unwrap();
        Fixed(a, sparsity, ctx)
    }

    /// The solutions of A x = b by `TreeLu` and by faer, for a fixed b,
    /// and A^T x = b by both; and whether `TreeLu` used faer.
    fn both_solves(op: &Fixed) -> ([Vec<f64>; 4], bool) {
        let mut tree = TreeLu::default();
        tree.set_sparsity(op);
        tree.set_linearisation(op);
        let mut faer = FaerSparseLU::<f64>::default();
        faer.set_sparsity(op);
        faer.set_linearisation(op);
        // a fixed b with entries of both signs and several magnitudes
        let b: Vec<f64> = (0..N)
            .map(|k| ((k * 7919) % 101) as f64 / 7.0 - 7.0)
            .collect();
        let solve = |transpose: bool, by_tree: bool| {
            let mut x = V::from_slice(&b, op.2);
            match (by_tree, transpose) {
                (true, false) => tree.solve_in_place(&mut x),
                (true, true) => tree.solve_transpose_in_place(&mut x),
                (false, false) => faer.solve_in_place(&mut x),
                (false, true) => faer.solve_transpose_in_place(&mut x),
            }
            .unwrap();
            x.clone_as_vec()
        };
        let xs = [
            solve(false, true),
            solve(false, false),
            solve(true, true),
            solve(true, false),
        ];
        (xs, tree.use_faer)
    }

    /// max |a - b| / max |b|.
    fn rel_diff(a: &[f64], b: &[f64]) -> f64 {
        let scale = b.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        a.iter()
            .zip(b)
            .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
            / scale
    }

    #[test]
    fn the_tree_solve_matches_faer_on_real_newton_matrices() {
        // the resting chain, and the hyperpolarised-soma one mid-step, from
        // a short step's c to a long one's
        let mut active = state();
        active[..45]
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v = -40.0 + i as f64);
        for (rhs, x, t) in [(model(), state(), 1.0), (step_model(-5000.0), active, 6.0)] {
            for c in [1e-4, 1e-2, 0.3] {
                let op = newton_matrix(rhs.clone(), &x, t, c, &[]);
                let ([tree, faer, tree_t, faer_t], used_faer) = both_solves(&op);
                assert!(!used_faer, "fell back at c = {c}");
                assert!(
                    rel_diff(&tree, &faer) < 1e-12,
                    "c = {c}: {}",
                    rel_diff(&tree, &faer)
                );
                assert!(rel_diff(&tree_t, &faer_t) < 1e-12, "c = {c}: transposed");
            }
        }
    }

    #[test]
    fn an_off_pattern_matrix_goes_to_faer() {
        // a gate coupled to another compartment's V, and a V cycle
        for extra in [(45, 1), (0, 2)] {
            let op = newton_matrix(model(), &state(), 1.0, 1e-2, &[extra, (extra.1, extra.0)]);
            let before = tree_lu_fallbacks();
            let ([tree, faer, tree_t, faer_t], used_faer) = both_solves(&op);
            assert!(used_faer, "{extra:?} was taken as the pattern");
            assert!(tree_lu_fallbacks() > before);
            assert_eq!((tree, tree_t), (faer, faer_t));
        }
    }

    #[test]
    fn a_vanishing_tree_pivot_goes_to_faer() {
        // the leaf compartment 44: no V diagonal and no gate feedback, so
        // its pivot is 0 while faer's pivoting still solves
        let mut op = newton_matrix(model(), &state(), 1.0, 1e-2, &[]);
        let (idx, vals) = op.0.triplet_iter();
        let pairs: Vec<(usize, usize)> = idx.collect();
        let values: Vec<f64> = pairs
            .iter()
            .zip(vals)
            .map(|(&(r, c), v)| {
                if r == 44 && (c == 44 || c >= 45) {
                    0.0
                } else {
                    v
                }
            })
            .collect();
        op.0 = M::try_from_triplets(N, N, pairs, values, op.2).unwrap();
        let ([tree, faer, ..], used_faer) = both_solves(&op);
        assert!(used_faer);
        assert_eq!(tree, faer);
    }
}
