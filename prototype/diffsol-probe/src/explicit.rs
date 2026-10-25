//! A hand-written `OdeEquations` for the toy model whose Jacobian sparsity
//! is supplied, not detected: the path lane D would take to hand diffsol
//! `core.rs::tree_jac_sparsity` (the single source of truth shared with
//! Python's `multi._tree_jac_sparsity`).
//!
//! diffsol 0.18.0 needs, for `OdeBuilder::build_from_eqn` plus
//! `OdeSolverProblem::bdf_solver` (bound `OdeEquationsImplicit`):
//!
//! - `Op` on the equations and on every operator (sizes, context);
//! - `OdeEquationsRef<'a>` naming the operator types (`Rhs`, `Init`,
//!   `Root`, and `ParameterisedOp<UnitCallable>` for the unused `Mass`,
//!   `Out`, `Reset`) and `OdeEquations` returning them;
//! - `NonLinearOp` + `NonLinearOpJacobian` on the RHS, the latter with
//!   `jacobian_sparsity` returning the supplied pattern and
//!   `jacobian_inplace` filling it through a `JacobianColoring` built from
//!   the same pattern (without it diffsol falls back to N `Jv` calls);
//! - `ConstantOp` on the initial state, `NonLinearOp` on the root function.
//!
//! `OdeEquationsImplicit` itself is a blanket impl.

use diffsol::{
    ConstantOp, DiffsolError, JacobianColoring, Matrix, NonLinearOp, NonLinearOpJacobian,
    OdeEquations, OdeEquationsRef, Op, ParameterisedOp, UnitCallable, Vector,
};

use crate::model::{Jv, Model, N};
use crate::pattern::{to_pairs, to_sparsity};

/// The toy model with a supplied Jacobian pattern.
pub struct ExplicitEqn<'m, M: Matrix<T = f64>> {
    model: &'m Model,
    jv: Jv,
    y0: M::V,
    sparsity: M::Sparsity,
    coloring: JacobianColoring<M>,
    ctx: M::C,
}

impl<'m, M: Matrix<T = f64>> ExplicitEqn<'m, M> {
    /// Equations for `model` from `y0`, with the pattern `(rows, cols)` in
    /// the format `core.rs::tree_jac_sparsity` returns.
    pub fn new(
        model: &'m Model,
        jv: Jv,
        y0: &[f64],
        rows: &[usize],
        cols: &[usize],
    ) -> Result<Self, DiffsolError> {
        let ctx = M::C::default();
        let sparsity = to_sparsity::<M>(N, rows, cols)?;
        let coloring = JacobianColoring::new(&sparsity, &to_pairs(rows, cols), ctx.clone());
        Ok(Self {
            model,
            jv,
            y0: M::V::from_vec(y0.to_vec(), ctx.clone()),
            sparsity,
            coloring,
            ctx,
        })
    }
}

/// The right-hand side operator.
pub struct Rhs<'a, 'm, M: Matrix<T = f64>>(&'a ExplicitEqn<'m, M>);
/// The initial-state operator.
pub struct Init<'a, 'm, M: Matrix<T = f64>>(&'a ExplicitEqn<'m, M>);
/// The two root functions.
pub struct Root<'a, 'm, M: Matrix<T = f64>>(&'a ExplicitEqn<'m, M>);

macro_rules! impl_op {
    ($ty:ident, $nout:expr) => {
        impl<M: Matrix<T = f64>> Op for $ty<'_, '_, M> {
            type T = f64;
            type V = M::V;
            type M = M;
            type C = M::C;
            fn context(&self) -> &M::C {
                &self.0.ctx
            }
            fn nstates(&self) -> usize {
                N
            }
            fn nout(&self) -> usize {
                $nout
            }
            fn nparams(&self) -> usize {
                0
            }
        }
    };
}

impl_op!(Rhs, N);
impl_op!(Init, N);
impl_op!(Root, 2);

impl<M: Matrix<T = f64>> NonLinearOp for Rhs<'_, '_, M> {
    fn call_inplace(&self, x: &M::V, t: f64, y: &mut M::V) {
        let m = self.0.model;
        y.for_each_batch_host([x], |y, [x], _| m.rhs(x, t, y));
    }
}

impl<M: Matrix<T = f64>> NonLinearOpJacobian for Rhs<'_, '_, M> {
    fn jac_mul_inplace(&self, x: &M::V, t: f64, v: &M::V, y: &mut M::V) {
        let (m, kind) = (self.0.model, self.0.jv);
        y.for_each_batch_host([x, v], |y, [x, v], _| m.jv_by(kind, x, t, v, y));
    }
    fn jacobian_inplace(&self, x: &M::V, t: f64, y: &mut M) {
        self.0.coloring.jacobian_inplace(self, x, t, y);
    }
    fn jacobian_sparsity(&self) -> Option<M::Sparsity> {
        Some(self.0.sparsity.clone())
    }
}

impl<M: Matrix<T = f64>> ConstantOp for Init<'_, '_, M> {
    fn call_inplace(&self, _t: f64, y: &mut M::V) {
        y.copy_from(&self.0.y0);
    }
}

impl<M: Matrix<T = f64>> NonLinearOp for Root<'_, '_, M> {
    fn call_inplace(&self, x: &M::V, t: f64, y: &mut M::V) {
        let m = self.0.model;
        y.for_each_batch_host([x], |g, [x], _| m.root(x, t, g));
    }
}

impl<M: Matrix<T = f64>> Op for ExplicitEqn<'_, M> {
    type T = f64;
    type V = M::V;
    type M = M;
    type C = M::C;
    fn context(&self) -> &M::C {
        &self.ctx
    }
    fn nstates(&self) -> usize {
        N
    }
    fn nout(&self) -> usize {
        N
    }
    fn nparams(&self) -> usize {
        0
    }
}

impl<'a, M: Matrix<T = f64>> OdeEquationsRef<'a> for ExplicitEqn<'_, M> {
    type Rhs = Rhs<'a, 'a, M>;
    type Init = Init<'a, 'a, M>;
    type Root = Root<'a, 'a, M>;
    type Mass = ParameterisedOp<'a, UnitCallable<M>>;
    type Out = ParameterisedOp<'a, UnitCallable<M>>;
    type Reset = ParameterisedOp<'a, UnitCallable<M>>;
}

impl<M: Matrix<T = f64>> OdeEquations for ExplicitEqn<'_, M> {
    fn rhs(&self) -> Rhs<'_, '_, M> {
        Rhs(self)
    }
    fn mass(&self) -> Option<ParameterisedOp<'_, UnitCallable<M>>> {
        None
    }
    fn root(&self) -> Option<Root<'_, '_, M>> {
        Some(Root(self))
    }
    fn init(&self) -> Init<'_, '_, M> {
        Init(self)
    }
    fn set_params(&mut self, _p: &M::V) {}
    fn get_params(&self, _p: &mut M::V) {}
}
