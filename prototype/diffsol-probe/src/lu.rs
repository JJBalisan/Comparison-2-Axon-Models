//! Direct timing of `FaerSparseLU` on the BDF iteration matrix, to see
//! whether a branched pattern causes fill-in trouble.

use std::time::Instant;

use diffsol::{FaerSparseLU, FaerSparseMat, LinearSolver, Matrix, MatrixCommon, MatrixOp, Vector};

type M = FaerSparseMat<f64>;

/// Median timings (us) of `FaerSparseLU` on `J - I / c`, the BDF matrix
/// `-(I - c J) / c` for `c = h gamma`.
pub struct LuTiming {
    /// `set_problem`: COLAMD ordering, elimination tree, allocation.
    pub symbolic_us: f64,
    /// `set_linearisation`: copy of the matrix + numeric factorisation.
    pub numeric_us: f64,
    /// One triangular solve.
    pub solve_us: f64,
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn solve_in_place<LS: LinearSolver<M>>(ls: &LS, b: &mut <M as MatrixCommon>::V) {
    ls.solve_in_place(b).unwrap();
}

/// Time the LU of `jac - I / c` (`jac` must hold its diagonal).
pub fn time_lu(jac: &M, c: f64, reps: usize) -> LuTiming {
    let (idx, val) = jac.triplet_iter();
    let idx: Vec<(usize, usize)> = idx.collect();
    let val: Vec<f64> = idx
        .iter()
        .zip(val)
        .map(|(&(i, j), v)| if i == j { v - 1.0 / c } else { v })
        .collect();
    let n = jac.nrows();
    let ctx = *jac.context();
    let a = M::try_from_triplets(n, n, idx, val, ctx).unwrap();
    let op = MatrixOp::new(a);
    let x = <M as MatrixCommon>::V::zeros(n, ctx);
    let ones = <M as MatrixCommon>::V::from_vec(vec![1.0; n], ctx);
    let (mut sym, mut num, mut sol) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..reps {
        let mut ls = FaerSparseLU::<f64>::default();
        let t0 = Instant::now();
        ls.set_problem(&op);
        sym.push(t0.elapsed().as_secs_f64());
        // the first factorisation sizes faer's buffers; time the reuse too
        for _ in 0..2 {
            let t1 = Instant::now();
            ls.set_linearisation(&op, &x, 0.0);
            num.push(t1.elapsed().as_secs_f64());
        }
        let mut b = ones.clone();
        let t2 = Instant::now();
        solve_in_place(&ls, &mut b);
        sol.push(t2.elapsed().as_secs_f64());
        std::hint::black_box(b);
    }
    LuTiming {
        symbolic_us: median(sym) * 1e6,
        numeric_us: median(num) * 1e6,
        solve_us: median(sol) * 1e6,
    }
}
