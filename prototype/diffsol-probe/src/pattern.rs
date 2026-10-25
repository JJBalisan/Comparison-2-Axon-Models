//! The Jacobian pattern the real port computes, copied verbatim from
//! `rust/src/core.rs` on origin/main (`chain_jac_sparsity`,
//! `tree_jac_sparsity`, the Rust twin of `multi._tree_jac_sparsity`), plus
//! the conversion into diffsol's sparsity type.

use diffsol::matrix::sparsity::MatrixSparsity;
use diffsol::{DiffsolError, Matrix};

/// The Jacobian pattern of an unbranched chain of `n` compartments, as
/// `multi._chain_jac_sparsity` (`n` = 45 there): see [`tree_jac_sparsity`].
pub fn chain_jac_sparsity(n: usize) -> (Vec<usize>, Vec<usize>) {
    let parent: Vec<usize> = (0..n.saturating_sub(1)).collect();
    let child: Vec<usize> = (1..n).collect();
    tree_jac_sparsity(n, &parent, &child)
}

/// The Jacobian pattern of a tree, as `multi._tree_jac_sparsity`.
///
/// The state is V then the six gates, 7n long. dV couples to its own V, to
/// each V it shares an edge (`parent[e]`, `child[e]`) with, and to its own
/// gates; each gate to its own V and itself. Returns `(rows, cols)` in
/// column-major (CSC) order: columns ascending, rows ascending within a
/// column, each entry once.
pub fn tree_jac_sparsity(n: usize, parent: &[usize], child: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let mut nbrs: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    for (&p, &c) in parent.iter().zip(child) {
        nbrs[p].push(c);
        nbrs[c].push(p);
    }
    let (mut rows, mut cols) = (Vec::new(), Vec::new());
    for (j, nb) in nbrs.iter_mut().enumerate() {
        nb.sort_unstable();
        nb.dedup();
        // column j of V: the V rows it couples to, then its six gates
        for &r in nb.iter() {
            rows.push(r);
            cols.push(j);
        }
        for g in 1..=6 {
            rows.push(g * n + j);
            cols.push(j);
        }
    }
    for g in 1..=6 {
        for j in 0..n {
            // column of gate g: dV of the same compartment, then itself
            rows.extend([j, g * n + j]);
            cols.extend([g * n + j; 2]);
        }
    }
    (rows, cols)
}

/// `(rows, cols)` as the `(row, col)` pairs diffsol's sparsity and
/// colouring constructors take (same order, so CSC order is kept).
pub fn to_pairs(rows: &[usize], cols: &[usize]) -> Vec<(usize, usize)> {
    rows.iter().copied().zip(cols.iter().copied()).collect()
}

/// The conversion lane D needs: `(rows, cols)` into `M::Sparsity` (for
/// `FaerSparseMat` that is faer's `SymbolicSparseColMat<usize>`), through
/// the generic `MatrixSparsity::try_from_indices`.
pub fn to_sparsity<M: Matrix>(
    n: usize,
    rows: &[usize],
    cols: &[usize],
) -> Result<M::Sparsity, DiffsolError> {
    Ok(<M::Sparsity as MatrixSparsity<M>>::try_from_indices(
        n,
        n,
        to_pairs(rows, cols),
    )?)
}
