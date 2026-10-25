//! The Python binding of the Rust right-hand side of msoaxon's
//! 45-compartment model.
//!
//! The model itself is in [`core`] (plain Rust, no pyo3, usable from a Rust
//! integrator). This file only converts the numpy arguments, calls the core
//! constructor and wraps the result as the `MultiRhs` class that scipy's BDF
//! solver calls as `f(t, x)`; msoaxon builds one per `mso_axon` call.

pub mod core;

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;

impl From<core::CoreError> for PyErr {
    fn from(e: core::CoreError) -> PyErr {
        PyValueError::new_err(e.0)
    }
}

fn to_vec(a: &PyReadonlyArray1<'_, f64>) -> Vec<f64> {
    a.as_array().to_vec()
}

/// The right-hand side of one `mso_axon` run: f(t, x) -> dx/dt.
#[pyclass(frozen, module = "msoaxon_fast")]
struct MultiRhs {
    inner: core::MultiRhs,
}

#[pymethods]
impl MultiRhs {
    /// Builds the right-hand side of one run.
    ///
    /// Arrays are per compartment (length n), except g_ax, parent and child
    /// (one entry per axial edge, 0-based). `stim_type` is one of
    /// `multi.MSO_STIM_TYPES` or "none"; the stimulus keywords are the
    /// fields of `_common.stimulus`, with `k`/`k2` the 0-based input
    /// compartments and `epsg` = (1 / peak, decay, rise). `gating` maps
    /// every name declared in `Gating` to its value (`_fast._gating()`).
    #[new]
    #[pyo3(signature = (
        *, g_na, g_kht, g_klt, g_h, g_lk, vlk, na_vhalf, cap, sa, g_ax,
        parent, child, active, gating, stim_type, start=0.0,
        stop=0.0, amp=0.0, k=0, k2=0, f=0.0, epsg=(1.0, 1.0, 1.0),
        t_syn=None, g_syn=None, vsyn=0.0, diff=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        g_na: PyReadonlyArray1<'_, f64>,
        g_kht: PyReadonlyArray1<'_, f64>,
        g_klt: PyReadonlyArray1<'_, f64>,
        g_h: PyReadonlyArray1<'_, f64>,
        g_lk: PyReadonlyArray1<'_, f64>,
        vlk: PyReadonlyArray1<'_, f64>,
        na_vhalf: PyReadonlyArray1<'_, f64>,
        cap: PyReadonlyArray1<'_, f64>,
        sa: PyReadonlyArray1<'_, f64>,
        g_ax: PyReadonlyArray1<'_, f64>,
        parent: PyReadonlyArray1<'_, i64>,
        child: PyReadonlyArray1<'_, i64>,
        active: [bool; 6],
        gating: HashMap<String, f64>,
        stim_type: &str,
        start: f64,
        stop: f64,
        amp: f64,
        k: i64,
        k2: i64,
        f: f64,
        epsg: (f64, f64, f64),
        t_syn: Option<PyReadonlyArray1<'_, f64>>,
        g_syn: Option<PyReadonlyArray1<'_, f64>>,
        vsyn: f64,
        diff: Option<f64>,
    ) -> PyResult<Self> {
        let params = core::Params {
            g_na: to_vec(&g_na),
            g_kht: to_vec(&g_kht),
            g_klt: to_vec(&g_klt),
            g_h: to_vec(&g_h),
            g_lk: to_vec(&g_lk),
            vlk: to_vec(&vlk),
            na_vhalf: to_vec(&na_vhalf),
            cap: to_vec(&cap),
            sa: to_vec(&sa),
            g_ax: to_vec(&g_ax),
            parent: parent.as_array().to_vec(),
            child: child.as_array().to_vec(),
            active,
            gating: core::Gating::from_map(&gating)?,
        };
        let spec = core::StimSpec {
            stim_type: stim_type.to_string(),
            start,
            stop,
            amp,
            k,
            k2,
            f,
            epsg,
            t_syn: t_syn.as_ref().map(to_vec),
            g_syn: g_syn.as_ref().map(to_vec),
            vsyn,
            diff,
        };
        Ok(MultiRhs {
            inner: core::MultiRhs::new(params, &spec)?,
        })
    }

    /// Returns dx/dt at time t for the state x (length 7n) as a new array.
    ///
    /// A contiguous x is read in place; a strided view (scipy's num_jac
    /// passes columns of a 2-D array) is copied first.
    fn __call__<'py>(
        &self,
        py: Python<'py>,
        t: f64,
        x: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let len = x.len()?;
        if len != 7 * self.inner.n() {
            return Err(PyValueError::new_err(format!(
                "x has {len} entries, expected {}",
                7 * self.inner.n()
            )));
        }
        let mut out = vec![0.0f64; len];
        match x.as_slice() {
            Ok(xs) => self.inner.eval(t, xs, &mut out),
            Err(_) => self.inner.eval(t, &x.as_array().to_vec(), &mut out),
        }
        Ok(PyArray1::from_vec(py, out))
    }

    /// The number of compartments.
    #[getter]
    fn n(&self) -> usize {
        self.inner.n()
    }
}

/// `Table::interp` at one point: `np.interp(x, xs, ys, 0.0, 0.0)`, the
/// contract of `synaptic.interp_g`. Exposed so the tests can hold the port
/// of numpy's interpolation to numpy (tests/test_fast.py).
#[pyfunction]
fn interp(xs: PyReadonlyArray1<'_, f64>, ys: PyReadonlyArray1<'_, f64>, x: f64) -> PyResult<f64> {
    let table = core::Table {
        x: to_vec(&xs),
        y: to_vec(&ys),
    };
    if table.x.is_empty() || table.x.len() != table.y.len() {
        return Err(PyValueError::new_err(
            "xs and ys must be non-empty and the same length",
        ));
    }
    Ok(table.interp(x))
}

/// Experiment: runs the diffsol sparse Faer BDF solver on a stiff 45-state
/// chain (x_i' = -k_i x_i + c (x_{i-1} - x_i)) to `t_end` and returns
/// (x[0] at t_end, number of steps taken). Only here so the linker keeps
/// diffsol's code while the packaging cost is measured.
#[pyfunction]
fn diffsol_probe(t_end: f64) -> PyResult<(f64, usize)> {
    use diffsol::{DenseMatrix, FaerSparseLU, FaerSparseMat, OdeBuilder, OdeSolverMethod};
    const N: usize = 45;
    let err = |e: diffsol::DiffsolError| PyValueError::new_err(e.to_string());
    let problem = OdeBuilder::<FaerSparseMat<f64>>::new()
        .rhs_implicit(
            |x, _p, _t, y| {
                for i in 0..N {
                    let k = 1.0 + 100.0 * i as f64;
                    let left = if i > 0 { x[i - 1] - x[i] } else { 0.0 };
                    y[i] = -k * x[i] + left;
                }
            },
            |_x, _p, _t, v, y| {
                for i in 0..N {
                    let k = 1.0 + 100.0 * i as f64;
                    let left = if i > 0 { v[i - 1] - v[i] } else { 0.0 };
                    y[i] = -k * v[i] + left;
                }
            },
        )
        .init(|_p, _t, y| y.fill(1.0), N)
        .build()
        .map_err(err)?;
    let mut solver = problem.bdf::<FaerSparseLU<f64>>().map_err(err)?;
    let (ys, ts, _) = solver.solve(t_end).map_err(err)?;
    let last = ts.len().saturating_sub(1);
    Ok((ys.get_index(0, last), ts.len()))
}

/// The Rust right-hand side of msoaxon's 45-compartment model.
#[pymodule]
fn msoaxon_fast(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MultiRhs>()?;
    m.add_function(wrap_pyfunction!(interp, m)?)?;
    m.add_function(wrap_pyfunction!(diffsol_probe, m)?)?;
    Ok(())
}
