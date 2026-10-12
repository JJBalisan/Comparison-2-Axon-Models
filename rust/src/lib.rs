//! The right-hand side of msoaxon's 45-compartment model, in Rust.
//!
//! A port of `msoaxon.multi._rhs` (with `external_current`,
//! `bilateral_current`, `axial_current` and the gating functions of
//! `msoaxon.constants`). The Python code is the reference. Every expression
//! keeps numpy's floating-point evaluation order, operation for operation, so
//! that on a platform where numpy and Rust call the same libm the result is
//! bit-identical. Two details matter:
//!
//! * `m**4` is `powf(4.0)` (numpy calls libm `pow`); `powi(4)` drifts.
//! * Nothing is reordered, and nothing is fused except where numpy's own
//!   compiled code is (`np.interp`'s line on arm64, see `line`).
//!
//! msoaxon builds one `MultiRhs` per `mso_axon` call, holding every per-model
//! constant, and scipy's BDF solver calls it as `f(t, x)`.

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::cmp::Ordering;

const V_NA: f64 = 69.0;
const V_K: f64 = -90.0;
const V_H: f64 = -35.0;

/// The unitary EPSG waveform (`synaptic.epsg_unitary`).
///
/// `inv_peak` is `1 / peak` as Python computes it: `1 / 0.21317` for the
/// default kinetics, else from the exact peak.
#[derive(Clone, Copy, Debug)]
struct Epsg {
    inv_peak: f64,
    decay: f64,
    rise: f64,
}

impl Epsg {
    #[inline]
    fn unitary(&self, t: f64) -> f64 {
        let pos = if t > 0.0 { 1.0 } else { 0.0 };
        (self.inv_peak * pos) * ((-t / self.decay).exp() - (-t / self.rise).exp())
    }
}

/// A conductance table, interpolated like `np.interp(t, x, y, 0.0, 0.0)`.
#[derive(Clone, Debug)]
struct Table {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl Table {
    /// numpy's `arr_interp` for one point, branch for branch.
    fn interp(&self, t: f64) -> f64 {
        let (dx, dy) = (&self.x, &self.y);
        let len = dx.len();
        if len == 1 {
            // NaN compares neither way and gives dy[0], as in numpy
            return match t.partial_cmp(&dx[0]) {
                Some(Ordering::Less | Ordering::Greater) => 0.0,
                _ => dy[0],
            };
        }
        if t.is_nan() {
            return t;
        }
        if t > dx[len - 1] || t < dx[0] {
            return 0.0;
        }
        // the largest j with dx[j] <= t (dx is increasing)
        let j = dx.partition_point(|&v| v <= t) - 1;
        if j == len - 1 || dx[j] == t {
            return dy[j];
        }
        let slope = (dy[j + 1] - dy[j]) / (dx[j + 1] - dx[j]);
        let mut res = line(slope, t - dx[j], dy[j]);
        if res.is_nan() {
            res = line(slope, t - dx[j + 1], dy[j + 1]);
            if res.is_nan() && dy[j] == dy[j + 1] {
                res = dy[j];
            }
        }
        res
    }
}

/// `slope * dx + y` as numpy's compiled `arr_interp` evaluates it.
///
/// The C compiler contracts that expression into one fused multiply-add
/// where the target has one in its baseline: numpy's arm64 wheels (clang)
/// do, its x86_64 wheels (no FMA before x86-64-v3) don't. Measured on
/// darwin-arm64 with numpy 2.5.3: the fused form matched np.interp on
/// 20,000 of 20,000 points, the unfused one differed on 14%.
#[inline]
fn line(slope: f64, dx: f64, y: f64) -> f64 {
    if cfg!(target_arch = "aarch64") {
        slope.mul_add(dx, y)
    } else {
        slope * dx + y
    }
}

/// One stimulus, with what `_common.stimulus` stores for it.
#[derive(Clone, Debug)]
enum Stim {
    None,
    Ramp {
        start: f64,
        amp: f64,
        k: usize,
        slope: f64,
        t_top: f64,
        t_end_local: f64,
    },
    Ramp2 {
        start: f64,
        amp: f64,
        k: usize,
        slope: f64,
    },
    Step {
        start: f64,
        stop: f64,
        amp: f64,
    },
    Sine {
        start: f64,
        stop: f64,
        amp: f64,
        k: usize,
        two_pi_f: f64,
    },
    Synaptic {
        start: f64,
        stop: f64,
        k: usize,
        table: Table,
        vsyn: f64,
        diff: Option<f64>,
    },
    Epsg {
        start: f64,
        stop: f64,
        amp: f64,
        k: usize,
        w: Epsg,
    },
    EpsgPair {
        start: f64,
        td: f64,
        amp: f64,
        k: usize,
        w: Epsg,
    },
    Bilateral {
        start: f64,
        td: f64,
        amp: f64,
        a: usize,
        b: usize,
        w: Epsg,
    },
}

impl Stim {
    /// Adds the input current density to `total`, as `_rhs` does.
    ///
    /// Like the Python, a stimulus that is off still adds an exact 0.0 to
    /// the soma.
    fn add_to(&self, t: f64, v: &[f64], sa: &[f64], total: &mut [f64]) {
        let (k, iext) = match *self {
            Stim::None => (0, 0.0),
            Stim::Ramp {
                start,
                amp,
                k,
                slope,
                t_top,
                t_end_local,
            } => {
                let on = if start <= t && t <= t_top { 1.0 } else { 0.0 };
                let i0 = on * amp * (t - start) * slope / 1000.0;
                if 5.0 < t && t <= t_end_local {
                    (k, -i0 * 1e3 / sa[0])
                } else {
                    (0, 0.0)
                }
            }
            Stim::Ramp2 {
                start,
                amp,
                k,
                slope,
            } => {
                let on = if t >= start { 1.0 } else { 0.0 };
                let a = on * amp * (t - start) * slope / 1000.0;
                let b = amp / 1000.0;
                let i0 = if b < a { b } else { a }; // Python's min(a, b)
                (k, -i0 * 1e3 / sa[0])
            }
            Stim::Step { start, stop, amp } => {
                if start < t && t <= stop {
                    (0, -amp / sa[0])
                } else {
                    (0, 0.0)
                }
            }
            Stim::Sine {
                start,
                stop,
                amp,
                k,
                two_pi_f,
            } => {
                let wave = (two_pi_f * (t - start) / 1000.0).sin();
                if start < t && t <= stop {
                    let pos = if wave > 0.0 { 1.0 } else { 0.0 };
                    (k, -(amp * wave * pos) / sa[0])
                } else {
                    (0, 0.0)
                }
            }
            Stim::Synaptic {
                start,
                stop,
                k,
                ref table,
                vsyn,
                diff,
            } => {
                if start < t && t <= stop {
                    let g = match diff {
                        None => table.interp(t),
                        Some(d) => table.interp(t) + table.interp(t + d),
                    };
                    (k, g * (v[0] - vsyn) / sa[0])
                } else {
                    (0, 0.0)
                }
            }
            Stim::Epsg {
                start,
                stop,
                amp,
                k,
                w,
            } => {
                let te = t - start;
                if start < t && t <= stop {
                    let on = if te >= 0.0 { 1.0 } else { 0.0 };
                    (k, amp * (0.0 - v[k]) * on * w.unitary(te) / -sa[k])
                } else {
                    (0, 0.0)
                }
            }
            Stim::EpsgPair {
                start,
                td,
                amp,
                k,
                w,
            } => {
                let te = t - start;
                let on0 = if te >= 0.0 { 1.0 } else { 0.0 };
                let on1 = if te >= td { 1.0 } else { 0.0 };
                let wave = on0 * w.unitary(te) + on1 * w.unitary(te - td);
                (k, amp * (0.0 - v[k]) * wave / -sa[k])
            }
            Stim::Bilateral {
                start,
                td,
                amp,
                a,
                b,
                w,
            } => {
                let te = t - start;
                let on0 = if te >= 0.0 { 1.0 } else { 0.0 };
                let on1 = if te >= td { 1.0 } else { 0.0 };
                let ia = amp * (0.0 - v[a]) * on0 * w.unitary(te) / -sa[a];
                let ib = amp * (0.0 - v[b]) * on1 * w.unitary(te - td) / -sa[b];
                total[a] += ia;
                total[b] += ib;
                return;
            }
        };
        total[k] += iext;
    }
}

/// The right-hand side of one `mso_axon` run: f(t, x) -> dx/dt.
#[pyclass(frozen, module = "msoaxon_fast")]
struct MultiRhs {
    n: usize,
    g_na: Vec<f64>,
    g_kht: Vec<f64>,
    g_klt: Vec<f64>,
    g_h: Vec<f64>,
    g_lk: Vec<f64>,
    vlk: Vec<f64>,
    na_vhalf: Vec<f64>,
    cap: Vec<f64>,
    sa: Vec<f64>,
    g_ax: Vec<f64>,
    parent: Vec<usize>,
    child: Vec<usize>,
    active: [bool; 6],
    p_temp: f64,
    a_temp: f64,
    stim: Stim,
}

impl MultiRhs {
    /// Writes dx/dt for state `x` at time `t` into `out` (both 7n long).
    fn eval(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let n = self.n;
        let (v, rest) = x.split_at(n);
        let (m, rest) = rest.split_at(n);
        let (h, rest) = rest.split_at(n);
        let (p, rest) = rest.split_at(n);
        let (w, rest) = rest.split_at(n);
        let (z, a) = rest.split_at(n);

        let mut total = vec![0.0f64; n];
        for i in 0..n {
            let vi = v[i];
            let ina = self.g_na[i] * m[i].powf(4.0) * (0.993 * h[i] + 0.007) * (vi - V_NA);
            let ikht = self.g_kht[i] * p[i] * (vi - V_K);
            let iklt = self.g_klt[i] * w[i].powf(4.0) * z[i] * (vi - V_K);
            let ih = self.g_h[i] * a[i] * (vi - V_H);
            let ilk = self.g_lk[i] * (vi - self.vlk[i]);
            total[i] = ina + ikht + iklt + ih + ilk;
        }
        self.stim.add_to(t, v, &self.sa, &mut total);

        // axial_current: np.subtract.at, then np.add.at, both in edge order
        // (for the chain this is msoAxon.m's I[:-1] -= flow; I[1:] += flow)
        let mut axial = vec![0.0f64; n];
        for (e, (&pa, &ch)) in self.parent.iter().zip(&self.child).enumerate() {
            axial[pa] -= self.g_ax[e] * (v[pa] - v[ch]);
        }
        for (e, (&pa, &ch)) in self.parent.iter().zip(&self.child).enumerate() {
            axial[ch] += self.g_ax[e] * (v[pa] - v[ch]);
        }

        let (dv, rest) = out.split_at_mut(n);
        let (dm, rest) = rest.split_at_mut(n);
        let (dh, rest) = rest.split_at_mut(n);
        let (dp, rest) = rest.split_at_mut(n);
        let (dw, rest) = rest.split_at_mut(n);
        let (dz, da) = rest.split_at_mut(n);
        for i in 0..n {
            let ax = -axial[i] * 1e9 / self.sa[i];
            dv[i] = -(total[i] + ax) / self.cap[i];
        }
        let [act_m, act_h, act_p, act_w, act_z, act_a] = self.active;
        for i in 0..n {
            let vi = v[i];
            dm[i] = if act_m {
                let minf = 1.0 / (1.0 + ((vi + 46.0) / -11.0).exp());
                let taum = (0.141 + (-0.0826 / (1.0 + ((-20.5 - vi) / 10.8).exp()))) / 3.0;
                (minf - m[i]) / taum
            } else {
                0.0
            };
            dh[i] = if act_h {
                let hinf = 1.0 / (1.0 + ((vi + self.na_vhalf[i]) / 7.77).exp());
                let tauh = (4.0 + (-3.74 / (1.0 + ((-40.6 - vi) / 5.05).exp()))) / 3.0;
                (hinf - h[i]) / tauh
            } else {
                0.0
            };
            dp[i] = if act_p {
                let pinf = 1.0 / (1.0 + (-(vi + 23.0) / 6.0).exp());
                let taup = self.p_temp
                    * (100.0
                        / (4.0 * ((vi + 60.0) / 32.0).exp() + 5.0 * (-(vi + 60.0) / 22.0).exp())
                        + 5.0);
                (pinf - p[i]) / taup
            } else {
                0.0
            };
            dw[i] = if act_w {
                let winf = 1.0 / (1.0 + ((vi + 57.34) / -11.7).exp());
                let tauw = 21.5
                    / (6.0 * ((vi + 60.0) / 7.0).exp() + 24.0 * (-(vi + 60.0) / 50.6).exp())
                    + 0.35;
                (winf - w[i]) / tauw
            } else {
                0.0
            };
            dz[i] = if act_z {
                let zinf = (1.0 - 0.27) / (1.0 + ((vi + 67.0) / 6.16).exp()) + 0.27;
                let tauz =
                    170.0 / (5.0 * ((vi + 60.0) / 10.0).exp() + (-(vi + 70.0) / 8.0).exp()) + 10.7;
                (zinf - z[i]) / tauz
            } else {
                0.0
            };
            da[i] = if act_a {
                let ainf = 1.0 / (1.0 + (0.1 * (vi + 80.4)).exp());
                let s = vi + 61.5;
                let taua = self.a_temp * (79.0 + 417.0 * (-(s * s) / 800.0).exp());
                (ainf - a[i]) / taua
            } else {
                0.0
            };
        }
    }
}

fn to_vec(a: &PyReadonlyArray1<'_, f64>) -> Vec<f64> {
    a.as_array().to_vec()
}

fn check_len(name: &str, v: &[f64], n: usize) -> PyResult<()> {
    if v.len() != n {
        return Err(PyValueError::new_err(format!(
            "{name} has {} entries, expected {n}",
            v.len()
        )));
    }
    Ok(())
}

fn index(name: &str, i: i64, n: usize) -> PyResult<usize> {
    usize::try_from(i)
        .ok()
        .filter(|&u| u < n)
        .ok_or_else(|| PyValueError::new_err(format!("{name}={i} is not in 0..{n}")))
}

#[pymethods]
impl MultiRhs {
    /// Builds the right-hand side of one run.
    ///
    /// Arrays are per compartment (length n), except g_ax, parent and child
    /// (one entry per axial edge, 0-based). `stim_type` is one of
    /// `multi.MSO_STIM_TYPES` or "none"; the stimulus keywords are the
    /// fields of `_common.stimulus`, with `k`/`k2` the 0-based input
    /// compartments and `epsg` = (1 / peak, decay, rise).
    #[new]
    #[pyo3(signature = (
        *, g_na, g_kht, g_klt, g_h, g_lk, vlk, na_vhalf, cap, sa, g_ax,
        parent, child, active, p_temp, a_temp, stim_type, start=0.0,
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
        p_temp: f64,
        a_temp: f64,
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
        let n = g_na.len()?;
        let arrays = [
            ("g_kht", to_vec(&g_kht)),
            ("g_klt", to_vec(&g_klt)),
            ("g_h", to_vec(&g_h)),
            ("g_lk", to_vec(&g_lk)),
            ("vlk", to_vec(&vlk)),
            ("na_vhalf", to_vec(&na_vhalf)),
            ("cap", to_vec(&cap)),
            ("sa", to_vec(&sa)),
        ];
        for (name, v) in &arrays {
            check_len(name, v, n)?;
        }
        let [g_kht, g_klt, g_h, g_lk, vlk, na_vhalf, cap, sa] = arrays.map(|(_, v)| v);
        let g_ax = to_vec(&g_ax);
        let parent = parent.as_array().to_vec();
        let child = child.as_array().to_vec();
        if parent.len() != g_ax.len() || child.len() != g_ax.len() {
            return Err(PyValueError::new_err(
                "g_ax, parent and child must have one entry per edge",
            ));
        }
        let parent = parent
            .into_iter()
            .map(|i| index("parent", i, n))
            .collect::<PyResult<Vec<_>>>()?;
        let child = child
            .into_iter()
            .map(|i| index("child", i, n))
            .collect::<PyResult<Vec<_>>>()?;

        let w = Epsg {
            inv_peak: epsg.0,
            decay: epsg.1,
            rise: epsg.2,
        };
        let slope = 1.0 / (stop - start);
        let stim = match stim_type {
            "none" => Stim::None,
            "ramp" => Stim::Ramp {
                start,
                amp,
                k: index("k", k, n)?,
                slope,
                t_top: start + 1.0 / slope,
                t_end_local: stop + 5.0,
            },
            "ramp2" => Stim::Ramp2 {
                start,
                amp,
                k: index("k", k, n)?,
                slope,
            },
            "step" => Stim::Step { start, stop, amp },
            "sine" => Stim::Sine {
                start,
                stop,
                amp,
                k: index("k", k, n)?,
                two_pi_f: 2.0 * std::f64::consts::PI * f,
            },
            "Synaptic" | "SynapticPair" => {
                let (Some(x), Some(y)) = (t_syn, g_syn) else {
                    return Err(PyValueError::new_err(format!(
                        "{stim_type} needs t_syn and g_syn"
                    )));
                };
                let table = Table {
                    x: to_vec(&x),
                    y: to_vec(&y),
                };
                if table.x.is_empty() || table.x.len() != table.y.len() {
                    return Err(PyValueError::new_err(
                        "t_syn and g_syn must be non-empty and the same length",
                    ));
                }
                let diff = if stim_type == "SynapticPair" {
                    Some(diff.ok_or_else(|| PyValueError::new_err("SynapticPair needs diff"))?)
                } else {
                    None
                };
                Stim::Synaptic {
                    start,
                    stop,
                    k: index("k", k, n)?,
                    table,
                    vsyn,
                    diff,
                }
            }
            "EPSG" => Stim::Epsg {
                start,
                stop,
                amp,
                k: index("k", k, n)?,
                w,
            },
            "EPSGpair" => Stim::EpsgPair {
                start,
                td: stop - start,
                amp,
                k: index("k", k, n)?,
                w,
            },
            "EPSGbilateral" => Stim::Bilateral {
                start,
                td: stop - start,
                amp,
                a: index("k", k, n)?,
                b: index("k2", k2, n)?,
                w,
            },
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown stim_type {other:?}"
                )))
            }
        };
        Ok(MultiRhs {
            n,
            g_na: to_vec(&g_na),
            g_kht,
            g_klt,
            g_h,
            g_lk,
            vlk,
            na_vhalf,
            cap,
            sa,
            g_ax,
            parent,
            child,
            active,
            p_temp,
            a_temp,
            stim,
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
        if len != 7 * self.n {
            return Err(PyValueError::new_err(format!(
                "x has {len} entries, expected {}",
                7 * self.n
            )));
        }
        let mut out = vec![0.0f64; len];
        match x.as_slice() {
            Ok(xs) => self.eval(t, xs, &mut out),
            Err(_) => self.eval(t, &x.as_array().to_vec(), &mut out),
        }
        Ok(PyArray1::from_vec(py, out))
    }

    /// The number of compartments.
    #[getter]
    fn n(&self) -> usize {
        self.n
    }
}

/// The Rust right-hand side of msoaxon's 45-compartment model.
#[pymodule]
fn msoaxon_fast(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MultiRhs>()?;
    Ok(())
}
