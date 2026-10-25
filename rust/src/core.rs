//! The right-hand side of msoaxon's 45-compartment model, in plain Rust.
//!
//! No pyo3 here: `lib.rs` is the thin Python binding. A Rust integrator
//! builds a [`MultiRhs`] from [`Params`] and a [`StimSpec`], calls
//! [`MultiRhs::eval`] and asks [`chain_jac_sparsity`] for the Jacobian
//! pattern.
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
//! This file has no model number of its own: the reversal potentials and the
//! gating kinetics' constants come from the caller (see `Gating`; Python
//! passes them from `_fast._gating()`), so a change to `multi.py` or
//! `constants.py` reaches both implementations. What stays here is
//! structure: the exponents (`m**4`), the `1.0`s and the unit conversions of
//! `_rhs`.

use std::cmp::Ordering;
use std::collections::HashMap;

/// Declares `Gating`, one `f64` per name, and its constructor from the dict
/// that `_fast._gating()` passes. Each field is named for the
/// `msoaxon.constants` function (or `multi` constant) its number belongs to;
/// the Python side lists the same names with the values and their sources.
macro_rules! gating {
    ($($(#[$doc:meta])* $field:ident),* $(,)?) => {
        /// The constants of the currents and the gating kinetics.
        #[derive(Clone, Debug)]
        pub struct Gating {
            $($(#[$doc])* pub $field: f64,)*
        }

        impl Gating {
            /// Reads every field from `map`; a missing or unknown name is an
            /// error, so a typo on either side cannot go unnoticed.
            pub fn from_map(map: &HashMap<String, f64>) -> Result<Self, CoreError> {
                const NAMES: &[&str] = &[$(stringify!($field)),*];
                if let Some(name) = map.keys().find(|k| !NAMES.contains(&k.as_str())) {
                    return Err(CoreError::new(format!(
                        "unknown gating constant {name:?}"
                    )));
                }
                let get = |name: &str| {
                    map.get(name).copied().ok_or_else(|| {
                        CoreError::new(format!("missing gating constant {name:?}"))
                    })
                };
                Ok(Gating {
                    $($field: get(stringify!($field))?,)*
                })
            }
        }
    };
}

gating! {
    /// `multi.V_NA`, `V_K`, `V_H`: the reversal potentials [mV].
    v_na, v_k, v_h,
    /// `na_open`: the inactivating fraction and the floor, 0.993 h + 0.007.
    na_open_h, na_open_floor,
    /// `minf`: (V + v) / k.
    minf_v, minf_k,
    /// `taum`: (a + b / (1 + exp((v - V) / k))) / div.
    taum_a, taum_b, taum_v, taum_k, taum_div,
    /// `hinf`: (V + vhalf) / k, with vhalf per compartment.
    hinf_k,
    /// `tauh`: (a + b / (1 + exp((v - V) / k))) / div.
    tauh_a, tauh_b, tauh_v, tauh_k, tauh_div,
    /// `pinf`: 1 / (1 + exp(-(V + v) / k)).
    pinf_v, pinf_k,
    /// `taup`: temp * (num / (a exp((V + v) / ka) + b exp(-(V + v) / kb)) + c).
    taup_temp, taup_num, taup_a, taup_ka, taup_b, taup_kb, taup_v, taup_c,
    /// `winf`: (V + v) / k.
    winf_v, winf_k,
    /// `tauw`: num / (a exp((V + v) / ka) + b exp(-(V + v) / kb)) + c.
    tauw_num, tauw_a, tauw_ka, tauw_b, tauw_kb, tauw_v, tauw_c,
    /// `zinf`: (1 - floor) / (1 + exp((V + v) / k)) + floor.
    zinf_floor, zinf_v, zinf_k,
    /// `tauz`: num / (a exp((V + va) / ka) + exp(-(V + vb) / kb)) + c.
    tauz_num, tauz_a, tauz_va, tauz_ka, tauz_vb, tauz_kb, tauz_c,
    /// `ainf`: 1 / (1 + exp(k (V + v))).
    ainf_k, ainf_v,
    /// `taua`: temp * (c + a exp(-(V + v)^2 / k)).
    taua_temp, taua_c, taua_a, taua_v, taua_k,
}

/// The unitary EPSG waveform (`synaptic.epsg_unitary`).
///
/// `inv_peak` is `1 / peak` as Python computes it: `1 / 0.21317` for the
/// default kinetics, else from the exact peak.
#[derive(Clone, Copy, Debug)]
pub struct Epsg {
    pub inv_peak: f64,
    pub decay: f64,
    pub rise: f64,
}

impl Epsg {
    #[inline]
    fn unitary(&self, t: f64) -> f64 {
        let pos = if t > 0.0 { 1.0 } else { 0.0 };
        // Far before onset the exponentials overflow and 0 * inf is NaN;
        // like the Python, they see t no lower than -1 (the product is
        // already 0 there, and t > -1 keeps its bits). NaN passes through.
        let te = if t < -1.0 { -1.0 } else { t };
        (self.inv_peak * pos) * ((-te / self.decay).exp() - (-te / self.rise).exp())
    }
}

/// A conductance table, interpolated like `np.interp(t, x, y, 0.0, 0.0)`.
#[derive(Clone, Debug)]
pub struct Table {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Table {
    /// numpy's `arr_interp` for one point, branch for branch.
    pub fn interp(&self, t: f64) -> f64 {
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
                    (k, g * (v[k] - vsyn) / sa[0])
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

/// A validation failure of the constructors; `lib.rs` turns it into a Python
/// `ValueError` carrying the same text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoreError(pub String);

impl CoreError {
    fn new(msg: impl Into<String>) -> Self {
        CoreError(msg.into())
    }
}

impl std::fmt::Display for CoreError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CoreError {}

/// What `MultiRhs` needs per model, the arrays `#[pymethods] new` takes.
///
/// The arrays are per compartment (length n = `g_na.len()`), except `g_ax`,
/// `parent` and `child` (one entry per axial edge, 0-based).
#[derive(Clone, Debug)]
pub struct Params {
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
    /// Which of the six gates (m, h, p, w, z, a) move.
    pub active: [bool; 6],
    pub gating: Gating,
}

/// The stimulus, with the keywords of `_common.stimulus`.
///
/// `stim_type` is one of `multi.MSO_STIM_TYPES` or "none"; `k` and `k2` are
/// the 0-based input compartments and `epsg` is (1 / peak, decay, rise).
/// `Default` is "none" with the defaults of the Python signature.
#[derive(Clone, Debug)]
pub struct StimSpec {
    pub stim_type: String,
    pub start: f64,
    pub stop: f64,
    pub amp: f64,
    pub k: i64,
    pub k2: i64,
    pub f: f64,
    pub epsg: (f64, f64, f64),
    pub t_syn: Option<Vec<f64>>,
    pub g_syn: Option<Vec<f64>>,
    pub vsyn: f64,
    pub diff: Option<f64>,
}

impl Default for StimSpec {
    fn default() -> Self {
        StimSpec {
            stim_type: "none".to_string(),
            start: 0.0,
            stop: 0.0,
            amp: 0.0,
            k: 0,
            k2: 0,
            f: 0.0,
            epsg: (1.0, 1.0, 1.0),
            t_syn: None,
            g_syn: None,
            vsyn: 0.0,
            diff: None,
        }
    }
}

fn check_len(name: &str, v: &[f64], n: usize) -> Result<(), CoreError> {
    if v.len() != n {
        return Err(CoreError::new(format!(
            "{name} has {} entries, expected {n}",
            v.len()
        )));
    }
    Ok(())
}

fn index(name: &str, i: i64, n: usize) -> Result<usize, CoreError> {
    usize::try_from(i)
        .ok()
        .filter(|&u| u < n)
        .ok_or_else(|| CoreError::new(format!("{name}={i} is not in 0..{n}")))
}

/// The right-hand side of one `mso_axon` run: f(t, x) -> dx/dt.
#[derive(Clone, Debug)]
pub struct MultiRhs {
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
    gating: Gating,
    stim: Stim,
}

impl MultiRhs {
    /// Builds the right-hand side of one run; the checks (and their texts)
    /// are those of the Python constructor.
    pub fn new(params: Params, spec: &StimSpec) -> Result<Self, CoreError> {
        let Params {
            g_na,
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
            gating,
        } = params;
        let n = g_na.len();
        let arrays = [
            ("g_kht", &g_kht),
            ("g_klt", &g_klt),
            ("g_h", &g_h),
            ("g_lk", &g_lk),
            ("vlk", &vlk),
            ("na_vhalf", &na_vhalf),
            ("cap", &cap),
            ("sa", &sa),
        ];
        for (name, v) in arrays {
            check_len(name, v, n)?;
        }
        if parent.len() != g_ax.len() || child.len() != g_ax.len() {
            return Err(CoreError::new(
                "g_ax, parent and child must have one entry per edge",
            ));
        }
        let parent = parent
            .into_iter()
            .map(|i| index("parent", i, n))
            .collect::<Result<Vec<_>, _>>()?;
        let child = child
            .into_iter()
            .map(|i| index("child", i, n))
            .collect::<Result<Vec<_>, _>>()?;

        let StimSpec {
            ref stim_type,
            start,
            stop,
            amp,
            k,
            k2,
            f,
            epsg,
            ref t_syn,
            ref g_syn,
            vsyn,
            diff,
        } = *spec;
        let w = Epsg {
            inv_peak: epsg.0,
            decay: epsg.1,
            rise: epsg.2,
        };
        let slope = 1.0 / (stop - start);
        let stim = match stim_type.as_str() {
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
                    return Err(CoreError::new(format!("{stim_type} needs t_syn and g_syn")));
                };
                let table = Table {
                    x: x.clone(),
                    y: y.clone(),
                };
                if table.x.is_empty() || table.x.len() != table.y.len() {
                    return Err(CoreError::new(
                        "t_syn and g_syn must be non-empty and the same length",
                    ));
                }
                let diff = if stim_type == "SynapticPair" {
                    Some(diff.ok_or_else(|| CoreError::new("SynapticPair needs diff"))?)
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
            other => return Err(CoreError::new(format!("unknown stim_type {other:?}"))),
        };
        Ok(MultiRhs {
            n,
            g_na,
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
            gating,
            stim,
        })
    }

    /// The number of compartments; the state has 7n entries.
    pub fn n(&self) -> usize {
        self.n
    }

    /// The axial edges, 0-based: (parent, child), one entry per edge.
    pub fn edges(&self) -> (&[usize], &[usize]) {
        (&self.parent, &self.child)
    }

    /// Writes dx/dt for state `x` at time `t` into `out` (both 7n long).
    pub fn eval(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let n = self.n;
        let c = &self.gating;
        let (v, rest) = x.split_at(n);
        let (m, rest) = rest.split_at(n);
        let (h, rest) = rest.split_at(n);
        let (p, rest) = rest.split_at(n);
        let (w, rest) = rest.split_at(n);
        let (z, a) = rest.split_at(n);

        let mut total = vec![0.0f64; n];
        for i in 0..n {
            let vi = v[i];
            let ina = self.g_na[i]
                * m[i].powf(4.0)
                * (c.na_open_h * h[i] + c.na_open_floor)
                * (vi - c.v_na);
            let ikht = self.g_kht[i] * p[i] * (vi - c.v_k);
            let iklt = self.g_klt[i] * w[i].powf(4.0) * z[i] * (vi - c.v_k);
            let ih = self.g_h[i] * a[i] * (vi - c.v_h);
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
                let minf = 1.0 / (1.0 + ((vi + c.minf_v) / c.minf_k).exp());
                let taum = (c.taum_a + (c.taum_b / (1.0 + ((c.taum_v - vi) / c.taum_k).exp())))
                    / c.taum_div;
                (minf - m[i]) / taum
            } else {
                0.0
            };
            dh[i] = if act_h {
                let hinf = 1.0 / (1.0 + ((vi + self.na_vhalf[i]) / c.hinf_k).exp());
                let tauh = (c.tauh_a + (c.tauh_b / (1.0 + ((c.tauh_v - vi) / c.tauh_k).exp())))
                    / c.tauh_div;
                (hinf - h[i]) / tauh
            } else {
                0.0
            };
            dp[i] = if act_p {
                let pinf = 1.0 / (1.0 + (-(vi + c.pinf_v) / c.pinf_k).exp());
                let taup = c.taup_temp
                    * (c.taup_num
                        / (c.taup_a * ((vi + c.taup_v) / c.taup_ka).exp()
                            + c.taup_b * (-(vi + c.taup_v) / c.taup_kb).exp())
                        + c.taup_c);
                (pinf - p[i]) / taup
            } else {
                0.0
            };
            dw[i] = if act_w {
                let winf = 1.0 / (1.0 + ((vi + c.winf_v) / c.winf_k).exp());
                let tauw = c.tauw_num
                    / (c.tauw_a * ((vi + c.tauw_v) / c.tauw_ka).exp()
                        + c.tauw_b * (-(vi + c.tauw_v) / c.tauw_kb).exp())
                    + c.tauw_c;
                (winf - w[i]) / tauw
            } else {
                0.0
            };
            dz[i] = if act_z {
                let zinf = (1.0 - c.zinf_floor) / (1.0 + ((vi + c.zinf_v) / c.zinf_k).exp())
                    + c.zinf_floor;
                let tauz = c.tauz_num
                    / (c.tauz_a * ((vi + c.tauz_va) / c.tauz_ka).exp()
                        + (-(vi + c.tauz_vb) / c.tauz_kb).exp())
                    + c.tauz_c;
                (zinf - z[i]) / tauz
            } else {
                0.0
            };
            da[i] = if act_a {
                let ainf = 1.0 / (1.0 + (c.ainf_k * (vi + c.ainf_v)).exp());
                let s = vi + c.taua_v;
                let taua = c.taua_temp * (c.taua_c + c.taua_a * (-(s * s) / c.taua_k).exp());
                (ainf - a[i]) / taua
            } else {
                0.0
            };
        }
    }
}

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

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Numpy and Rust call the same libm on the machine the port was written
    /// on (darwin-arm64): bit-identical there, close elsewhere.
    const BITWISE: bool = cfg!(all(target_os = "macos", target_arch = "aarch64"));

    fn assert_same(got: f64, want: f64, rel: f64) {
        if BITWISE {
            assert_eq!(got.to_bits(), want.to_bits(), "{got:e} vs {want:e}");
        } else {
            assert!(
                (got - want).abs() <= rel * want.abs(),
                "{got:e} vs {want:e}"
            );
        }
    }

    fn default_epsg() -> Epsg {
        Epsg {
            inv_peak: 1.0 / 0.21317,
            decay: 0.18,
            rise: 0.1,
        }
    }

    fn table(x: &[f64], y: &[f64]) -> Table {
        Table {
            x: x.to_vec(),
            y: y.to_vec(),
        }
    }

    #[test]
    fn interp_matches_numpy_on_the_knots_and_between() {
        // values from np.interp(x, xs, ys, left=0.0, right=0.0)
        let t = table(&[0.0, 0.3, 1.0, 2.5], &[0.0, 2.0, 1.5, -1.0]);
        let cases = [
            (0.1, 0.6666666666666667),
            (0.65, 1.75),
            (1.7, 0.33333333333333337),
            (2.0, -0.16666666666666674),
        ];
        for (x, want) in cases {
            assert_same(t.interp(x), want, 1e-15);
        }
        // a knot gives its own value exactly, whatever the slopes
        for (x, want) in [(0.0f64, 0.0f64), (0.3, 2.0), (1.0, 1.5), (2.5, -1.0)] {
            assert_eq!(t.interp(x).to_bits(), want.to_bits(), "knot {x}");
        }
    }

    #[test]
    fn interp_is_zero_outside_the_table_not_clamped() {
        // synaptic.interp_g passes left=right=0.0; np.interp's default
        // would return the end values instead
        let t = table(&[0.0, 0.3, 1.0], &[4.0, 2.0, 1.5]);
        for x in [
            -1e-300,
            -0.1,
            f64::NEG_INFINITY,
            1.0 + 1e-15,
            2.6,
            f64::INFINITY,
        ] {
            assert_eq!(t.interp(x).to_bits(), 0.0f64.to_bits(), "x={x}");
        }
    }

    #[test]
    fn interp_on_a_one_point_table() {
        // numpy: the point itself gives y, anything else 0, and NaN (which
        // compares neither way) gives y as well
        let t = table(&[1.5], &[7.0]);
        assert_eq!(t.interp(1.5), 7.0);
        assert_eq!(t.interp(1.4), 0.0);
        assert_eq!(t.interp(1.6), 0.0);
        assert_eq!(t.interp(f64::NAN), 7.0);
    }

    #[test]
    fn interp_passes_nan_through_a_longer_table() {
        assert!(table(&[0.0, 1.0], &[0.0, 1.0]).interp(f64::NAN).is_nan());
    }

    #[test]
    fn interp_on_a_flat_segment_is_that_value() {
        let t = table(&[0.0, 1.0, 2.0], &[3.0, 3.0, 5.0]);
        assert_eq!(t.interp(0.37), 3.0);
    }

    /// The values `msoaxon._fast._gating()` passes (multi.V_NA, constants.py's
    /// literals and Q10 factors), for the fixed model below.
    fn gating() -> Gating {
        let q10 = |a: f64, b: f64| 3f64.powf((a - b) / 10.0);
        let values = [
            ("v_na", 69.0),
            ("v_k", -90.0),
            ("v_h", -35.0),
            ("na_open_h", 0.993),
            ("na_open_floor", 0.007),
            ("minf_v", 46.0),
            ("minf_k", -11.0),
            ("taum_a", 0.141),
            ("taum_b", -0.0826),
            ("taum_v", -20.5),
            ("taum_k", 10.8),
            ("taum_div", 3.0),
            ("hinf_k", 7.77),
            ("tauh_a", 4.0),
            ("tauh_b", -3.74),
            ("tauh_v", -40.6),
            ("tauh_k", 5.05),
            ("tauh_div", 3.0),
            ("pinf_v", 23.0),
            ("pinf_k", 6.0),
            ("taup_temp", q10(22.0, 35.0)),
            ("taup_num", 100.0),
            ("taup_a", 4.0),
            ("taup_ka", 32.0),
            ("taup_b", 5.0),
            ("taup_kb", 22.0),
            ("taup_v", 60.0),
            ("taup_c", 5.0),
            ("winf_v", 57.34),
            ("winf_k", -11.7),
            ("tauw_num", 21.5),
            ("tauw_a", 6.0),
            ("tauw_ka", 7.0),
            ("tauw_b", 24.0),
            ("tauw_kb", 50.6),
            ("tauw_v", 60.0),
            ("tauw_c", 0.35),
            ("zinf_floor", 0.27),
            ("zinf_v", 67.0),
            ("zinf_k", 6.16),
            ("tauz_num", 170.0),
            ("tauz_a", 5.0),
            ("tauz_va", 60.0),
            ("tauz_ka", 10.0),
            ("tauz_vb", 70.0),
            ("tauz_kb", 8.0),
            ("tauz_c", 10.7),
            ("ainf_k", 0.1),
            ("ainf_v", 80.4),
            ("taua_temp", q10(32.0, 35.0)),
            ("taua_c", 79.0),
            ("taua_a", 417.0),
            ("taua_v", 61.5),
            ("taua_k", 800.0),
        ];
        let map = values
            .map(|(k, v)| (k.to_string(), v))
            .into_iter()
            .collect();
        Gating::from_map(&map).unwrap()
    }

    /// A 45-compartment chain in the pattern of `constants.py` (soma, two
    /// AIS, 21 internode/node pairs) with an EPSG pair at the soma; the
    /// values are fixed, not the real morphology.
    pub(crate) fn model() -> MultiRhs {
        let n = 45;
        let axon = |head: [f64; 3], internode: f64, node: f64| -> Vec<f64> {
            (0..n)
                .map(|i| match i {
                    0..=2 => head[i],
                    _ if (i - 3) % 2 == 0 => internode,
                    _ => node,
                })
                .collect()
        };
        let edges = n - 1;
        MultiRhs {
            n,
            g_na: axon([0.2, 4.0, 4.0], 0.0, 4.0),
            g_kht: (0..n).map(|i| if i == 0 { 0.1 } else { 0.0 }).collect(),
            g_klt: axon([1.55; 3], 0.0, 1.55),
            g_h: (0..n).map(|i| if i < 3 { 0.02 } else { 0.0 }).collect(),
            g_lk: axon([0.0005; 3], 0.0002, 0.05),
            vlk: (0..n).map(|i| -65.0 - 0.05 * i as f64).collect(),
            na_vhalf: (0..n).map(|i| if i == 0 { 77.0 } else { 62.5 }).collect(),
            cap: (0..n).map(|i| 1.0 + 0.01 * i as f64).collect(),
            sa: axon([8750.0, 200.0, 20.7], 207.3, 2.07),
            g_ax: (0..edges)
                .map(|e| 2.0e-8 * (1.0 + 0.05 * e as f64))
                .collect(),
            parent: (0..edges).collect(),
            child: (1..n).collect(),
            active: [true; 6],
            gating: gating(),
            stim: Stim::EpsgPair {
                start: 5.0,
                td: 0.2,
                amp: 100.0,
                k: 0,
                w: default_epsg(),
            },
        }
    }

    /// [`model`] with a current step into the soma, on for 0 < t <= 50 ms
    /// (a negative `amp` hyperpolarises it).
    pub(crate) fn step_model(amp: f64) -> MultiRhs {
        MultiRhs {
            stim: Stim::Step {
                start: 0.0,
                stop: 50.0,
                amp,
            },
            ..model()
        }
    }

    /// [`model`] with a synaptic conductance into compartment `k`.
    pub(crate) fn synaptic_model(k: usize) -> MultiRhs {
        MultiRhs {
            stim: Stim::Synaptic {
                start: 0.0,
                stop: 50.0,
                k,
                table: table(&[0.0, 1.0, 50.0], &[0.0, 5.0, 0.0]),
                vsyn: 0.0,
                diff: None,
            },
            ..model()
        }
    }

    /// V from -75 mV in 0.7 mV steps, the six gates a fixed pattern in
    /// 0.05..0.95.
    pub(crate) fn state() -> Vec<f64> {
        (0..7 * 45)
            .map(|j| {
                if j < 45 {
                    -75.0 + 0.7 * j as f64
                } else {
                    0.05 + 0.9 * ((j * 37 % 101) as f64 / 100.0)
                }
            })
            .collect()
    }

    #[test]
    fn eval_keeps_the_values_captured_before_the_gating_moved_to_python() {
        let mut out = vec![0.0; 7 * 45];
        model().eval(5.1, &state(), &mut out);
        // every 23rd of the 315 outputs, from the code before this change
        let sampled: [f64; 14] = [
            -4.413877638942935e0,
            1.5664623361151598e-3,
            -1.6087343526101417e1,
            -1.1683193512611474e0,
            8.837880264411313e-2,
            1.374628095108558e-1,
            -1.701359809791315e-1,
            -2.3789557480257964e-1,
            -1.9384294020535559e-1,
            -2.774157034887882e-1,
            6.681458243998471e-3,
            -9.598634051855456e-3,
            3.922629037735645e-4,
            -1.367808281399423e-3,
        ];
        for (i, want) in sampled.into_iter().enumerate() {
            let got = out[23 * i];
            if BITWISE {
                assert_eq!(got.to_bits(), want.to_bits(), "out[{}]", 23 * i);
            } else {
                // elementwise, so the cancellation in (x_inf - x) / tau shows
                assert!(
                    (got - want).abs() <= 1e-9 * want.abs() + 1e-12,
                    "out[{}]: {got:e} vs {want:e}",
                    23 * i
                );
            }
        }
        if BITWISE {
            // FNV-1a over the bits of all 315 outputs
            let mut hash: u64 = 0xcbf29ce484222325;
            for v in &out {
                for byte in v.to_bits().to_le_bytes() {
                    hash = (hash ^ byte as u64).wrapping_mul(0x100000001b3);
                }
            }
            assert_eq!(hash, 0x4c8fdc374248577f);
        }
    }

    #[test]
    fn eval_follows_a_changed_constant() {
        let (x, mut out, mut moved) = (state(), vec![0.0; 315], vec![0.0; 315]);
        let mut rhs = model();
        rhs.eval(5.1, &x, &mut out);
        rhs.gating.v_na += 1.0;
        rhs.eval(5.1, &x, &mut moved);
        assert_ne!(out[0], moved[0]); // dV of the soma
        assert_eq!(out[45..], moved[45..]); // the gates don't see V_Na
    }

    #[test]
    fn unitary_keeps_its_value_at_a_reachable_time() {
        // captured from the code before the clamp was added
        assert_same(default_epsg().unitary(0.1), 9.657736997044162e-1, 1e-14);
    }

    #[test]
    fn unitary_is_negative_zero_far_before_onset() {
        // exp(200 / 0.1) overflows: without the clamp this was 0 * inf = NaN
        let far = default_epsg().unitary(-200.0);
        assert_eq!(far, 0.0);
        assert!(far.is_sign_negative(), "expected -0.0, got {far:?}");
        // the same -0.0 just either side of the clamp, and at onset
        for t in [-0.5, -1.0, -1.5, -70.0, -71.0] {
            let u = default_epsg().unitary(t);
            assert!(u == 0.0 && u.is_sign_negative(), "t={t}: {u:?}");
        }
        assert!(default_epsg().unitary(f64::NAN).is_nan());
    }

    #[test]
    fn chain_sparsity_has_the_shape_of_multi_chain_jac_sparsity() {
        let n = 45;
        let (rows, cols) = chain_jac_sparsity(n);
        // tridiagonal V block + 6 n (dV, gates) + 6 n (gates, V) + 6 n diagonal
        assert_eq!(rows.len(), (3 * n - 2) + 18 * n);
        assert_eq!(rows.len(), cols.len());
        // CSC order: columns ascending, rows strictly ascending within one
        for w in 0..rows.len() - 1 {
            assert!((cols[w], rows[w]) < (cols[w + 1], rows[w + 1]), "entry {w}");
        }
        assert!(rows.iter().chain(&cols).all(|&i| i < 7 * n));
        let column = |j: usize| -> Vec<usize> {
            (0..rows.len())
                .filter(|&e| cols[e] == j)
                .map(|e| rows[e])
                .collect()
        };
        // first V column: itself, the next V, its six gates
        assert_eq!(column(0), [0, 1, 45, 90, 135, 180, 225, 270]);
        // an interior V column also couples to both neighbours
        assert_eq!(column(10), [9, 10, 11, 55, 100, 145, 190, 235, 280]);
        assert_eq!(column(44), [43, 44, 89, 134, 179, 224, 269, 314]);
        // a gate column: dV of the same compartment and itself
        assert_eq!(column(45 + 3), [3, 48]);
        assert_eq!(column(6 * 45 + 44), [44, 314]);
    }

    #[test]
    fn tree_sparsity_of_a_chain_is_the_chain_pattern_and_a_branch_widens_it() {
        let n = 45;
        let chain = chain_jac_sparsity(n);
        let parent: Vec<usize> = (0..n - 1).collect();
        let child: Vec<usize> = (1..n).collect();
        assert_eq!(tree_jac_sparsity(n, &parent, &child), chain);
        // compartment 2 also feeds a branch compartment 44 (edge 2 -> 44)
        let (rows, cols) = tree_jac_sparsity(n, &[0, 1, 2, 2], &[1, 2, 3, 44]);
        let col2: Vec<usize> = (0..rows.len())
            .filter(|&e| cols[e] == 2)
            .map(|e| rows[e])
            .collect();
        assert_eq!(col2[..4], [1, 2, 3, 44]);
    }
}
