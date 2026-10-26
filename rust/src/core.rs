//! The right-hand side of msoaxon's 45-compartment model, in plain Rust.
//!
//! No pyo3 here: `lib.rs` is the thin Python binding. A Rust integrator
//! builds a [`MultiRhs`] from [`Params`] and a [`StimSpec`], calls
//! [`MultiRhs::eval`], asks [`tree_jac_sparsity`] for the Jacobian pattern
//! and [`MultiRhs::jacobian`] for its values.
//!
//! A port of `msoaxon.multi._rhs` (with `external_current`,
//! `bilateral_current`, `events_current`, `axial_current` and the gating
//! functions of `msoaxon.constants`). The Python code is the reference. Every expression
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
    /// The "Events" stimulus: one EPSG-shaped conductance per event, each at
    /// its own compartment and voltage (`multi.events_current`).
    Events {
        onset: Vec<f64>,
        k: Vec<usize>,
        /// `scale * g`: numpy multiplies the run's `I` into each peak first,
        /// so folding it in once keeps the bits and saves a multiply a call.
        scale_g: Vec<f64>,
        e_rev: Vec<f64>,
        w: Vec<Epsg>,
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
            Stim::Events {
                ref onset,
                ref k,
                ref scale_g,
                ref e_rev,
                ref w,
            } => {
                // np.bincount(k, weights=-i / sa[k], minlength=n): a fresh
                // zero per compartment, each event added in list order, then
                // the sums added to every compartment (an exact 0.0 where no
                // event enters, which also turns an event's -0.0 into +0.0)
                let mut density = vec![0.0f64; total.len()];
                for e in 0..onset.len() {
                    let ke = k[e];
                    let i = scale_g[e] * (e_rev[e] - v[ke]) * w[e].unitary(t - onset[e]);
                    density[ke] += -i / sa[ke];
                }
                for (sum, d) in total.iter_mut().zip(&density) {
                    *sum += d;
                }
                return;
            }
        };
        total[k] += iext;
    }

    /// Calls `f(k, d)` with d the derivative, with respect to V[k], of the
    /// current density [`add_to`](Self::add_to) adds to compartment k, under
    /// the same on/off conditions; once per conductance input (so twice for
    /// a compartment two inputs share), never for a current that does not
    /// depend on V (ramp, ramp2, step, sine).
    fn slopes(&self, t: f64, sa: &[f64], mut f: impl FnMut(usize, f64)) {
        match *self {
            Stim::None
            | Stim::Ramp { .. }
            | Stim::Ramp2 { .. }
            | Stim::Step { .. }
            | Stim::Sine { .. } => {}
            Stim::Synaptic {
                start,
                stop,
                k,
                ref table,
                diff,
                ..
            } => {
                if start < t && t <= stop {
                    let g = match diff {
                        None => table.interp(t),
                        Some(d) => table.interp(t) + table.interp(t + d),
                    };
                    // g (V - vsyn) / sa[0]
                    f(k, g / sa[0]);
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
                    // amp (0 - V) on g(te) / -sa[k]
                    f(k, amp * on * w.unitary(te) / sa[k]);
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
                f(k, amp * wave / sa[k]);
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
                f(a, amp * on0 * w.unitary(te) / sa[a]);
                f(b, amp * on1 * w.unitary(te - td) / sa[b]);
            }
            Stim::Events {
                ref onset,
                ref k,
                ref scale_g,
                ref w,
                ..
            } => {
                // -(scale g (e_rev - V) g(t - onset)) / sa[k], per event
                for e in 0..onset.len() {
                    let ke = k[e];
                    f(ke, scale_g[e] * w[e].unitary(t - onset[e]) / sa[ke]);
                }
            }
        }
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
    /// The event table of "Events"; `amp` is the run's `I`, scaling every peak.
    pub events: Option<EventTable>,
}

/// The arrays of `synaptic.event_arrays`, one entry per event.
///
/// `k` is the 0-based compartment, `inv_peak` is `1 / peak` as Python
/// computes it (`synaptic.epsg_inv_peak`).
#[derive(Clone, Debug, Default)]
pub struct EventTable {
    pub t: Vec<f64>,
    pub k: Vec<i64>,
    pub g: Vec<f64>,
    pub rise: Vec<f64>,
    pub decay: Vec<f64>,
    pub inv_peak: Vec<f64>,
    pub e_rev: Vec<f64>,
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
            events: None,
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
            ref events,
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
            "Events" => {
                let Some(ev) = events else {
                    return Err(CoreError::new("Events needs an event table"));
                };
                let count = ev.t.len();
                if count == 0 {
                    return Err(CoreError::new("Events needs at least one event"));
                }
                for (name, v) in [
                    ("g", &ev.g),
                    ("rise", &ev.rise),
                    ("decay", &ev.decay),
                    ("inv_peak", &ev.inv_peak),
                    ("e_rev", &ev.e_rev),
                ] {
                    check_len(&format!("events {name}"), v, count)?;
                }
                if ev.k.len() != count {
                    return Err(CoreError::new(format!(
                        "events k has {} entries, expected {count}",
                        ev.k.len()
                    )));
                }
                Stim::Events {
                    onset: ev.t.clone(),
                    k: ev
                        .k
                        .iter()
                        .map(|&i| index("events k", i, n))
                        .collect::<Result<Vec<_>, _>>()?,
                    scale_g: ev.g.iter().map(|g| amp * g).collect(),
                    e_rev: ev.e_rev.clone(),
                    w: (0..count)
                        .map(|e| Epsg {
                            inv_peak: ev.inv_peak[e],
                            decay: ev.decay[e],
                            rise: ev.rise[e],
                        })
                        .collect(),
                }
            }
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

    /// Writes the exact Jacobian of [`eval`](Self::eval) at (t, x) into
    /// `out`: one value per entry of [`tree_jac_sparsity`] of this model's
    /// edges, in that order (so `out` has that pattern's length).
    ///
    /// dV_i = (-I_i + axial_i 1e9 / sa_i) / cap_i, so dV_i/dV_i holds every
    /// conductance (leak, the channels' g times open fraction, the
    /// stimulus's dI/dV) and the axial diagonal, dV_i/dV_j of an edge its
    /// g_ax (scaled by row i's sa and cap), and dV_i/d(gate) the open
    /// fraction's derivative times the driving force. A gate row is
    /// (x_inf'(V) - (x_inf - x) tau'(V) / tau) / tau and -1 / tau, or zero
    /// where that gate is frozen. The values are those of `eval`'s
    /// expressions, not its float order: only `eval` is bit-pinned.
    pub fn jacobian(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let n = self.n;
        let c = &self.gating;
        let (v, rest) = x.split_at(n);
        let (m, rest) = rest.split_at(n);
        let (h, rest) = rest.split_at(n);
        let (p, rest) = rest.split_at(n);
        let (w, rest) = rest.split_at(n);
        let (z, a) = rest.split_at(n);

        // where tree_jac_sparsity puts each entry: V column j starts at
        // start[j] with its V rows (nbrs[j]) then its six gate rows; gate
        // g's column of compartment j is the pair (dV_j, gate) at
        // gates + 2 (g n + j), g counted from 0
        let nbrs = tree_neighbours(n, &self.parent, &self.child);
        let mut start = Vec::with_capacity(n);
        let mut gates = 0;
        for nb in &nbrs {
            start.push(gates);
            gates += nb.len() + 6;
        }
        assert_eq!(out.len(), gates + 12 * n, "out is not the pattern's length");
        let at = |r: usize, j: usize| start[j] + nbrs[j].binary_search(&r).unwrap();
        out.fill(0.0);

        let [act_m, act_h, act_p, act_w, act_z, act_a] = self.active;
        for i in 0..n {
            let (vi, cap) = (v[i], self.cap[i]);
            let m3 = m[i] * m[i] * m[i];
            let w3 = w[i] * w[i] * w[i];
            let na_h = c.na_open_h * h[i] + c.na_open_floor;
            let g_na = self.g_na[i] * m3 * m[i] * na_h;
            let g_klt = self.g_klt[i] * w3 * w[i] * z[i];
            let g = g_na + self.g_kht[i] * p[i] + g_klt + self.g_h[i] * a[i] + self.g_lk[i];
            out[at(i, i)] -= g / cap;

            // d I_i / d(gate), and d(gate)/dV, d(gate)/d(gate) of the
            // gate's relaxation, gate by gate in state order
            let (e_na, e_k, e_h) = (vi - c.v_na, vi - c.v_k, vi - c.v_h);
            let di = [
                self.g_na[i] * 4.0 * m3 * na_h * e_na,
                self.g_na[i] * m3 * m[i] * c.na_open_h * e_na,
                self.g_kht[i] * e_k,
                self.g_klt[i] * 4.0 * w3 * z[i] * e_k,
                self.g_klt[i] * w3 * w[i] * e_k,
                self.g_h[i] * e_h,
            ];
            let rows = [
                act_m.then(|| {
                    let (inf, ds) = logistic((vi + c.minf_v) / c.minf_k);
                    let (q, dq) = logistic((c.taum_v - vi) / c.taum_k);
                    let tau = (c.taum_a + c.taum_b * q) / c.taum_div;
                    let d_tau = c.taum_b * (-dq / c.taum_k) / c.taum_div;
                    relax(inf, ds / c.minf_k, m[i], tau, d_tau)
                }),
                act_h.then(|| {
                    let (inf, ds) = logistic((vi + self.na_vhalf[i]) / c.hinf_k);
                    let (q, dq) = logistic((c.tauh_v - vi) / c.tauh_k);
                    let tau = (c.tauh_a + c.tauh_b * q) / c.tauh_div;
                    let d_tau = c.tauh_b * (-dq / c.tauh_k) / c.tauh_div;
                    relax(inf, ds / c.hinf_k, h[i], tau, d_tau)
                }),
                act_p.then(|| {
                    let (inf, ds) = logistic(-(vi + c.pinf_v) / c.pinf_k);
                    let s = vi + c.taup_v;
                    let (e1, e2) = ((s / c.taup_ka).exp(), (-s / c.taup_kb).exp());
                    let d = c.taup_a * e1 + c.taup_b * e2;
                    let d_d = c.taup_a * e1 / c.taup_ka - c.taup_b * e2 / c.taup_kb;
                    let tau = c.taup_temp * (c.taup_num / d + c.taup_c);
                    let d_tau = finite_over(-c.taup_temp * c.taup_num * d_d, d);
                    relax(inf, -ds / c.pinf_k, p[i], tau, d_tau)
                }),
                act_w.then(|| {
                    let (inf, ds) = logistic((vi + c.winf_v) / c.winf_k);
                    let s = vi + c.tauw_v;
                    let (e1, e2) = ((s / c.tauw_ka).exp(), (-s / c.tauw_kb).exp());
                    let d = c.tauw_a * e1 + c.tauw_b * e2;
                    let d_d = c.tauw_a * e1 / c.tauw_ka - c.tauw_b * e2 / c.tauw_kb;
                    let tau = c.tauw_num / d + c.tauw_c;
                    relax(
                        inf,
                        ds / c.winf_k,
                        w[i],
                        tau,
                        finite_over(-c.tauw_num * d_d, d),
                    )
                }),
                act_z.then(|| {
                    let (s, ds) = logistic((vi + c.zinf_v) / c.zinf_k);
                    let inf = (1.0 - c.zinf_floor) * s + c.zinf_floor;
                    let d_inf = (1.0 - c.zinf_floor) * ds / c.zinf_k;
                    let e1 = ((vi + c.tauz_va) / c.tauz_ka).exp();
                    let e2 = (-(vi + c.tauz_vb) / c.tauz_kb).exp();
                    let d = c.tauz_a * e1 + e2;
                    let d_d = c.tauz_a * e1 / c.tauz_ka - e2 / c.tauz_kb;
                    let tau = c.tauz_num / d + c.tauz_c;
                    relax(inf, d_inf, z[i], tau, finite_over(-c.tauz_num * d_d, d))
                }),
                act_a.then(|| {
                    let (inf, ds) = logistic(c.ainf_k * (vi + c.ainf_v));
                    let s = vi + c.taua_v;
                    let bump = c.taua_a * (-(s * s) / c.taua_k).exp();
                    let tau = c.taua_temp * (c.taua_c + bump);
                    let d_tau = c.taua_temp * bump * (-2.0 * s / c.taua_k);
                    relax(inf, ds * c.ainf_k, a[i], tau, d_tau)
                }),
            ];
            let gate_rows = start[i] + nbrs[i].len();
            for (g, (di, row)) in di.into_iter().zip(rows).enumerate() {
                let col = gates + 2 * (g * n + i);
                out[col] = -di / cap;
                if let Some((dv, dx)) = row {
                    out[gate_rows + g] = dv;
                    out[col + 1] = dx;
                }
            }
        }
        self.stim
            .slopes(t, &self.sa, |k, s| out[at(k, k)] -= s / self.cap[k]);

        // axial_i = sum over i's edges of +-g_ax (V_parent - V_child): an
        // edge adds -g_ax to d axial_r / dV_r and +g_ax to d axial_r / dV_o
        // at both its ends r (o the other end)
        for (e, (&pa, &ch)) in self.parent.iter().zip(&self.child).enumerate() {
            let g = self.g_ax[e] * 1e9;
            for (r, o) in [(pa, ch), (ch, pa)] {
                let gr = g / self.sa[r] / self.cap[r];
                out[at(r, r)] -= gr;
                out[at(r, o)] += gr;
            }
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
    let (mut rows, mut cols) = (Vec::new(), Vec::new());
    for (j, nb) in tree_neighbours(n, parent, child).iter().enumerate() {
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

/// Per compartment, itself and every compartment it shares an edge with,
/// ascending and each once: the V rows of its V column in
/// [`tree_jac_sparsity`].
fn tree_neighbours(n: usize, parent: &[usize], child: &[usize]) -> Vec<Vec<usize>> {
    let mut nbrs: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    for (&p, &c) in parent.iter().zip(child) {
        nbrs[p].push(c);
        nbrs[c].push(p);
    }
    for nb in &mut nbrs {
        nb.sort_unstable();
        nb.dedup();
    }
    nbrs
}

/// A logistic 1 / (1 + exp(u)) and its derivative with respect to u.
///
/// The derivative is -1 / (2 + e + 1 / e) with e = exp(u): no 1 - s to
/// cancel and no inf * 0, so it is accurate (and -0.0, not NaN) at both
/// saturated ends.
#[inline]
fn logistic(u: f64) -> (f64, f64) {
    let e = u.exp();
    (1.0 / (1.0 + e), -1.0 / (2.0 + e + 1.0 / e))
}

/// `num / d^2`, the form of d tau / dV for a tau = k / d(V) + c (num being
/// -k d'), or 0 where d has overflowed: `eval`'s tau is then flat at c,
/// while the quotient would be inf / inf = NaN.
#[inline]
fn finite_over(num: f64, d: f64) -> f64 {
    if d.is_finite() {
        num / (d * d)
    } else {
        0.0
    }
}

/// d/dV and d/dx of a gate's (x_inf(V) - x) / tau(V), from x_inf, tau and
/// their derivatives with respect to V.
#[inline]
fn relax(inf: f64, d_inf: f64, x: f64, tau: f64, d_tau: f64) -> (f64, f64) {
    ((d_inf - (inf - x) / tau * d_tau) / tau, -1.0 / tau)
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

    /// The `Params` of [`model`], to build a model through `MultiRhs::new`.
    pub(crate) fn model_params() -> Params {
        let base = model();
        Params {
            g_na: base.g_na,
            g_kht: base.g_kht,
            g_klt: base.g_klt,
            g_h: base.g_h,
            g_lk: base.g_lk,
            vlk: base.vlk,
            na_vhalf: base.na_vhalf,
            cap: base.cap,
            sa: base.sa,
            g_ax: base.g_ax,
            parent: base.parent.iter().map(|&i| i as i64).collect(),
            child: base.child.iter().map(|&i| i as i64).collect(),
            active: base.active,
            gating: base.gating,
        }
    }

    /// An EPSG at the soma, a slower IPSG-like event at compartment 10 that
    /// starts first, and a second soma event.
    fn events_table() -> EventTable {
        EventTable {
            t: vec![5.0, 4.8, 5.2],
            k: vec![0, 10, 0],
            g: vec![60.0, 30.0, 40.0],
            rise: vec![0.1, 0.28, 0.1],
            decay: vec![0.18, 1.5, 0.18],
            inv_peak: vec![1.0 / 0.21317, 2.5, 1.0 / 0.21317],
            e_rev: vec![0.0, -80.0, 0.0],
        }
    }

    /// [`model`] with [`events_table`], `scale` being the run's I.
    pub(crate) fn events_model(scale: f64) -> MultiRhs {
        let spec = StimSpec {
            stim_type: "Events".to_string(),
            amp: scale,
            events: Some(events_table()),
            ..StimSpec::default()
        };
        MultiRhs::new(model_params(), &spec).unwrap()
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

    /// One model per stimulus type (two for the bilateral input, the second
    /// on one compartment), and an EPSG pair with three gates frozen, each
    /// with a time inside its stimulus and away from its switches.
    pub(crate) fn stim_models() -> Vec<(&'static str, MultiRhs, f64)> {
        let spec = |stim_type: &str, start: f64, stop: f64, amp: f64, k: i64, k2: i64| StimSpec {
            stim_type: stim_type.to_string(),
            start,
            stop,
            amp,
            k,
            k2,
            f: 100.0,
            epsg: (1.0 / 0.21317, 0.18, 0.1),
            t_syn: Some(vec![0.0, 1.0, 50.0]),
            g_syn: Some(vec![0.0, 5.0, 0.0]),
            vsyn: 0.0,
            diff: Some(0.4),
            events: Some(events_table()),
        };
        let build = |s: StimSpec| MultiRhs::new(model_params(), &s).unwrap();
        let mut frozen = model_params();
        frozen.active = [true, false, true, false, true, false];
        vec![
            ("none", build(spec("none", 0.0, 0.0, 0.0, 0, 0)), 1.0),
            ("ramp", build(spec("ramp", 1.0, 20.0, 500.0, 0, 0)), 10.0),
            ("ramp2", build(spec("ramp2", 1.0, 20.0, 500.0, 0, 0)), 10.0),
            ("step", build(spec("step", 0.0, 50.0, 1000.0, 0, 0)), 10.0),
            ("sine", build(spec("sine", 0.0, 50.0, 1000.0, 0, 0)), 10.3),
            (
                "Synaptic",
                build(spec("Synaptic", 0.0, 50.0, 0.0, 2, 0)),
                10.3,
            ),
            (
                "SynapticPair",
                build(spec("SynapticPair", 0.0, 50.0, 0.0, 2, 0)),
                10.3,
            ),
            ("EPSG", build(spec("EPSG", 5.0, 25.0, 100.0, 0, 0)), 5.13),
            (
                "EPSGpair",
                build(spec("EPSGpair", 5.0, 5.2, 100.0, 3, 0)),
                5.27,
            ),
            (
                "EPSGbilateral",
                build(spec("EPSGbilateral", 5.0, 5.1, 100.0, 0, 3)),
                5.17,
            ),
            (
                "EPSGbilateral, one compartment",
                build(spec("EPSGbilateral", 5.0, 5.1, 100.0, 2, 2)),
                5.17,
            ),
            ("Events", build(spec("Events", 0.0, 0.0, 1.5, 0, 0)), 5.25),
            (
                "frozen gates",
                MultiRhs::new(frozen, &spec("EPSGpair", 5.0, 5.2, 100.0, 0, 0)).unwrap(),
                5.27,
            ),
        ]
    }

    /// Three states in the regimes the Jacobian must cover: every V near
    /// rest, an action potential part-way up the axon, and the soma far
    /// below rest; the gates as in [`state`].
    fn jac_states() -> Vec<(&'static str, Vec<f64>)> {
        let with_v = |v: &dyn Fn(usize) -> f64| -> Vec<f64> {
            let mut x = state();
            x[..45].iter_mut().enumerate().for_each(|(i, x)| *x = v(i));
            x
        };
        vec![
            ("rest", with_v(&|i| -68.0 + 0.01 * i as f64)),
            ("upstroke", with_v(&|i| -60.0 + 1.5 * i as f64)),
            (
                "hyperpolarised soma",
                with_v(&|i| if i == 0 { -120.0 } else { -66.0 }),
            ),
        ]
    }

    #[test]
    fn logistic_is_finite_and_exact_at_both_saturated_ends() {
        for u in [-800.0, -40.0, 0.0, 40.0, 800.0] {
            let (s, ds) = logistic(u);
            assert!(
                s.is_finite() && ds.is_finite() && ds <= 0.0,
                "u={u}: {s} {ds}"
            );
            // ds/du = -s (1 - s), by its own formula where neither saturates
            if u.abs() < 100.0 {
                assert!((ds + s * (1.0 - s)).abs() <= 1e-15, "u={u}");
            }
        }
        assert_eq!(logistic(0.0), (0.5, -0.25));
    }

    #[test]
    fn the_jacobian_stays_finite_where_a_time_constant_overflows() {
        // a diverging Newton predictor can reach any V; the denominators of
        // tau_p, tau_w and tau_z overflow to inf beyond ~|V| 5 V while eval
        // stays finite, and so must the Jacobian
        for v in [-1e4, -6e3, 8e3, 3e4] {
            let mut x = state();
            x[..45].iter_mut().for_each(|x| *x = v);
            let (mut f, mut jac) = (vec![0.0; 315], vec![0.0; 943]);
            model().eval(5.1, &x, &mut f);
            assert!(f.iter().all(|y| y.is_finite()), "eval at V = {v}");
            model().jacobian(5.1, &x, &mut jac);
            assert!(jac.iter().all(|y| y.is_finite()), "jacobian at V = {v}");
        }
    }

    #[test]
    fn the_jacobian_is_the_tree_pattern_and_central_differences_of_eval() {
        // every V column perturbed by a central difference: rows off the
        // pattern must not move at all (bit for bit), rows on it must agree
        // with the analytic value. h = 1e-5 (1 + |x|): the truncation error
        // h^2 f''' / 6 is ~1e-9 of an entry for the steepest kinetics (a
        // 5 mV slope factor), the round-off ~1e-11 of the row's scale.
        let mut worst = (0.0_f64, "", "");
        for (name, rhs, t) in stim_models() {
            let n = rhs.n();
            let ns = 7 * n;
            let (rows, cols) = tree_jac_sparsity(n, &rhs.parent, &rhs.child);
            let mut slot = vec![None; ns * ns];
            for (e, (&r, &c)) in rows.iter().zip(&cols).enumerate() {
                slot[r * ns + c] = Some(e);
            }
            for (label, x) in jac_states() {
                let mut jac = vec![f64::NAN; rows.len()];
                rhs.jacobian(t, &x, &mut jac);
                assert!(jac.iter().all(|v| v.is_finite()), "{name}, {label}");
                // each row's scale: its largest analytic entry
                let mut scale = vec![0.0_f64; ns];
                for (e, &r) in rows.iter().enumerate() {
                    scale[r] = scale[r].max(jac[e].abs());
                }
                let (mut fp, mut fm) = (vec![0.0; ns], vec![0.0; ns]);
                for c in 0..ns {
                    let h = 1e-5 * (1.0 + x[c].abs());
                    let mut y = x.clone();
                    y[c] = x[c] + h;
                    rhs.eval(t, &y, &mut fp);
                    y[c] = x[c] - h;
                    rhs.eval(t, &y, &mut fm);
                    for r in 0..ns {
                        match slot[r * ns + c] {
                            None => assert_eq!(
                                fp[r], fm[r],
                                "{name}, {label}: ({r}, {c}) is off the pattern"
                            ),
                            Some(e) => {
                                let fd = (fp[r] - fm[r]) / (2.0 * h);
                                let err = (jac[e] - fd).abs() / scale[r].max(f64::MIN_POSITIVE);
                                if err > worst.0 {
                                    worst = (err, name, label);
                                }
                                assert!(
                                    err < 1e-7,
                                    "{name}, {label}: ({r}, {c}) {} vs {fd} ({err:e})",
                                    jac[e]
                                );
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "central differences: worst {:e} ({}, {})",
            worst.0, worst.1, worst.2
        );
    }

    #[test]
    fn the_jacobian_of_a_branched_tree_sums_its_edges() {
        // compartment 2 also feeds 44 (the chain's 43 -> 44 edge moved), and
        // a doubled edge 0 -> 1 couples twice
        let mut rhs = model();
        rhs.parent = vec![0, 0];
        rhs.child = vec![1, 1];
        rhs.parent.extend(1..43);
        rhs.child.extend(2..44);
        rhs.parent.push(2);
        rhs.child.push(44);
        rhs.g_ax = (0..rhs.parent.len())
            .map(|e| 1e-8 * (1.0 + 0.1 * e as f64))
            .collect();
        let (rows, cols) = tree_jac_sparsity(45, &rhs.parent, &rhs.child);
        let x = state();
        let mut jac = vec![0.0; rows.len()];
        rhs.jacobian(5.1, &x, &mut jac);
        let entry = |r: usize, c: usize| {
            let e = (0..rows.len())
                .find(|&e| (rows[e], cols[e]) == (r, c))
                .unwrap();
            jac[e]
        };
        let coupling = |r: usize, g: f64| g * 1e9 / rhs.sa[r] / rhs.cap[r];
        let (g01, g244) = (rhs.g_ax[0] + rhs.g_ax[1], rhs.g_ax[44]);
        assert!((entry(0, 1) - coupling(0, g01)).abs() <= 1e-12 * coupling(0, g01));
        assert!((entry(44, 2) - coupling(44, g244)).abs() <= 1e-12 * coupling(44, g244));
        assert!((entry(2, 44) - coupling(2, g244)).abs() <= 1e-12 * coupling(2, g244));
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

    #[test]
    fn events_leave_the_jacobian_pattern_unchanged() {
        let rhs = events_model(1.0);
        let n = rhs.n();
        let (parent, child) = rhs.edges();
        assert_eq!(tree_jac_sparsity(n, parent, child), chain_jac_sparsity(n));
        // and the claim behind it: an event reads and drives only its own
        // compartment, so dV of the soma and of compartment 10 does not move
        // when the voltage of a compartment that is not their neighbour does
        let (x, mut base, mut moved) = (state(), vec![0.0; 315], vec![0.0; 315]);
        rhs.eval(5.1, &x, &mut base);
        for j in [3, 20, 30] {
            let mut y = x.clone();
            y[j] += 1.0;
            rhs.eval(5.1, &y, &mut moved);
            for r in [0, 10] {
                assert_eq!(base[r].to_bits(), moved[r].to_bits(), "V{j} -> dV{r}");
            }
        }
        // control: the event at 10 does read V[10]
        let mut y = x.clone();
        y[10] += 1.0;
        rhs.eval(5.1, &y, &mut moved);
        assert_ne!(base[10].to_bits(), moved[10].to_bits());
    }

    #[test]
    fn one_event_is_the_epsg_pair_half_bit_for_bit() {
        // Python: test_multi_single_event_is_the_epsgpair_half_bit_for_bit.
        // One soma event with g = 1 and scale = the amplitude against an
        // EPSG pair whose second waveform lies beyond the horizon.
        let single = |stim| MultiRhs { stim, ..model() };
        let w = default_epsg();
        let event = single(Stim::Events {
            onset: vec![5.0],
            k: vec![0],
            scale_g: vec![100.0 * 1.0],
            e_rev: vec![0.0],
            w: vec![w],
        });
        let pair = single(Stim::EpsgPair {
            start: 5.0,
            td: 14.0,
            amp: 100.0,
            k: 0,
            w,
        });
        let x = state();
        let (mut a, mut b) = (vec![0.0; 315], vec![0.0; 315]);
        for t in [5.0 + 1e-9, 5.05, 5.3, 7.0, 15.0] {
            event.eval(t, &x, &mut a);
            pair.eval(t, &x, &mut b);
            for (i, (u, v)) in a.iter().zip(&b).enumerate() {
                assert_eq!(u.to_bits(), v.to_bits(), "t={t} out[{i}]");
            }
        }
    }

    #[test]
    fn events_touch_only_their_compartments_and_nothing_before_onset() {
        let (x, mut one, mut two) = (state(), vec![0.0; 315], vec![0.0; 315]);
        events_model(1.0).eval(5.3, &x, &mut one);
        events_model(2.0).eval(5.3, &x, &mut two);
        // the events act on compartments 0 and 10 alone
        let mut none = vec![0.0; 315];
        MultiRhs {
            stim: Stim::None,
            ..events_model(1.0)
        }
        .eval(5.3, &x, &mut none);
        for i in 0..45 {
            let touched = i == 0 || i == 10;
            assert_eq!(one[i] != none[i], touched, "dV[{i}]");
            assert_eq!(two[i] != none[i], touched, "dV[{i}]");
        }
        // before every onset the current is an exact zero
        let (mut early, mut quiet) = (vec![0.0; 315], vec![0.0; 315]);
        events_model(1.0).eval(4.0, &x, &mut early);
        MultiRhs {
            stim: Stim::None,
            ..events_model(1.0)
        }
        .eval(4.0, &x, &mut quiet);
        assert_eq!(early, quiet);
    }

    #[test]
    fn an_event_table_is_validated() {
        let build = |edit: &dyn Fn(&mut StimSpec)| {
            let mut spec = StimSpec {
                stim_type: "Events".to_string(),
                events: Some(events_table()),
                ..StimSpec::default()
            };
            edit(&mut spec);
            MultiRhs::new(model_params(), &spec).err().map(|e| e.0)
        };
        assert_eq!(build(&|_| {}), None);
        assert_eq!(
            build(&|s| s.events = None).as_deref(),
            Some("Events needs an event table")
        );
        assert_eq!(
            build(&|s| s.events = Some(EventTable::default())).as_deref(),
            Some("Events needs at least one event")
        );
        let msg = build(&|s| s.events.as_mut().unwrap().k[2] = 45).unwrap();
        assert!(msg.contains("events k=45 is not in 0..45"), "{msg}");
        let msg = build(&|s| s.events.as_mut().unwrap().g.push(1.0)).unwrap();
        assert!(msg.contains("events g has 4 entries, expected 3"), "{msg}");
        let msg = build(&|s| s.events.as_mut().unwrap().k.push(0)).unwrap();
        assert!(msg.contains("events k has 4 entries, expected 3"), "{msg}");
    }
}
