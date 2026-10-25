//! Toy stiff cable model with the sparsity of msoaxon's active-full system.
//!
//! State layout (315 = 7 * 45): voltage block `V[0..45]` followed by six
//! gate blocks `g_k[0..45]`.  The voltage Jacobian block is tridiagonal
//! (a sealed cable) and each gate block is diagonal, with each gate driven
//! only by the voltage of its own compartment, exactly like
//! `multi._chain_jac_sparsity`.  [`Model::tree`] swaps the chain for a
//! branched morphology given as a parent/child edge list (the pattern of
//! `multi._tree_jac_sparsity`).

use std::cell::Cell;

use crate::pattern::tree_jac_sparsity;

/// Which Jacobian-vector product the solver is handed.
///
/// All of them return the same `J v` for finite `v`; they differ only in
/// what they do with the NaN-tainted probe vectors diffsol uses to detect
/// the sparsity pattern.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Jv {
    /// Branch-free analytic product ([`Model::jv`]).
    Analytic,
    /// Forward difference of the RHS ([`Model::jv_fd`]).
    Fd,
    /// Analytic, skipping inputs with `v == 0.0` ([`Model::jv_skip_zero`]).
    SkipZero,
    /// Analytic, gate inputs used only if `|v| > 0` ([`Model::jv_guarded`]).
    Guarded,
    /// Analytic, skipping zero coefficients ([`Model::jv_coef_skip`]).
    CoefSkip,
}

/// Number of compartments.
pub const NC: usize = 45;
/// Number of gate variables per compartment.
pub const NG: usize = 6;
/// Total number of states.
pub const N: usize = NC * (1 + NG);

/// Parameters and the piecewise-constant stimulus of the toy model.
pub struct Model {
    /// Axial coupling (1/ms); the largest cable eigenvalue is about 4 d.
    pub d: f64,
    /// Leak conductance (1/ms).
    pub gl: f64,
    /// Leak reversal (mV).
    pub el: f64,
    /// Gate conductances.
    pub gbar: [f64; NG],
    /// Gate reversal potentials (mV).
    pub ek: [f64; NG],
    /// Half-activation voltages (mV).
    pub vh: [f64; NG],
    /// Activation slopes (mV).
    pub sl: [f64; NG],
    /// Gate time constants (ms); spans the stiff ratio.
    pub tau: [f64; NG],
    /// Current injected into compartment 0 for the running segment.
    pub stim: Cell<f64>,
    /// Amplitude of the stimulus pulse.
    pub amp: f64,
    /// Pulse start and end (ms).
    pub pulse: (f64, f64),
    /// Voltage threshold of root 0 (first compartment).
    pub thr0: f64,
    /// Voltage threshold of root 1 (distal compartment).
    pub thr1: f64,
    /// Axial edges `(parent[e], child[e])`; the chain is `e -> e + 1`.
    pub parent: Vec<usize>,
    pub child: Vec<usize>,
    /// Neighbour lists built from the edges; `None` keeps the chain's
    /// original arithmetic (bit-identical to the probe before trees).
    nbrs: Option<Vec<Vec<usize>>>,
}

impl Default for Model {
    fn default() -> Self {
        Self {
            d: 100.0,
            gl: 0.3,
            el: -65.0,
            gbar: [10.0, 8.0, 4.0, 0.3, 0.2, 0.5],
            ek: [50.0, -80.0, -80.0, -70.0, 20.0, -75.0],
            vh: [-40.0, -45.0, -50.0, -60.0, -30.0, -55.0],
            sl: [6.0, 8.0, 8.0, 10.0, 7.0, 9.0],
            // 0.01 ms .. 20 ms; with d = 100 the stiffness ratio is ~ 4e4.
            tau: [0.01, 0.1, 0.5, 2.0, 8.0, 20.0],
            stim: Cell::new(0.0),
            amp: 2000.0,
            pulse: (2.0, 2.3),
            thr0: -20.0,
            thr1: 0.0,
            parent: (0..NC - 1).collect(),
            child: (1..NC).collect(),
            nbrs: None,
        }
    }
}

fn ginf(m: &Model, k: usize, v: f64) -> f64 {
    1.0 / (1.0 + (-(v - m.vh[k]) / m.sl[k]).exp())
}

fn dginf(m: &Model, k: usize, v: f64) -> f64 {
    let s = ginf(m, k, v);
    s * (1.0 - s) / m.sl[k]
}

impl Model {
    /// The same model with the stimulus pulse switched off.
    pub fn quiet() -> Self {
        Self {
            amp: 0.0,
            ..Self::default()
        }
    }

    /// The model on a tree given by its edge list (`NC` compartments).
    pub fn with_edges(parent: Vec<usize>, child: Vec<usize>) -> Self {
        let mut nbrs = vec![Vec::new(); NC];
        for (&p, &c) in parent.iter().zip(&child) {
            nbrs[p].push(c);
            nbrs[c].push(p);
        }
        Self {
            parent,
            child,
            nbrs: Some(nbrs),
            ..Self::default()
        }
    }

    /// A 45-compartment tree: a 25-compartment trunk `0..25` with two
    /// 10-compartment branches `25..35` and `35..45`, both off compartment 5.
    pub fn tree() -> Self {
        let (mut parent, mut child) = (Vec::new(), Vec::new());
        let mut edge = |p: usize, c: usize| {
            parent.push(p);
            child.push(c);
        };
        for i in 0..24 {
            edge(i, i + 1);
        }
        for start in [25, 35] {
            edge(5, start);
            for i in start..start + 9 {
                edge(i, i + 1);
            }
        }
        Self::with_edges(parent, child)
    }

    /// The explicit Jacobian pattern: `core.rs::tree_jac_sparsity` on the
    /// model's edges (`chain_jac_sparsity` for the default chain).
    pub fn pattern(&self) -> (Vec<usize>, Vec<usize>) {
        tree_jac_sparsity(NC, &self.parent, &self.child)
    }

    /// The axial term of compartment `i` before the factor `d`.
    #[inline]
    fn axial(&self, x: &[f64], i: usize) -> f64 {
        match &self.nbrs {
            None => {
                let vm = if i > 0 { x[i - 1] } else { x[i] };
                let vp = if i + 1 < NC { x[i + 1] } else { x[i] };
                vm - 2.0 * x[i] + vp
            }
            Some(nb) => nb[i].iter().fold(0.0, |a, &j| a + (x[j] - x[i])),
        }
    }

    /// Resting state: V at the leak reversal, gates at steady state.
    pub fn rest(&self) -> Vec<f64> {
        let mut y = vec![0.0; N];
        for i in 0..NC {
            y[i] = self.el;
            for k in 0..NG {
                y[NC * (1 + k) + i] = ginf(self, k, self.el);
            }
        }
        y
    }

    /// Right-hand side `y = f(x)`.
    pub fn rhs(&self, x: &[f64], _t: f64, y: &mut [f64]) {
        for i in 0..NC {
            let mut dv = self.d * self.axial(x, i);
            dv += self.gl * (self.el - x[i]);
            for k in 0..NG {
                let g = x[NC * (1 + k) + i];
                dv += self.gbar[k] * g * (self.ek[k] - x[i]);
            }
            if i == 0 {
                dv += self.stim.get();
            }
            y[i] = dv;
        }
        for k in 0..NG {
            for i in 0..NC {
                let g = x[NC * (1 + k) + i];
                y[NC * (1 + k) + i] = (ginf(self, k, x[i]) - g) / self.tau[k];
            }
        }
    }

    /// Analytic Jacobian-vector product `y = J(x) v`.
    pub fn jv(&self, x: &[f64], _t: f64, v: &[f64], y: &mut [f64]) {
        for i in 0..NC {
            let mut acc = self.d * self.axial(v, i) - self.gl * v[i];
            for k in 0..NG {
                let g = x[NC * (1 + k) + i];
                acc -= self.gbar[k] * g * v[i];
                acc += self.gbar[k] * (self.ek[k] - x[i]) * v[NC * (1 + k) + i];
            }
            y[i] = acc;
        }
        for k in 0..NG {
            for i in 0..NC {
                let j = NC * (1 + k) + i;
                y[j] = (dginf(self, k, x[i]) * v[i] - v[j]) / self.tau[k];
            }
        }
    }

    /// `J v` with the input-sparsity shortcut `if v == 0.0 { skip }` on
    /// every gate term.  NaN == 0.0 is false, so NaN probes are NOT
    /// skipped: this shortcut is harmless for detection.
    pub fn jv_skip_zero(&self, x: &[f64], _t: f64, v: &[f64], y: &mut [f64]) {
        for i in 0..NC {
            let mut acc = self.d * self.axial(v, i) - self.gl * v[i];
            for k in 0..NG {
                acc -= self.gbar[k] * x[NC * (1 + k) + i] * v[i];
                let vg = v[NC * (1 + k) + i];
                if vg == 0.0 {
                    continue;
                }
                acc += self.gbar[k] * (self.ek[k] - x[i]) * vg;
            }
            y[i] = acc;
        }
        for k in 0..NG {
            for i in 0..NC {
                let j = NC * (1 + k) + i;
                let mut acc = -v[j];
                if v[i] != 0.0 {
                    acc += dginf(self, k, x[i]) * v[i];
                }
                y[j] = acc / self.tau[k];
            }
        }
    }

    /// The same shortcut written the other way round: a gate input is used
    /// only `if v.abs() > 0.0`.  Exact for every finite `v`, but NaN > 0 is
    /// false, so the probe silently drops all 6 * NC dV/dgate entries.
    pub fn jv_guarded(&self, x: &[f64], _t: f64, v: &[f64], y: &mut [f64]) {
        for i in 0..NC {
            let mut acc = self.d * self.axial(v, i) - self.gl * v[i];
            for k in 0..NG {
                acc -= self.gbar[k] * x[NC * (1 + k) + i] * v[i];
                let vg = v[NC * (1 + k) + i];
                if vg.abs() > 0.0 {
                    acc += self.gbar[k] * (self.ek[k] - x[i]) * vg;
                }
            }
            y[i] = acc;
        }
        for k in 0..NG {
            for i in 0..NC {
                let j = NC * (1 + k) + i;
                y[j] = (dginf(self, k, x[i]) * v[i] - v[j]) / self.tau[k];
            }
        }
    }

    /// `J v` skipping terms whose coefficient is exactly zero at `x`: what
    /// the probe sees then depends on the state it is evaluated at (y0).
    /// At V = E_k the dV/dg_k coefficient `gbar_k (E_k - V)` vanishes.
    pub fn jv_coef_skip(&self, x: &[f64], _t: f64, v: &[f64], y: &mut [f64]) {
        for i in 0..NC {
            let mut acc = self.d * self.axial(v, i) - self.gl * v[i];
            for k in 0..NG {
                acc -= self.gbar[k] * x[NC * (1 + k) + i] * v[i];
                let c = self.gbar[k] * (self.ek[k] - x[i]);
                if c != 0.0 {
                    acc += c * v[NC * (1 + k) + i];
                }
            }
            y[i] = acc;
        }
        for k in 0..NG {
            for i in 0..NC {
                let j = NC * (1 + k) + i;
                let c = dginf(self, k, x[i]);
                let mut acc = -v[j];
                if c != 0.0 {
                    acc += c * v[i];
                }
                y[j] = acc / self.tau[k];
            }
        }
    }

    /// Dispatch on [`Jv`].
    pub fn jv_by(&self, kind: Jv, x: &[f64], t: f64, v: &[f64], y: &mut [f64]) {
        match kind {
            Jv::Analytic => self.jv(x, t, v, y),
            Jv::Fd => self.jv_fd(x, t, v, y),
            Jv::SkipZero => self.jv_skip_zero(x, t, v, y),
            Jv::Guarded => self.jv_guarded(x, t, v, y),
            Jv::CoefSkip => self.jv_coef_skip(x, t, v, y),
        }
    }

    /// A start state on the edge of the domain: every V at E_K = -80 mV
    /// (the reversal of gates 1 and 2) and gates alternately exactly 0 and
    /// exactly 1.
    pub fn edge_state(&self) -> Vec<f64> {
        let mut y = vec![0.0; N];
        for i in 0..NC {
            y[i] = self.ek[1];
            for k in 0..NG {
                y[NC * (1 + k) + i] = ((i + k) % 2) as f64;
            }
        }
        y
    }

    /// Forward-difference Jacobian-vector product built from `rhs` only.
    ///
    /// The step is a fixed `1e-7 * (1 + max|x|)`; probe vectors are 0/1
    /// columns (or NaN-tainted columns during sparsity detection), so no
    /// norm of `v` is taken and structural zeros stay exact.
    pub fn jv_fd(&self, x: &[f64], t: f64, v: &[f64], y: &mut [f64]) {
        let xmax = x.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        let eps = 1e-7 * (1.0 + xmax);
        let xp: Vec<f64> = x.iter().zip(v).map(|(a, b)| a + eps * b).collect();
        let mut f0 = vec![0.0; N];
        self.rhs(x, t, &mut f0);
        self.rhs(&xp, t, y);
        for (a, b) in y.iter_mut().zip(&f0) {
            *a = (*a - b) / eps;
        }
    }

    /// Terminal events; both start negative and cross zero upwards.
    pub fn root(&self, x: &[f64], _t: f64, g: &mut [f64]) {
        g[0] = x[0] - self.thr0;
        g[1] = x[NC / 2] - self.thr1;
    }

    /// Stimulus current for a segment starting at `t`: a 0.3 ms pulse.
    pub fn stim_at(&self, t: f64) -> f64 {
        if (self.pulse.0..self.pulse.1).contains(&t) {
            self.amp
        } else {
            0.0
        }
    }
}
