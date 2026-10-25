//! Toy stiff cable model with the sparsity of msoaxon's active-full system.
//!
//! State layout (315 = 7 * 45): voltage block `V[0..45]` followed by six
//! gate blocks `g_k[0..45]`.  The voltage Jacobian block is tridiagonal
//! (a sealed cable) and each gate block is diagonal, with each gate driven
//! only by the voltage of its own compartment, exactly like
//! `multi._chain_jac_sparsity`.

use std::cell::Cell;

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
            let vm = if i > 0 { x[i - 1] } else { x[i] };
            let vp = if i + 1 < NC { x[i + 1] } else { x[i] };
            let mut dv = self.d * (vm - 2.0 * x[i] + vp);
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
            let vm = if i > 0 { v[i - 1] } else { v[i] };
            let vp = if i + 1 < NC { v[i + 1] } else { v[i] };
            let mut acc = self.d * (vm - 2.0 * v[i] + vp) - self.gl * v[i];
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
