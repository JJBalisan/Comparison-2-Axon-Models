//! One test per probe item; see REPORT.md for the narrative.

use diffsol_probe::model::{Model, N, NC};
use diffsol_probe::solve::{dense, faer_dense, sparse, Cfg, Trace};

fn run(cfg: &Cfg) -> (Model, Trace) {
    let m = Model::default();
    let y0 = m.rest();
    let tr = sparse::solve(&m, &y0, cfg).unwrap();
    (m, tr)
}

fn staged() -> Cfg {
    Cfg {
        breaks: vec![0.0, 2.0, 2.3, 12.0],
        ..Cfg::default()
    }
}

fn reference() -> Cfg {
    Cfg {
        rtol: 1e-12,
        atol: 1e-14,
        ..staged()
    }
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |m, (x, y)| m.max((x - y).abs()))
}

/// Item 1: RootFound stops at the root and interpolate places it.
#[test]
fn item1_terminal_root_is_found_and_placed() {
    let (m, tr) = run(&staged());
    let r = tr.root.expect("root expected");
    assert_eq!(r.index, 0, "V[0] crosses its threshold first");
    assert!(r.t > 2.0 && r.t <= r.step_end);
    // The solver stopped inside the step: nothing past the root was kept.
    assert!(tr.final_y.is_empty());
    let rel = ((r.y[0] - m.thr0) / m.thr0).abs();
    assert!(rel < 1e-10, "g(interpolate(root)) relative {rel:e}");
    // Both events start negative; the other one has not fired.
    assert!(r.y[NC / 2] < m.thr1);
}

/// Item 1: the root time agrees with a 1e-12 reference to solver accuracy.
#[test]
fn item1_root_time_matches_reference() {
    let (_, tr) = run(&staged());
    let (_, rf) = run(&reference());
    let (a, b) = (tr.root.unwrap(), rf.root.unwrap());
    let err = (a.t - b.t).abs() / b.t;
    assert!(err < 1e-5, "root time error {err:e}");
}

/// Item 2: every accepted step is available from the step() loop.
#[test]
fn item2_every_accepted_step_is_collected() {
    let cfg = Cfg {
        record: true,
        terminal: false,
        ..staged()
    };
    let (_, tr) = run(&cfg);
    assert_eq!(tr.t.len(), tr.steps);
    assert_eq!(tr.y.len(), tr.steps);
    assert!(tr.y.iter().all(|y| y.len() == N));
    assert!(tr.t.windows(2).all(|w| w[1] > w[0]));
    assert!((tr.t.last().unwrap() - 12.0).abs() < 1e-12);
}

/// Item 3: re-initialising at extra breakpoints leaves the answer unchanged.
#[test]
fn item3_reinit_matches_single_segment() {
    let one = Cfg {
        breaks: vec![0.0, 12.0],
        terminal: false,
        ..Cfg::default()
    };
    let many = Cfg {
        breaks: (0..=24).map(|i| 0.5 * i as f64).collect(),
        ..one.clone()
    };
    let quiet = Model::quiet();
    let y0 = quiet.rest();
    let a = sparse::solve(&quiet, &y0, &one).unwrap();
    let b = sparse::solve(&quiet, &y0, &many).unwrap();
    assert_eq!(b.segments, 24);
    assert!(max_abs_diff(&a.final_y, &b.final_y) < 1e-6);
}

/// Item 3: a stimulus change at a breakpoint is honoured.
#[test]
fn item3_stimulus_breakpoint_changes_dynamics() {
    let cfg = Cfg {
        terminal: false,
        ..staged()
    };
    let (_, with) = run(&cfg);
    let quiet = Model::quiet();
    let q = sparse::solve(&quiet, &quiet.rest(), &cfg).unwrap();
    assert!(max_abs_diff(&with.final_y, &q.final_y) > 1.0);
}

/// Item 4: sparse/dense and analytic/FD Jacobians solve the same problem.
#[test]
fn item4_backends_agree() {
    let cfg = Cfg {
        breaks: vec![0.0, 2.0, 2.3, 3.0],
        terminal: false,
        ..Cfg::default()
    };
    let m = Model::default();
    let y0 = m.rest();
    let base = sparse::solve(&m, &y0, &cfg).unwrap().final_y;
    let fd = Cfg {
        fd: true,
        ..cfg.clone()
    };
    for (name, y) in [
        ("sparse fd", sparse::solve(&m, &y0, &fd).unwrap().final_y),
        ("dense", dense::solve(&m, &y0, &cfg).unwrap().final_y),
        ("dense fd", dense::solve(&m, &y0, &fd).unwrap().final_y),
        ("faer", faer_dense::solve(&m, &y0, &cfg).unwrap().final_y),
    ] {
        let d = max_abs_diff(&base, &y);
        assert!(d < 1e-3, "{name} differs by {d:e}");
    }
}

/// Item 5: there is no max-step setting; hopping set_stop_time bounds h.
#[test]
fn item5_stop_time_hops_bound_the_step() {
    let cfg = Cfg {
        max_step: Some(0.005),
        record: true,
        terminal: false,
        ..staged()
    };
    let (_, tr) = run(&cfg);
    let mut prev = 0.0;
    for t in &tr.t {
        assert!(t - prev <= 0.005 + 1e-9, "step {} too long", t - prev);
        prev = *t;
    }
}

/// Item 5: dense output on a 0.005 grid recovers the spike peak.
#[test]
fn item5_dense_output_grid_finds_peak() {
    let hop = Cfg {
        max_step: Some(0.005),
        record: true,
        terminal: false,
        rtol: 1e-9,
        atol: 1e-11,
        ..staged()
    };
    let (_, rf) = run(&hop);
    let truth = rf.y.iter().fold(f64::MIN, |m, y| m.max(y[0]));
    let grid = Cfg {
        grid: Some(0.005),
        terminal: false,
        ..staged()
    };
    let (_, g) = run(&grid);
    let peak = g.gy.iter().fold(f64::MIN, |m, y| m.max(y[0]));
    assert!((peak - truth).abs() < 0.05, "peak {peak} vs {truth}");
}
