//! One test per probe item; see REPORT.md for the narrative.

use diffsol::matrix::sparsity::MatrixSparsity;
use diffsol::FaerSparseMat;
use diffsol_probe::model::{Jv, Model, N, NC};
use diffsol_probe::pattern::{chain_jac_sparsity, to_pairs, to_sparsity};
use diffsol_probe::solve::{
    dense, detected_pattern, faer_dense, sparse, sparse_explicit, Cfg, Trace,
};

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
        jv: Jv::Fd,
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

fn sorted_pattern(m: &Model) -> Vec<(usize, usize)> {
    let (rows, cols) = m.pattern();
    let mut nz = to_pairs(&rows, &cols);
    nz.sort_unstable();
    nz
}

fn full_run(jv: Jv) -> Cfg {
    Cfg {
        jv,
        terminal: false,
        record: true,
        ..staged()
    }
}

/// Item 4a: the (rows, cols) pattern converts to faer's symbolic matrix.
#[test]
fn item4a_pattern_converts_to_the_diffsol_sparsity() {
    let m = Model::default();
    let (rows, cols) = m.pattern();
    assert_eq!((rows.clone(), cols.clone()), chain_jac_sparsity(NC));
    assert_eq!(rows.len(), (3 * NC - 2) + 18 * NC);
    let s = to_sparsity::<FaerSparseMat<f64>>(N, &rows, &cols).unwrap();
    assert_eq!(s.compute_nnz(), rows.len());
    let back = MatrixSparsity::<FaerSparseMat<f64>>::indices(&s);
    assert_eq!(back, to_pairs(&rows, &cols), "CSC order is kept");
}

/// Item 4a: build_from_eqn with the supplied pattern reproduces the
/// auto-detected run bit for bit (same pattern, same colouring).
#[test]
fn item4a_explicit_pattern_matches_auto_detection() {
    let m = Model::default();
    let y0 = m.rest();
    assert_eq!(detected_pattern(&m, &y0, Jv::Fd), sorted_pattern(&m));
    for jv in [Jv::Analytic, Jv::Fd] {
        let cfg = full_run(jv);
        let a = sparse::solve(&m, &y0, &cfg).unwrap();
        let e = sparse_explicit::solve(&m, &y0, &cfg).unwrap();
        assert_eq!(a.jac_nnz, Some(943));
        assert_eq!(e.jac_nnz, Some(943));
        assert_eq!(a.steps, e.steps, "{jv:?}");
        assert_eq!(a.lu_setups, e.lu_setups, "{jv:?}");
        assert_eq!(a.t, e.t, "{jv:?}: step times");
        assert_eq!(a.y, e.y, "{jv:?}: states");
    }
}

/// Item 4b: `if v == 0.0 { skip }` is safe (NaN == 0.0 is false), but the
/// same shortcut written `if v.abs() > 0.0` hides every dV/dgate entry from
/// the NaN probe and the auto-detected run fails; the explicit pattern
/// makes the same Jv work.
#[test]
fn item4b_guarded_jv_undercounts_and_breaks_the_auto_path() {
    let m = Model::default();
    let y0 = m.rest();
    assert_eq!(detected_pattern(&m, &y0, Jv::SkipZero).len(), 943);
    assert_eq!(detected_pattern(&m, &y0, Jv::Guarded).len(), 943 - 6 * NC);
    let cfg = full_run(Jv::Guarded);
    assert!(sparse::solve(&m, &y0, &cfg).is_err());
    let e = sparse_explicit::solve(&m, &y0, &cfg).unwrap();
    let base = sparse_explicit::solve(&m, &y0, &full_run(Jv::Analytic)).unwrap();
    assert_eq!(e.steps, base.steps);
    assert!(max_abs_diff(&e.final_y, &base.final_y) < 1e-9);
}

/// Item 4b: branch-free Jv (analytic or FD) is detected the same at any y0
/// (0 * NaN is NaN); a zero-coefficient skip makes detection depend on y0.
#[test]
fn item4b_detection_depends_on_y0_only_through_branches() {
    let m = Model::default();
    let (rest, edge) = (m.rest(), m.edge_state());
    let want = sorted_pattern(&m);
    for jv in [Jv::Analytic, Jv::Fd] {
        assert_eq!(detected_pattern(&m, &rest, jv), want, "{jv:?} rest");
        assert_eq!(detected_pattern(&m, &edge, jv), want, "{jv:?} edge");
    }
    assert_eq!(detected_pattern(&m, &rest, Jv::CoefSkip), want);
    // V = E_K kills dV/dg for the two gates reversing at E_K
    assert_eq!(
        detected_pattern(&m, &edge, Jv::CoefSkip).len(),
        943 - 2 * NC
    );
    let cfg = full_run(Jv::CoefSkip);
    assert!(sparse::solve(&m, &edge, &cfg).is_err());
    assert!(sparse_explicit::solve(&m, &edge, &cfg).is_ok());
}

/// Item 4c: on a branched tree the detected pattern is
/// `tree_jac_sparsity`, and both paths complete with the same steps.
#[test]
fn item4c_tree_patterns_agree_and_the_solve_completes() {
    let m = Model::tree();
    let y0 = m.rest();
    let want = sorted_pattern(&m);
    assert_eq!(want.len(), 943, "a tree on n nodes has n - 1 edges");
    assert_ne!(want, sorted_pattern(&Model::default()));
    // compartment 5 couples to 4, 6 and both branch roots
    assert!(want.contains(&(25, 5)) && want.contains(&(35, 5)));
    for jv in [Jv::Analytic, Jv::Fd] {
        assert_eq!(detected_pattern(&m, &y0, jv), want, "{jv:?}");
        let cfg = full_run(jv);
        let a = sparse::solve(&m, &y0, &cfg).unwrap();
        let e = sparse_explicit::solve(&m, &y0, &cfg).unwrap();
        assert!((e.final_t - 12.0).abs() < 1e-9);
        assert_eq!(a.steps, e.steps);
        assert_eq!(a.final_y, e.final_y);
    }
    // the branches change the dynamics
    let c = sparse::solve(&Model::default(), &y0, &full_run(Jv::Analytic)).unwrap();
    let t = sparse::solve(&m, &y0, &full_run(Jv::Analytic)).unwrap();
    assert!(max_abs_diff(&c.final_y, &t.final_y) > 1e-3);
}
