//! Prints the measurements quoted in REPORT.md.  Run with
//! `cargo run --release`.

use std::time::{Duration, Instant};

use diffsol_probe::lu::time_lu;
use diffsol_probe::model::{Jv, Model, N};
use diffsol_probe::solve::{
    dense, detected_pattern, faer_dense, jacobian_at, sparse, sparse_explicit, Cfg, Trace,
};

const RUNS: usize = 9;

type Solve = fn(&Model, &[f64], &Cfg) -> Result<Trace, diffsol::DiffsolError>;

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn staged() -> Cfg {
    Cfg {
        breaks: vec![0.0, 2.0, 2.3, 12.0],
        terminal: false,
        ..Cfg::default()
    }
}

fn timed(f: Solve, m: &Model, cfg: &Cfg) -> (f64, Trace) {
    let y0 = m.rest();
    let mut walls = Vec::new();
    let mut last = None;
    for _ in 0..RUNS {
        let tr = f(m, &y0, cfg).unwrap();
        walls.push(tr.wall.as_secs_f64());
        last = Some(tr);
    }
    (median(walls), last.unwrap())
}

fn item1(m: &Model) {
    println!("== item 1: terminal root");
    let y0 = m.rest();
    let cfg = Cfg {
        terminal: true,
        ..staged()
    };
    let tr = sparse::solve(m, &y0, &cfg).unwrap();
    let r = tr.root.unwrap();
    let rel = ((r.y[0] - m.thr0) / m.thr0).abs();
    let tight = Cfg {
        rtol: 1e-12,
        atol: 1e-14,
        ..cfg
    };
    let rf = sparse::solve(m, &y0, &tight).unwrap().root.unwrap();
    println!(
        "root index {} t {:.12} step_end {:.12}",
        r.index, r.t, r.step_end
    );
    println!("|g(interpolate(t_root))|/|thr| = {rel:.2e}");
    println!(
        "reference (rtol 1e-12) t {:.12}; rel time error {:.2e}",
        rf.t,
        (r.t - rf.t).abs() / rf.t
    );
}

fn item2(m: &Model) {
    println!("== item 2: every accepted step");
    let cfg = Cfg {
        record: true,
        ..staged()
    };
    let tr = sparse::solve(m, &m.rest(), &cfg).unwrap();
    println!(
        "{} steps recorded over t in [0, 12], each with {} states",
        tr.t.len(),
        tr.y[0].len()
    );
}

fn item3(m: &Model) {
    println!("== item 3: re-initialisation cost (quiet model, t in [0, 12])");
    let quiet = Model::quiet();
    let y0 = quiet.rest();
    for (name, f, nseg) in [
        ("sparse", sparse::solve as Solve, 240_usize),
        ("dense nalgebra", dense::solve as Solve, 24),
    ] {
        let one = Cfg {
            breaks: vec![0.0, 12.0],
            terminal: false,
            ..Cfg::default()
        };
        let many = Cfg {
            breaks: (0..=nseg).map(|i| 12.0 * i as f64 / nseg as f64).collect(),
            ..one.clone()
        };
        let mut a = Vec::new();
        let mut b = Vec::new();
        let (mut sa, mut sb) = (0, 0);
        for _ in 0..5 {
            let t1 = f(&quiet, &y0, &one).unwrap();
            let t2 = f(&quiet, &y0, &many).unwrap();
            a.push(t1.wall.as_secs_f64());
            b.push(t2.wall.as_secs_f64());
            (sa, sb) = (t1.steps, t2.steps);
        }
        let (a, b) = (median(a), median(b));
        println!(
            "{name}: 1 segment {:.0} us / {sa} steps; {nseg} segments {:.0} us \
             / {sb} steps; extra per re-init {:.1} us",
            a * 1e6,
            b * 1e6,
            (b - a) / (nseg - 1) as f64 * 1e6
        );
    }
    let _ = m;
}

fn item4(m: &Model) {
    println!("== item 4: per-step cost, N = {N}, staged run to t = 12");
    let y0 = m.rest();
    let mut buf = vec![0.0; N];
    let reps = 20_000;
    let t0 = Instant::now();
    for _ in 0..reps {
        m.rhs(&y0, 0.0, &mut buf);
        std::hint::black_box(&buf);
    }
    let rhs_us = t0.elapsed().as_secs_f64() / reps as f64 * 1e6;
    println!("model.rhs: {rhs_us:.2} us per call");
    println!("config, steps, setups, median wall (us), us/step");
    let probes: [(&str, Solve, bool); 6] = [
        ("sparse+FaerSparseLU, analytic Jv", sparse::solve, false),
        ("sparse+FaerSparseLU, FD Jv", sparse::solve, true),
        (
            "dense NalgebraMat+NalgebraLU, analytic Jv",
            dense::solve,
            false,
        ),
        ("dense NalgebraMat+NalgebraLU, FD Jv", dense::solve, true),
        (
            "dense FaerMat+FaerLU, analytic Jv",
            faer_dense::solve,
            false,
        ),
        ("dense FaerMat+FaerLU, FD Jv", faer_dense::solve, true),
    ];
    for (name, f, fd) in probes {
        let jv = if fd { Jv::Fd } else { Jv::Analytic };
        let cfg = Cfg { jv, ..staged() };
        let (w, tr) = timed(f, m, &cfg);
        println!(
            "{name}: {} steps, {} setups, {:.0} us, {:.2} us/step",
            tr.steps,
            tr.lu_setups,
            w * 1e6,
            w * 1e6 / tr.steps as f64
        );
        if let Some(nnz) = tr.jac_nnz {
            println!("  detected Jacobian nnz {nnz} of {}", N * N);
        }
    }
}

fn peak(tr: &Trace, grid: bool) -> f64 {
    let ys = if grid { &tr.gy } else { &tr.y };
    ys.iter().fold(f64::MIN, |a, y| a.max(y[0]))
}

fn item5(m: &Model) {
    println!("== item 5: max step");
    let y0 = m.rest();
    let truth = sparse::solve(
        m,
        &y0,
        &Cfg {
            max_step: Some(0.005),
            record: true,
            rtol: 1e-10,
            atol: 1e-12,
            ..staged()
        },
    )
    .unwrap();
    println!("reference peak V[0] {:.6}", peak(&truth, false));
    let runs: [(&str, Cfg, bool); 4] = [
        (
            "default (no cap), steps only",
            Cfg {
                record: true,
                ..staged()
            },
            false,
        ),
        (
            "cap 0.005 via set_stop_time, re-armed every step",
            Cfg {
                record: true,
                max_step: Some(0.005),
                hop_always: true,
                ..staged()
            },
            false,
        ),
        (
            "cap 0.005 via set_stop_time, only when h > cap",
            Cfg {
                record: true,
                max_step: Some(0.005),
                ..staged()
            },
            false,
        ),
        (
            "default + interpolate grid 0.005",
            Cfg {
                grid: Some(0.005),
                ..staged()
            },
            true,
        ),
    ];
    for (name, cfg, grid) in runs {
        let mut w = Vec::new();
        let mut last = None;
        for _ in 0..RUNS {
            let t = sparse::solve(m, &y0, &cfg).unwrap();
            w.push(t.wall.as_secs_f64());
            last = Some(t);
        }
        let tr = last.unwrap();
        println!(
            "{name}: peak {:.6}, {} steps, {:.0} us",
            peak(&tr, grid),
            tr.steps,
            median(w) * 1e6
        );
    }
    let _: Option<Duration> = None;
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |m, (x, y)| m.max((x - y).abs()))
}

/// One line per run: steps, set-ups, nnz, build and median solve time.
fn report(name: &str, f: Solve, m: &Model, y0: &[f64], cfg: &Cfg) -> Option<Trace> {
    let mut walls = Vec::new();
    let mut builds = Vec::new();
    let mut last = None;
    for _ in 0..RUNS {
        match f(m, y0, cfg) {
            Ok(tr) => {
                walls.push(tr.wall.as_secs_f64());
                builds.push(tr.build.as_secs_f64());
                last = Some(tr);
            }
            Err(e) => {
                println!("{name}: FAILED: {e}");
                return None;
            }
        }
    }
    let tr = last.unwrap();
    let w = median(walls) * 1e6;
    println!(
        "{name}: nnz {:?}, {} steps, {} setups, build {:.0} us, solve {:.0} us, \
         {:.2} us/step",
        tr.jac_nnz,
        tr.steps,
        tr.lu_setups,
        median(builds) * 1e6,
        w,
        w / tr.steps as f64
    );
    Some(tr)
}

fn item4a(m: &Model) {
    println!("== item 4a: explicit pattern (build_from_eqn) vs auto-detected");
    let y0 = m.rest();
    for jv in [Jv::Analytic, Jv::Fd] {
        let cfg = Cfg {
            jv,
            record: true,
            ..staged()
        };
        let a = report(&format!("auto     {jv:?}"), sparse::solve, m, &y0, &cfg).unwrap();
        let e = report(
            &format!("explicit {jv:?}"),
            sparse_explicit::solve,
            m,
            &y0,
            &cfg,
        )
        .unwrap();
        let same_t = a.t == e.t;
        let same_y = a.y == e.y;
        println!(
            "  same step times bitwise {same_t}, same states bitwise {same_y}, \
             max |final diff| {:.1e}",
            max_abs_diff(&a.final_y, &e.final_y)
        );
    }
}

fn item4b(m: &Model) {
    println!("== item 4b: what NaN-probe detection sees (explicit pattern 943)");
    let (rows, _) = m.pattern();
    println!("explicit pattern length {}", rows.len());
    for (yname, y0) in [
        ("rest", m.rest()),
        ("edge (V = E_K, gates 0/1)", m.edge_state()),
    ] {
        for jv in [
            Jv::Analytic,
            Jv::Fd,
            Jv::SkipZero,
            Jv::Guarded,
            Jv::CoefSkip,
        ] {
            let nz = detected_pattern(m, &y0, jv);
            println!("  y0 {yname}, Jv {jv:?}: detected nnz {}", nz.len());
        }
    }
    println!("-- consequence, staged run t = 0..12 from rest");
    let y0 = m.rest();
    let tight = Cfg {
        rtol: 1e-10,
        atol: 1e-12,
        ..staged()
    };
    let truth = sparse_explicit::solve(m, &y0, &tight).unwrap().final_y;
    for (name, f, jv) in [
        ("auto     Analytic", sparse::solve as Solve, Jv::Analytic),
        ("auto     Guarded ", sparse::solve as Solve, Jv::Guarded),
        (
            "explicit Guarded ",
            sparse_explicit::solve as Solve,
            Jv::Guarded,
        ),
    ] {
        let cfg = Cfg { jv, ..staged() };
        if let Some(tr) = report(name, f, m, &y0, &cfg) {
            println!(
                "  max |final - rtol1e-10 reference| {:.2e} mV",
                max_abs_diff(&tr.final_y, &truth)
            );
        }
    }
    println!("-- consequence, staged run t = 0..12 from the edge state");
    let y0 = m.edge_state();
    let truth = sparse_explicit::solve(m, &y0, &tight).unwrap().final_y;
    for (name, f, jv) in [
        ("auto     Analytic", sparse::solve as Solve, Jv::Analytic),
        ("auto     CoefSkip", sparse::solve as Solve, Jv::CoefSkip),
        (
            "explicit CoefSkip",
            sparse_explicit::solve as Solve,
            Jv::CoefSkip,
        ),
    ] {
        let cfg = Cfg { jv, ..staged() };
        if let Some(tr) = report(name, f, m, &y0, &cfg) {
            println!(
                "  max |final - rtol1e-10 reference| {:.2e} mV",
                max_abs_diff(&tr.final_y, &truth)
            );
        }
    }
}

fn item4c() {
    println!("== item 4c: dendrite tree (25 trunk + 2 x 10 branches off 5)");
    let chain = Model::default();
    let tree = Model::tree();
    for (name, m) in [("chain", &chain), ("tree ", &tree)] {
        let y0 = m.rest();
        let (rows, cols) = m.pattern();
        let mut want: Vec<(usize, usize)> = rows.iter().copied().zip(cols).collect();
        want.sort_unstable();
        let got = detected_pattern(m, &y0, Jv::Fd);
        println!(
            "{name}: explicit pattern {} entries, auto-detected (FD Jv) {}, identical {}",
            want.len(),
            got.len(),
            want == got
        );
        for jv in [Jv::Analytic, Jv::Fd] {
            let cfg = Cfg { jv, ..staged() };
            let a = report(
                &format!("  {name} auto     {jv:?}"),
                sparse::solve,
                m,
                &y0,
                &cfg,
            );
            let e = report(
                &format!("  {name} explicit {jv:?}"),
                sparse_explicit::solve,
                m,
                &y0,
                &cfg,
            );
            if let (Some(a), Some(e)) = (a, e) {
                println!(
                    "    auto vs explicit: steps {} / {}, max |final diff| {:.1e}",
                    a.steps,
                    e.steps,
                    max_abs_diff(&a.final_y, &e.final_y)
                );
            }
        }
        let jac = jacobian_at(m, &y0, Jv::Analytic);
        for c in [1e-3, 1e-1] {
            let lu = time_lu(&jac, c, 201);
            println!(
                "  {name} FaerSparseLU on J - I/{c}: symbolic {:.1} us, numeric {:.1} us, \
                 solve {:.1} us",
                lu.symbolic_us, lu.numeric_us, lu.solve_us
            );
        }
    }
}

fn main() {
    let m = Model::default();
    let which = std::env::args().nth(1).unwrap_or_default();
    let all = which.is_empty();
    if all || which == "1" {
        item1(&m);
    }
    if all || which == "2" {
        item2(&m);
    }
    if all || which == "3" {
        item3(&m);
    }
    if all || which == "4" {
        item4(&m);
    }
    if all || which == "5" {
        item5(&m);
    }
    if all || which == "4a" {
        item4a(&m);
    }
    if all || which == "4b" {
        item4b(&m);
    }
    if all || which == "4c" {
        item4c();
    }
}
