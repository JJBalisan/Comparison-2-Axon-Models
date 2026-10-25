//! Prints the measurements quoted in REPORT.md.  Run with
//! `cargo run --release`.

use std::time::{Duration, Instant};

use diffsol_probe::model::{Model, N};
use diffsol_probe::solve::{dense, faer_dense, sparse, Cfg, Trace};

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
        let cfg = Cfg { fd, ..staged() };
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
}
