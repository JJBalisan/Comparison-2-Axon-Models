//! diffsol-real command line.
//!
//! ```text
//! diffsol-real <case.json> <out.json> [--fd probe|scaled]
//!     [--max-step case|none|<ms>] [--extra-cuts <count>,<from>,<to>]
//!     [--no-stops] [--repeat <n>]
//! diffsol-real --bench <case.json> [<t>]
//! ```
//!
//! A run writes the trace of the median-wall repeat (every accepted step:
//! t, soma, axon column) with the step, evaluation and solver counts.
//! `--bench` times `MultiRhs::eval` (10k calls on y0 at time t, default
//! 5.1 ms, median of 9 batches).

use std::error::Error;
use std::fs;

use diffsol_real::{bench_eval, median_wall, solve, Case, Cfg, Fd};

fn load(path: &str) -> Result<Case, Box<dyn Error>> {
    Ok(serde_json::from_str(&fs::read_to_string(path)?)?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("--bench") {
        let case = load(&args[1])?;
        let t = args.get(2).map_or(Ok(5.1), |s| s.parse())?;
        let us = bench_eval(&case, t, &case.y0, 10_000, 9);
        println!("{}: {us:.3} us per eval (t = {t})", case.name);
        return Ok(());
    }
    let case = load(&args[0])?;
    let out = &args[1];
    let mut cfg = Cfg::for_case(&case);
    let mut repeat = 1;
    let mut it = args[2..].iter();
    while let Some(flag) = it.next() {
        let val = || {
            it.clone()
                .next()
                .cloned()
                .ok_or(format!("{flag} needs a value"))
        };
        match flag.as_str() {
            "--fd" => {
                cfg.fd = match val()?.as_str() {
                    "probe" => Fd::Probe,
                    "scaled" => Fd::Scaled,
                    other => return Err(format!("unknown --fd {other}").into()),
                };
                it.next();
            }
            "--max-step" => {
                cfg.max_step = match val()?.as_str() {
                    "case" => Some(case.max_step),
                    "none" => None,
                    x => Some(x.parse()?),
                };
                it.next();
            }
            "--extra-cuts" => {
                let v: Vec<f64> = val()?
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<_, _>>()?;
                let (k, a, b) = (v[0] as usize, v[1], v[2]);
                cfg.extra_cuts = (1..=k)
                    .map(|i| a + (b - a) * i as f64 / (k + 1) as f64)
                    .collect();
                it.next();
            }
            "--no-stops" => cfg.stops = false,
            "--repeat" => {
                repeat = val()?.parse()?;
                it.next();
            }
            other => return Err(format!("unknown flag {other}").into()),
        }
    }
    let (tr, wall) = median_wall(repeat, || solve(&case, &cfg));
    let mut tr = tr?;
    tr.wall_us = wall.as_secs_f64() * 1e6;
    fs::write(out, serde_json::to_string(&tr)?)?;
    println!(
        "{}: {} steps ({} accepted), {} evals ({} Jv), {} LU, nnz {}, \
         wall {:.0} us (build {:.0}), t_last {} (raw {}), inexact ends {}, root {:?}",
        tr.name,
        tr.steps,
        tr.stats.accepted_steps,
        tr.evals,
        tr.jv_calls,
        tr.stats.lu_setups,
        tr.jac_nnz,
        tr.wall_us,
        tr.build_us,
        tr.t.last().copied().unwrap_or(f64::NAN),
        tr.t_last_raw,
        tr.inexact_segment_ends,
        tr.root.as_ref().map(|r| (r.t, r.index)),
    );
    Ok(())
}
