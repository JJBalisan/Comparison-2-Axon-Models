//! The real model under diffsol: the golden EPSGpair run (fixture dumped by
//! `dump_cases.py --fixture`) reaches t_end with scipy's final state.

use diffsol_real::{core, solve, Case, Cfg, Fd};

fn fixture() -> Case {
    let text = include_str!("../fixtures/epsgpair.json");
    serde_json::from_str(text).expect("fixture parses")
}

fn scipy_final(case: &Case) -> Vec<f64> {
    case.scipy["final"]
        .as_array()
        .expect("final state")
        .iter()
        .map(|v| v.as_f64().expect("number"))
        .collect()
}

#[test]
fn epsgpair_reaches_t_end_with_scipys_final_state() {
    let case = fixture();
    for fd in [Fd::Probe, Fd::Scaled] {
        let cfg = Cfg {
            fd,
            ..Cfg::for_case(&case)
        };
        let tr = solve(&case, &cfg).expect("solve");
        assert_eq!(tr.t.last().copied(), Some(case.t_end));
        assert_eq!(tr.t_last_raw, case.t_end);
        assert!(tr.root.is_none());
        let want = scipy_final(&case);
        let n = want.len() / 7;
        for (i, (got, w)) in tr.final_y.iter().zip(&want).take(n).enumerate() {
            let d = (got - w).abs();
            assert!(d < 1e-3, "V[{i}] off by {d} mV ({fd:?})");
        }
        // every accepted step is output, strictly increasing
        assert!(tr.t.windows(2).all(|w| w[1] > w[0]));
        assert_eq!(tr.steps, tr.stats.accepted_steps + 1);
    }
}

#[test]
fn detected_sparsity_is_the_chain_pattern() {
    let case = fixture();
    let tr = solve(&case, &Cfg::for_case(&case)).expect("solve");
    assert_eq!(tr.jac_nnz, core::chain_jac_sparsity(45).0.len());
    assert_eq!(tr.jac_nnz, 943);
}

#[test]
fn spike_event_stops_at_the_criterion() {
    let mut case = fixture();
    case.stops = Some(diffsol_real::Stops {
        axon_col: 2,
        factor: 10.0,
        settled: false,
        t_off: 8.8,
        settled_fraction: 0.5,
    });
    let tr = solve(&case, &Cfg::for_case(&case)).expect("solve");
    let root = tr.root.expect("the pair spikes");
    assert_eq!(root.index, 0);
    assert!((root.y[2] - root.y[0] - 10.0).abs() < 1e-9);
    assert!(root.t < root.step_end + 1e-12);
    assert_eq!(tr.t.last().copied(), Some(root.t));
    assert!(root.t > 5.0 && root.t < 6.0);
}
