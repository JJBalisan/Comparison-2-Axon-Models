"""End-to-end runs of the analysis scripts in their quick modes.

--quick coarsens grids and tolerances without changing code paths, so each run
exercises the same code as a full one in seconds. The assertions are invariants
that survive the coarse settings (orderings, broad ranges), not exact values;
the full-mode numbers are documented in PYTHON.md.
"""

import json
import math
import pathlib
import re
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]

# The scripts don't depend on the Python version, so CI runs these on 3.12
# only (deselect locally with -m "not slow"), in exact mode (no Rust
# extension). na_availability's --quick (pooled; no searches where the 45-
# compartment rest is unsteady, which cost minutes) took 47-56 s wall on 4
# workers, 90 CPU-s, with this machine at a load average of ~215.
pytestmark = pytest.mark.slow


def run_script(name, *args):
    r = subprocess.run(
        [sys.executable, str(REPO / "scripts" / name), *map(str, args)],
        capture_output=True,
        text=True,
        env={"MPLBACKEND": "Agg", "PATH": ""},
        check=False,
    )
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.mark.parametrize(
    "flags, rin_range",
    [
        ((), (2.0, 3.5)),  # totals conserved (default): 2.6 MOhm
        (("--klt-from-soma", "--dendrite-ra", 200), (4.0, 6.0)),  # 5.0 MOhm
    ],
    ids=["default", "klt-from-soma-ra200"],
)
def test_dendrites_script(tmp_path, flags, rin_range):
    run_script("dendrites.py", "--quick", "--out-dir", tmp_path, *flags)
    res = json.loads((tmp_path / "dendrites.json").read_text())
    assert (tmp_path / "dendrites.png").stat().st_size > 0

    assert (
        rin_range[0]
        < res["A_rin_tau"]["dendritic"]["rin_steady"]
        < rin_range[1]
    )
    B = res[
        "B_epsp"
    ]  # attenuation along the cable, sharpening by dendritic KLT
    assert (
        B["dendritic, EPSG distal (180 um)"]["amp"]
        < B["dendritic, EPSG mid-dendrite (100 um)"]["amp"]
        < B["dendritic, EPSG at soma"]["amp"]
    )
    assert (
        B["mid-dendrite, dendritic KLT removed"]["half_width_us"]
        > B["dendritic, EPSG mid-dendrite (100 um)"]["half_width_us"]
    )
    C = res["C_threshold"]  # Scott et al 2010: bilateral beats unilateral
    assert (
        C["bilateral (one EPSG per dendrite)"]
        < C["unilateral (both on the lateral dendrite)"]
    )
    assert math.isclose(C["both at the soma, lumped model"], 48.0, rel_tol=0.02)
    assert 0.9 < res["D_summation"]["ratio"] < 1.1  # near-linear summation
    for widths in res["E_window_bilateral_dendritic_us"].values():
        assert 50 < widths["margin_0.03"] < 400
    G = res["G_somatic_spike"]
    assert math.isfinite(G["rheobase_pA"])
    assert G["amplitudes_mV"]["3.0"] > G["amplitudes_mV"]["1.5"] > 5


def test_coincidence_window_script(tmp_path):
    run_script(
        "coincidence_window.py", "--quick", "--trials", 8, "--out-dir", tmp_path
    )
    res = json.loads((tmp_path / "coincidence_window.json").read_text())
    assert (tmp_path / "coincidence_window.png").stat().st_size > 0

    rows = {(r["model"], r["epsg"]): r for r in res["half_widths_us"]}
    assert len(rows) == 6  # 3 model configurations x 2 EPSG kinetics
    lumped = rows["45-compartment", "model EPSG (decay 0.18 ms)"]
    assert math.isclose(lumped["threshold0"], 48.0, rel_tol=0.02)
    for r in rows.values():
        assert 50 < r["margin_0.03"] < 400
        assert (
            r["margin_0.005"] < r["margin_0.03"]
        )  # larger inputs stay above threshold longer
    noisy = res["noisy_check"]
    assert 0 < noisy["peak_probability"] <= 1
    # nan when the few quick trials never fall below half: must still be written
    assert isinstance(noisy["half_width_us"], float)


def test_somatic_spike_script(tmp_path):
    run_script("somatic_spike.py", "--quick", "--out-dir", tmp_path)
    res = json.loads((tmp_path / "somatic_spike.json").read_text())
    assert (tmp_path / "somatic_spike.png").stat().st_size > 0

    assert len(res) == 3
    control = res["45-compartment (msoAxon.m)"]
    dtx = res["45-compartment, KLT blocked at soma + AIS (dendrotoxin)"]
    for case in res.values():
        amps = [s["amplitude"] for s in case["spikes"]]
        assert [s["multiple"] for s in case["spikes"]] == [1.5, 3.0]
        assert amps[1] > amps[0] > 0  # graded with the stimulus
    # Scott et al 2005: blocking KLT lowers rheobase and enlarges the spike
    assert dtx["rheobase_pA"] < control["rheobase_pA"]
    assert dtx["spikes"][0]["amplitude"] > control["spikes"][0]["amplitude"]


def test_na_availability_script(tmp_path):
    run_script("na_availability.py", "--quick", "--out-dir", tmp_path)
    res = json.loads((tmp_path / "na_availability.json").read_text())

    assert set(res) == {"45-compartment", "2-compartment"}
    for cells in res.values():
        assert [(c["rest_mV"], c["vhalf_mV"]) for c in cells] == [
            (-68.0, -62.5),
            (-68.0, -77.0),
            (-58.0, -62.5),
            (-58.0, -77.0),
        ]
        # the issue's table: 1 - hinf at rest with the model's own slope
        got = [round(c["inactivated"], 2) for c in cells]
        assert got == [0.33, 0.76, 0.64, 0.92]
        for c in cells:
            # Rin is read from a step only where the rest is steady
            # (RIPPLE_MAX = 0.01 mV in the script)
            assert (c["rin_steady_MOhm"] is None) == (
                c["rest_ripple_mV"] >= 0.01
            )
            if c["rheobase_pA"] is None:
                continue  # --quick skips the searches at an unsteady rest
            assert 50 < c["window_3pct_us"] < 400
            assert 100 < c["epsp_half_width_us"] < 2000
        # less available Na at rest takes more current to fire
        for lo, hi in ((0, 1), (2, 3)):
            if cells[lo]["rheobase_pA"] and cells[hi]["rheobase_pA"]:
                assert cells[hi]["rheobase_pA"] > cells[lo]["rheobase_pA"]
    # the default cells are the published ones (PYTHON.md: 3.8 nA, 2.4 MOhm)
    for model in res.values():
        assert math.isclose(model[0]["rheobase_pA"], 3800, rel_tol=0.03)
        assert math.isclose(model[0]["rin_steady_MOhm"], 2.4, rel_tol=0.05)
    # the 2-compartment model is steady at both rests
    assert all(c["rin_steady_MOhm"] for c in res["2-compartment"])


def test_making_threshold_graphs_script(tmp_path):
    out = run_script(
        "making_threshold_graphs.py",
        "--only",
        "sine",
        "--points",
        1,
        "--out-dir",
        tmp_path,
    )
    # one sweep point, BinarySearch.m's rounding to 10 (unrounded: 11054.08 /
    # 11541.14, see test_analysis.py). One flipped search decision (~1.8) can
    # carry 11054.08 past 11055 and round it to 11060, so compare within one
    # rounding step (rel 1e-3) rather than the printed text
    multi = re.search(r"sine: multi=\[([0-9.]+)\]", out)
    two = re.search(r"two  =\[([0-9.]+)\]", out)
    assert float(multi[1]) == pytest.approx(11054.08, rel=1e-3)
    assert float(two[1]) == pytest.approx(11541.14, rel=1e-3)
    assert (tmp_path / "sine_thresholds_cpt3.png").stat().st_size > 0


def test_phase0_harness_script(tmp_path):
    out_json = tmp_path / "phase0.json"
    out = run_script("phase0_harness.py", "--mode", "exact", "--json", out_json)
    assert "| exact |" in out
    res = json.loads(out_json.read_text())["modes"]["exact"]
    # the exact path produces the golden traces: well inside the tolerance
    assert res["max_margin"] < 1
    assert res["t_last_all_equal"]
    assert res["threshold_s"] > 0 < res["window_s"]
    assert res["steps_epsgpair_20ms"] > 100
    assert res["steps_step_110ms"] > 100
