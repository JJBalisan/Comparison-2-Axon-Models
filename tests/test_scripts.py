"""End-to-end runs of the analysis scripts in their quick modes.

--quick coarsens grids and tolerances without changing code paths, so each run
exercises the same code as a full one in seconds. The assertions are invariants
that survive the coarse settings (orderings, broad ranges), not exact values;
the full-mode numbers are documented in PYTHON.md.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# about 60 s locally, ~4 min on a CI runner; the scripts don't depend on the Python
# version, so CI runs these on 3.12 only (deselect locally with -m "not slow")
pytestmark = pytest.mark.slow


def run_script(name, *args):
    r = subprocess.run([sys.executable, str(REPO / "scripts" / name), *map(str, args)],
                       capture_output=True, text=True, env={"MPLBACKEND": "Agg", "PATH": ""})
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.mark.parametrize("flags, rin_range", [
    ((), (2.0, 3.5)),  # totals conserved (default): 2.6 MOhm
    (("--klt-from-soma", "--dendrite-ra", 200), (4.0, 6.0)),  # 5.0 MOhm
], ids=["default", "klt-from-soma-ra200"])
def test_dendrites_script(tmp_path, flags, rin_range):
    run_script("dendrites.py", "--quick", "--out-dir", tmp_path, *flags)
    res = json.loads((tmp_path / "dendrites.json").read_text())
    assert (tmp_path / "dendrites.png").stat().st_size > 0

    assert rin_range[0] < res["A_rin_tau"]["dendritic"]["rin_steady"] < rin_range[1]
    B = res["B_epsp"]  # attenuation along the cable, sharpening by dendritic KLT
    assert B["dendritic, EPSG distal (180 um)"]["amp"] < B["dendritic, EPSG mid-dendrite (100 um)"]["amp"] \
        < B["dendritic, EPSG at soma"]["amp"]
    assert B["mid-dendrite, dendritic KLT removed"]["half_width_us"] \
        > B["dendritic, EPSG mid-dendrite (100 um)"]["half_width_us"]
    C = res["C_threshold"]  # Scott et al 2010: bilateral beats unilateral
    assert C["bilateral (one EPSG per dendrite)"] < C["unilateral (both on the lateral dendrite)"]
    assert math.isclose(C["both at the soma, lumped model"], 48.0, rel_tol=0.02)
    assert 0.9 < res["D_summation"]["ratio"] < 1.1  # near-linear summation
    for widths in res["E_window_bilateral_dendritic_us"].values():
        assert 50 < widths["margin_0.03"] < 400
    G = res["G_somatic_spike"]
    assert math.isfinite(G["rheobase_pA"])
    assert G["amplitudes_mV"]["3.0"] > G["amplitudes_mV"]["1.5"] > 5


def test_coincidence_window_script(tmp_path):
    run_script("coincidence_window.py", "--quick", "--trials", 8, "--out-dir", tmp_path)
    res = json.loads((tmp_path / "coincidence_window.json").read_text())
    assert (tmp_path / "coincidence_window.png").stat().st_size > 0

    rows = {(r["model"], r["epsg"]): r for r in res["half_widths_us"]}
    assert len(rows) == 6  # 3 model configurations x 2 EPSG kinetics
    lumped = rows["45-compartment", "model EPSG (decay 0.18 ms)"]
    assert math.isclose(lumped["threshold0"], 48.0, rel_tol=0.02)
    for r in rows.values():
        assert 50 < r["margin_0.03"] < 400
        assert r["margin_0.005"] < r["margin_0.03"]  # larger inputs stay above threshold longer
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


def test_making_threshold_graphs_script(tmp_path):
    out = run_script("making_threshold_graphs.py", "--only", "sine", "--points", 1,
                     "--out-dir", tmp_path)
    # one sweep point, BinarySearch.m's rounding to 10 (unrounded: 11054.08 / 11541.14)
    assert "sine: multi=[11050.0]" in out and "two  =[11540.0]" in out
    assert (tmp_path / "sine_thresholds_cpt3.png").stat().st_size > 0
