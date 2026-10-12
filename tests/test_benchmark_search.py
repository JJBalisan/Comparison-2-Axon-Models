"""benchmark_search.py --check must compare every yes/no the searches make.

The EPSG-pair, bilateral-dendrite and window searches decide through
_dispatch.response (imported by name into coincidence.py), not spikes(), so
--check has to wrap response too or it reports "mismatches 0" without having
compared anything. A fake model whose spike decision depends on the settled
stop must be counted as a mismatch on both paths.
"""

import importlib.util
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from msoaxon import _dispatch, coincidence, somatic, threshold as threshold_mod

REPO = Path(__file__).resolve().parents[1]
FACTOR = coincidence.FACTOR
T_END = 10.0


def _load_benchmark():
    """Imports scripts/benchmark_search.py as a module.

    Returns:
        The loaded module.
    """
    spec = importlib.util.spec_from_file_location(
        "benchmark_search", REPO / "scripts" / "benchmark_search.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_run_model(disagree, peak_differs=False):
    """Builds a run_model stand-in that never simulates anything.

    Args:
        disagree: If True, the run spikes with the settled stop on and does
            not spike with it off, so the two settings give different
            yes/no answers. If False, both settings spike.
        peak_differs: If True, the peak axon - soma differs between the two
            settings while the spike decision stays the same.

    Returns:
        A function with run_model's signature returning (t, x) for a model
        whose soma is column 0 and axon is column 1 or node - 1.
    """

    def run_model(
        model,
        stim,
        start,
        stop,
        amp,
        node,
        t_end,
        v0,
        input_node,
        syn=None,
        *,
        stop_on_spike,
        stop_when_settled,
        **kw,
    ):
        spiked = stop_when_settled or not disagree
        t = np.array([0.0, 1.0 if spiked else t_end])
        x = np.zeros((2, max(node, 2)))
        axon = node - 1 if model == "multi" else 1
        x[-1, axon] = stop_on_spike if spiked else 0.0
        if peak_differs and not stop_when_settled:
            x[0, axon] = 2 * stop_on_spike
        return t, x

    return run_model


@pytest.fixture
def instrumented(monkeypatch):
    """Runs benchmark_search.instrument on a fake model, undone afterwards.

    Yields:
        A function (disagree, peak_differs=False) -> stats that installs the
        fake model and the --check wrappers and returns the live stats dict.
    """
    bench = _load_benchmark()
    monkeypatch.setattr(_dispatch, "STOP_WHEN_SETTLED", True)
    monkeypatch.setattr(_dispatch, "run_model", _dispatch.run_model)
    monkeypatch.setattr(_dispatch, "response", _dispatch.response)
    monkeypatch.setattr(coincidence, "response", coincidence.response)
    for mod in (_dispatch, coincidence, somatic, threshold_mod):
        monkeypatch.setattr(
            mod, "spikes", getattr(mod, "spikes", None), raising=False
        )

    def install(disagree, peak_differs=False):
        monkeypatch.setattr(
            _dispatch, "run_model", _fake_run_model(disagree, peak_differs)
        )
        stats = defaultdict(float)
        bench.instrument(stats, True)
        return stats

    return install


def _epsg_pair_decision():
    """One yes/no of the EPSG-pair threshold search, through coincidence."""
    return coincidence._spikes(
        "multi",
        100.0,
        0.1,
        3,
        -68.0,
        coincidence.EPSG_TAU,
        ("EPSGpair", 1, None),
        {},
    )


def _spikes_decision():
    """One yes/no through _dispatch.spikes, as rheobase and BinarySearch."""
    return _dispatch.spikes(
        "multi", "step", 0.0, 1.0, 100.0, 3, T_END, -68.0, 1, factor=FACTOR
    )


def test_check_flags_disagreeing_response(instrumented):
    stats = instrumented(disagree=True)
    assert _epsg_pair_decision() is True  # the settled-stop answer is kept
    assert stats["mismatches"] == 1
    assert stats["runs"] == 1  # the check's own run is not counted


def test_check_flags_disagreeing_spikes(instrumented):
    stats = instrumented(disagree=True)
    assert _spikes_decision() is True
    assert stats["mismatches"] == 1
    assert stats["runs"] == 1


def test_check_passes_agreeing_model(instrumented):
    stats = instrumented(disagree=False)
    assert _epsg_pair_decision() is True
    assert _spikes_decision() is True
    assert stats["mismatches"] == 0
    assert stats["runs"] == 2


def test_check_compares_spike_flag_not_peak(instrumented):
    stats = instrumented(disagree=False, peak_differs=True)
    assert _epsg_pair_decision() is True
    assert stats["mismatches"] == 0
    assert _dispatch.STOP_WHEN_SETTLED is True  # restored after the re-run
