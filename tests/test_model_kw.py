"""Model keywords through somatic.* and threshold.binary_search (issue #40).

somatic and the MATLAB-port threshold search used to take only mem=, so
mode=, r1=/tau_est= and the rest reached the models only through the
environment or a direct _dispatch call. Most tests here stand a fake in for
_dispatch and check which keywords arrive.
"""

import numpy as np
import pytest

from msoaxon import (
    _bisect,
    _dispatch,
    coincidence,
    multi,
    somatic,
    threshold,
    two,
)

GOLDWYN_KW = {k: v for k, v in two.GOLDWYN_2019.items() if k != "v0"}
GOLDWYN_V0 = two.GOLDWYN_2019["v0"]
# the Goldwyn 2019 rheobase at the default 10 mV criterion (PYTHON.md, #29)
GOLDWYN_RHEOBASE = 2070.0
KW = {"mode": "exact", "max_step": 0.05, "model_type": "active-full"}


def test_goldwyn_rheobase_through_the_public_api():
    rb = somatic.rheobase(
        "two", v0=GOLDWYN_V0, model_kw=GOLDWYN_KW, rel_tol=1e-2
    )
    assert rb == pytest.approx(GOLDWYN_RHEOBASE, rel=1e-2)


def test_goldwyn_rheobase_matches_a_direct_dispatch_search():
    def response(amp):
        return _dispatch.response(
            "two",
            "step",
            somatic.START,
            somatic.STOP,
            amp,
            3,
            somatic.T_END,
            GOLDWYN_V0,
            1,
            factor=somatic.FACTOR,
            **GOLDWYN_KW,
        )

    direct = _bisect.smallest_crossing(
        response, 2000.0, 1e6, 1e-2, somatic.FACTOR, descent=2.0
    )
    public = somatic.rheobase(
        "two", v0=GOLDWYN_V0, model_kw=GOLDWYN_KW, rel_tol=1e-2
    )
    assert public == pytest.approx(direct, rel=1e-9)


@pytest.fixture
def seen(monkeypatch):
    """Fakes the runs; returns the model keywords each one was given."""
    got = []

    def spikes(model, stim, start, stop, amp, *args, factor, **kw):
        got.append(kw)
        return False

    def response(model, stim, start, stop, amp, *args, factor, **kw):
        got.append(kw)
        return False, 0.0

    def run_model(model, *args, **kw):
        got.append(kw)
        return np.array([0.0, somatic.T_END]), np.zeros((2, 2))

    monkeypatch.setattr(_dispatch, "spikes", spikes)
    monkeypatch.setattr(_dispatch, "response", response)
    monkeypatch.setattr(_dispatch, "run_model", run_model)
    return got


@pytest.mark.parametrize(
    "call",
    [
        lambda **kw: somatic.fires("two", 100.0, **kw),
        lambda **kw: somatic.response("two", 100.0, **kw),
        lambda **kw: somatic.spike_amplitude("two", 100.0, **kw),
    ],
    ids=["fires", "response", "spike_amplitude"],
)
def test_model_kw_reaches_the_model(seen, call):
    call(model_kw=KW)
    assert seen and all(KW.items() <= kw.items() for kw in seen)


def test_rheobase_model_kw_reaches_every_run(seen):
    # nothing fires, so the search climbs to its ceiling and gives up
    assert somatic.rheobase("two", model_kw=KW, rel_tol=1e-2) == np.inf
    assert len(seen) > 1
    assert all(KW.items() <= kw.items() for kw in seen)


def test_no_model_kw_passes_nothing_extra(seen):
    somatic.fires("two", 100.0)
    somatic.spike_amplitude("two", 100.0)
    assert seen == [{}, {"stop_on_spike": somatic.FACTOR}]


def test_spike_amplitude_keeps_its_fine_max_step(monkeypatch):
    runs = []

    def run_model(model, *args, **kw):
        runs.append(kw)
        if "stop_on_spike" in kw:  # the spike is at 6 ms
            return np.array([0.0, 5.0, 6.0]), np.zeros((3, 2))
        t = np.linspace(0.0, 8.0, 801)
        return t, np.zeros((len(t), 2))

    monkeypatch.setattr(_dispatch, "run_model", run_model)
    somatic.spike_amplitude("two", 100.0, model_kw={"max_step": 0.5})
    assert [kw["max_step"] for kw in runs] == [0.5, 0.005]


def test_mem_is_a_shorthand_for_model_kw_mem(seen):
    mem = multi.membrane(-68.0)
    somatic.fires("multi", 100.0, mem=mem)
    somatic.fires("multi", 100.0, model_kw={"mem": mem})
    assert seen[0]["mem"] is seen[1]["mem"] is mem


def test_mem_and_other_keywords_combine(seen):
    mem = multi.membrane(-68.0)
    somatic.fires("multi", 100.0, mem=mem, model_kw={"mode": "exact"})
    assert seen[0]["mem"] is mem and seen[0]["mode"] == "exact"


def test_model_kw_is_not_modified(seen):
    mem = multi.membrane(-68.0)
    kw = {"mode": "exact"}
    somatic.fires("multi", 100.0, mem=mem, model_kw=kw)
    assert kw == {"mode": "exact"}


@pytest.mark.parametrize(
    "call",
    [
        somatic.fires,
        somatic.response,
        somatic.rheobase,
        somatic.spike_amplitude,
    ],
)
def test_mem_given_both_ways_raises(seen, call):
    mem = multi.membrane(-68.0)
    args = ("multi",) if call is somatic.rheobase else ("multi", 100.0)
    with pytest.raises(ValueError, match="mem"):
        call(*args, mem=mem, model_kw={"mem": mem})
    assert not seen


def test_relative_criterion_measures_height_with_the_model_kw(
    monkeypatch, seen
):
    asked = []

    def spike_height(model, node, v0, epsg_tau=None, model_kw=None, **site):
        asked.append(model_kw)
        return 30.0

    monkeypatch.setattr(coincidence, "spike_height", spike_height)
    somatic.fires(
        "two", 1.0, criterion=coincidence.Relative(0.5), model_kw=GOLDWYN_KW
    )
    assert asked == [GOLDWYN_KW]


def test_model_kw_distinguishes_spike_height_cache_entries():
    def key(**model_kw):
        return coincidence._height_key("two", 3, -68.0, (0.1, 0.18), model_kw)

    variants = [
        {},
        {"r1": 8.5},
        {"tau_est": 0.34},
        {"model_type": "passive"},
        {"max_step": 0.05},
    ]
    keys = [key(**v) for v in variants]
    assert len(set(keys)) == len(keys)
    multi_keys = {
        coincidence._height_key("multi", 3, -68.0, (0.1, 0.18), v)
        for v in ({"mode": "exact"}, {"morph": multi.LUMPED})
    }
    assert len(multi_keys) == 2


def test_binary_search_model_kw_reaches_the_model(seen):
    threshold.binary_search(
        "EPSGpair", 1, 3, 10.0, 10.0, zoom=1.0, workers=1, model_kw=KW
    )
    assert len(seen) >= 2  # both models
    assert all(KW.items() <= kw.items() for kw in seen)


def test_binary_search_defaults_are_unchanged(seen):
    threshold.binary_search("EPSGpair", 1, 3, 10.0, 10.0, workers=1)
    assert seen and all(
        kw == {"model_type": threshold.MODEL_TYPE} for kw in seen
    )


def test_binary_search_model_kw_overrides_model_type(seen):
    threshold.binary_search(
        "EPSGpair",
        1,
        3,
        10.0,
        10.0,
        workers=1,
        model_kw={"model_type": "passive"},
    )
    assert {kw["model_type"] for kw in seen} == {"passive"}
