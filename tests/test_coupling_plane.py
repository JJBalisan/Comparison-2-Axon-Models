"""scripts/coupling_plane.py: the grid, the locating and the point statuses.

None of these runs a simulation: the grid and the locating are arithmetic, and
the statuses are checked with the searches replaced (the real ones run in
test_scripts.py's --quick test and in test_two_knobs.py).
"""

import importlib.util
import json
import math
import pathlib

import pytest

from msoaxon import _bisect, coincidence, constants

REPO = pathlib.Path(__file__).resolve().parents[1]


def _script():
    spec = importlib.util.spec_from_file_location(
        "coupling_plane", REPO / "scripts" / "coupling_plane.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


script = _script()


def _pt(c12, c21, window_us, threshold, status="fires"):
    return {
        "couple12": c12,
        "couple21": c21,
        "window_us": window_us,
        "threshold": threshold,
        "spike_height": 38.0 if status == "fires" else math.nan,
        "status": status,
    }


# --- the grid ----------------------------------------------------------------


def test_axis_is_regular_and_free_of_float_drift():
    assert script.axis(0.1) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    assert script.axis(0.4) == [0.1, 0.5, 0.9]
    assert script.axis(0.25, 0.25, 0.75) == [0.25, 0.5, 0.75]
    assert script.axis(1.0, 0.3, 0.3) == [0.3]


@pytest.mark.parametrize(
    "step, lo, hi",
    [
        (0.0, 0.1, 0.9),
        (-0.1, 0.1, 0.9),
        (math.nan, 0.1, 0.9),
        (0.1, 0.0, 0.9),
        (0.1, 0.1, 1.0),
        (0.1, 0.5, 0.4),
    ],
)
def test_axis_rejects_what_is_no_axis(step, lo, hi):
    with pytest.raises(ValueError):
        script.axis(step, lo, hi)


def test_grid_is_goldwyns_region():
    points = script.grid(0.1)
    # the lower right half of a 9 x 9 grid, diagonal included
    assert len(points) == 45
    assert len(set(points)) == 45
    for c12, c21 in points:
        assert c12 >= c21  # "the upper left half ... is empty"
        assert c12 / c21 <= 10
    assert (0.1, 0.1) in points and (0.9, 0.9) in points
    assert (0.1, 0.2) not in points
    # ascending in the forward constant, then the backward
    assert points == sorted(points)


def test_grid_drops_ratios_above_the_limit():
    # 0.5 / 0.05 is exactly 10 and stays; 0.55 / 0.05 is 11 and goes
    points = script.grid(0.05, 0.05, 0.55)
    assert (0.5, 0.05) in points
    assert (0.55, 0.05) not in points
    assert (0.55, 0.1) in points
    assert (0.55, 0.05) in script.grid(0.05, 0.05, 0.55, max_ratio=11)


def test_grid_rejects_a_ratio_below_one():
    with pytest.raises(ValueError, match="max_ratio"):
        script.grid(0.1, max_ratio=0.5)


def test_quick_grid_has_six_points():
    assert script.grid(script.QUICK_STEP) == [
        (0.1, 0.1),
        (0.5, 0.1),
        (0.5, 0.5),
        (0.9, 0.1),
        (0.9, 0.5),
        (0.9, 0.9),
    ]


def test_landmarks_are_the_lookup_and_goldwyns_stars():
    marks = script.landmarks()
    default = marks[f"node {script.NODE} (TwoCpt.m)"]
    assert default == (
        pytest.approx(constants.COUPLING1[1]),
        pytest.approx(constants.COUPLING2[1]),
    )
    assert [v for k, v in marks.items() if "Goldwyn" in k] == [
        (0.3, 0.2),
        (0.8, 0.2),
        (0.8, 0.7),
    ]
    # the stars and the lookup's own point are all inside the plane swept
    for c12, c21 in marks.values():
        assert c12 >= c21 and c12 / c21 <= script.MAX_RATIO


# --- locating the 45-compartment model ---------------------------------------


def test_distance_is_relative_in_both_quantities():
    assert script.distance(200.0, 48.0, 200.0, 48.0) == 0.0
    # 10% off in the window alone, and 10% in the threshold alone
    assert script.distance(220.0, 48.0, 200.0, 48.0) == pytest.approx(0.1)
    assert script.distance(200.0, 43.2, 200.0, 48.0) == pytest.approx(0.1)
    # both at once add in quadrature
    assert script.distance(220.0, 52.8, 200.0, 48.0) == pytest.approx(
        0.1 * math.sqrt(2)
    )


@pytest.mark.parametrize(
    "w, th",
    [(math.nan, 48.0), (200.0, math.inf), (math.inf, math.nan)],
)
def test_distance_is_nan_without_a_value(w, th):
    assert math.isnan(script.distance(w, th))


def test_locate_finds_the_closest_point_of_a_synthetic_plane():
    ref_w, ref_th = 200.0, 48.0
    points = [
        _pt(0.1, 0.1, 80.0, 40.0),  # far in the window
        _pt(0.5, 0.2, 205.0, 55.0),  # 2.5% and 14.6% off
        _pt(0.7, 0.3, 190.0, 49.0),  # 5% and 2% off: the closest
        _pt(0.9, 0.5, 200.0, math.inf, "silent"),  # a perfect window, no th
        _pt(0.8, 0.8, math.nan, math.nan, "spontaneous"),
        _pt(0.9, 0.1, None, None, "wide"),
        # the 10 mV criterion cleared by a graded response: no window, but
        # a threshold and a window that would be closest if they counted
        _pt(0.2, 0.1, 200.0, 48.0, "graded"),
    ]
    best, d = script.locate(points, ref_w, ref_th)
    assert (best["couple12"], best["couple21"]) == (0.7, 0.3)
    assert d == pytest.approx(math.hypot(0.05, 1 / 48))


def test_locate_takes_the_first_of_a_tie():
    points = [_pt(0.4, 0.2, 190.0, 48.0), _pt(0.6, 0.2, 210.0, 48.0)]
    assert script.locate(points, 200.0, 48.0)[0]["couple12"] == 0.4


def test_locate_with_an_exact_match_has_distance_zero():
    best, d = script.locate([_pt(0.9, 0.1, 200.0, 48.0)], 200.0, 48.0)
    assert d == 0.0 and best["couple12"] == 0.9


def test_locate_refuses_a_plane_with_nothing_to_compare():
    with pytest.raises(ValueError, match="finite"):
        script.locate(
            [
                _pt(0.5, 0.1, math.nan, math.inf, "silent"),
                _pt(0.9, 0.9, None, None),
            ]
        )
    with pytest.raises(ValueError):
        script.locate([])


def test_locate_defaults_to_the_published_45_compartment_row():
    best, d = script.locate([_pt(0.9, 0.1, 201.1, 48.01)])
    assert d == 0.0
    assert script.MULTI_WINDOW_US == 201.1
    assert script.MULTI_THRESHOLD_NS == 48.01


def test_extremes_ignore_points_without_a_value():
    points = [
        _pt(0.1, 0.1, 170.0, 78.0),
        _pt(0.8, 0.2, 141.0, 80.0),
        _pt(0.9, 0.1, 201.0, 34.0),
        _pt(0.5, 0.5, math.nan, math.inf, "silent"),
        _pt(0.9, 0.9, None, None, "spontaneous"),
        _pt(0.3, 0.1, 100.0, 20.0, "graded"),  # does not spike: not counted
    ]
    ex = script.extremes(points)
    assert ex["window_min"]["couple12"] == 0.8
    assert ex["window_max"]["couple12"] == 0.9
    assert ex["threshold_min"]["couple12"] == 0.9
    assert ex["threshold_max"]["couple12"] == 0.8
    assert script.extremes([_pt(0.5, 0.5, math.nan, math.inf, "silent")]) == {}


# --- the statuses of measure_point -------------------------------------------


def _fake(monkeypatch, threshold, window=None, height=38.0):
    """Replaces the searches and the spike height; returns what they saw."""
    seen = []

    def fake_threshold(model, delay, model_kw=None, **kw):
        seen.append((model, delay, model_kw, kw))
        if isinstance(threshold, Exception):
            raise threshold
        return threshold

    def fake_window(model, margin, model_kw=None, **kw):
        seen.append((model, margin, model_kw, kw))
        return window, kw["threshold0"]

    def fake_height(model_kw, th):
        seen.append(("two", "height", model_kw, {"threshold": th}))
        return height

    monkeypatch.setattr(script, "_spike_height", fake_height)
    monkeypatch.setattr(coincidence, "threshold", fake_threshold)
    monkeypatch.setattr(coincidence, "window", fake_window)
    return seen


TASK = (0.8, 0.2, "exact", 500.0, 1e-3, 1e-3)


def test_a_point_that_fires_has_a_threshold_and_a_window_in_us(monkeypatch):
    seen = _fake(monkeypatch, 80.0, window=0.141)
    p = script.measure_point(TASK)
    assert p["status"] == "fires"
    assert p["threshold"] == 80.0
    assert p["window_us"] == pytest.approx(141.0)
    (m, d, kw, search), (_, _, hkw, hsearch), (_, margin, wkw, wsearch) = seen
    # the coupling goes in through model_kw, the ceiling and tolerances to
    # the searches, and the window reuses the threshold found
    assert kw == hkw == wkw
    assert kw == {"couple12": 0.8, "couple21": 0.2, "mode": "exact"}
    assert p["spike_height"] == 38.0 and hsearch == {"threshold": 80.0}
    assert (m, d, margin) == ("two", 0.0, script.MARGIN)
    assert search == {"rel_tol": 1e-3, "ceiling": 500.0}
    assert wsearch["threshold0"] == 80.0 and wsearch["delay_tol"] == 1e-3


def test_a_cell_that_never_spikes_is_silent_and_has_no_window(monkeypatch):
    seen = _fake(monkeypatch, math.inf)
    p = script.measure_point(TASK)
    assert p["status"] == "silent" and p["threshold"] == math.inf
    assert math.isnan(p["window_us"])
    assert len(seen) == 1  # no spike height or window without a threshold
    assert math.isnan(p["spike_height"])


def test_a_cell_that_fires_unprompted_is_spontaneous(monkeypatch):
    _fake(monkeypatch, _bisect.FiresWithoutInputError("fires at rest"))
    p = script.measure_point(TASK)
    assert p["status"] == "spontaneous"
    assert math.isnan(p["threshold"]) and math.isnan(p["window_us"])


def test_a_graded_response_has_a_threshold_but_no_window(monkeypatch):
    # the 10 mV criterion is cleared, but the cell's spike is only 18 mV tall
    seen = _fake(monkeypatch, 80.3, window=0.141, height=18.0)
    p = script.measure_point(TASK)
    assert p["status"] == "graded"
    assert p["threshold"] == 80.3 and p["spike_height"] == 18.0
    assert math.isnan(p["window_us"])
    assert len(seen) == 2  # the threshold and the height, no window search
    assert not script.spikes(p)


def test_the_spike_limit_is_above_the_criterion():
    assert script.MIN_SPIKE_MV > coincidence.FACTOR


def test_a_window_beyond_the_search_is_wide(monkeypatch):
    _fake(monkeypatch, 33.9, window=math.nan)
    p = script.measure_point(TASK)
    assert p["status"] == "wide" and p["threshold"] == 33.9
    assert math.isnan(p["window_us"])


def test_every_status_is_listed():
    for status in ("fires", "graded", "silent", "spontaneous", "wide"):
        assert status in script.STATUSES


# --- the table and the JSON --------------------------------------------------


def test_markdown_has_one_row_per_point_and_dashes_for_missing_values():
    points = [
        _pt(0.9, 0.1, 201.3, 33.91) | {"label": "node 3"},
        _pt(0.3, 0.2, math.nan, math.inf, "silent") | {"label": ""},
    ]
    rows = script.markdown(points).splitlines()
    assert len(rows) == 4  # header, rule, two points
    assert rows[2].startswith("| 0.3 | 0.2 |") and "| - | - | - |" in rows[2]
    assert "| node 3 | fires | 33.91 | 38.0 | 201.3 |" in rows[3]


def test_clean_makes_the_json_strict():
    assert script._clean(
        {"a": [math.nan, math.inf, 1.5], "b": (2, -math.inf)}
    ) == {
        "a": [None, None, 1.5],
        "b": [2, None],
    }


def test_a_sweep_where_nothing_spikes_still_writes_its_results(
    monkeypatch, tmp_path, capsys
):
    def silent(points, *_):
        return [
            {
                "couple12": c12,
                "couple21": c21,
                "status": "silent",
                "threshold": math.inf,
                "window_us": math.nan,
                "spike_height": math.nan,
            }
            for c12, c21 in points
        ]

    plotted = []
    monkeypatch.setattr(script, "sweep", silent)
    monkeypatch.setattr(script, "plot", lambda *a: plotted.append(a))
    monkeypatch.setattr(
        "sys.argv", ["coupling_plane.py", "--quick", "--out-dir", str(tmp_path)]
    )
    script.main()
    assert "none, no point spikes" in capsys.readouterr().out
    res = json.loads((tmp_path / "coupling_plane.json").read_text())
    assert res["located"] is None and res["extremes"] == {}
    assert res["default"]["status"] == "silent"
    assert plotted and plotted[0][2] is None  # the plot leaves it unmarked
