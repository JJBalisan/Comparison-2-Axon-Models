"""Search arguments outside their range raise instead of hanging (issue #69).

Every case here used to loop for ever (or for ~1000 runs). The model entry
points are patched to fail the test if they are reached, so a missing guard
shows up as a failure rather than as a hung test run, and a passing test
shows that the error arrives before any simulation (a Relative criterion's
spike height included). The repo has no pytest-timeout, so none is used.

Relative(f) with f outside (0, 1) is already refused when it is built:
tests/test_criterion.py::test_relative_needs_a_fraction_strictly_between_0_and_1.
"""

import numpy as np
import pytest

from msoaxon import _bisect, _dispatch, coincidence, somatic, threshold

BAD = [0.0, -1.0, float("nan")]


@pytest.fixture(autouse=True)
def no_simulation(monkeypatch):
    """Fails the test if a search reaches the model."""

    def reached(*args, **kw):
        raise AssertionError("a model run happened before the check")

    for name in ("run_model", "response", "spikes"):
        monkeypatch.setattr(_dispatch, name, reached)


def never(*args):
    raise AssertionError("the predicate was called before the check")


# --- _bisect: the two searches, with a fake predicate ---


@pytest.mark.parametrize("bad", BAD)
def test_smallest_firing_rejects_a_guess_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"guess must be > 0"):
        _bisect.smallest_firing(never, bad, 100.0, 1e-3)


@pytest.mark.parametrize("bad", BAD)
def test_smallest_firing_rejects_a_tolerance_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"rel_tol must be > 0"):
        _bisect.smallest_firing(never, 1.0, 100.0, bad)


@pytest.mark.parametrize("bad", BAD)
def test_smallest_crossing_rejects_a_guess_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"guess must be > 0"):
        _bisect.smallest_crossing(never, bad, 100.0, 1e-3, 10.0)


@pytest.mark.parametrize("bad", BAD)
def test_smallest_crossing_rejects_a_tolerance_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"rel_tol must be > 0"):
        _bisect.smallest_crossing(never, 1.0, 100.0, bad, 10.0)


@pytest.mark.parametrize("bad", [1.0, 0.5, 0.0, -2.0, float("nan")])
def test_smallest_crossing_rejects_a_descent_that_is_not_above_1(bad):
    # descent has no public surface (coincidence.threshold uses the
    # default, somatic.rheobase passes 2.0), so only this entry point
    with pytest.raises(ValueError, match=r"descent must be > 1"):
        _bisect.smallest_crossing(never, 1.0, 100.0, 1e-3, 10.0, descent=bad)


def test_the_message_names_the_value():
    with pytest.raises(ValueError, match=r"guess must be > 0, got -5\.0"):
        _bisect.smallest_firing(never, -5.0, 100.0, 1e-3)


def test_valid_arguments_still_search():
    """The smallest legal values are accepted: tiny, but positive."""
    got = _bisect.smallest_firing(lambda x: x >= 3.0, 1.0, 100.0, 1e-3)
    assert 3.0 <= got <= 3.0 * (1 + 1e-3) * 2
    got = _bisect.smallest_crossing(
        lambda x: (x >= 3.0, 10.0 * x / 3.0 if x < 3.0 else 10.0),
        1.0,
        100.0,
        1e-3,
        10.0,
        descent=1.0001,
    )
    assert got >= 3.0


# --- coincidence ---


@pytest.mark.parametrize("bad", BAD)
@pytest.mark.parametrize("name", ["guess", "rel_tol"])
def test_threshold_rejects_guess_and_tolerance_not_positive(name, bad):
    with pytest.raises(ValueError, match=f"{name} must be > 0"):
        coincidence.threshold("two", 0.0, **{name: bad})


@pytest.mark.parametrize("name", ["guess", "rel_tol"])
def test_threshold_checks_before_a_relative_criterion_runs_anything(name):
    # Relative(0.5) would first search for the spike height: a full
    # threshold search, which no_simulation would trip on
    with pytest.raises(ValueError, match=f"{name} must be > 0"):
        coincidence.threshold(
            "two",
            0.0,
            criterion=coincidence.Relative(0.5),
            **{name: 0.0},
        )


@pytest.mark.parametrize("bad", BAD)
def test_threshold_curve_rejects_a_tolerance_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"rel_tol must be > 0"):
        coincidence.threshold_curve(
            "two",
            [0.0, 0.1],
            rel_tol=bad,
            criterion=coincidence.Relative(0.5),
            workers=1,
        )


@pytest.mark.parametrize("bad", BAD)
def test_window_rejects_a_boundary_tolerance_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"delay_tol must be > 0"):
        coincidence.window(
            "two",
            0.03,
            delay_tol=bad,
            threshold0=50.0,
            criterion=coincidence.Relative(0.5),
        )


@pytest.mark.parametrize("bad", BAD)
def test_window_rejects_a_threshold_tolerance_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"rel_tol must be > 0"):
        coincidence.window(
            "two", 0.03, rel_tol=bad, criterion=coincidence.Relative(0.5)
        )


def test_window_ignores_rel_tol_when_threshold0_is_given():
    # rel_tol only tolerances the search for th0, which is skipped: the
    # run that follows is the first model call, and it is the patched one
    with pytest.raises(AssertionError, match="before the check"):
        coincidence.window("two", 0.03, rel_tol=0.0, threshold0=50.0)


# --- somatic ---


@pytest.mark.parametrize("bad", BAD)
@pytest.mark.parametrize("name", ["guess", "rel_tol"])
def test_rheobase_rejects_guess_and_tolerance_not_positive(name, bad):
    with pytest.raises(ValueError, match=f"{name} must be > 0"):
        somatic.rheobase("two", **{name: bad})


def test_rheobase_checks_before_a_relative_criterion_runs_anything():
    with pytest.raises(ValueError, match=r"guess must be > 0"):
        somatic.rheobase("two", guess=0.0, criterion=coincidence.Relative(0.5))


# --- threshold ---


@pytest.mark.parametrize("bad", BAD)
def test_binary_search_rejects_a_zoom_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"zoom must be > 0"):
        threshold.binary_search(
            "EPSGpair", 2, 3, 10.0, 150.0, zoom=bad, workers=1
        )


@pytest.mark.parametrize("bad", BAD)
def test_the_halving_search_rejects_a_zoom_that_is_not_positive(bad):
    with pytest.raises(ValueError, match=r"zoom must be > 0"):
        threshold._search(never, "step", 1, 10.0, 150.0, bad, 1 / 25)


# --- what did not need a guard ---


def test_a_ceiling_below_the_guess_ends_the_search_at_once():
    """Checked while writing this: `ceiling` cannot make a search hang."""
    assert coincidence.threshold("two", 0.0, guess=50.0, ceiling=1.0) == np.inf
