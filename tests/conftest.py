"""Fixtures shared by the test modules."""

import pytest

from msoaxon import _solve


@pytest.fixture
def fresh_quiet_cache():
    """Starts the test, and leaves it, with an empty quiet-prefix cache."""
    _solve._QUIET_CACHE.clear()
    yield
    _solve._QUIET_CACHE.clear()
