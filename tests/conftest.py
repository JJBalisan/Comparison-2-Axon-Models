"""Fixtures shared by the test modules."""

import importlib.util
import pathlib
import sys

import pytest

from msoaxon import _solve

# Tests that build a forkserver context directly. Windows has no forkserver
# (the pool there spawns), so they cannot run; the other tests cover the
# spawn path.
NO_FORKSERVER_ON_WINDOWS = {
    "tests/test_analysis.py::test_pools_never_fork",
    "tests/test_fast.py::test_pool_workers_inherit_the_mode_env",
}


def pytest_collection_modifyitems(config, items):
    """Skips the tests that need a forkserver context on Windows."""
    if sys.platform != "win32":
        return
    skip = pytest.mark.skip(reason="Windows has no forkserver start method")
    for item in items:
        if item.nodeid in NO_FORKSERVER_ON_WINDOWS:
            item.add_marker(skip)


@pytest.fixture
def fresh_quiet_cache():
    """Starts the test, and leaves it, with an empty quiet-prefix cache."""
    _solve._QUIET_CACHE.clear()
    yield
    _solve._QUIET_CACHE.clear()


@pytest.fixture
def load_script():
    """Returns a function that imports scripts/<name>.py as a module."""

    def load(name):
        path = pathlib.Path(__file__).parents[1] / "scripts" / f"{name}.py"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    return load
