"""The analysis layer: parallel helpers, threshold sweeps, coincidence trials,
somatic measurements and passive measurements."""

import numpy as np
import pytest

from msoaxon import _solve
from msoaxon._parallel import map_tasks
from msoaxon.coincidence import threshold_curve


def test_map_tasks_handles_no_tasks():
    assert map_tasks(abs, []) == []
    assert threshold_curve("two", []).shape == (0,)
