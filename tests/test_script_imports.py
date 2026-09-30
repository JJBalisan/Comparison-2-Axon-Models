"""Pool workers re-import the running script as __mp_main__ (spawn and forkserver
both do), so the scripts that open worker pools must stay cheap to import: in
particular they import matplotlib only where they draw, in the main process."""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PROBE = ("import runpy, sys; runpy.run_path(sys.argv[1], run_name='__mp_main__'); "
         "print('matplotlib' in sys.modules)")


@pytest.mark.parametrize("script", ["dendrites.py", "coincidence_window.py",
                                    "making_threshold_graphs.py"])
def test_pool_scripts_import_without_matplotlib(script):
    r = subprocess.run([sys.executable, "-c", PROBE, str(REPO / "scripts" / script)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "False"
