"""The package version matches pyproject.toml, the changelog and the crate."""

import pathlib
import tomllib

import msoaxon

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_version_matches_pyproject():
    with (ROOT / "pyproject.toml").open("rb") as f:
        assert msoaxon.__version__ == tomllib.load(f)["project"]["version"]


def test_changelog_has_an_entry_for_this_version():
    heading = f"## {msoaxon.__version__} ("
    assert heading in (ROOT / "CHANGELOG.md").read_text()


def test_crate_version_matches_the_package():
    # msoaxon-fast carries msoaxon's version, bumped in the same release PR
    # (#57); the crate's pyproject.toml, Cargo.toml and Cargo.lock all say so.
    with (ROOT / "rust" / "pyproject.toml").open("rb") as f:
        assert tomllib.load(f)["project"]["version"] == msoaxon.__version__
    with (ROOT / "rust" / "Cargo.toml").open("rb") as f:
        assert tomllib.load(f)["package"]["version"] == msoaxon.__version__
    with (ROOT / "rust" / "Cargo.lock").open("rb") as f:
        crate = [
            p for p in tomllib.load(f)["package"] if p["name"] == "msoaxon-fast"
        ]
    assert [p["version"] for p in crate] == [msoaxon.__version__]
