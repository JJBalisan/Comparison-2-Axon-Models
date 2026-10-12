"""The package version matches pyproject.toml and the changelog."""

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
