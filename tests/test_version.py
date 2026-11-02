import tomllib
from importlib.metadata import version
from pathlib import Path

import lattice_mc

PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def test_pyproject_reads_version_from_package():
    project = tomllib.loads(PYPROJECT.read_text())
    assert "version" not in project["project"]
    assert "version" in project["project"]["dynamic"]
    assert project["tool"]["setuptools"]["dynamic"]["version"] == {"attr": "lattice_mc.__version__"}


def test_version_matches_installed_package_metadata():
    assert version("lattice_mc") == lattice_mc.__version__
