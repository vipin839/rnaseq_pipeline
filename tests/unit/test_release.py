"""Release consistency: one authoritative version, reflected everywhere a user or packager reads it."""
import re

from conftest import ROOT
from rnaseq_pipeline import __version__


def test_version_is_semver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)


def test_changelog_top_entry_matches_version():
    first = next(line for line in (ROOT / "CHANGELOG.md").read_text().splitlines() if line.startswith("## "))
    assert first.split()[1] == __version__, f"CHANGELOG top entry {first!r} != package version {__version__}"


def test_bioconda_recipe_version_matches():
    recipe = (ROOT / "packaging" / "bioconda" / "rnaseq-pipeline" / "meta.yaml").read_text()
    assert f'{{% set version = "{__version__}" %}}' in recipe


def test_pyproject_reads_version_from_package():
    text = (ROOT / "pyproject.toml").read_text()
    assert 'dynamic = ["version"]' in text and 'attr = "rnaseq_pipeline.__version__"' in text
    # a literal `version = "x.y.z"` would be a second source of truth (the dynamic attr line is fine)
    assert not re.search(r'^version\s*=\s*["\']', text, re.M), "a second, hard-coded version would drift"


def test_no_stale_version_file():
    assert not (ROOT / "VERSION").exists()


def test_readme_does_not_hardcode_an_old_version():
    head = (ROOT / "README.md").read_text().splitlines()[0]
    found = re.findall(r"\d+\.\d+\.\d+", head)
    assert not found or found == [__version__], head
