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


def test_bioconda_recipe_requires_what_the_environment_files_install():
    """The one-command Bioconda install must not allow tool versions older than the runtime environment files the
    pipeline is tested with (found in P3: the recipe allowed e.g. fastqc 0.11.9, samtools 1.15, rseqc 4.0)."""
    import yaml
    envs = ROOT / "src" / "rnaseq_pipeline" / "envs"
    want = {}
    for f in ("rnaseq-tools.yml", "rnaseq-r.yml"):
        for d in yaml.safe_load((envs / f).read_text())["dependencies"]:
            name = re.split(r"[=<>]", d)[0]
            if name not in ("python", "pytest", "pyyaml"):
                want[name] = d[len(name):]
    run = (ROOT / "packaging" / "bioconda" / "rnaseq-pipeline" / "meta.yaml").read_text().split("run:")[1].split("test:")[0]
    have = {}
    for line in run.splitlines():
        spec = line.strip().lstrip("- ").split("#")[0].strip()
        if spec:
            name, _, ver = spec.partition(" ")
            have[name] = ver.replace(" ", "")
    mismatch = {n: (v, have.get(n)) for n, v in want.items() if have.get(n) != v}
    assert not mismatch, f"environment file spec vs recipe: {mismatch}"
