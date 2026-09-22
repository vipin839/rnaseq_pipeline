"""Installation tests: build the wheel, install it in clean environments, run the global command.

These tests do not trust the developer machine: the installed command is run with an empty HOME, a minimal
PATH and no conda, from several working directories (including one with spaces). Slow (a few minutes) and
needs network access to PyPI for PyYAML/setuptools; enable with RNASEQ_PACKAGING_TESTS=1 (CI does).
"""
import os
import shutil
import subprocess
import sys
import venv
from pathlib import Path

import pytest

from conftest import ROOT

pytestmark = pytest.mark.skipif(os.environ.get("RNASEQ_PACKAGING_TESTS") != "1",
                                reason="set RNASEQ_PACKAGING_TESTS=1 to run installation tests")

VERSION = next(line.split('"')[1] for line in (ROOT / "src" / "rnaseq_pipeline" / "__init__.py").read_text()
               .splitlines() if line.startswith("__version__"))


@pytest.fixture(scope="module")
def wheel(tmp_path_factory):
    out = tmp_path_factory.mktemp("dist")
    subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "-q", "-w", str(out), str(ROOT)], check=True)
    whl = next(out.glob("rnaseq_pipeline-*.whl"))
    assert VERSION in whl.name
    return whl


def clean_env(home, bin_dir):
    """No developer PATH, no conda, empty HOME."""
    return {"HOME": str(home), "PATH": f"{bin_dir}:/usr/bin:/bin", "LANG": "C.UTF-8"}


def run(cmd, env, cwd, stdin=""):
    return subprocess.run([str(c) for c in cmd], env=env, cwd=cwd, input=stdin, capture_output=True, text=True,
                          timeout=600)


@pytest.fixture(scope="module")
def venv_install(wheel, tmp_path_factory):
    d = tmp_path_factory.mktemp("venv")
    venv.EnvBuilder(with_pip=True, clear=True).create(d / "env")
    py = d / "env" / "bin" / "python"
    subprocess.run([py, "-m", "pip", "install", "-q", str(wheel)], check=True)
    home = d / "home"
    home.mkdir()
    return d / "env" / "bin", home


@pytest.mark.parametrize("where", ["home", "tmp", "dir with spaces"])
def test_installed_command_from_any_directory(venv_install, tmp_path, where):
    bin_dir, home = venv_install
    cwd = {"home": home, "tmp": Path("/tmp"), "dir with spaces": tmp_path / "my data dir"}[where]
    cwd.mkdir(exist_ok=True)
    env = clean_env(home, bin_dir)
    for exe in ("rnaseq-pipeline", "rnaseq_pipeline"):
        r = run([bin_dir / exe, "--version"], env, cwd)
        assert r.returncode == 0 and f"Pipeline Version: {VERSION}" in r.stdout, r.stderr
    r = run([bin_dir / "rnaseq-pipeline", "--help"], env, cwd)
    assert r.returncode == 0 and "--validate-project" in r.stdout
    r = run([bin_dir / "python", "-m", "rnaseq_pipeline", "--version"], env, cwd)
    assert r.returncode == 0 and VERSION in r.stdout


def test_installed_package_data_and_runtime_env_files(venv_install):
    bin_dir, home = venv_install
    env = clean_env(home, bin_dir)
    for kind in ("tools", "r"):
        r = run([bin_dir / "rnaseq-pipeline", "--runtime-env", kind], env, home)
        path = Path(r.stdout.strip())
        assert r.returncode == 0 and path.is_file() and "dependencies:" in path.read_text()
        assert str(ROOT) not in str(path), "must come from the installed package, not the checkout"


def test_check_without_tools_fails_cleanly(venv_install):
    """No tools, no R, no conda: the doctor must FAIL with remedies — never a traceback."""
    bin_dir, home = venv_install
    r = run([bin_dir / "rnaseq-pipeline", "--check", "--quick"], clean_env(home, bin_dir), home)
    assert r.returncode == 1
    assert "Traceback" not in r.stdout + r.stderr
    assert "HEALTH: FAIL" in r.stdout and "hisat2" in r.stdout and "main menu 4" in r.stdout


def test_menu_starts_and_exits_without_tools(venv_install):
    bin_dir, home = venv_install
    r = run([bin_dir / "rnaseq-pipeline"], clean_env(home, bin_dir), home, stdin="8\n")
    assert r.returncode == 0 and "MAIN MENU" in r.stdout and "Goodbye" in r.stdout
    assert "Traceback" not in r.stdout + r.stderr


def test_projects_and_state_live_in_home_not_package(venv_install, tmp_path):
    """Creating a project from the installed command must not write into site-packages."""
    bin_dir, home = venv_install
    site = next((bin_dir.parent / "lib").glob("python*/site-packages/rnaseq_pipeline"))
    before = {p: p.stat().st_mtime_ns for p in site.rglob("*") if p.is_file()}
    answers = "1\ny\nInstallTest\n\n8\n"          # continue despite missing tools, name, default dir, then quit
    run([bin_dir / "rnaseq-pipeline"], clean_env(home, bin_dir), home, stdin=answers)
    assert (home / "rnaseq_projects" / "RNAseq_InstallTest" / "project.json").exists()
    after = {p: p.stat().st_mtime_ns for p in site.rglob("*") if p.is_file() and "__pycache__" not in str(p)}
    changed = [str(p) for p in after if p in before and after[p] != before[p]]
    assert not changed, f"package files modified: {changed}"


def test_pipx_install(wheel, tmp_path):
    pipx = shutil.which("pipx") or os.environ.get("RNASEQ_PIPX")
    if not pipx:
        pytest.skip("pipx not available")
    env = {"HOME": str(tmp_path / "home"), "PIPX_HOME": str(tmp_path / "pipx"),
           "PIPX_BIN_DIR": str(tmp_path / "bin"), "PATH": f"{tmp_path / 'bin'}:/usr/bin:/bin", "LANG": "C.UTF-8"}
    (tmp_path / "home").mkdir()
    r = run([pipx, "install", str(wheel), "--python", sys.executable], env, tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    exe = tmp_path / "bin" / "rnaseq-pipeline"
    assert exe.exists()
    r = run(["rnaseq-pipeline", "--version"], env, tmp_path / "home")
    assert r.returncode == 0 and VERSION in r.stdout
