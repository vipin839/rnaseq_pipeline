"""Conda/Mamba environment discovery, creation and reproducibility exports."""
import os
import platform
import shutil
from pathlib import Path

import yaml

from . import PACKAGE_DIR, USER_STATE_DIR, runner, system_check, ui

ENV_FILES = {"tools": PACKAGE_DIR / "envs" / "rnaseq-tools.yml",
             "r": PACKAGE_DIR / "envs" / "rnaseq-r.yml"}


def _package_name(spec):
    for i, ch in enumerate(spec):
        if ch in "=<>! ":
            return spec[:i]
    return spec


def env_file(key, missing_cpu_features=None):
    """(path, changes): the environment file to create `key` from, adapted to this processor.

    Some tools' current Bioconda builds need x86-64-v3 instructions (AVX2/BMI2) and are killed with SIGILL on older
    CPUs, and conda cannot tell (the packages do not declare it). On such a CPU the tools are pinned to the newest
    version whose build runs there (dependency_manager.OLD_CPU_VERSIONS). The adapted copy is written to the user
    state directory; the installed package directory is never written to. `changes` lists "old -> new" specs.
    """
    from .dependency_manager import OLD_CPU_VERSIONS
    shipped = ENV_FILES[key]
    missing = system_check.missing_x86_64_v3() if missing_cpu_features is None else missing_cpu_features
    if not missing:
        return shipped, []
    spec = yaml.safe_load(shipped.read_text())
    changes, deps = [], []
    for d in spec.get("dependencies", []):
        name = _package_name(d) if isinstance(d, str) else None
        if name in OLD_CPU_VERSIONS:
            new = f"{name}={OLD_CPU_VERSIONS[name]}"
            if d != new:
                changes.append(f"{d} -> {new}")
            d = new
        deps.append(d)
    if not changes:
        return shipped, []
    spec["dependencies"] = deps
    out = USER_STATE_DIR / "envs" / f"{shipped.stem}.cpu-compatible.yml"
    out.parent.mkdir(parents=True, exist_ok=True)
    header = (f"# Generated from {shipped.name} for a CPU without {', '.join(missing)}.\n"
              f"# Changed: {'; '.join(changes)} (the newer builds need AVX2/BMI2 and would crash with SIGILL).\n")
    tmp = out.with_suffix(".yml.part")
    tmp.write_text(header + yaml.safe_dump(spec, sort_keys=False))
    os.replace(tmp, out)
    return out, changes


def rscript_on_path():
    """The first Rscript on PATH whose R library has DESeq2, else the first Rscript on PATH.

    Several R installations can be on PATH: the tools environment gets a bare R as a dependency of RSeQC, and it
    must not be mistaken for the R environment just because it comes first.
    """
    found = []
    for d in os.environ.get("PATH", "").split(os.pathsep):
        f = Path(d or ".") / "Rscript"
        if f.is_file() and os.access(f, os.X_OK) and f not in found:
            found.append(f)
    for f in found:
        if (f.parent.parent / "lib" / "R" / "library" / "DESeq2").is_dir():
            return f
    return found[0] if found else None


class Environments:
    def __init__(self, cfg):
        env_cfg = cfg.get("environment", {})
        self.root = Path(os.path.expanduser(env_cfg.get("conda_root", "~/miniforge3")))
        self.names = {"tools": env_cfg.get("tools_env", "rnaseq-tools"), "r": env_cfg.get("r_env", "rnaseq-r")}
        self.path_fallback = env_cfg.get("use_path_fallback", True)

    # ---------- discovery ----------
    def conda_bin(self):
        for name in ("mamba", "conda"):
            p = self.root / "bin" / name
            if p.exists():
                return p
        for name in ("mamba", "micromamba", "conda"):
            w = shutil.which(name)
            if w:
                return Path(w)
        return None

    def prefix(self, key):
        name = self.names[key]
        cands = [self.root / "envs" / name]
        # micromamba (and setup-micromamba in CI) keeps environments under $MAMBA_ROOT_PREFIX, not beside its binary
        for var in ("MAMBA_ROOT_PREFIX", "CONDA_ROOT"):
            if os.environ.get(var):
                cands.append(Path(os.environ[var]) / "envs" / name)
        exe = self.conda_bin()
        if exe:
            cands.append(exe.resolve().parent.parent / "envs" / name)
        for c in cands:
            if (c / "conda-meta").is_dir():
                return c
        return None

    def bin_dir(self, key):
        p = self.prefix(key)
        return p / "bin" if p else None

    def rscript(self):
        b = self.bin_dir("r")
        if b and (b / "Rscript").exists():
            return b / "Rscript"
        if self.path_fallback:
            return rscript_on_path()
        return None

    def activate(self):
        """Make tools-env binaries visible to all child processes."""
        dirs = [self.bin_dir("tools")]
        runner.configure_tool_path([d for d in dirs if d])
        if not self.path_fallback:
            os.environ["PATH"] = os.pathsep.join(str(d) for d in dirs if d)

    # ---------- creation / installation ----------
    def create_commands(self, key):
        exe = self.conda_bin()
        return [str(exe), "env", "create", "-y", "-n", self.names[key], "-f", str(env_file(key)[0])]

    def install_commands(self, key, packages):
        exe = self.conda_bin()
        return [str(exe), "install", "-y", "-n", self.names[key], "-c", "conda-forge", "-c", "bioconda", *packages]

    def run_install(self, cmd, log_file):
        ui.info("This can take 5–30 minutes depending on network speed.")
        runner.run(cmd, stage="environment", log_file=log_file, description="Installing (see log for details)")

    # ---------- reproducibility exports ----------
    def export(self, out_dir, sysinfo=None):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        exe = self.conda_bin()
        written = []
        env_yml, versions = [], []
        for key in ("tools", "r"):
            if not exe or not self.prefix(key):
                continue
            y = runner.tool_output([exe, "env", "export", "-p", self.prefix(key), "--no-builds"], timeout=300)
            if y:
                env_yml.append(f"# ===== {self.names[key]} ({self.prefix(key)}) =====\n{y}")
            lst = runner.tool_output([exe, "list", "-p", self.prefix(key)], timeout=300)
            if lst:
                versions.append(f"# ===== {self.names[key]} =====\n{lst}")
        if env_yml:
            (out_dir / "environment.yml").write_text("\n".join(env_yml))
            written.append(out_dir / "environment.yml")
        if versions:
            (out_dir / "package_versions.txt").write_text("\n".join(versions))
            written.append(out_dir / "package_versions.txt")
        info = sysinfo or {}
        lines = [f"{k}: {v}" for k, v in info.items()]
        lines.append(f"platform: {platform.platform()}")
        (out_dir / "system_information.txt").write_text("\n".join(lines) + "\n")
        written.append(out_dir / "system_information.txt")
        rs = self.rscript()
        if rs:
            si = runner.tool_output([rs, "-e", "suppressMessages(library(DESeq2)); sessionInfo()"], timeout=300)
            if si:
                (out_dir / "R_sessionInfo.txt").write_text(si)
                written.append(out_dir / "R_sessionInfo.txt")
        return written


def from_config(cfg):
    return Environments(cfg)
