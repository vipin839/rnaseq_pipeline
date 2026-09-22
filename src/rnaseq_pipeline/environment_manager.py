"""Conda/Mamba environment discovery, creation and reproducibility exports."""
import os
import platform
import shutil
from pathlib import Path

from . import PACKAGE_DIR, runner, ui
from . import config as C

ENV_FILES = {"tools": PACKAGE_DIR / "envs" / "rnaseq-tools.yml",
             "r": PACKAGE_DIR / "envs" / "rnaseq-r.yml"}


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
            w = shutil.which("Rscript")
            return Path(w) if w else None
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
        return [str(exe), "env", "create", "-y", "-n", self.names[key], "-f", str(ENV_FILES[key])]

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
