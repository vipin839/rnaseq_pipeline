"""Detection of OS, CPU, RAM, disk, filesystem, interpreters and temp space."""
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from . import ui

GB = 1024 ** 3

# Filesystems with restrictive per-file limits or unsafe semantics for pipelines.
FS_LIMITS = {
    "vfat": ("4 GB max file size, no symlinks/permissions", 4 * GB),
    "msdos": ("4 GB max file size", 4 * GB),
    "fat": ("4 GB max file size", 4 * GB),
    "exfat": ("no symlinks/permissions (large files OK)", None),
    "9p": ("WSL Windows-drive mount: slow, unreliable locking/permissions", None),
    "drvfs": ("WSL Windows-drive mount: slow, unreliable locking/permissions", None),
    "fuseblk": ("FUSE (often NTFS): slow, limited permission support", None),
    "tmpfs": ("RAM-backed: contents lost on reboot, consumes memory", None),
}

MIN_RAM_GB = 8
MIN_CORES = 2


def mount_for(path):
    """Return (mountpoint, fstype) for the filesystem containing path."""
    path = os.path.realpath(path)
    best = ("/", "unknown")
    try:
        with open("/proc/mounts") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                mnt = parts[1].replace("\\040", " ")
                if (path == mnt or path.startswith(mnt.rstrip("/") + "/")) and len(mnt) >= len(best[0]):
                    best = (mnt, parts[2])
    except OSError:
        pass
    return best


def mem_gb():
    try:
        with open("/proc/meminfo") as f:
            info = {l.split(":")[0]: int(l.split()[1]) for l in f if ":" in l}
        return info["MemTotal"] / 1024 ** 2, info.get("MemAvailable", info["MemFree"]) / 1024 ** 2
    except (OSError, KeyError, ValueError):
        return 0.0, 0.0


def distro():
    try:
        with open("/etc/os-release") as f:
            d = dict(l.rstrip().split("=", 1) for l in f if "=" in l)
        return d.get("PRETTY_NAME", "").strip('"')
    except OSError:
        return platform.platform()


def writable(path):
    try:
        Path(path).mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path, prefix=".wtest"):
            pass
        return True
    except OSError:
        return False


def _version(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return (r.stdout + r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr) else None
    except (OSError, subprocess.TimeoutExpired):
        return None


def collect(project_dir=None, temp_dir=None, r_bin=None, conda_bin=None):
    project_dir = Path(project_dir or os.getcwd())
    probe = project_dir if project_dir.exists() else project_dir.parent
    while not probe.exists():
        probe = probe.parent
    temp_dir = Path(temp_dir or tempfile.gettempdir())
    total_ram, avail_ram = mem_gb()
    du = shutil.disk_usage(probe)
    mnt, fstype = mount_for(probe)
    tmnt, tfstype = mount_for(temp_dir if temp_dir.exists() else tempfile.gettempdir())
    tdu = shutil.disk_usage(temp_dir if temp_dir.exists() else tempfile.gettempdir())
    wsl = "microsoft" in platform.release().lower()
    return {
        "os": platform.system(), "distribution": distro(), "kernel": platform.release(),
        "wsl": wsl, "architecture": platform.machine(),
        "cpu_cores": os.cpu_count() or 1,
        "ram_total_gb": round(total_ram, 1), "ram_available_gb": round(avail_ram, 1),
        "project_path": str(probe), "disk_total_gb": round(du.total / GB, 1),
        "disk_free_gb": round(du.free / GB, 1),
        "mountpoint": mnt, "filesystem": fstype,
        "max_file_size": (FS_LIMITS.get(fstype, (None, None))[1] or "no restrictive limit known"),
        "fs_note": FS_LIMITS.get(fstype, (None, None))[0],
        "project_writable": writable(probe),
        "temp_dir": str(temp_dir), "temp_free_gb": round(tdu.free / GB, 1), "temp_filesystem": tfstype,
        "temp_writable": writable(temp_dir),
        "python": f"{platform.python_version()} ({sys.executable})",
        "python_ok": sys.version_info >= (3, 10),
        "r": _version([str(r_bin), "--version"]) if r_bin else None,
        "conda": _version([str(conda_bin), "--version"]) if conda_bin else None,
    }


def assess(info, needed_disk_gb=None, cfg=None):
    """Return (critical_problems, warnings)."""
    crit, warns = [], []
    if info["os"] != "Linux":
        crit.append(f"unsupported OS {info['os']} (Linux required)")
    if not info["python_ok"]:
        crit.append("Python 3.10+ required")
    if not info["project_writable"]:
        crit.append(f"project location not writable: {info['project_path']}")
    if not info["temp_writable"]:
        crit.append(f"temporary directory not writable: {info['temp_dir']}")
    fs = info["filesystem"]
    if fs in ("vfat", "msdos", "fat"):
        crit.append(f"project filesystem is {fs}: {FS_LIMITS[fs][0]}. BAM/FASTQ/genome files routinely exceed "
                    "4 GB — place the project on an ext4/xfs filesystem")
    elif fs in FS_LIMITS:
        warns.append(f"project filesystem is {fs}: {FS_LIMITS[fs][0]}")
    if info["ram_total_gb"] < MIN_RAM_GB:
        warns.append(f"only {info['ram_total_gb']} GB RAM; HISAT2 alignment to a mammalian genome needs ~8 GB")
    if info["cpu_cores"] < MIN_CORES:
        warns.append(f"only {info['cpu_cores']} CPU core(s); processing will be slow")
    min_free = (cfg or {}).get("storage", {}).get("min_free_gb_warn", 50)
    if info["disk_free_gb"] < min_free:
        warns.append(f"only {info['disk_free_gb']} GB free disk (warning threshold {min_free} GB)")
    if needed_disk_gb and info["disk_free_gb"] < needed_disk_gb:
        crit.append(f"insufficient disk: need ~{needed_disk_gb:.1f} GB, {info['disk_free_gb']} GB free")
    if info["temp_filesystem"] == "tmpfs" and info["temp_free_gb"] < 20:
        warns.append(f"temp dir {info['temp_dir']} is RAM-backed with {info['temp_free_gb']} GB; "
                     "the pipeline uses <project>/temp instead")
    if not info.get("r"):
        warns.append("R not found (needed for DESeq2 stage)")
    if not info.get("conda"):
        warns.append("conda/mamba not found (needed for automatic tool installation)")
    return crit, warns


def display(info):
    ui.section("SYSTEM CHECK")
    ui.kv([
        ("OS", f"{info['distribution']} ({info['os']} {info['kernel']}{', WSL2' if info['wsl'] else ''})"),
        ("Architecture", info["architecture"]),
        ("CPU cores", info["cpu_cores"]),
        ("RAM", f"{info['ram_total_gb']} GB total, {info['ram_available_gb']} GB available"),
        ("Available disk", f"{info['disk_free_gb']} GB free of {info['disk_total_gb']} GB ({info['project_path']})"),
        ("Filesystem", f"{info['filesystem']} on {info['mountpoint']}"
                       + (f"  — {info['fs_note']}" if info['fs_note'] else "")),
        ("Max file size", info["max_file_size"] if isinstance(info["max_file_size"], str)
         else f"{info['max_file_size'] / GB:.0f} GB"),
        ("Writable", "yes" if info["project_writable"] else "NO"),
        ("Python", info["python"]),
        ("R", info["r"] or "not found"),
        ("Conda", info["conda"] or "not found"),
        ("Temporary space", f"{info['temp_dir']} ({info['temp_filesystem']}, {info['temp_free_gb']} GB free)"),
    ])


def mounted_storage():
    """List real mounted filesystems with free space (for choosing a project location)."""
    out = []
    seen = set()
    try:
        with open("/proc/mounts") as f:
            for line in f:
                dev, mnt, fs = line.split()[:3]
                if fs in ("proc", "sysfs", "cgroup", "cgroup2", "devpts", "mqueue", "overlay", "tmpfs",
                          "devtmpfs", "securityfs", "debugfs", "tracefs", "bpf", "fusectl", "configfs",
                          "binfmt_misc", "pstore", "hugetlbfs", "autofs", "nsfs", "squashfs"):
                    continue
                mnt = mnt.replace("\\040", " ")
                if mnt in seen or mnt.startswith(("/mnt/wslg", "/snap", "/boot", "/usr/lib/wsl", "/init")):
                    continue
                seen.add(mnt)
                try:
                    du = shutil.disk_usage(mnt)
                except OSError:
                    continue
                out.append((mnt, fs, round(du.free / GB, 1), writable_quick(mnt)))
    except OSError:
        pass
    return out


def writable_quick(path):
    return os.access(path, os.W_OK)
