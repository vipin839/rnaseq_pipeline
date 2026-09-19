"""Disk-space estimation, pre-flight checks and optional (validated) cleanup.

Estimation factors (rules of thumb, documented in PIPELINE_METHODS.md):
  trimmed FASTQ.gz      ~ 0.9 x raw FASTQ.gz
  sorted BAM            ~ 0.8 x FASTQ.gz of the sample (+ ~1x of the largest sample as sort temp)
  fasterq-dump temp     ~ 8 x .sra size (uncompressed FASTQ) for the sample being converted
  StringTie per sample  ~ 0.05 GB (mammalian annotation)
"""
import shutil
from pathlib import Path

from . import ui

GB = 1024 ** 3


def fastq_bytes(project, samples, trimmed=None):
    total, largest = 0, 0
    for sid in samples:
        rec = project.samples[sid]
        size = 0
        if rec.get("r1"):
            for p in project.fastqs(sid, trimmed=trimmed):
                if p and Path(p).exists():
                    size += Path(p).stat().st_size
        else:
            size = sum(f[2] or 0 for f in (rec.get("download") or {}).get("files", []))
        total += size
        largest = max(largest, size)
    return total, largest


def estimate(project, stage, samples, ref_pkg=None):
    """Return estimated new disk usage in GB for a stage."""
    fq, largest = fastq_bytes(project, samples)
    if stage == "data_acquired":
        max_reads = int((project.config().get("download") or {}).get("max_reads") or 0)
        rec_sizes = []
        for s in samples:
            rec = project.samples[s]
            if rec.get("r1"):
                continue
            size = sum(f[2] or 0 for f in (rec.get("download") or {}).get("files", []))
            n = str((rec.get("metadata") or {}).get("read_count", ""))
            if max_reads and n.isdigit() and int(n) > 0:
                size *= min(1.0, max_reads / int(n))  # pilot subset: only the first N reads are fetched
            rec_sizes.append(size)
        return sum(rec_sizes) / GB * 1.05 + (max(rec_sizes, default=0) * 8 / GB
                                               if project.state.get("data_source") == "sra" else 0)
    if stage == "trimming_completed":
        return fq * 0.9 / GB
    if stage == "alignment_completed":
        return (fq * 0.8 + largest) / GB
    if stage == "stringtie_completed":
        return 0.05 * len(samples)
    if stage == "reference_ready" and ref_pkg:
        return (ref_pkg.get("genome", {}).get("approx_size_gb", 1) * 4.5
                + ref_pkg.get("annotation", {}).get("approx_size_gb", 0.05) * 12
                + ref_pkg.get("index_approx_size_gb", 5))
    return 0.1


def preflight(path, needed_gb, cfg):
    """Returns True if OK to proceed. Asks the user when borderline; refuses when insufficient."""
    free = shutil.disk_usage(path).free / GB
    factor = cfg.get("storage", {}).get("borderline_factor", 1.5)
    ui.kv([("Estimated new storage", f"{needed_gb:.1f} GB"), ("Free space", f"{free:.1f} GB ({path})")])
    if needed_gb > free:
        ui.error(f"insufficient disk space: need ~{needed_gb:.1f} GB, only {free:.1f} GB free")
        ui.info("Free space or move the project to a larger filesystem, then resume.")
        return False
    if needed_gb * factor > free:
        ui.warn(f"storage is borderline (less than {factor}x the estimate is free)")
        return ui.ask_yes_no("Proceed anyway?", default=False)
    return True


def cleanup_candidates(project, checkpoints):
    """Intermediate files that may be removed ONLY because a validated downstream output exists."""
    cands = []
    if checkpoints.exists("alignment_completed") and not checkpoints.verify_files("alignment_completed"):
        for sid, rec in project.samples.items():
            t = rec.get("trimmed")
            if t:
                for m in ("r1", "r2"):
                    if t.get(m) and project.abs(t[m]).exists():
                        cands.append((project.abs(t[m]), "trimmed FASTQ (BAM validated; raw FASTQ kept)"))
    if checkpoints.exists("fastq_verified") and not checkpoints.verify_files("fastq_verified"):
        for d in project.path("data", "raw").glob("*/*.sra"):
            cands.append((d, ".sra cache (FASTQ extracted and validated)"))
    for d in project.path("temp").iterdir() if project.path("temp").exists() else []:
        cands.append((d, "temporary working file"))
    return cands
