"""HISAT2 alignment streamed into samtools sort (no SAM on disk), then indexed and validated."""
import os
import re
from pathlib import Path

from . import PipelineError, runner


def sort_memory(cfg, threads, mem_gb):
    v = cfg.get("samtools_sort_memory_per_thread", "auto")
    if v not in (None, "auto"):
        if not re.match(r"^\d+[KMG]$", str(v)):
            raise PipelineError("samtools_sort_memory_per_thread must look like 768M or 2G")
        return str(v)
    # leave half the memory budget for HISAT2 (~8 GB for mammals); split the rest among sort threads
    per = int(max(256, min(2048, (mem_gb * 1024 * 0.4) / max(1, threads))))
    return f"{per}M"


def hisat2_command(project, sid, ref, cfg, threads, summary_file):
    r1, r2 = project.fastqs(sid)
    hp = cfg.get("hisat2_parameters", {})
    cmd = ["hisat2", "-p", str(threads), "-x", ref["index_prefix"], "--new-summary",
           "--summary-file", summary_file, "--rg-id", sid, "--rg", f"SM:{sid}"]
    if hp.get("dta", True):
        cmd.append("--dta")
    ss = ref.get("splice_sites")
    if (hp.get("use_known_splice_sites", True) and not ref.get("index_has_splice_sites") and ss
            and Path(ss).exists() and Path(ss).stat().st_size > 0):
        cmd += ["--known-splicesite-infile", ss]
    elif ref.get("no_splice_sites"):
        cmd.append("--no-spliced-alignment")  # genome without introns (e.g. bacteria)
    if r2:
        cmd += ["-1", r1, "-2", r2]
        if hp.get("no_mixed"):
            cmd.append("--no-mixed")
        if hp.get("no_discordant"):
            cmd.append("--no-discordant")
    else:
        cmd += ["-U", r1]
    cmd += [str(x) for x in hp.get("extra_args", [])]
    return cmd


def parse_hisat2_summary(path):
    text = Path(path).read_text()
    out = {}
    m = re.search(r"Total (pairs|reads):\s+(\d+)", text)
    if m:
        out["unit"] = m.group(1)
        out["total"] = int(m.group(2))
    m = re.search(r"Overall alignment rate:\s+([\d.]+)%", text)
    if m:
        out["overall_alignment_rate"] = float(m.group(1))
    m = re.search(r"Aligned concordantly 1 time:\s+\d+\s+\(([\d.]+)%\)", text) or \
        re.search(r"Aligned 1 time:\s+\d+\s+\(([\d.]+)%\)", text)
    if m:
        out["uniquely_aligned_pct"] = float(m.group(1))
    m = re.search(r"Aligned concordantly >1 times:\s+\d+\s+\(([\d.]+)%\)", text) or \
        re.search(r"Aligned >1 times:\s+\d+\s+\(([\d.]+)%\)", text)
    if m:
        out["multi_aligned_pct"] = float(m.group(1))
    return out


def align_sample(project, sid, ref, cfg, threads, mem_gb, log_file):
    bam_dir = project.path("alignment", "bam")
    rep_dir = project.path("alignment", "reports")
    tmp_dir = project.path("temp", f"sort_{sid}")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    partial = bam_dir / f"{sid}.partial.bam"
    final = bam_dir / f"{sid}.sorted.bam"
    summary = rep_dir / f"{sid}.hisat2.summary"
    for stale in (partial, Path(str(partial) + ".bai")):
        if stale.exists():
            stale.unlink()  # our own incomplete output from an interrupted run
    for t in tmp_dir.glob("*"):
        t.unlink()
    h = hisat2_command(project, sid, ref, cfg, threads, summary)
    s = ["samtools", "sort", "-@", str(max(1, threads - 1)), "-m", sort_memory(cfg, threads, mem_gb),
         "-T", tmp_dir / "chunk", "-o", partial, "-"]
    runner.run_pipeline([h, s], stage="alignment", sample=sid, log_file=log_file,
                        description=f"{sid}: HISAT2 | samtools sort")
    if runner.DRY_RUN:
        return final
    runner.run(["samtools", "quickcheck", "-v", partial], stage="alignment", sample=sid, log_file=log_file)
    os.replace(partial, final)
    runner.run(["samtools", "index", "-@", str(threads), final], stage="alignment", sample=sid,
               log_file=log_file, description=f"{sid}: indexing BAM")
    try:
        tmp_dir.rmdir()
    except OSError:
        pass
    return final
