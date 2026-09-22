"""Conditional adapter/quality trimming with fastp (raw reads are never modified)."""
import json
import os
from pathlib import Path

from . import PipelineError, runner


def fastp_command(r1, r2, o1, o2, json_path, html_path, sample, params, threads):
    p = params
    cmd = ["fastp", "-i", r1, "-o", o1]
    if r2:
        cmd += ["-I", r2, "-O", o2]
        if p.get("detect_adapter_for_pe", True):
            cmd.append("--detect_adapter_for_pe")
    cmd += ["--qualified_quality_phred", str(p["qualified_quality_phred"]),
            "--unqualified_percent_limit", str(p["unqualified_percent_limit"]),
            "--length_required", str(p["length_required"])]
    if p.get("cut_right"):
        cmd += ["--cut_right", "--cut_right_window_size", str(p["cut_right_window_size"]),
                "--cut_right_mean_quality", str(p["cut_right_mean_quality"])]
    cmd += ["--thread", str(max(1, min(threads, 16))), "--json", json_path, "--html", html_path,
            "--report_title", f"fastp: {sample}"]
    cmd += [str(x) for x in p.get("extra_args", [])]
    return cmd


def trim_sample(project, sid, params, threads, log_file):
    r1, r2 = project.fastqs(sid, trimmed=False)
    out = project.path("data", "trimmed")
    rep = project.path("qc", "fastp")
    rep.mkdir(parents=True, exist_ok=True)
    final1 = out / f"{sid}_R1.trimmed.fastq.gz"
    final2 = out / f"{sid}_R2.trimmed.fastq.gz" if r2 else None
    tmp1 = out / f"{sid}_R1.partial.fastq.gz"
    tmp2 = out / f"{sid}_R2.partial.fastq.gz" if r2 else None
    j, h = rep / f"{sid}.fastp.json", rep / f"{sid}.fastp.html"
    for raw in (r1, r2):
        if raw and Path(raw).resolve() in (final1.resolve(), final2.resolve() if final2 else None):
            raise PipelineError("refusing to overwrite raw reads", sample=sid, stage="trimming")
    cmd = fastp_command(r1, r2, tmp1, tmp2, j, h, sid, params, threads)
    runner.run(cmd, stage="trimming", sample=sid, log_file=log_file, description=f"{sid}: fastp")
    if runner.DRY_RUN:
        return None
    try:
        rep_json = json.loads(Path(j).read_text())
        after = rep_json["summary"]["after_filtering"]["total_reads"]
        before = rep_json["summary"]["before_filtering"]["total_reads"]
    except (OSError, KeyError, json.JSONDecodeError) as e:
        raise PipelineError(f"fastp JSON report unreadable: {e}", sample=sid, stage="trimming")
    per_file = after // 2 if r2 else after
    os.replace(tmp1, final1)
    if r2:
        os.replace(tmp2, final2)
    project.samples[sid]["trimmed"] = {"r1": project.rel(final1), "r2": project.rel(final2) if final2 else None,
                                       "fastp_json": project.rel(j), "fastp_html": project.rel(h),
                                       "reads_before": before, "reads_after": after,
                                       "expected_reads_per_file": per_file, "command": runner.fmt(cmd)}
    project.save()
    return per_file
