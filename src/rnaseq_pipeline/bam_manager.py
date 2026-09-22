"""BAM integrity validation and samtools-based QC."""
import json
from pathlib import Path

from . import runner


def header_sorted(bam):
    out = runner.tool_output(["samtools", "view", "-H", str(bam)]) or ""
    for line in out.splitlines():
        if line.startswith("@HD"):
            return "SO:coordinate" in line
    return False


def flagstat(bam, threads=1):
    out = runner.tool_output(["samtools", "flagstat", "-@", str(threads), "-O", "json", str(bam)], timeout=7200)
    if not out:
        return None
    try:
        start = out.index("{")
        return json.loads(out[start:])["QC-passed reads"]
    except (ValueError, KeyError, json.JSONDecodeError):
        return None


def validate_bam(bam, expected_input_reads, paired, threads=1):
    """Return (ok, problems, metrics). expected_input_reads = FASTQ reads (pairs for PE)."""
    bam = Path(bam)
    problems, m = [], {}
    if not bam.exists() or bam.stat().st_size == 0:
        return False, ["BAM missing or empty"], m
    q = runner.tool_output(["samtools", "quickcheck", "-v", str(bam)])
    if q is None or q.strip():
        problems.append(f"samtools quickcheck failed (truncated/corrupt BAM, missing EOF block): {q.strip()[:200] if q else ''}")
        return False, problems, m
    if not header_sorted(bam):
        problems.append("BAM header is not SO:coordinate (not sorted)")
    bai = Path(str(bam) + ".bai")
    if not bai.exists():
        problems.append("BAM index (.bai) missing")
    elif bai.stat().st_mtime < bam.stat().st_mtime:
        problems.append("BAM index is older than the BAM")
    fs = flagstat(bam, threads)
    if fs is None:
        problems.append("samtools flagstat failed")
        return False, problems, m
    primary = fs.get("primary", fs["total"] - fs.get("secondary", 0) - fs.get("supplementary", 0))
    m.update({
        "total_records": fs["total"], "primary": primary, "secondary": fs.get("secondary", 0),
        "supplementary": fs.get("supplementary", 0), "mapped": fs.get("mapped", 0),
        "primary_mapped": fs.get("primary mapped", fs.get("mapped", 0)),
        "properly_paired": fs.get("properly paired", 0), "singletons": fs.get("singletons", 0),
        "duplicates": fs.get("duplicates", 0),
    })
    m["unmapped"] = primary - m["primary_mapped"]
    m["mapped_pct"] = round(100 * m["primary_mapped"] / primary, 2) if primary else 0.0
    if paired and primary:
        m["properly_paired_pct"] = round(100 * m["properly_paired"] / primary, 2)
    expected_primary = expected_input_reads * (2 if paired else 1) if expected_input_reads is not None else None
    m["expected_primary"] = expected_primary
    if expected_primary is not None and primary != expected_primary:
        problems.append(f"primary alignments ({primary:,}) != input reads ({expected_primary:,}); "
                        "BAM is incomplete or inputs changed")
    if paired and fs.get("paired in sequencing", 0) != primary:
        problems.append("not all primary records are flagged as paired")
    return not problems, problems, m


def qc_reports(bam, out_dir, sid, threads, log_file):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    runner.run(["samtools", "idxstats", bam], stage="bam_qc", sample=sid, stdout_file=out_dir / f"{sid}.idxstats",
               log_file=log_file)
    runner.run(["samtools", "flagstat", "-@", str(threads), bam], stage="bam_qc", sample=sid,
               stdout_file=out_dir / f"{sid}.flagstat", log_file=log_file)
    runner.run(["samtools", "stats", "-@", str(threads), bam], stage="bam_qc", sample=sid,
               stdout_file=out_dir / f"{sid}.samtools_stats", log_file=log_file,
               description=f"{sid}: samtools stats")
    return [out_dir / f"{sid}.{s}" for s in ("idxstats", "flagstat", "samtools_stats")]


def parse_stats_summary(path):
    out = {}
    try:
        for line in Path(path).read_text().splitlines():
            if line.startswith("SN\t"):
                _, key, val = line.split("\t")[:3]
                out[key.rstrip(":")] = val
    except OSError:
        pass
    return out


SUMMARY_COLUMNS = ["sample", "status", "input_reads", "unit", "overall_alignment_rate", "uniquely_aligned_pct",
                   "multi_aligned_pct", "total_records", "primary", "primary_mapped", "mapped_pct", "unmapped",
                   "properly_paired", "properly_paired_pct", "secondary", "supplementary", "problems"]


def write_summary(rows, path):
    with open(path, "w") as f:
        f.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in SUMMARY_COLUMNS) + "\n")
    return path
