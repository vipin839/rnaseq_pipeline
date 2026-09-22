"""featureCounts gene-level counting and output validation."""
import os
from pathlib import Path

from . import PipelineError, runner
from .strandedness import FEATURECOUNTS_FLAG

LOW_ASSIGNED_WARN = 30.0  # % assigned reads below which a warning is shown (documented in PIPELINE_METHODS)


def command(bams, gtf, out_txt, strand, paired, params, threads, tmp_dir):
    p = params
    cmd = ["featureCounts", "-T", str(threads), "-a", gtf, "-o", out_txt, "-t", p["feature_type"],
           "-g", p["attribute"], "-s", FEATURECOUNTS_FLAG[strand], "-Q", str(p.get("min_mapping_quality", 0)),
           "--tmpDir", tmp_dir]
    if paired:
        cmd += ["-p", "--countReadPairs"]
        if p.get("require_both_ends_mapped", True):
            cmd.append("-B")
        if not p.get("count_chimeric", False):
            cmd.append("-C")
    if p.get("count_multimapping"):
        cmd.append("-M")
    if p.get("fractional"):
        cmd.append("--fraction")
    if p.get("primary_only"):
        cmd.append("--primary")
    cmd += [str(x) for x in p.get("extra_args", [])]
    cmd += [str(b) for b in bams]
    return cmd


def run(project, samples, bams, gtf, strand, paired, params, threads, log_file):
    out_dir = project.path("featurecounts")
    tmp_txt = out_dir / "featurecounts.partial.txt"
    final_txt = out_dir / "featurecounts.txt"
    final_sum = out_dir / "featurecounts.summary"
    tmp_dir = project.path("temp", "featurecounts")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = command(bams, gtf, tmp_txt, strand, paired, params, threads, tmp_dir)
    runner.run(cmd, stage="featurecounts", log_file=log_file, description=f"featureCounts on {len(bams)} BAM(s)")
    if runner.DRY_RUN:
        return final_txt, final_sum, cmd
    os.replace(tmp_txt, final_txt)
    os.replace(Path(str(tmp_txt) + ".summary"), final_sum)
    return final_txt, final_sum, cmd


def parse(txt_path, bams):
    """Parse featureCounts table. Returns (gene_ids, {bam: [counts]}) with strict validation."""
    problems = []
    with open(txt_path) as f:
        first = f.readline()
        if not first.startswith("# Program:featureCounts"):
            problems.append("missing featureCounts program header")
        head = f.readline().rstrip("\n").split("\t")
        fixed = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        if head[:6] != fixed:
            raise PipelineError(f"unexpected featureCounts columns: {head[:6]}", stage="featurecounts")
        cols = head[6:]
        if [os.path.abspath(c) for c in cols] != [os.path.abspath(str(b)) for b in bams]:
            raise PipelineError("featureCounts sample columns do not match the BAM list/order", stage="featurecounts")
        genes, counts, seen = [], {b: [] for b in cols}, set()
        for lineno, line in enumerate(f, 3):
            c = line.rstrip("\n").split("\t")
            if len(c) != len(head):
                raise PipelineError(f"featureCounts line {lineno}: {len(c)} columns, expected {len(head)}",
                                    stage="featurecounts")
            gid = c[0]
            if not gid or gid in seen:
                raise PipelineError(f"empty or duplicate gene ID at line {lineno}: {gid!r}", stage="featurecounts")
            seen.add(gid)
            genes.append(gid)
            for b, v in zip(cols, c[6:]):
                try:
                    x = float(v)
                except ValueError:
                    raise PipelineError(f"non-numeric count {v!r} (gene {gid})", stage="featurecounts") from None
                if x < 0:
                    raise PipelineError(f"negative count for {gid}", stage="featurecounts")
                counts[b].append(x)
    if problems:
        raise PipelineError("; ".join(problems), stage="featurecounts")
    if not genes:
        raise PipelineError("featureCounts table has no genes", stage="featurecounts")
    return genes, counts


def reconcile(stats, input_fragments, alignment_rows):
    """Independent check of featureCounts against the BAM/FASTQ it was given.

    stats: parse_summary() output; input_fragments: {sample: reads (SE) or read pairs (PE)} from FASTQ
    validation; alignment_rows: {sample: alignment_summary row}. Returns (errors, warnings).
    Invariants (measured on real and synthetic data, see docs/VERIFICATION_MATRIX.md):
      * featureCounts sees every fragment in the BAM: total >= input fragments
      * it cannot assign more fragments than exist: assigned <= input fragments
      * the only surplus comes from secondary/supplementary alignments (warning otherwise)
    """
    errors, warns = [], []
    for sid, st in stats.items():
        n = input_fragments.get(sid)
        if not n:
            continue
        if st["total"] < n:
            errors.append(f"{sid}: featureCounts processed {st['total']:,} fragments but the BAM holds {n:,} — "
                          "it did not read the whole BAM")
        if st["assigned"] > n:
            errors.append(f"{sid}: {st['assigned']:,} fragments assigned but only {n:,} exist")
        row = alignment_rows.get(sid) or {}
        extra = int(row.get("secondary") or 0) + int(row.get("supplementary") or 0)
        if row and st["total"] > n + extra:
            warns.append(f"{sid}: featureCounts total {st['total']:,} exceeds fragments + secondary/supplementary "
                         f"alignments ({n + extra:,})")
    return errors, warns


def parse_summary(sum_path, bams, samples):
    rows = {}
    with open(sum_path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = line.rstrip("\n").split("\t")
            rows[c[0]] = [int(x) for x in c[1:]]
    if len(head) - 1 != len(bams):
        raise PipelineError("featureCounts summary sample columns mismatch", stage="featurecounts")
    out = {}
    for i, sid in enumerate(samples):
        total = sum(v[i] for v in rows.values())
        assigned = rows.get("Assigned", [0] * len(samples))[i]
        out[sid] = {"assigned": assigned, "total": total,
                    "assigned_pct": round(100 * assigned / total, 2) if total else 0.0,
                    **{k: v[i] for k, v in rows.items() if v[i]}}
    return out
