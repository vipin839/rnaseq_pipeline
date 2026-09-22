"""Clean gene count matrix from featureCounts + validation."""
import csv
from pathlib import Path

from . import PipelineError


def build(genes, counts_by_bam, bam_to_sample, samples):
    """Return {sample: [int counts]} in `samples` order."""
    matrix = {}
    for bam, vals in counts_by_bam.items():
        sid = bam_to_sample[bam]
        if any(v != int(v) for v in vals):
            raise PipelineError(f"non-integer counts for {sid} (fractional counting?). DESeq2 requires integers",
                                stage="count_matrix",
                                remedy="disable featurecounts_parameters.fractional, or round deliberately")
        matrix[sid] = [int(v) for v in vals]
    missing = [s for s in samples if s not in matrix]
    if missing:
        raise PipelineError(f"samples missing from featureCounts output: {missing}", stage="count_matrix")
    return {s: matrix[s] for s in samples}


def check_against_assigned(matrix, assigned):
    """Column sums must equal featureCounts' Assigned per sample. Returns a list of mismatches."""
    out = []
    for sid, vals in matrix.items():
        total = sum(vals)
        if sid not in assigned:
            out.append(f"{sid}: no featureCounts summary")
        elif total != assigned[sid]:
            out.append(f"{sid}: column sum {total:,} != featureCounts Assigned {assigned[sid]:,}")
    return out


def collapse(genes, matrix, groups):
    """Sum technical replicates. groups: {bio_unit: [run ids]}. Returns new matrix keyed by bio_unit."""
    out = {}
    for unit, runs in groups.items():
        out[unit] = [sum(matrix[r][i] for r in runs) for i in range(len(genes))]
    return out


def write(genes, matrix, tsv, csv_path):
    samples = list(matrix)
    tsv, csv_path = Path(tsv), Path(csv_path)
    with open(tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene_ID"] + samples)
        for i, g in enumerate(genes):
            w.writerow([g] + [matrix[s][i] for s in samples])
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["Gene_ID"] + samples)
        for i, g in enumerate(genes):
            w.writerow([g] + [matrix[s][i] for s in samples])
    return tsv, csv_path


def validate_file(path, expected_samples=None, expected_genes=None):
    """Validate a written count matrix TSV. Returns summary dict; raises PipelineError on any problem."""
    path = Path(path)
    problems = []
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        head = next(r, None)
        if not head or head[0] != "Gene_ID":
            raise PipelineError("count matrix must start with a Gene_ID column", stage="count_matrix")
        samples = head[1:]
        if not samples:
            raise PipelineError("count matrix has no sample columns", stage="count_matrix")
        if len(set(samples)) != len(samples):
            problems.append("duplicate sample names")
        if any(s in ("Chr", "Start", "End", "Strand", "Length") for s in samples):
            problems.append("annotation columns present in the biological matrix")
        genes, seen = 0, set()
        totals = [0] * len(samples)
        zero_rows = 0
        for lineno, row in enumerate(r, 2):
            if len(row) != len(head):
                problems.append(f"line {lineno}: wrong number of columns")
                break
            gid = row[0]
            if not gid:
                problems.append(f"line {lineno}: empty gene ID")
            elif gid in seen:
                problems.append(f"duplicate gene ID {gid}")
            seen.add(gid)
            rowsum = 0
            for i, v in enumerate(row[1:]):
                if v in ("", "NA", "NaN", "nan"):
                    problems.append(f"missing value for {gid}/{samples[i]}")
                    continue
                try:
                    x = int(v)
                except ValueError:
                    problems.append(f"non-integer value {v!r} for {gid}/{samples[i]}")
                    continue
                if x < 0:
                    problems.append(f"negative value for {gid}/{samples[i]}")
                totals[i] += x
                rowsum += x
            if rowsum == 0:
                zero_rows += 1
            genes += 1
            if len(problems) > 20:
                break
    if genes == 0:
        problems.append("no genes")
    if expected_samples is not None and list(samples) != list(expected_samples):
        problems.append(f"sample columns {samples} do not match expected {list(expected_samples)}")
    if expected_genes is not None and genes != expected_genes:
        problems.append(f"{genes} genes, expected {expected_genes} (annotation gene count)")
    zero_libs = [s for s, t in zip(samples, totals) if t == 0]
    if zero_libs:
        problems.append(f"samples with zero total counts: {zero_libs}")
    if problems:
        raise PipelineError("count matrix validation failed:\n  - " + "\n  - ".join(problems[:20]),
                            stage="count_matrix")
    return {"genes": genes, "samples": samples, "library_sizes": dict(zip(samples, totals)),
            "all_zero_genes": zero_rows}


def check_metadata_match(matrix_samples, metadata_samples):
    ms, md = list(matrix_samples), list(metadata_samples)
    problems = []
    if set(ms) - set(md):
        problems.append(f"in count matrix but not in metadata: {sorted(set(ms) - set(md))}")
    if set(md) - set(ms):
        problems.append(f"in metadata but not in count matrix: {sorted(set(md) - set(ms))}")
    if not problems and ms != md:
        problems.append("sample order differs between count matrix and metadata")
    return problems


def write_summary(path, info, extra_lines=()):
    lines = [f"Genes: {info['genes']}", f"Samples: {len(info['samples'])}",
             f"Genes with zero counts in all samples: {info['all_zero_genes']}", "", "Library sizes (assigned reads):"]
    lines += [f"  {s}\t{n:,}" for s, n in info["library_sizes"].items()]
    lines += ["", *extra_lines]
    Path(path).write_text("\n".join(lines) + "\n")
    return path
