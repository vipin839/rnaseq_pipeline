"""StringTie2 reference-guided transcript quantification (independent of the DESeq2 count matrix).

Mode: `-e -G annotation.gtf` — abundances are estimated ONLY for annotated
transcripts (no novel transcript assembly), so results are comparable across
samples without a merge step. Output TPM/FPKM are NOT used for DESeq2.
"""
import os
import re

from . import PipelineError, runner
from .strandedness import STRINGTIE_FLAG

ABUND_HEADER = ["Gene ID", "Gene Name", "Reference", "Strand", "Start", "End", "Coverage", "FPKM", "TPM"]


def run_sample(project, sid, bam, gtf, strand, params, threads, log_file):
    sdir = project.path("stringtie", sid)
    sdir.mkdir(parents=True, exist_ok=True)
    out_gtf = sdir / f"{sid}.gtf"
    abund = project.path("stringtie", "abundance", f"{sid}.gene_abund.tab")
    tmp_gtf = sdir / f"{sid}.partial.gtf"
    tmp_ab = abund.with_name(f"{sid}.partial.gene_abund.tab")
    cmd = ["stringtie", bam, "-G", gtf, "-o", tmp_gtf, "-A", tmp_ab, "-p", str(threads)]
    if params.get("reference_guided_only", True):
        cmd.append("-e")
    if params.get("ballgown_tables"):
        cmd.append("-B")
    flag = STRINGTIE_FLAG[strand]
    if flag:
        cmd.append(flag)
    cmd += [str(x) for x in params.get("extra_args", [])]
    runner.run(cmd, stage="stringtie", sample=sid, log_file=log_file, description=f"{sid}: StringTie2")
    if runner.DRY_RUN:
        return out_gtf, abund
    problems = validate_outputs(tmp_gtf, tmp_ab)
    if problems:
        raise PipelineError(f"StringTie2 output invalid for {sid}: {'; '.join(problems)}", sample=sid,
                            stage="stringtie", remedy=f"see {log_file}")
    os.replace(tmp_gtf, out_gtf)
    os.replace(tmp_ab, abund)
    return out_gtf, abund


TPM_RE = re.compile(r'TPM "([^"]+)"')


def validate_outputs(gtf, abund):
    problems = []
    n_tx = 0
    try:
        with open(gtf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t")
                if len(c) != 9:
                    problems.append("malformed GTF line")
                    break
                if c[2] == "transcript":
                    n_tx += 1
                    if not TPM_RE.search(c[8]):
                        problems.append("transcript line without TPM")
                        break
    except OSError as e:
        problems.append(f"GTF unreadable: {e}")
    if n_tx == 0:
        problems.append("no transcripts in output GTF")
    try:
        with open(abund) as f:
            head = f.readline().rstrip("\n").split("\t")
            rows = sum(1 for _ in f)
        if head != ABUND_HEADER:
            problems.append(f"unexpected gene abundance header {head}")
        if rows == 0:
            problems.append("gene abundance table is empty")
    except OSError as e:
        problems.append(f"abundance table unreadable: {e}")
    return problems


def merge_tables(project, samples):
    """Combine per-sample outputs into gene- and transcript-level TPM matrices (stringtie/merged/)."""
    merged = project.path("stringtie", "merged")
    merged.mkdir(parents=True, exist_ok=True)
    gene = {}
    names = {}
    for sid in samples:
        with open(project.path("stringtie", "abundance", f"{sid}.gene_abund.tab")) as f:
            next(f)
            for line in f:
                c = line.rstrip("\n").split("\t")
                gene.setdefault(c[0], {}).setdefault(sid, 0.0)
                gene[c[0]][sid] += float(c[8])  # same gene on several loci (e.g. PAR) -> sum TPM
                names.setdefault(c[0], c[1])
    tx = {}
    for sid in samples:
        with open(project.path("stringtie", sid, f"{sid}.gtf")) as f:
            for line in f:
                c = line.split("\t")
                if len(c) == 9 and c[2] == "transcript":
                    tid = re.search(r'transcript_id "([^"]+)"', c[8]).group(1)
                    gid = re.search(r'gene_id "([^"]+)"', c[8]).group(1)
                    tpm = float(TPM_RE.search(c[8]).group(1))
                    tx.setdefault(tid, {"gene": gid})[sid] = tpm
    gpath = merged / "gene_tpm_matrix.tsv"
    with open(gpath, "w") as f:
        f.write("Gene_ID\tGene_Name\t" + "\t".join(samples) + "\n")
        for g in sorted(gene):
            f.write(f"{g}\t{names[g]}\t" + "\t".join(f"{gene[g].get(s, 0.0):.4f}" for s in samples) + "\n")
    tpath = merged / "transcript_tpm_matrix.tsv"
    with open(tpath, "w") as f:
        f.write("Transcript_ID\tGene_ID\t" + "\t".join(samples) + "\n")
        for t in sorted(tx):
            f.write(f"{t}\t{tx[t]['gene']}\t" + "\t".join(f"{tx[t].get(s, 0.0):.4f}" for s in samples) + "\n")
    readme = merged / "README.txt"
    readme.write_text(
        "StringTie2 reference-guided (-e) transcript-level quantification.\n"
        "TPM values are NOT the DESeq2 input; the DESeq2 gene count matrix comes from featureCounts\n"
        "(counts/gene_count_matrix.tsv).\n")
    return [gpath, tpath, readme], len(gene), len(tx)
