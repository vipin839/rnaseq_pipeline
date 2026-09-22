"""End-to-end test through the real interactive CLI on a controlled synthetic dataset with known truth.

Covers: project creation, local paired-end import, FASTQ validation, FastQC/MultiQC, quality gate,
conditional fastp + post-trim QC, custom reference validation + HISAT2 index, alignment + BAM validation,
BAM QC, strandedness inference, StringTie2, featureCounts, count matrix, design confirmation,
Python->R transition, DESeq2, plots, report, manifest, checkpoints and resume.
Run:  ~/miniforge3/envs/rnaseq-tools/bin/python -m pytest tests/integration -q
"""
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import RSCRIPT, ROOT, have

sys.path.insert(0, str(ROOT / "tests" / "data"))
import make_synthetic  # noqa: E402

pytestmark = pytest.mark.skipif(
    not have("fastqc", "hisat2", "samtools", "featureCounts", "stringtie", "fastp", "multiqc")
    or RSCRIPT is None,
    reason="bioinformatics tools / R env not installed")

LAUNCHER = ROOT / "rnaseq_pipeline"


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    d = tmp_path_factory.mktemp("synth")
    make_synthetic.make(d, pairs_per_sample=40000)
    return d


def run_cli(args, answers, timeout=1800):
    env = dict(os.environ, PYTHONNOUSERSITE="1")
    r = subprocess.run([str(LAUNCHER), *args], input="\n".join(answers) + "\n", capture_output=True, text=True,
                       timeout=timeout, env=env)
    return r


def new_project_answers(name, parent, fastq_dir, fasta, gtf):
    return ["1", name, str(parent), "5", str(fastq_dir), "", "1", "y",
            "5", str(fasta), str(gtf), "Synthetic", "SYNTH1", "", "synth-v1", "y", "n", "1", "y",
            "1", "1", "1", "1", "1",           # data, fastq, raw qc, gate screen, accept recommendation
            "1", "1", "1", "1",                # trimming, reference, alignment, bam qc
            "1", "y",                          # strandedness + confirm inferred value
            "1", "1", "1",                     # stringtie, featurecounts, count matrix
            "1", "2", "Ctrl", "Control", "y", "2", "Treat", "Treatment", "y", "7", "1", "y", "8", "y",
            "1", "y",                          # DESeq2 + proceed
            "1",                               # report
            "8"]                               # exit


@pytest.fixture(scope="module")
def project(synth, tmp_path_factory):
    work = tmp_path_factory.mktemp("work")
    cfg = work / "override.yaml"
    cfg.write_text(f"reference_store: {work / 'refstore'}\n")
    r = run_cli(["--config", str(cfg), "--projects-dir", str(work / "projects")],
                new_project_answers("E2E", work / "projects", synth / "fastq", synth / "genome.fa",
                                    synth / "annotation.gtf"))
    (work / "cli_stdout.txt").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-2000:]
    return work / "projects" / "RNAseq_E2E", r.stdout, cfg


def test_all_stages_checkpointed(project):
    p, out, _ = project
    for st in ["project_initialized", "data_acquired", "fastq_verified", "raw_qc_completed", "quality_assessed",
               "trimming_completed", "reference_ready", "alignment_completed", "bam_qc_completed",
               "strandedness_determined", "stringtie_completed", "featurecounts_completed",
               "count_matrix_completed", "design_confirmed", "deseq2_completed", "report_generated"]:
        assert (p / "checkpoints" / f"{st}.json").exists(), st
    assert "FAILED" not in out


def test_trimming_and_qc(project):
    p, out, _ = project
    assert "fastp trimming recommended" in out
    state = json.loads((p / "project.json").read_text())
    assert state["qc_decision"]["decision"] == "trim"
    for s in ("Ctrl1", "Treat3"):
        assert (p / "data" / "trimmed" / f"{s}_R1.trimmed.fastq.gz").exists()
        assert (p / "data" / "fastq" / f"{s}_R1.fastq.gz").is_symlink()  # raw reads untouched (linked)
    assert (p / "qc" / "multiqc_raw" / "multiqc_report.html").exists()
    assert (p / "qc" / "multiqc_trimmed" / "multiqc_report.html").exists()
    for f in ("quality_assessment.json", "quality_assessment.tsv", "quality_assessment.html"):
        assert (p / "qc" / "assessment" / "raw" / f).exists()


def test_alignment_integrity(project):
    p, _, _ = project
    rows = list(csv.DictReader(open(p / "alignment" / "reports" / "alignment_summary.tsv"), delimiter="\t"))
    assert len(rows) == 6 and all(r["status"] == "PASS" for r in rows)
    for r in rows:
        assert int(r["primary"]) == 2 * int(r["input_reads"])
        assert float(r["overall_alignment_rate"]) > 95
    assert not list((p / "alignment" / "bam").glob("*.partial.bam"))
    assert not list((p / "alignment" / "sam").glob("*.sam"))  # SAM never written to disk


def test_strandedness_inferred_correctly(project):
    p, _, _ = project
    s = json.loads((p / "project.json").read_text())["strandedness"]
    assert s["value"] == "reverse" and s["featurecounts_flag"] == "2" and s["source"].startswith("RSeQC")


def test_count_matrix(project):
    p, _, _ = project
    lines = (p / "counts" / "gene_count_matrix.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["Gene_ID", "Ctrl1", "Ctrl2", "Ctrl3", "Treat1", "Treat2", "Treat3"]
    assert len(lines) == 61
    assert (p / "stringtie" / "merged" / "gene_tpm_matrix.tsv").exists()


def test_deseq2_recovers_truth(project, synth):
    p, _, _ = project
    truth = {r["gene_id"]: r["direction"] for r in csv.DictReader(open(synth / "truth.tsv"), delimiter="\t")}
    d = p / "results" / "deseq2" / "condition_Treatment_vs_Control"
    res = {r["gene_id"]: r for r in csv.DictReader(open(d / "full_results.tsv"), delimiter="\t")}
    true_de = [g for g, x in truth.items() if x != "none"]
    recovered = sum(1 for g in true_de if res.get(g, {}).get("regulation") == truth[g])
    false_pos = sum(1 for g, x in truth.items() if x == "none" and res.get(g, {}).get("regulation", "ns") != "ns")
    assert recovered >= 14, f"only {recovered}/16 true DE genes recovered"
    assert false_pos <= 1
    for f in ("significant_results.tsv", "upregulated.tsv", "downregulated.tsv", "top_upregulated.tsv",
              "top_downregulated.tsv", "statistics.txt"):
        assert (d / f).exists()
    stats = (d / "statistics.txt").read_text()
    assert "padj < 0.05 AND log2FoldChange >= 1" in stats
    for plot in ("pca", "sample_distance_heatmap", "sample_correlation_heatmap", "library_sizes"):
        assert (p / "results" / "plots" / f"{plot}.png").exists() and (p / "results" / "plots" / f"{plot}.pdf").exists()
    for plot in ("ma_plot", "volcano_plot", "significant_genes_heatmap", "top_genes_heatmap"):
        assert (p / "results" / "plots" / "condition_Treatment_vs_Control" / f"{plot}.png").exists()
    assert (d / "plot_data" / "volcano_data.tsv").exists()


def test_report_manifest_logs(project):
    p, _, _ = project
    html = (p / "reports" / "final_pipeline_report.html").read_text()
    for section in ("Project overview", "Alignment (HISAT2)", "Library strandedness", "Differential genes",
                    "Reproducibility"):
        assert section in html
    m = json.loads((p / "pipeline_manifest" / "manifest.json").read_text())
    assert m["tool_versions"]["hisat2"] and m["r_version"] and m["strandedness"]["value"] == "reverse"
    for f in ("environment.yml", "package_versions.txt", "system_information.txt", "R_sessionInfo.txt"):
        assert (p / "pipeline_manifest" / f).exists()
    for f in ("pipeline.log", "command_history.log", "commands.jsonl", "software_versions.tsv"):
        assert (p / "logs" / f).exists()
    cmds = [json.loads(l) for l in (p / "logs" / "commands.jsonl").read_text().splitlines()]
    assert any(c["command"][0][0] == "hisat2" for c in cmds if c.get("command"))
    assert list((p / "config").glob("config_used_*.yaml"))


def test_resume_revalidates_without_rerun(project):
    p, _, cfg = project
    before = len((p / "logs" / "commands.jsonl").read_text().splitlines())
    r = run_cli(["--project", str(p)], ["1", "10"], timeout=600)
    assert r.returncode == 0, r.stdout[-3000:]
    assert r.stdout.count("already complete (checkpoint revalidated)") == 15
    after_cmds = (p / "logs" / "commands.jsonl").read_text().splitlines()[before:]
    assert not any("hisat2" in l or "featureCounts" in l for l in after_cmds)


def test_dry_run_executes_nothing(project):
    p, _, _ = project
    before = len((p / "logs" / "commands.jsonl").read_text().splitlines())
    r = run_cli(["--dry-run", "--project", str(p)], [], timeout=600)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "DRY RUN" in r.stdout
    assert len((p / "logs" / "commands.jsonl").read_text().splitlines()) == before


def test_corrupted_fastq_stops_pipeline(synth, tmp_path):
    """A truncated FASTQ must be detected before alignment; the user declines to continue; logs are kept."""
    fq = tmp_path / "fq"
    fq.mkdir()
    for f in (synth / "fastq").glob("*.fastq.gz"):
        shutil.copy(f, fq / f.name)
    bad = fq / "Treat2_R2.fastq.gz"
    data = bad.read_bytes()
    bad.write_bytes(data[: len(data) // 3])
    cfg = tmp_path / "o.yaml"
    cfg.write_text(f"reference_store: {tmp_path / 'refstore'}\n")
    answers = ["1", "Broken", str(tmp_path / "projects"), "5", str(fq), "", "1", "y",
               "5", str(synth / "genome.fa"), str(synth / "annotation.gtf"), "Synthetic", "SYNTH1", "", "v1", "y",
               "n", "1", "y", "1", "1", "n", "8"]
    r = run_cli(["--config", str(cfg), "--projects-dir", str(tmp_path / "projects")], answers, timeout=900)
    p = tmp_path / "projects" / "RNAseq_Broken"
    assert "STATUS: FAILED" in r.stdout and "Treat2" in r.stdout
    rep = (p / "data" / "metadata" / "fastq_validation_report.tsv").read_text()
    assert "FAILED" in rep and "Treat2_R2" in rep
    assert not (p / "checkpoints" / "fastq_verified.json").exists()
    assert (p / "logs" / "pipeline_errors.log").read_text().strip()
    st = json.loads((p / "project.json").read_text())
    assert st["samples"]["Treat2"]["status"] == "FAILED"
    assert not list((p / "alignment" / "bam").glob("*.bam"))
