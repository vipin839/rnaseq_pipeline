"""End-to-end test through the real interactive CLI on a controlled synthetic dataset with known truth.

Covers: project creation, local paired-end import, FASTQ validation, FastQC/MultiQC, quality gate,
conditional fastp + post-trim QC, custom reference validation + HISAT2 index, alignment + BAM validation,
BAM QC, strandedness inference, StringTie2, featureCounts, count matrix, design confirmation,
Python->R transition, DESeq2, plots, report, manifest, checkpoints and resume.
Run:  ~/miniforge3/envs/rnaseq-tools/bin/python -m pytest tests/integration -q
"""
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path
import sys

import pytest

from conftest import RSCRIPT, ROOT, have

sys.path.insert(0, str(ROOT / "tests" / "data"))
import make_synthetic  # noqa: E402

pytestmark = pytest.mark.skipif(
    not have("fastqc", "hisat2", "samtools", "featureCounts", "stringtie", "fastp", "multiqc")
    or RSCRIPT is None,
    reason="bioinformatics tools / R env not installed")

LAUNCHER = ROOT / "rnaseq_pipeline"
FAKE_KEY = "fakeNCBIkey0123456789abcdef"  # must never appear in any project file


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
    cfg.write_text(f"reference_store: {work / 'refstore'}\n"
                   f"ncbi:\n  email: \"e2e@example.org\"\n  api_key: \"{FAKE_KEY}\"\n")
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
    # F9: inputs, reference and commands are identifiable from the manifest alone
    assert len(m["input_files"]) == 6
    assert all(len(v["r1"]["sha256"]) == 64 and len(v["r2"]["sha256"]) == 64 for v in m["input_files"].values())
    assert len(m["reference_checksums"]["genome"]["sha256"]) == 64
    assert len(m["reference_checksums"]["annotation"]["sha256"]) == 64
    assert m["commands"]["count"] > 20 and m["commands"]["log"] == "logs/commands.jsonl"
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


def test_no_credentials_anywhere_in_project(project):
    """F0 regression: the API key/email given via --config must not be written to any project file."""
    p, out, _ = project
    assert FAKE_KEY not in out
    for f in p.rglob("*"):
        if f.is_file() and not f.is_symlink() and f.stat().st_size < 50_000_000:
            data = f.read_bytes()
            assert FAKE_KEY.encode() not in data, f"API key leaked into {f.relative_to(p)}"
            assert b"e2e@example.org" not in data, f"email leaked into {f.relative_to(p)}"


# ---------------------------------------------------------------- resume / invalidation matrix
def _status_after(project_dir, tmp_path, mutate):
    """Copy the finished project, apply `mutate(copy)`, return {stage_key: status} as resume would see it."""
    from rnaseq_pipeline import environment_manager as EM, system_check as SC, workflow
    from rnaseq_pipeline.project import Project
    dst = tmp_path / "copy"
    shutil.copytree(project_dir, dst, symlinks=True)
    (dst / ".lock").unlink(missing_ok=True)
    p = Project.open(dst)
    mutate(p)
    p = Project.open(dst)
    cfg = p.config()
    envs = EM.Environments(cfg)
    envs.activate()
    ctx = workflow.Context(p, cfg, envs, SC.collect(p.root))
    memo = {}
    return {st.key: workflow.stage_status(ctx, st, memo=memo)[0] for st in workflow.STAGES}, memo


def _set(p, key, value):
    from rnaseq_pipeline import config as C
    cfg = C.load_yaml(p.config_path)
    C.set_value(cfg, key, value)
    C.save_yaml(cfg, p.config_path)


STAGE_ORDER = ["data_acquired", "fastq_verified", "raw_qc_completed", "quality_assessed", "trimming_completed",
               "reference_ready", "alignment_completed", "bam_qc_completed", "strandedness_determined",
               "stringtie_completed", "featurecounts_completed", "count_matrix_completed", "design_confirmed",
               "deseq2_completed", "report_generated"]


@pytest.mark.parametrize("name,mutate,first_invalid", [
    ("unchanged", lambda p: None, None),
    ("alpha changed", lambda p: _set(p, "alpha", 0.01), "deseq2_completed"),
    ("log2FC threshold changed", lambda p: _set(p, "log2fc_threshold", 2), "deseq2_completed"),
    ("featureCounts MAPQ changed", lambda p: _set(p, "featurecounts_parameters.min_mapping_quality", 10),
     "featurecounts_completed"),
    ("fastp min length changed", lambda p: _set(p, "fastp_parameters.length_required", 50), "trimming_completed"),
    ("hisat2 option changed", lambda p: _set(p, "hisat2_parameters.no_mixed", True), "alignment_completed"),
    ("non-result setting changed", lambda p: _set(p, "threads", 2), None),
    ("BAM truncated", lambda p: (lambda f: f.write_bytes(f.read_bytes()[:-100]))(
        p.path("alignment", "bam", "Ctrl1.sorted.bam")), "alignment_completed"),
    ("BAM index deleted", lambda p: p.path("alignment", "bam", "Ctrl1.sorted.bam.bai").unlink(),
     "alignment_completed"),
    ("count matrix edited", lambda p: (lambda f: f.write_text(f.read_text().replace("\t", "\t1", 1)))(
        p.path("counts", "gene_count_matrix.tsv")), "count_matrix_completed"),
    ("metadata edited after confirmation", lambda p: (lambda f: f.write_text(f.read_text() + "\n"))(
        p.path("counts", "sample_metadata.tsv")), "design_confirmed"),
    ("result table edited", lambda p: (lambda f: f.write_text(f.read_text() + "FAKE\t1\t5\t0\t1\t0\t0\tup\n"))(
        next(p.path("results", "deseq2").glob("*/upregulated.tsv"))), "deseq2_completed"),
    ("report deleted", lambda p: p.path("reports", "final_pipeline_report.html").unlink(), "report_generated"),
])
def test_resume_invalidation_matrix(project, tmp_path, name, mutate, first_invalid):
    p, _, _ = project
    status, memo = _status_after(p, tmp_path, mutate)
    invalid = [k for k in STAGE_ORDER if status[k] != "VALID"]
    if first_invalid is None:
        assert not invalid, f"{name}: unexpectedly invalid {invalid}"
        return
    assert invalid, f"{name}: nothing was invalidated"
    assert invalid[0] == first_invalid, f"{name}: first invalid stage {invalid[0]}, expected {first_invalid}"
    # cascade: the stages that consume this stage's output must be invalid too
    must_cascade = {"deseq2_completed": ["report_generated"],
                    "featurecounts_completed": ["count_matrix_completed", "design_confirmed", "deseq2_completed",
                                                "report_generated"],
                    "trimming_completed": ["alignment_completed", "featurecounts_completed", "deseq2_completed",
                                           "report_generated"],
                    "alignment_completed": ["featurecounts_completed", "count_matrix_completed", "deseq2_completed",
                                            "report_generated"],
                    "count_matrix_completed": ["design_confirmed", "deseq2_completed", "report_generated"],
                    "design_confirmed": ["deseq2_completed", "report_generated"]}.get(first_invalid, [])
    not_cascaded = [k for k in must_cascade if status[k] == "VALID"]
    assert not not_cascaded, f"{name}: downstream stages still VALID: {not_cascaded}"
    if name.endswith("changed"):  # a settings change must be reported as such, naming the setting
        assert "setting" in (memo[first_invalid][1] or ""), memo[first_invalid][1]


# ---------------------------------------------------------------- F8: report validated independently
def _ctx(project_dir):
    from rnaseq_pipeline import environment_manager as EM, system_check as SC, workflow
    from rnaseq_pipeline.project import Project
    p = Project.open(project_dir)
    cfg = p.config()
    envs = EM.Environments(cfg)
    envs.activate()
    return workflow.Context(p, cfg, envs, SC.collect(p.root))


def test_report_validates_against_result_files(project):
    from rnaseq_pipeline import report_manager
    p, _, _ = project
    assert report_manager.validate(_ctx(p)) == []


@pytest.mark.parametrize("tamper,expect", [
    (lambda t: t.replace("data-check='up:condition_Treatment_vs_Control'>", "data-check='up:condition_Treatment_vs_Control'>9", 1),
     "report shows up:"),
    (lambda t: t.replace("pca.png", "pca_missing.png", 1), "broken link"),
    (lambda t: t.replace("</html>", ""), "truncated"),
])
def test_report_tampering_detected(project, tmp_path, tamper, expect):
    from rnaseq_pipeline import report_manager
    p, _, _ = project
    dst = tmp_path / "copy"
    shutil.copytree(p, dst, symlinks=True)
    rep = dst / "reports" / "final_pipeline_report.html"
    rep.write_text(tamper(rep.read_text()))
    problems = report_manager.validate(_ctx(dst))
    assert any(expect in x for x in problems), problems


# ---------------------------------------------------------------- F15: non-interactive project health and doctor
def _copy_project(project_dir, tmp_path):
    dst = tmp_path / "copy"
    shutil.copytree(project_dir, dst, symlinks=True)
    (dst / ".lock").unlink(missing_ok=True)
    return dst


def test_validate_project_passes_on_finished_project(project, tmp_path):
    p, _, cfg = project
    r = run_cli(["--config", str(cfg), "--validate-project", str(_copy_project(p, tmp_path))], [], timeout=600)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "PROJECT HEALTH: PASS" in r.stdout
    assert FAKE_KEY not in r.stdout + r.stderr


def test_validate_project_fails_on_edited_count_matrix(project, tmp_path):
    p, _, cfg = project
    dst = _copy_project(p, tmp_path)
    m = dst / "counts" / "gene_count_matrix.tsv"
    lines = m.read_text().splitlines()
    cells = lines[1].split("\t")
    cells[1] = str(int(cells[1]) + 1)       # one count changed by one read
    lines[1] = "\t".join(cells)
    m.write_text("\n".join(lines) + "\n")
    r = run_cli(["--config", str(cfg), "--validate-project", str(dst)], [], timeout=600)
    assert r.returncode == 1, r.stdout[-3000:]
    assert "PROJECT HEALTH: FAIL" in r.stdout


def test_health_check_runs_real_mini_job(tmp_path):
    # --check runs HISAT2 -> samtools -> featureCounts (and the others) on a 2 kb genome, plus DESeq2 in R;
    # network is not required for PASS/WARNING, so the result must not be FAIL on a machine with the tools.
    env = dict(os.environ, PYTHONNOUSERSITE="1")
    r = subprocess.run([str(LAUNCHER), "--check"], capture_output=True, text=True, timeout=900, env=env,
                       cwd=tmp_path)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-2000:]
    assert "HEALTH: FAIL" not in r.stdout
    assert "200/200 read pairs assigned" in r.stdout, r.stdout[-4000:]


# ---------------------------------------------------------------- P1/H1: DESeq2 output checked independently of R
def _deseq2_copy(project_dir, tmp_path):
    """Copy the finished project. The parameter file keeps the ORIGINAL absolute paths on purpose: validation must
    re-point them to the copy (a moved project), otherwise it would check the untouched original and every tamper
    test below would fail. Only the summary is re-pointed, so the tamper helpers edit the copy."""
    dst = _copy_project(project_dir, tmp_path)
    path = dst / "results" / "deseq2" / "deseq2_summary.json"
    path.write_text(path.read_text().replace(str(project_dir), str(dst)))
    params = json.loads((dst / "results" / "deseq2" / "deseq2_params.json").read_text())
    assert str(project_dir) in params["counts"]
    return dst, params


def _read_tsv(path):
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return head, [dict(zip(head, l.split("\t"))) for l in lines[1:] if l.strip()]


def _write_tsv(path, head, rows):
    path.write_text("\t".join(head) + "\n" + "".join("\t".join(r[h] for h in head) + "\n" for r in rows))


def _flip_direction(dst):
    """What a swapped baseline looks like: every fold change negated, up and down exchanged, all tables and the
    summary rewritten consistently (the existing threshold checks alone cannot notice)."""
    sfile = dst / "results" / "deseq2" / "deseq2_summary.json"
    summary = json.loads(sfile.read_text())
    params = json.loads((dst / "results" / "deseq2" / "deseq2_params.json").read_text())
    a, t = params["alpha"], params["log2fc_threshold"]
    for c in summary["contrasts"].values():
        cdir = Path(c["dir"])
        head, rows = _read_tsv(cdir / "full_results.tsv")
        for r in rows:
            for k in ("log2FoldChange", "log2FoldChange_shrunk", "stat"):
                if k in r and r[k] not in ("NA", ""):
                    r[k] = repr(-float(r[k]))
            if r["regulation"] in ("up", "down"):
                r["regulation"] = {"up": "down", "down": "up"}[r["regulation"]]
        _write_tsv(cdir / "full_results.tsv", head, rows)
        sig = [r for r in rows if r["regulation"] != "ns"]
        up = [r for r in sig if r["regulation"] == "up"]
        down = [r for r in sig if r["regulation"] == "down"]
        for name, sel in (("significant_results", sig), ("upregulated", up), ("downregulated", down),
                          ("top_upregulated", up[:params["top_n_genes"]]),
                          ("top_downregulated", down[:params["top_n_genes"]])):
            _write_tsv(cdir / f"{name}.tsv", head, sel)
        c["up"], c["down"] = len(up), len(down)
    sfile.write_text(json.dumps(summary))
    assert a and t is not None


def _drop_one_gene(dst):
    sfile = dst / "results" / "deseq2" / "deseq2_summary.json"
    summary = json.loads(sfile.read_text())
    gone = None
    for c in summary["contrasts"].values():
        head, rows = _read_tsv(Path(c["dir"]) / "full_results.tsv")
        gone = gone or next(r["gene_id"] for r in reversed(rows) if r["regulation"] == "ns")
        _write_tsv(Path(c["dir"]) / "full_results.tsv", head, [r for r in rows if r["gene_id"] != gone])
        c["genes_tested"] -= 1
    nc = dst / "results" / "deseq2" / "normalized_counts.tsv"
    head, rows = _read_tsv(nc)
    _write_tsv(nc, head, [r for r in rows if r["gene_id"] != gone])
    summary["genes_tested"] -= 1
    sfile.write_text(json.dumps(summary))


def _edit_first_contrast(dst, fn):
    summary = json.loads((dst / "results" / "deseq2" / "deseq2_summary.json").read_text())
    path = Path(next(iter(summary["contrasts"].values()))["dir"]) / "full_results.tsv"
    head, rows = _read_tsv(path)
    fn(rows)
    _write_tsv(path, head, rows)


def _scale_first_size_factor(dst):
    f = dst / "results" / "deseq2" / "size_factors.tsv"
    head, rows = _read_tsv(f)
    rows[0]["size_factor"] = repr(float(rows[0]["size_factor"]) * 1.1)
    _write_tsv(f, head, rows)


def _double_a_basemean(rows):
    rows[-1]["baseMean"] = repr(float(rows[-1]["baseMean"]) * 2)


def _padj_below_pvalue(rows):
    r = next(r for r in reversed(rows) if r["regulation"] == "ns" and r["pvalue"] not in ("NA", "")
             and float(r["pvalue"]) > 0.2)
    r["padj"] = repr(float(r["pvalue"]) / 2)      # still not significant, so the threshold checks cannot notice


@pytest.mark.parametrize("tamper,expect", [
    (_flip_direction, "direction"),
    (_drop_one_gene, "low-count filter"),
    (_scale_first_size_factor, "size factor"),
    (lambda d: _edit_first_contrast(d, _double_a_basemean), "baseMean"),
    (lambda d: _edit_first_contrast(d, _padj_below_pvalue), "padj"),
], ids=["direction-swapped", "gene-dropped", "size-factor", "baseMean", "padj-below-pvalue"])
def test_deseq2_output_rederived_independently(project, tmp_path, tamper, expect):
    from rnaseq_pipeline import PipelineError, r_bridge
    from rnaseq_pipeline.project import Project
    p, _, _ = project
    dst, params = _deseq2_copy(p, tmp_path)
    assert r_bridge.validate_outputs(Project.open(dst), params)      # untouched copy passes
    tamper(dst)
    with pytest.raises(PipelineError) as e:
        r_bridge.validate_outputs(Project.open(dst), params)
    assert expect in str(e.value), str(e.value)
