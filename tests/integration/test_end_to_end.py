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


def _edit_json(path, **changes):
    d = json.loads(path.read_text())
    d.update(changes)
    path.write_text(json.dumps(d, indent=2))


def _edit_state_strand(p, value):
    s = json.loads(p.path("project.json").read_text())
    s["strandedness"]["value"] = value
    p.path("project.json").write_text(json.dumps(s, indent=2))


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
    # P1/H4
    ("strandedness setting changed", lambda p: _set(p, "strandedness", "forward"), "strandedness_determined"),
    ("strandedness decision edited", lambda p: _edit_json(p.path("alignment", "reports", "strandedness.json"),
                                                          value="forward", featurecounts_flag="1"),
     "strandedness_determined"),
    ("strandedness in project state edited", lambda p: _edit_state_strand(p, "forward"), "strandedness_determined"),
    ("trimming decision setting changed", lambda p: _set(p, "trim_adapters", "never"), "quality_assessed"),
    ("trimming decision edited", lambda p: _edit_json(p.path("qc", "assessment", "qc_decision.json"),
                                                      decision="skip"), "quality_assessed"),
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
                    "design_confirmed": ["deseq2_completed", "report_generated"],
                    "strandedness_determined": ["stringtie_completed", "featurecounts_completed",
                                                "count_matrix_completed", "deseq2_completed", "report_generated"],
                    "quality_assessed": ["trimming_completed", "alignment_completed", "featurecounts_completed",
                                         "deseq2_completed", "report_generated"]}.get(first_invalid, [])
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


# ---------------------------------------------------------------- P1/H7: tools not run from a conda environment
def test_report_without_conda_environment(project, tmp_path, monkeypatch):
    """Tools on PATH (modules, containers, system packages): there is nothing to export, and the finished analysis
    must still get a valid report that says so — not fail at the last step over a 'missing' file."""
    from rnaseq_pipeline import environment_manager, manifest, report_manager
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    for f in ("environment.yml", "package_versions.txt"):
        (dst / "pipeline_manifest" / f).unlink()
    monkeypatch.setattr(environment_manager.Environments, "conda_bin", lambda self: None)
    ctx = _ctx(dst)
    manifest.write(ctx)
    rep = report_manager.generate(ctx)
    assert report_manager.validate(ctx, rep) == []
    body = rep.read_text()
    assert "not produced" in body and "software_versions.tsv" in body
    assert not (dst / "pipeline_manifest" / "environment.yml").exists()


def test_missing_conda_export_still_detected(project, tmp_path):
    """The other side: conda WAS used, and its export later disappears -> the report must not validate."""
    from rnaseq_pipeline import manifest, report_manager
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    manifest.write(ctx)
    rep = report_manager.generate(ctx)
    assert report_manager.validate(ctx, rep) == []
    (dst / "pipeline_manifest" / "environment.yml").unlink()
    rep = report_manager.generate(ctx)
    assert any("missing" in x for x in report_manager.validate(ctx, rep))


# ---------------------------------------------------------------- P1/H2: plot data reconciled with the results
def _contrast_dir(dst):
    summary = json.loads((dst / "results" / "deseq2" / "deseq2_summary.json").read_text())
    return Path(next(iter(summary["contrasts"].values()))["dir"])


def _edit_tsv(path, fn):
    head, rows = _read_tsv(path)
    rows = fn(rows) or rows
    _write_tsv(path, head, rows)


def _volcano_hides_an_up_gene(dst):
    def fn(rows):
        next(r for r in rows if r["regulation"] == "up")["regulation"] = "ns"
    _edit_tsv(_contrast_dir(dst) / "plot_data" / "volcano_data.tsv", fn)


def _ma_drops_a_gene(dst):
    _edit_tsv(_contrast_dir(dst) / "plot_data" / "ma_plot_data.tsv", lambda rows: rows[:-1])


def _library_size_wrong(dst):
    def fn(rows):
        rows[0]["library_size"] = str(int(float(rows[0]["library_size"])) + 1000)
    _edit_tsv(dst / "results" / "tables" / "library_sizes.tsv", fn)


def _heatmap_shows_wrong_gene(dst):
    full = _read_tsv(_contrast_dir(dst) / "full_results.tsv")[1]
    outsider = next(r["gene_id"] for r in full if r["regulation"] == "ns")
    def fn(rows):
        rows[0]["gene_id"] = outsider
    _edit_tsv(_contrast_dir(dst) / "plot_data" / "significant_heatmap_zscores.tsv", fn)


def _correlation_not_symmetric(dst):
    def fn(rows):
        k = [c for c in rows[0] if c != "sample"][1]
        rows[0][k] = "0.1"
    _edit_tsv(dst / "results" / "tables" / "sample_correlation.tsv", fn)


@pytest.mark.parametrize("tamper,expect", [
    (_volcano_hides_an_up_gene, "volcano"),
    (_ma_drops_a_gene, "MA"),
    (_library_size_wrong, "library size"),
    (_heatmap_shows_wrong_gene, "heatmap"),
    (_correlation_not_symmetric, "correlation"),
], ids=["volcano-regulation", "ma-gene-missing", "library-size", "heatmap-genes", "correlation-matrix"])
def test_plot_data_reconciled_with_results(project, tmp_path, tamper, expect):
    from rnaseq_pipeline import PipelineError, r_bridge
    from rnaseq_pipeline.project import Project
    p, _, _ = project
    dst, params = _deseq2_copy(p, tmp_path)
    assert r_bridge.validate_outputs(Project.open(dst), params)
    tamper(dst)
    with pytest.raises(PipelineError) as e:
        r_bridge.validate_outputs(Project.open(dst), params)
    assert expect in str(e.value), str(e.value)


# ---------------------------------------------------------------- P1/H3: every per-sample report number re-checked
@pytest.mark.parametrize("key", ["reads:Ctrl1", "input:Ctrl1", "rate:Ctrl1", "primary_mapped:Ctrl1",
                                 "mapped_pct:Ctrl1", "assigned:Ctrl1", "counted:Ctrl1", "version:hisat2",
                                 "version:R"])
def test_report_per_sample_numbers_checked(project, tmp_path, key):
    """Each number is re-derived from the file the TOOL wrote (FASTQ validation report, HISAT2 summary, samtools
    flagstat, featureCounts summary, count matrix, manifest record) — not from the table the report was built from."""
    import re
    from rnaseq_pipeline import report_manager
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    rep = dst / "reports" / "final_pipeline_report.html"
    text = rep.read_text()
    m = re.search(rf"data-check='{re.escape(key)}'>([^<]*)<", text)
    assert m, f"report does not mark {key}"
    assert report_manager.validate(_ctx(dst)) == []
    rep.write_text(text.replace(m.group(0), f"data-check='{key}'>{m.group(1)}9<", 1))
    problems = report_manager.validate(_ctx(dst))
    assert any(key in x for x in problems), problems


def test_report_validation_after_project_moved(project, tmp_path):
    """Recorded absolute paths must not tie the report check to the project's old location."""
    from rnaseq_pipeline import report_manager
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    s = dst / "results" / "deseq2" / "deseq2_summary.json"
    s.write_text(s.read_text().replace(str(p), "/nonexistent/old/location"))
    assert report_manager.validate(_ctx(dst)) == []


def test_reference_file_change_invalidates_downstream(project, tmp_path):
    """P1/H4: the annotation in the shared reference store changes after the analysis (edited or replaced).
    The change is undone (bytes and time) before the test ends, so no other test sees it."""
    import os
    p, _, _ = project
    cp = json.loads((p / "checkpoints" / "reference_ready.json").read_text())
    gtf = Path(next(o["path"] for o in cp["outputs"] if o["path"].endswith("annotation.gtf")))
    original, st = gtf.read_bytes(), gtf.stat()
    try:
        status, memo = _status_after(p, tmp_path, lambda _: gtf.write_bytes(original + b"# edited\n"))
    finally:
        gtf.write_bytes(original)
        os.utime(gtf, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert gtf.read_bytes() == original
    invalid = [k for k in STAGE_ORDER if status[k] != "VALID"]
    assert invalid and invalid[0] == "reference_ready", (invalid, memo.get("reference_ready"))
    for k in ("alignment_completed", "featurecounts_completed", "deseq2_completed", "report_generated"):
        assert status[k] != "VALID", k


# ---------------------------------------------------------------- P1/H4: killed during alignment, then resumed
def _drive_cli(args, answer, env, timeout=900, kill_when=None):
    """Run the real CLI and answer each prompt by what it says (answer(text_since_last_prompt) -> reply).
    kill_when(): if given and true, send SIGTERM to the CLI (what closing the terminal / `kill` does)."""
    import signal
    import threading
    import time
    proc = subprocess.Popen([str(LAUNCHER), *args], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, env=env)
    out, pending = [], []

    def reader():
        while True:
            ch = proc.stdout.read(1)
            if not ch:
                break
            out.append(ch)
            pending.append(ch)
    t = threading.Thread(target=reader, daemon=True)
    t.start()
    deadline, killed = time.time() + timeout, False
    while proc.poll() is None and time.time() < deadline:
        time.sleep(0.05)
        if kill_when and not killed and kill_when():
            proc.send_signal(signal.SIGTERM)
            killed = True
        text = b"".join(pending).decode(errors="replace")
        if text.endswith((": ", "]: ")) and not killed:
            reply = answer(text)
            pending.clear()
            try:
                proc.stdin.write((reply + "\n").encode())
                proc.stdin.flush()
            except BrokenPipeError:
                break
    if proc.poll() is None:
        proc.kill()
    proc.wait(timeout=30)
    t.join(timeout=5)
    return proc.returncode, b"".join(out).decode(errors="replace")


def test_killed_during_alignment_then_resumed(project, tmp_path):
    import shutil as sh
    from rnaseq_pipeline import config as C
    p, _, cfg = project
    dst = _copy_project(p, tmp_path)
    # a tools environment whose hisat2 waits while a flag file exists, so the kill lands mid-alignment
    root = tmp_path / "conda"
    tools = Path(sh.which("hisat2")).parent
    (root / "envs" / "rnaseq-tools" / "conda-meta").mkdir(parents=True)
    bin_dir = root / "envs" / "rnaseq-tools" / "bin"
    bin_dir.mkdir()
    for f in tools.iterdir():
        (bin_dir / f.name).symlink_to(f)
    (bin_dir / "hisat2").unlink()
    slow, started = tmp_path / "slow", tmp_path / "started"
    (bin_dir / "hisat2").write_text(f'#!/bin/sh\ncase "$*" in *--version*) exec "{tools / "hisat2"}" "$@";; esac\n'
                                    f'if [ -e "{slow}" ]; then touch "{started}"; sleep 31.73; fi\n'
                                    f'exec "{tools / "hisat2"}" "$@"\n')
    (bin_dir / "hisat2").chmod(0o755)
    (root / "envs" / "rnaseq-r").symlink_to(RSCRIPT.parent.parent)
    pc = C.load_yaml(dst / "config" / "project_config.yaml")
    C.set_value(pc, "environment.conda_root", str(root))
    C.save_yaml(pc, dst / "config" / "project_config.yaml")
    for f in ("Ctrl1.sorted.bam", "Ctrl1.sorted.bam.bai"):
        (dst / "alignment" / "bam" / f).unlink()
    env = dict(os.environ, PYTHONNOUSERSITE="1")
    args = ["--config", str(cfg), "--project", str(dst)]
    seen = {"menu": 0}

    def answer(text):
        if "PROJECT:" in text and text.endswith("Select: "):
            seen["menu"] += 1
            return "1" if seen["menu"] == 1 else "10"
        return "y" if text.endswith(("[Y/n]: ", "[y/N]: ")) else "1"

    slow.touch()
    code, out = _drive_cli(args, answer, env, timeout=300, kill_when=started.exists)
    assert code == 143, out[-3000:]
    assert "stopped by SIGTERM during ALIGNMENT" in out
    assert not (dst / "alignment" / "bam" / "Ctrl1.sorted.bam").exists()      # nothing partial accepted
    assert subprocess.run(["pgrep", "-f", r"^sleep 31\.73$"], capture_output=True).returncode != 0
    status, _ = _status_after(dst, tmp_path / "s1", lambda _: None)
    assert status["alignment_completed"] != "VALID"

    slow.unlink()
    seen["menu"] = 0
    before = len((dst / "logs" / "commands.jsonl").read_text().splitlines())
    code, out = _drive_cli(args, answer, env, timeout=900)
    assert code == 0, out[-3000:]
    status, memo = _status_after(dst, tmp_path / "s2", lambda _: None)
    assert all(v == "VALID" for v in status.values()), {k: memo.get(k) for k, v in status.items() if v != "VALID"}
    new = (dst / "logs" / "commands.jsonl").read_text().splitlines()[before:]
    aligned = [json.loads(l)["sample"] for l in new if '"hisat2' in l]
    assert aligned == ["Ctrl1"], aligned                                        # the other samples were reused
    assert (dst / "counts" / "gene_count_matrix.tsv").read_bytes() == \
        (p / "counts" / "gene_count_matrix.tsv").read_bytes()                   # same result as an uninterrupted run


# ---------------------------------------------------------------- P1/H5: failure injection on a finished project
def _fake_tool(tmp_path, name, body):
    d = tmp_path / "fakebin"
    d.mkdir(exist_ok=True)
    f = d / name
    f.write_text("#!/bin/sh\n" + body + "\n")
    f.chmod(0o755)
    return f


def _stage(key):
    from rnaseq_pipeline import workflow
    return workflow.STAGES.index(workflow.BY_KEY[key]) + 1, workflow.BY_KEY[key]


def test_outdated_tool_refused_before_the_stage_starts(project, tmp_path, monkeypatch):
    import shutil as sh
    from rnaseq_pipeline import PipelineError, runner, workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    real = sh.which("featureCounts", path=runner.child_env()["PATH"])
    fake = _fake_tool(tmp_path, "featureCounts",
                      f'if [ "$1" = "-v" ]; then echo "featureCounts v1.5.0"; exit 0; fi\nexec "{real}" "$@"')
    monkeypatch.setitem(runner._tool_env, "path_prefix", [str(fake.parent)] + runner._tool_env["path_prefix"])
    before = (dst / "logs" / "commands.jsonl").read_text()
    with pytest.raises(PipelineError) as e:
        workflow.run_stage(ctx, *_stage("featurecounts_completed"))
    assert "featureCounts" in str(e.value) and "1.5.0" in str(e.value) and "2.0.0" in str(e.value), str(e.value)
    assert (dst / "logs" / "commands.jsonl").read_text() == before           # nothing was run


def test_permission_denied_is_explained(project, tmp_path):
    from rnaseq_pipeline import PipelineError, workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    counts = dst / "counts"                   # read-only like a project owned by another user / read-only copy
    files = [f for f in counts.iterdir() if f.is_file()]
    for f in files:
        f.chmod(0o444)
    counts.chmod(0o555)
    try:
        with pytest.raises(PipelineError) as e:
            workflow.run_stage(ctx, *_stage("count_matrix_completed"))
    finally:
        counts.chmod(0o755)
        for f in files:
            f.chmod(0o644)
    msg = str(e.value) + str(e.value.cause) + str(e.value.remedy)
    assert "permission" in msg.lower() and "counts" in msg, msg


def test_low_memory_asks_before_alignment(project, tmp_path, monkeypatch):
    from rnaseq_pipeline import UserAbort, ui, workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    ctx.sysinfo["ram_available_gb"] = 0.1
    asked = []
    monkeypatch.setattr(ui, "ask_yes_no", lambda q, default=None: asked.append(q) or False)
    with pytest.raises(UserAbort) as e:
        workflow.run_stage(ctx, *_stage("alignment_completed"))
    assert asked and "memory" in str(e.value)
    assert (dst / "alignment" / "bam" / "Ctrl1.sorted.bam").exists()        # existing results untouched


@pytest.mark.parametrize("mode,expect,remedy", [
    ("killed", "SIGKILL", "memory"),
    ("missing-package", "no package called", "--check"),
])
def test_r_failure_is_explained(project, tmp_path, monkeypatch, mode, expect, remedy):
    from rnaseq_pipeline import PipelineError, ui, workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    full = next((dst / "results" / "deseq2").glob("*/full_results.tsv"))
    body = {"killed": f'echo "[INFO] fitting"\nhead -c 200 "{full}" > "{full}.tmp" && mv "{full}.tmp" "{full}"\n'
                      'kill -KILL $$',
            "missing-package": "echo \"[ERROR] there is no package called ‘DESeq2’\"\nexit 1"}[mode]
    fake = _fake_tool(tmp_path, "Rscript", body)
    monkeypatch.setattr(type(ctx.envs), "rscript", lambda self: fake)
    monkeypatch.setattr(ui, "ask_yes_no", lambda q, default=None: True)
    with pytest.raises(PipelineError) as e:
        workflow.run_stage(ctx, *_stage("deseq2_completed"))
    assert expect in str(e.value), str(e.value)
    assert remedy in (e.value.remedy or ""), e.value.remedy
    status, _ = _status_after(dst, tmp_path / "s", lambda _: None)
    if mode == "killed":                          # a half-written result is never accepted
        assert status["deseq2_completed"] != "VALID" and status["report_generated"] != "VALID"


def test_count_matrix_rebuilds_in_a_moved_project(project, tmp_path):
    """P1/H9: featureCounts records absolute BAM paths; after moving/copying the project, re-running the count-matrix
    stage (e.g. after a technical-replicate decision) must still match the columns to the samples."""
    from rnaseq_pipeline import workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    assert str(p) in (dst / "featurecounts" / "featurecounts.txt").read_text().splitlines()[1]
    workflow.run_stage(_ctx(dst), *_stage("count_matrix_completed"))
    assert (dst / "counts" / "gene_count_matrix.tsv").read_bytes() == (p / "counts" / "gene_count_matrix.tsv").read_bytes()


def test_read_only_project_is_explained_not_a_traceback(project, tmp_path):
    p, _, cfg = project
    dst = _copy_project(p, tmp_path)
    for d in [dst, *[x for x in dst.rglob("*") if x.is_dir() and not x.is_symlink()]]:
        d.chmod(0o555)
    try:
        r = run_cli(["--config", str(cfg), "--project", str(dst)], ["10"], timeout=300)
        check = run_cli(["--config", str(cfg), "--validate-project", str(dst)], [], timeout=300)
    finally:
        for d in [dst, *[x for x in dst.rglob("*") if x.is_dir() and not x.is_symlink()]]:
            d.chmod(0o755)
    out = r.stdout + r.stderr
    assert "Traceback" not in out and "read-only" in out and "--validate-project" in out, out[-2000:]
    assert check.returncode == 0 and "PROJECT HEALTH: PASS" in check.stdout     # checking needs no writes


# ---------------------------------------------------------------- P1/M2: QC reports complete and for the right files
def test_fastqc_report_must_belong_to_its_file(project, tmp_path):
    """All synthetic samples have the same read count, so a report of another file passes a count check."""
    import shutil as sh
    from rnaseq_pipeline import PipelineError, qc_manager
    p, _, _ = project
    d = tmp_path / "fastqc"
    sh.copytree(p / "qc" / "fastqc_raw", d)
    fq = {s: p / "data" / "fastq" / f"{s}_R1.fastq.gz" for s in ("Ctrl1", "Ctrl2")}
    expected = {str(f): 40000 for f in fq.values()}
    qc_manager.validate_fastqc(list(fq.values()), d, expected, "raw_qc_completed")          # correct: passes
    for ext in (".zip", ".html"):
        sh.copy(d / f"Ctrl2_R1_fastqc{ext}", d / f"Ctrl1_R1_fastqc{ext}")                 # swapped report
    with pytest.raises(PipelineError) as e:
        qc_manager.validate_fastqc(list(fq.values()), d, expected, "raw_qc_completed")
    assert "Ctrl1_R1" in str(e.value) and "Ctrl2_R1" in str(e.value), str(e.value)


def test_multiqc_must_include_every_input(project, tmp_path):
    """MultiQC skips files it cannot parse (and overwrites duplicate sample names) and still exits 0."""
    import shutil as sh
    from rnaseq_pipeline import PipelineError, qc_manager
    p, _, _ = project
    d = tmp_path / "fastqc"
    sh.copytree(p / "qc" / "fastqc_raw", d)
    zips = sorted(d.glob("*_fastqc.zip"))
    qc_manager.run_multiqc([d], tmp_path / "mq_ok", "ok", tmp_path / "ok.log", "raw_qc_completed",
                           expected_sources=zips)
    (d / "Treat2_R1_fastqc.zip").write_bytes((d / "Treat2_R1_fastqc.zip").read_bytes()[:500])   # unreadable
    with pytest.raises(PipelineError) as e:
        qc_manager.run_multiqc([d], tmp_path / "mq", "t", tmp_path / "mq.log", "raw_qc_completed",
                               expected_sources=zips)
    assert "Treat2_R1_fastqc.zip" in str(e.value), str(e.value)


# ---------------------------------------------------------------- P1/M3: quality gate "Review" option
def test_quality_gate_review_shows_evidence_then_returns_to_decision(project, tmp_path, monkeypatch, capsys):
    from rnaseq_pipeline import ui, workflow
    p, _, _ = project
    dst = _copy_project(p, tmp_path)
    ctx = _ctx(dst)
    offered, picks = [], iter(["Review", "Run fastp"])

    def choose(title, options, default=None):
        offered.append(list(options))
        want = next(picks)
        return next(i for i, o in enumerate(options) if o.startswith(want))
    monkeypatch.setattr(ui, "choose", choose)
    workflow.run_stage(ctx, *_stage("quality_assessed"))
    out = capsys.readouterr().out
    assert len(offered) == 2 and offered[0] == offered[1]                     # back to the same decision
    assert [o.split(" (")[0] for o in offered[0]] == ["Accept recommendation", "Run fastp", "Skip trimming",
                                                      "Review the evidence per sample", "Stop pipeline"]
    assert "REVIEW — evidence per sample" in out and "Ctrl1" in out
    assert "adapter" in out.lower() and "normal for RNA-seq" in out           # technical vs biological explained
    assert "multiqc_report.html" in out
    assert json.loads((dst / "qc" / "assessment" / "qc_decision.json").read_text())["decision"] == "trim"


def test_aligner_status_shown_before_alignment(project):
    p, out, _ = project
    block = out[out.index("ALIGNER"):out.index("ALIGNER") + 600]
    assert "HISAT2" in block and "SUPPORTED AND VALIDATED" in block
    assert "STAR" in block and "NOT AVAILABLE" in block
    assert json.loads((p / "checkpoints" / "alignment_completed.json").read_text())["params"]["settings"]["aligner"] \
        == "hisat2"


def test_setting_added_in_a_later_version_does_not_invalidate(project, tmp_path):
    """A checkpoint written before a setting existed was made with that setting's default behaviour (e.g. older
    projects aligned with HISAT2 before 'aligner' existed); upgrading must not force recomputation."""
    def drop_aligner(p):
        f = p.path("checkpoints", "alignment_completed.json")
        d = json.loads(f.read_text())
        del d["params"]["settings"]["aligner"]
        f.write_text(json.dumps(d))
    p, _, _ = project
    status, memo = _status_after(p, tmp_path, drop_aligner)
    assert status["alignment_completed"] == "VALID", memo.get("alignment_completed")
    # but a recorded value that differs from the current one is still a change
    status, _ = _status_after(p, tmp_path / "b", lambda q: _set(q, "hisat2_parameters.no_mixed", True))
    assert status["alignment_completed"] != "VALID"


# ---------------------------------------------------------------- P1/M7: per-area project health
AREAS = ["SYSTEM", "INSTALLATION", "CONFIGURATION", "DATA", "REFERENCE", "QC", "ALIGNMENT", "QUANTIFICATION",
         "DESIGN", "DESEQ2 + RESULTS", "REPORT", "SECURITY", "RESUME", "SCIENTIFIC VALIDATION"]


def _area_summary(out):
    import re
    block = out[out.index("AREA SUMMARY"):]
    return {a: re.search(rf"^\s*{re.escape(a)}\s+(PASS|WARNING|FAIL)\b", block, re.M).group(1) for a in AREAS}


def test_project_health_reports_every_area(project, tmp_path):
    p, _, cfg = project
    r = run_cli(["--config", str(cfg), "--validate-project", str(_copy_project(p, tmp_path))], [], timeout=600)
    assert r.returncode == 0, r.stdout[-3000:]
    assert all(v == "PASS" for v in _area_summary(r.stdout).values()), _area_summary(r.stdout)
    assert "BLOCKING ISSUES: none" in r.stdout


def test_project_health_blocking_issue_named(project, tmp_path):
    p, _, cfg = project
    dst = _copy_project(p, tmp_path)
    m = dst / "counts" / "gene_count_matrix.tsv"
    m.write_text(m.read_text().replace("\t", "\t1", 1))
    r = run_cli(["--config", str(cfg), "--validate-project", str(dst)], [], timeout=600)
    s = _area_summary(r.stdout)
    assert r.returncode == 1 and s["QUANTIFICATION"] == "FAIL" and s["SCIENTIFIC VALIDATION"] == "FAIL", s
    assert s["DATA"] == "PASS" and s["ALIGNMENT"] == "PASS"                  # upstream areas unaffected
    blocking = r.stdout[r.stdout.index("BLOCKING ISSUES"):]
    assert "GENE COUNT MATRIX" in blocking.split("NON-BLOCKING")[0], blocking[:800]
