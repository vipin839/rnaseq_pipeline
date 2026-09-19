"""Python -> R transition: parameter file, Rscript execution, output validation."""
import json
from pathlib import Path

from . import PIPELINE_ROOT, PipelineError, runner, ui

R_MAIN = PIPELINE_ROOT / "R" / "deseq2_pipeline.R"


def build_params(project, cfg, design, orgdb=None):
    counts = project.path("counts", "gene_count_matrix.tsv")
    meta = project.abs(design["metadata_file"])
    return {
        "counts": str(counts), "metadata": str(meta),
        "formula": design["formula"], "variables": design["variables"],
        "variable_of_interest": design["variable_of_interest"], "reference_level": design["reference_level"],
        "contrasts": design["contrasts"],
        "alpha": float(cfg["alpha"]), "log2fc_threshold": float(cfg["log2fc_threshold"]),
        "lfc_test_threshold": float(cfg["lfc_test_threshold"]),
        "min_count_filter": cfg["min_count_filter"],
        "independent_filtering": bool(cfg["independent_filtering"]), "cooks_cutoff": bool(cfg["cooks_cutoff"]),
        "lfc_shrinkage": cfg["lfc_shrinkage"], "transformation": cfg["transformation"],
        "top_n_genes": int(cfg["top_n_genes"]), "heatmap_max_genes": int(cfg["heatmap_max_genes"]),
        "plot_formats": list(cfg["plot_formats"]), "plot_dpi": int(cfg["plot_dpi"]),
        "annotation": cfg.get("annotation", {}), "orgdb": orgdb,
        "out_dir": str(project.path("results", "deseq2")),
        "plot_dir": str(project.path("results", "plots")),
        "table_dir": str(project.path("results", "tables")),
        "reference": project.state.get("reference", {}).get("label"),
        "strandedness": (project.state.get("strandedness") or {}).get("value"),
    }


def write_params(project, params):
    path = project.path("results", "deseq2", "deseq2_params.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(params, indent=2))
    return path


def run(rscript, params_path, log_file, validate_only=False):
    cmd = [rscript, "--vanilla", R_MAIN, params_path] + (["--validate-only"] if validate_only else [])
    res = runner.run(cmd, stage="deseq2", log_file=log_file, capture=True, check=False,
                     description="R: validating design" if validate_only else "R: DESeq2 analysis")
    for line in (res.stdout or "").splitlines():
        if line.startswith("[") and "]" in line:
            tag, _, msg = line[1:].partition("] ")
            if tag in ("OK", "INFO", "WARNING", "ERROR", "SKIPPED"):
                ui.status(tag, "R: " + msg)
    if res.returncode != 0:
        err = [l for l in (res.stdout or "").splitlines() if l.startswith("[ERROR]")]
        raise PipelineError("R/DESeq2 step failed: " + (err[-1][8:] if err else "see log"), stage="deseq2",
                            cause=runner.tail_file(log_file) if log_file else None,
                            remedy=f"see {log_file}; fix the metadata/design and re-run")
    return res


EXPECTED = ["full_results.tsv", "significant_results.tsv", "upregulated.tsv", "downregulated.tsv",
            "top_upregulated.tsv", "top_downregulated.tsv", "statistics.txt"]


def validate_outputs(project, params):
    """Check R outputs structurally and for internal consistency. Returns summary dict."""
    out_dir = Path(params["out_dir"])
    sfile = out_dir / "deseq2_summary.json"
    if not sfile.exists():
        raise PipelineError("DESeq2 summary missing — the R step did not complete", stage="deseq2")
    summary = json.loads(sfile.read_text())
    problems = []
    for f in ("normalized_counts.tsv", "size_factors.tsv", "dds.rds", "R_sessionInfo.txt"):
        if not (out_dir / f).exists():
            problems.append(f"missing {f}")
    with open(out_dir / "normalized_counts.tsv") as fh:
        n_norm = sum(1 for _ in fh) - 1
    if n_norm != summary["genes_tested"]:
        problems.append(f"normalized counts rows {n_norm} != genes tested {summary['genes_tested']}")
    if len(summary["contrasts"]) != len(params["contrasts"]):
        problems.append("number of contrast results differs from requested contrasts")
    for name, c in summary["contrasts"].items():
        cdir = Path(c["dir"])
        for f in EXPECTED:
            if not (cdir / f).exists():
                problems.append(f"{name}: missing {f}")
        if problems:
            continue
        full = _rows(cdir / "full_results.tsv")
        up, down, sig = _rows(cdir / "upregulated.tsv"), _rows(cdir / "downregulated.tsv"), \
            _rows(cdir / "significant_results.tsv")
        if len(full) != c["genes_tested"]:
            problems.append(f"{name}: full_results has {len(full)} rows, expected {c['genes_tested']}")
        if len(up) != c["up"] or len(down) != c["down"] or len(sig) != c["significant"]:
            problems.append(f"{name}: up/down/significant table sizes inconsistent with summary")
        a, t = params["alpha"], params["log2fc_threshold"]
        for r in up:
            if not (float(r["padj"]) < a and float(r["log2FoldChange"]) >= t):
                problems.append(f"{name}: upregulated gene {r['gene_id']} violates thresholds")
                break
        for r in down:
            if not (float(r["padj"]) < a and float(r["log2FoldChange"]) <= -t):
                problems.append(f"{name}: downregulated gene {r['gene_id']} violates thresholds")
                break
        need = {"gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"}
        if full and not need <= set(full[0]):
            problems.append(f"{name}: results lack columns {need - set(full[0])}")
        for kind, files in (c.get("plots") or {}).items():
            for f in (files if isinstance(files, list) else [files]):
                if not Path(f).exists() or Path(f).stat().st_size == 0:
                    problems.append(f"{name}: plot missing {f}")
    for kind, files in (summary.get("qc_plots") or {}).items():
        for f in (files if isinstance(files, list) else [files]):
            if not Path(f).exists():
                problems.append(f"QC plot missing: {f}")
    if problems:
        raise PipelineError("DESeq2 output validation failed:\n  - " + "\n  - ".join(problems), stage="deseq2")
    return summary


def _rows(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f if l.strip()]
