"""Python -> R transition: parameter file, Rscript execution, output validation."""
import json
from pathlib import Path

from . import PACKAGE_DIR, PipelineError, runner, ui

R_MAIN = PACKAGE_DIR / "R" / "deseq2_pipeline.R"


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


def relocate(params, project):
    """The parameter file records absolute paths; re-point them to where the project is now (moved or copied)."""
    old_root = Path(params["out_dir"]).parents[1]
    new = dict(params)
    for k in ("counts", "metadata", "out_dir", "plot_dir", "table_dir"):
        if params.get(k):
            try:
                new[k] = str(project.root / Path(params[k]).relative_to(old_root))
            except ValueError:
                pass  # outside the project: keep as recorded
    return new


def validate_outputs(project, params):
    """Check R outputs structurally and for internal consistency. Returns summary dict."""
    params = relocate(params, project)
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
        cdir = out_dir / name
        c["dir"] = str(cdir)
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
        problems += [f"{name}: {x}" for x in
                     check_threshold_sets(full, up, down, sig, params["alpha"], params["log2fc_threshold"])]
        need = {"gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"}
        if full and not need <= set(full[0]):
            problems.append(f"{name}: results lack columns {need - set(full[0])}")
        for files in (c.get("plots") or {}).values():
            for f in (files if isinstance(files, list) else [files]):
                if not Path(f).exists() or Path(f).stat().st_size == 0:
                    problems.append(f"{name}: plot missing {f}")
    for files in (summary.get("qc_plots") or {}).values():
        for f in (files if isinstance(files, list) else [files]):
            if not Path(f).exists():
                problems.append(f"QC plot missing: {f}")
    if not problems:
        problems += independent_checks(params, summary)
    if not problems:
        problems += plot_data_checks(params, summary)
    if problems:
        raise PipelineError("DESeq2 output validation failed:\n  - " + "\n  - ".join(problems), stage="deseq2",
                            cause="the R results are inconsistent with the count matrix and design they were "
                                  "computed from (R/DESeq2 malfunction, or result files edited afterwards)",
                            remedy="re-run the DESeq2 stage; if it repeats, report it with logs/deseq2/")
    return summary


# ---------------------------------------------------------------- independent re-derivation (no trust in R)
SIZE_FACTOR_RTOL = 1e-6
DIRECTION_MIN_AGREEMENT = 0.9   # a correct contrast agrees on ~100% of strong genes; a swapped one on ~0%
DIRECTION_MIN_GENES = 5


def _matrix(path):
    with open(path) as f:
        samples = f.readline().rstrip("\n").split("\t")[1:]
        rows = {}
        for line in f:
            if line.strip():
                parts = line.rstrip("\n").split("\t")
                rows[parts[0]] = [float(x) for x in parts[1:]]
    return samples, rows


def _median(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def independent_checks(params, summary):
    """Recompute from the count matrix and the metadata what DESeq2 must have produced, and compare:

    * the genes tested = the genes passing the documented low-count filter (none dropped, none added),
    * size factors = median-of-ratios on those genes; normalized counts = counts / size factor,
    * baseMean = mean of the normalized counts; p-values and padj within [0, 1] and padj >= pvalue,
    * the direction of every contrast: log2FC of strong genes has the sign of log2(mean numerator group /
      mean denominator group) of the normalized counts (a swapped baseline flips it).
    """
    import math
    problems = []
    out_dir = Path(params["out_dir"])
    samples, counts = _matrix(params["counts"])
    meta = {r["sample"]: r for r in _rows(params["metadata"])}
    voi = params["variable_of_interest"]
    groups = {}
    for s in samples:
        groups.setdefault(meta[s][voi], []).append(s)
    mcf = params["min_count_filter"]
    min_samples = (min(len(v) for v in groups.values()) if mcf["min_samples"] == "smallest_group"
                   else int(mcf["min_samples"]))
    kept = {g for g, v in counts.items() if sum(x >= mcf["min_count"] for x in v) >= min_samples}

    nsamples, norm = _matrix(out_dir / "normalized_counts.tsv")
    if nsamples != samples:
        return [f"normalized counts samples {nsamples} differ from the count matrix {samples}"]
    for label, genes in [("normalized counts", set(norm))] + [
            (f"{n}: full_results", {r["gene_id"] for r in _rows(Path(c["dir"]) / "full_results.tsv")})
            for n, c in summary["contrasts"].items()]:
        if genes != kept:
            problems.append(f"{label} does not match the low-count filter (>= {mcf['min_count']} reads in >= "
                            f"{min_samples} samples): {len(kept - genes)} gene(s) missing, "
                            f"{len(genes - kept)} unexpected (e.g. {sorted((kept - genes) | (genes - kept))[:3]})")
    if problems:
        return problems

    # size factors: median of ratios over genes without zeros (DESeq2's default estimator)
    logs = {g: [math.log(x) for x in counts[g]] for g in kept if all(x > 0 for x in counts[g])}
    sf = {r["sample"]: float(r["size_factor"]) for r in _rows(out_dir / "size_factors.tsv")}
    if not logs:
        problems.append("no gene without zero counts: size factors cannot be verified")
    else:
        for j, s in enumerate(samples):
            expected = math.exp(_median([v[j] - sum(v) / len(v) for v in logs.values()]))
            if s not in sf or abs(sf[s] - expected) > SIZE_FACTOR_RTOL * expected:
                problems.append(f"size factor of {s} is {sf.get(s)}, median-of-ratios gives {expected:.6g}")
    if problems:
        return problems
    worst = max(abs(norm[g][j] - counts[g][j] / sf[s]) for g in kept for j, s in enumerate(samples))
    if worst > 1e-3:
        problems.append(f"normalized counts differ from counts / size factor (max difference {worst:.3g})")

    for name, c in summary["contrasts"].items():
        full = _rows(Path(c["dir"]) / "full_results.tsv")
        ids = [r["gene_id"] for r in full]
        if len(ids) != len(set(ids)):
            problems.append(f"{name}: gene IDs repeated in full_results")
        bad_bm, bad_p = [], []
        for r in full:
            g, bm = r["gene_id"], _num(r.get("baseMean"))
            mean = sum(norm[g]) / len(norm[g])
            if bm is None or abs(bm - mean) > 1e-3 + 1e-6 * mean:
                bad_bm.append(g)
            pv, pa = _num(r.get("pvalue")), _num(r.get("padj"))
            if any(x is not None and not 0 <= x <= 1 for x in (pv, pa)) or (
                    pv is not None and pa is not None and pa < pv * (1 - 1e-9)):
                bad_p.append(g)
        if bad_bm:
            problems.append(f"{name}: baseMean is not the mean of the normalized counts for {len(bad_bm)} "
                            f"gene(s) (e.g. {bad_bm[:3]})")
        if bad_p:
            problems.append(f"{name}: pvalue/padj outside [0, 1] or padj < pvalue for {len(bad_p)} gene(s) "
                            f"(e.g. {bad_p[:3]})")
        num, den = _contrast_levels(params, name)
        if num is None:
            problems.append(f"{name}: not one of the requested contrasts")
            continue
        agree = total = 0
        for r in full:
            lfc, pa = _num(r.get("log2FoldChange")), _num(r.get("padj"))
            if lfc is None or pa is None or pa >= params["alpha"] or abs(lfc) < 1:
                continue
            g = r["gene_id"]
            m_num = sum(norm[g][samples.index(s)] for s in groups[num]) / len(groups[num])
            m_den = sum(norm[g][samples.index(s)] for s in groups[den]) / len(groups[den])
            total += 1
            agree += (m_num > m_den) == (lfc > 0)
        if total >= DIRECTION_MIN_GENES and agree / total < DIRECTION_MIN_AGREEMENT:
            problems.append(f"{name}: direction check failed — only {agree} of {total} significant genes have a "
                            f"log2FC sign matching mean({num}) vs mean({den}) of the normalized counts "
                            "(numerator and baseline appear swapped)")
        c["direction_check"] = {"genes": total, "agreeing": agree}
    return problems


def _close(a, b, rel=1e-9):
    a, b = _num(a), _num(b)
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


def _zscored_rows(path, samples, genes, label):
    """A heatmap's data: exactly `genes` in order, one column per sample, every row z-scored (mean 0, sd 1, or a
    constant row set to 0). Values are rounded to 4 decimals by R, hence the tolerance."""
    rows = _rows(path)
    head = list(rows[0]) if rows else []
    problems = []
    if [r["gene_id"] for r in rows] != genes:
        problems.append(f"{label}: shows {len(rows)} gene(s) but should show the {len(genes)} gene(s) "
                        f"{genes[:3]}{'…' if len(genes) > 3 else ''} in this order")
    if head and head[1:] != samples:
        problems.append(f"{label}: columns {head[1:]} are not the samples {samples}")
    for r in rows:
        v = [_num(r.get(s)) for s in samples]
        if any(x is None for x in v):
            problems.append(f"{label}: non-numeric value for {r['gene_id']}")
            break
        mean = sum(v) / len(v)
        sd = (sum((x - mean) ** 2 for x in v) / (len(v) - 1)) ** 0.5 if len(v) > 1 else 0.0
        if abs(mean) > 1e-3 or not (abs(sd - 1) < 1e-3 or all(x == 0 for x in v)):
            problems.append(f"{label}: row {r['gene_id']} is not z-scored (mean {mean:.4f}, sd {sd:.4f})")
            break
    return problems


def _square_matrix(path, samples, label, diagonal, symmetric_tol, lo=None, hi=None):
    rows = _rows(path)
    if [r["sample"] for r in rows] != samples or (rows and [k for k in rows[0] if k != "sample"] != samples):
        return [f"{label}: rows/columns are not the samples {samples}"]
    for i, r in enumerate(rows):
        for j, s in enumerate(samples):
            v, w = _num(r[s]), _num(rows[j][samples[i]])
            if v is None or w is None or abs(v - w) > symmetric_tol:
                return [f"{label}: matrix is not symmetric ({samples[i]} vs {s}: {r[s]} / {rows[j][samples[i]]})"]
            if i == j and abs(v - diagonal) > symmetric_tol:
                return [f"{label}: diagonal value for {s} is {v}, expected {diagonal}"]
            if (lo is not None and v < lo - symmetric_tol) or (hi is not None and v > hi + symmetric_tol):
                return [f"{label}: value {v} for {samples[i]}/{s} outside [{lo}, {hi}]"]
    return []


def plot_data_checks(params, summary):
    """Every plot is drawn from a data table saved next to it. Reconcile those tables with the result files, so a
    plot cannot show different genes, values or up/down labels than the tables report:

    * MA / volcano data = full_results (same genes, baseMean, fold changes, padj, regulation labels),
    * heatmaps = exactly the significant genes (capped) / the top-N genes by padj, each row z-scored,
    * library sizes = count-matrix column sums; PCA/distance/correlation tables cover the samples and are
      well-formed (the VST values themselves come from DESeq2 and are not recomputed here).
    """
    problems = []
    tdir = Path(params["table_dir"])
    samples, counts = _matrix(params["counts"])
    meta = {r["sample"]: r for r in _rows(params["metadata"])}
    voi = params["variable_of_interest"]
    for name, c in summary["contrasts"].items():
        cdir = Path(c["dir"])
        full = _rows(cdir / "full_results.tsv")
        by_id = {r["gene_id"]: r for r in full}
        lfc_col = "log2FoldChange_shrunk" if full and "log2FoldChange_shrunk" in full[0] else "log2FoldChange"
        pdir = cdir / "plot_data"
        ma = _rows(pdir / "ma_plot_data.tsv") if (pdir / "ma_plot_data.tsv").exists() else None
        if ma is None:
            problems.append(f"{name}: MA plot data missing")
        elif [r["gene_id"] for r in ma] != [r["gene_id"] for r in full]:
            problems.append(f"{name}: MA plot data has {len(ma)} gene(s), full_results has {len(full)}")
        else:
            bad = [r["gene_id"] for r in ma if not (_close(r["baseMean"], by_id[r["gene_id"]]["baseMean"])
                   and _close(r["log2FC"], by_id[r["gene_id"]][lfc_col])
                   and r["regulation"] == by_id[r["gene_id"]]["regulation"])]
            if bad:
                problems.append(f"{name}: MA plot data differs from full_results for {len(bad)} gene(s) (e.g. {bad[:3]})")
        expected = [r for r in full if _num(r.get("padj")) is not None]
        vd = _rows(pdir / "volcano_data.tsv") if (pdir / "volcano_data.tsv").exists() else None
        if vd is None:
            problems.append(f"{name}: volcano plot data missing")
        elif [r["gene_id"] for r in vd] != [r["gene_id"] for r in expected]:
            problems.append(f"{name}: volcano plot data has {len(vd)} gene(s), expected the {len(expected)} with padj")
        else:
            bad = [r["gene_id"] for r in vd if not (_close(r["log2FoldChange"], by_id[r["gene_id"]]["log2FoldChange"])
                   and _close(r["padj"], by_id[r["gene_id"]]["padj"])
                   and r["regulation"] == by_id[r["gene_id"]]["regulation"])]
            if bad:
                problems.append(f"{name}: volcano plot data differs from full_results for {len(bad)} gene(s) "
                                f"(e.g. {bad[:3]}) — the plot would highlight different genes than the tables")
        sig = [r["gene_id"] for r in full if r["regulation"] != "ns"]
        hm = pdir / "significant_heatmap_zscores.tsv"
        if len(sig) >= 2:
            if not hm.exists():
                problems.append(f"{name}: significant-gene heatmap data missing")
            else:
                problems += [f"{name}: {x}" for x in
                             _zscored_rows(hm, samples, sig[:params["heatmap_max_genes"]], "significant-gene heatmap")]
        elif hm.exists():
            problems.append(f"{name}: significant-gene heatmap data exists but fewer than 2 genes are significant")
        top = [r["gene_id"] for r in expected][:params["top_n_genes"]]
        th = pdir / "top_genes_heatmap_zscores.tsv"
        if len(top) >= 2:
            if not th.exists():
                problems.append(f"{name}: top-gene heatmap data missing")
            else:
                problems += [f"{name}: {x}" for x in _zscored_rows(th, samples, top, "top-gene heatmap")]
    # QC plots
    lib = _rows(tdir / "library_sizes.tsv") if (tdir / "library_sizes.tsv").exists() else None
    if lib is None:
        problems.append("library size plot data missing")
    else:
        colsum = {s: sum(v[j] for v in counts.values()) for j, s in enumerate(samples)}
        if [r["sample"] for r in lib] != samples:
            problems.append(f"library size plot data samples {[r['sample'] for r in lib]} != {samples}")
        else:
            bad = [r["sample"] for r in lib if _num(r["library_size"]) != colsum[r["sample"]]
                   or r.get("group") != meta[r["sample"]][voi]]
            if bad:
                problems.append(f"library size plot data differs from the count-matrix column sums / groups for {bad}")
    pca = _rows(tdir / "pca_coordinates.tsv") if (tdir / "pca_coordinates.tsv").exists() else None
    if pca is None:
        problems.append("PCA plot data missing")
    elif [r["sample"] for r in pca] != samples or any(
            r.get(v) != meta[r["sample"]][v] for r in pca for v in params["variables"]):
        problems.append("PCA plot data: samples or their groups differ from the metadata")
    var = _rows(tdir / "pca_variance.tsv") if (tdir / "pca_variance.tsv").exists() else []
    total = sum(_num(r["variance_percent"]) or 0 for r in var)
    if var and abs(total - 100) > 0.05 * len(var) + 0.1:
        problems.append(f"PCA variance percentages sum to {total:.1f}, not 100")
    for f, label, diag, tol, lo, hi in (("sample_distances.tsv", "sample distance matrix", 0.0, 1e-6, 0.0, None),
                                        ("sample_correlation.tsv", "sample correlation matrix", 1.0, 1e-4, -1.0, 1.0)):
        if not (tdir / f).exists():
            problems.append(f"{label} data missing")
        else:
            problems += _square_matrix(tdir / f, samples, label, diag, tol, lo, hi)
    return problems


def _contrast_levels(params, name):
    voi = params["variable_of_interest"]
    for num, den in params["contrasts"]:
        if name == f"{voi}_{num}_vs_{den}":
            return num, den
    return None, None


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None  # "NA" (e.g. padj removed by independent filtering)


def check_threshold_sets(full, up, down, sig, alpha, t):
    """Independently re-derive the up/down/significant gene sets from the full table and the configured
    thresholds, and require the written tables to match them EXACTLY (nothing missing, nothing extra)."""
    exp_up, exp_down = set(), set()
    for r in full:
        padj, lfc = _num(r.get("padj")), _num(r.get("log2FoldChange"))
        if padj is None or lfc is None or not padj < alpha:
            continue
        if lfc >= t:
            exp_up.add(r["gene_id"])
        elif lfc <= -t:
            exp_down.add(r["gene_id"])
    problems = []
    for label, expected, table in (("upregulated", exp_up, up), ("downregulated", exp_down, down),
                                   ("significant", exp_up | exp_down, sig)):
        got = [r["gene_id"] for r in table]
        if len(got) != len(set(got)):
            problems.append(f"{label} table lists a gene more than once")
        missing, extra = expected - set(got), set(got) - expected
        if missing:
            problems.append(f"{len(missing)} gene(s) meet the {label} criteria but are missing from the table "
                            f"(e.g. {sorted(missing)[:3]})")
        if extra:
            problems.append(f"{len(extra)} gene(s) in the {label} table do not meet padj < {alpha} and "
                            f"|log2FC| >= {t} (e.g. {sorted(extra)[:3]})")
    labels = {r["gene_id"]: r.get("regulation") for r in full}
    wrong = [g for g in exp_up if labels.get(g) != "up"] + [g for g in exp_down if labels.get(g) != "down"]
    if wrong:
        problems.append(f"'regulation' column disagrees with the thresholds for {len(wrong)} gene(s)")
    return problems


def _rows(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f if l.strip()]
