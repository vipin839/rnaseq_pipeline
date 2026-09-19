"""Final self-contained HTML report (reports/final_pipeline_report.html) linking to underlying files."""
import csv
import html
import json
import os
from datetime import datetime
from pathlib import Path

from . import __version__

CSS = """
:root{--fg:#1d2330;--muted:#5b6475;--bg:#ffffff;--panel:#f5f7fa;--line:#dde2ea;--accent:#1f5fae;--ok:#2e7d32;--warn:#b26a00;--bad:#c62828}
@media (prefers-color-scheme: dark){:root{--fg:#e6e9ef;--muted:#9aa3b2;--bg:#12151b;--panel:#1b2029;--line:#2c3340;--accent:#6ea8ff;--ok:#66bb6a;--warn:#ffb74d;--bad:#ef5350}}
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;color:var(--fg);background:var(--bg);margin:0;line-height:1.5}
main{max-width:1150px;margin:0 auto;padding:24px 16px 80px}
nav{position:sticky;top:0;background:var(--panel);border-bottom:1px solid var(--line);padding:8px 16px;font-size:13px;overflow-x:auto;white-space:nowrap}
nav a{margin-right:12px;color:var(--accent);text-decoration:none}
h1{margin:8px 0 4px}h2{border-bottom:2px solid var(--line);padding-bottom:4px;margin-top:40px}
.muted{color:var(--muted)}.ok{color:var(--ok);font-weight:600}.warn{color:var(--warn);font-weight:600}.bad{color:var(--bad);font-weight:600}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0 16px;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:4px 8px;text-align:left;white-space:nowrap}th{background:var(--panel)}
.kv td:first-child{font-weight:600;width:230px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:16px}
figure{margin:0;background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:8px}
figure img{width:100%;height:auto;background:#fff}figcaption{font-size:12px;color:var(--muted)}
code,pre{background:var(--panel);border-radius:4px;padding:1px 4px;font-size:12px}pre{padding:8px;overflow-x:auto}
.banner{background:#fff3cd;color:#664d03;border:1px solid #ffe69c;padding:8px 12px;border-radius:6px}
"""


def esc(x):
    return html.escape(str(x)) if x is not None else ""


def link(report_dir, path, label=None):
    p = Path(path)
    if not p.exists():
        return f"<span class='muted'>{esc(label or p.name)} (missing)</span>"
    rel = os.path.relpath(p, report_dir)
    return f"<a href='{esc(rel)}'>{esc(label or p.name)}</a>"


def tsv_table(path, max_rows=50, cols=None):
    p = Path(path)
    if not p.exists():
        return "<p class='muted'>(not available)</p>"
    with open(p, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        head = next(r, [])
        idx = [head.index(c) for c in cols if c in head] if cols else list(range(len(head)))
        rows = []
        for i, row in enumerate(r):
            if i >= max_rows:
                break
            rows.append([row[j] if j < len(row) else "" for j in idx])
    th = "".join(f"<th>{esc(head[j])}</th>" for j in idx)
    tb = "".join("<tr>" + "".join(f"<td>{esc(fmt_cell(c))}</td>" for c in row) + "</tr>" for row in rows)
    return f"<table><tr>{th}</tr>{tb}</table>"


def fmt_cell(c):
    try:
        v = float(c)
        if v != v:
            return c
        if abs(v) >= 1e5 and v == int(v):
            return f"{int(v):,}"
        if v != int(v):
            return f"{v:.4g}"
    except ValueError:
        pass
    return c


def kv_table(pairs):
    return "<table class='kv'>" + "".join(f"<tr><td>{esc(k)}</td><td>{v}</td></tr>" for k, v in pairs) + "</table>"


def fig(report_dir, png, caption):
    p = Path(png)
    if not p.exists():
        return ""
    rel = os.path.relpath(p, report_dir)
    pdf = p.with_suffix(".pdf")
    extra = f" · {link(report_dir, pdf, 'PDF')}" if pdf.exists() else ""
    return f"<figure><a href='{esc(rel)}'><img src='{esc(rel)}' alt='{esc(caption)}' loading='lazy'></a>" \
           f"<figcaption>{esc(caption)}{extra}</figcaption></figure>"


def generate(ctx):
    p, cfg = ctx.project, ctx.cfg
    R = p.path("reports")
    s = p.state
    ref = s.get("reference") or {}
    strand = s.get("strandedness") or {}
    d = s.get("design") or {}
    sections = []

    def sec(id_, title, body):
        sections.append((id_, title, body))

    pilot = any(r.get("pilot_max_reads") for r in p.samples.values())
    ov = [("Project", esc(s["name"])), ("Project ID", esc(s["project_id"])), ("Created", esc(s["created"])),
          ("Report generated", datetime.now().strftime("%Y-%m-%d %H:%M")), ("Pipeline version", __version__),
          ("Data source", esc(s.get("data_source"))), ("Read type", esc(s.get("read_type"))),
          ("Samples (active / total)", f"{len(ctx.samples)} / {len(p.samples)}"),
          ("Reference", esc(ref.get("label"))), ("Strandedness", esc(f"{strand.get('value')} ({strand.get('source')})")),
          ("Design", esc(d.get("formula"))),
          ("Thresholds", f"padj &lt; {cfg['alpha']} and |log2FC| &ge; {cfg['log2fc_threshold']}")]
    body = (("<p class='banner'>PILOT MODE: only a subset of reads was downloaded. Results are for pipeline "
             "testing only and must not be interpreted biologically.</p>") if pilot else "") + kv_table(ov)
    sec("overview", "Project overview", body)

    rows = []
    for sid, r in p.samples.items():
        m = r.get("metadata") or {}
        rows.append((sid, r.get("status"), r.get("source"), r.get("accession") or "", r.get("layout"),
                     f"{(r.get('reads') or {}).get('R1', ''):,}" if (r.get("reads") or {}).get("R1") else "",
                     r.get("read_length") or "", m.get("geo_title") or m.get("sample_title") or "",
                     r.get("fail_reason") or ""))
    tb = "".join("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>" for row in rows)
    sec("input", "Input data",
        f"<table><tr><th>Sample</th><th>Status</th><th>Source</th><th>Accession</th><th>Layout</th><th>Reads (R1)</th>"
        f"<th>Max length</th><th>Description</th><th>Failure</th></tr>{tb}</table>"
        f"<p>FASTQ validation: {link(R, p.path('data', 'metadata', 'fastq_validation_report.tsv'))}</p>")

    sec("metadata", "Sample metadata / experimental design",
        tsv_table(p.path("counts", "sample_metadata.tsv")) +
        kv_table([("Formula", esc(d.get("formula"))), ("Variable of interest", esc(d.get("variable_of_interest"))),
                  ("Reference level", esc(d.get("reference_level"))),
                  ("Contrasts", esc(", ".join(f"{a} vs {b}" for a, b in d.get("contrasts", [])))),
                  ("Confirmed at", esc(d.get("confirmed_at"))),
                  ("Files", link(R, p.path("counts", "sample_metadata.tsv")) + " · " + link(R, p.path("counts", "design.json")))]))

    sec("reference", "Reference",
        kv_table([("Label", esc(ref.get("label"))), ("Organism", esc(ref.get("organism"))),
                  ("Assembly", esc(ref.get("assembly"))), ("Annotation", esc(ref.get("annotation_release"))),
                  ("Genes in annotation", esc(ref.get("genes"))),
                  ("Splice sites", "in index" if ref.get("index_has_splice_sites") else "supplied at alignment"),
                  ("Manifest", link(R, p.path("reference", "reference_manifest.yaml")))]))

    sec("versions", "Software versions", tsv_table(p.path("logs", "software_versions.tsv"), 100))

    qa_raw = p.path("qc", "assessment", "raw", "quality_assessment.tsv")
    sec("rawqc", "Raw read QC",
        f"<p>{link(R, p.path('qc', 'multiqc_raw', 'multiqc_report.html'), 'MultiQC report (raw)')} · "
        f"{link(R, p.path('qc', 'assessment', 'raw', 'quality_assessment.html'), 'Quality assessment')}</p>"
        + tsv_table(qa_raw, 100, ["sample", "mate", "status", "total_sequences", "gc_percent", "tail_mean_quality",
                                  "adapter_max_percent", "duplication_percent", "reasons"]))

    qd = s.get("qc_decision") or {}
    tr_rows = [(sid, r["trimmed"]["reads_before"], r["trimmed"]["reads_after"],
                f"{100 * r['trimmed']['reads_after'] / max(1, r['trimmed']['reads_before']):.1f}%",
                link(R, p.abs(r["trimmed"]["fastp_html"]), "fastp report"))
               for sid, r in p.samples.items() if r.get("trimmed")]
    tb = "".join("<tr>" + "".join(f"<td>{c if isinstance(c, str) and c.startswith('<') else esc(fmt_cell(str(c)))}</td>"
                                  for c in row) + "</tr>" for row in tr_rows)
    sec("trimming", "Trimming",
        kv_table([("Recommendation", esc(qd.get("recommended"))), ("Decision", esc(qd.get("decision"))),
                  ("fastp parameters", f"<code>{esc(json.dumps(cfg['fastp_parameters']))}</code>")]) +
        (f"<table><tr><th>Sample</th><th>Reads before</th><th>Reads after</th><th>Kept</th><th>Report</th></tr>{tb}</table>"
         if tr_rows else "<p class='muted'>Trimming was not performed.</p>"))
    if tr_rows:
        sec("posttrim", "Post-trimming QC",
            f"<p>{link(R, p.path('qc', 'multiqc_trimmed', 'multiqc_report.html'), 'MultiQC report (trimmed)')}</p>" +
            tsv_table(p.path("qc", "assessment", "trimmed", "quality_assessment.tsv"), 100,
                      ["sample", "mate", "status", "total_sequences", "tail_mean_quality", "adapter_max_percent", "reasons"]))

    sec("alignment", "Alignment (HISAT2)",
        tsv_table(p.path("alignment", "reports", "alignment_summary.tsv"), 200,
                  ["sample", "status", "input_reads", "overall_alignment_rate", "uniquely_aligned_pct", "mapped_pct",
                   "properly_paired_pct", "secondary", "supplementary", "unmapped", "problems"]) +
        f"<p>{link(R, p.path('alignment', 'reports', 'alignment_summary.tsv'))}</p>")
    sec("bamqc", "BAM QC",
        tsv_table(p.path("alignment", "reports", "bam_qc_summary.tsv"), 200) +
        f"<p>{link(R, p.path('qc', 'multiqc_alignment', 'multiqc_report.html'), 'MultiQC alignment report')}</p>")
    ev = strand.get("evidence") or {}
    ev_rows = "".join(f"<tr><td>{esc(k)}</td><td>{esc(v.get('forward'))}</td><td>{esc(v.get('reverse'))}</td>"
                      f"<td>{esc(v.get('undetermined'))}</td></tr>" for k, v in ev.items())
    sec("strand", "Library strandedness",
        kv_table([("Value", esc(strand.get("value"))), ("Source", esc(strand.get("source"))),
                  ("Confidence", esc(strand.get("confidence"))),
                  ("featureCounts -s", esc(strand.get("featurecounts_flag"))),
                  ("StringTie flag", esc(strand.get("stringtie_flag") or "(none)"))]) +
        (f"<table><tr><th>Sample</th><th>Forward fraction</th><th>Reverse fraction</th><th>Undetermined</th></tr>{ev_rows}</table>"
         if ev_rows else ""))
    sec("stringtie", "StringTie2 (transcript-level, reference-guided)",
        "<p>Reference-guided quantification (<code>-e -G annotation.gtf</code>) of annotated transcripts. "
        "These TPM values are an independent result and are <b>not</b> the DESeq2 input.</p>"
        f"<p>{link(R, p.path('stringtie', 'merged', 'gene_tpm_matrix.tsv'))} · "
        f"{link(R, p.path('stringtie', 'merged', 'transcript_tpm_matrix.tsv'))}</p>")
    fc_stats = p.path("featurecounts", "featurecounts_stats.json")
    fc_rows = ""
    if fc_stats.exists():
        st = json.loads(fc_stats.read_text())
        fc_rows = "".join(f"<tr><td>{esc(k)}</td><td>{v['assigned']:,}</td><td>{v['total']:,}</td>"
                          f"<td>{v['assigned_pct']}%</td></tr>" for k, v in st.items())
    sec("featurecounts", "featureCounts",
        f"<table><tr><th>Sample</th><th>Assigned</th><th>Total</th><th>Assigned %</th></tr>{fc_rows}</table>"
        f"<p>{link(R, p.path('featurecounts', 'featurecounts.txt'))} · {link(R, p.path('featurecounts', 'featurecounts.summary'))}</p>")
    cms = p.path("counts", "count_matrix_summary.txt")
    sec("counts", "Count matrix",
        f"<pre>{esc(cms.read_text() if cms.exists() else '(missing)')}</pre>"
        f"<p>{link(R, p.path('counts', 'gene_count_matrix.tsv'))} · {link(R, p.path('counts', 'gene_count_matrix.csv'))}</p>")

    ds = p.path("results", "deseq2", "deseq2_summary.json")
    summ = json.loads(ds.read_text()) if ds.exists() else {}
    plots = p.path("results", "plots")
    qc_figs = "".join(fig(R, plots / f"{n}.png", t) for n, t in (
        ("pca", "PCA"), ("sample_distance_heatmap", "Sample distance heatmap"),
        ("sample_correlation_heatmap", "Sample correlation heatmap"), ("library_sizes", "Library sizes"),
        ("normalized_count_distribution", "Normalized count distributions")))
    sec("deseq2", "DESeq2",
        kv_table([("Genes in matrix", esc(summ.get("genes_input"))),
                  ("Genes tested (after low-count filter)", esc(summ.get("genes_tested"))),
                  ("Low-count filter", esc(f">= {cfg['min_count_filter']['min_count']} counts in >= "
                                           f"{summ.get('min_samples_filter')} samples")),
                  ("Normalized counts", link(R, p.path("results", "deseq2", "normalized_counts.tsv"))),
                  ("Parameters", link(R, p.path("results", "deseq2", "deseq2_params.json")))]) +
        f"<div class='grid'>{qc_figs}</div>")
    body = ""
    for name, c in (summ.get("contrasts") or {}).items():
        cdir = Path(c["dir"])
        pdir = plots / name
        body += f"<h3>{esc(name)}</h3>" + kv_table([
            ("Comparison", esc(f"{c['numerator']} vs {c['denominator']} (denominator = baseline)")),
            ("Significance criteria", f"padj &lt; {summ.get('alpha')} AND |log2FC| &ge; {summ.get('log2fc_threshold')}"),
            ("Genes tested", esc(c["genes_tested"])),
            ("Significant", f"<b>{c['significant']}</b> (<span class='bad'>up {c['up']}</span>, "
                            f"<span style='color:var(--accent);font-weight:600'>down {c['down']}</span>)"),
            ("Files", " · ".join(link(R, cdir / f) for f in ("full_results.tsv", "significant_results.tsv",
                                                             "upregulated.tsv", "downregulated.tsv", "statistics.txt")))])
        body += "<div class='grid'>" + "".join(fig(R, pdir / f"{n}.png", t) for n, t in (
            ("volcano_plot", "Volcano plot"), ("ma_plot", "MA plot"),
            ("significant_genes_heatmap", "Significant genes heatmap"), ("top_genes_heatmap", "Top genes heatmap"))) + "</div>"
        cols = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]
        body += "<h4>Top upregulated</h4>" + tsv_table(cdir / "top_upregulated.tsv", 15, cols)
        body += "<h4>Top downregulated</h4>" + tsv_table(cdir / "top_downregulated.tsv", 15, cols)
    sec("de", "Differential genes (upregulated / downregulated)", body or "<p class='muted'>(no results)</p>")

    errs = p.path("logs", "pipeline_errors.log")
    lines = errs.read_text().splitlines()[-60:] if errs.exists() else []
    sec("errors", "Errors / warnings",
        (f"<pre>{esc(chr(10).join(lines))}</pre>" if lines else "<p class='ok'>No warnings or errors were logged.</p>")
        + f"<p>{link(R, errs)} · {link(R, p.path('logs', 'pipeline.log'))} · {link(R, p.path('logs', 'command_history.log'))}</p>")
    sec("params", "Parameters", f"<pre>{esc(json.dumps(cfg, indent=2, default=str))}</pre>")
    sec("repro", "Reproducibility",
        "<ul>" + "".join(f"<li>{link(R, p.path('pipeline_manifest', f))}</li>" for f in
                         ("manifest.json", "manifest.yaml", "environment.yml", "package_versions.txt",
                          "system_information.txt", "R_sessionInfo.txt")) + "</ul>"
        f"<p>Every executed command is recorded in {link(R, p.path('logs', 'commands.jsonl'))}.</p>")

    nav = "".join(f"<a href='#{i}'>{esc(t)}</a>" for i, t, _ in sections)
    content = "".join(f"<section id='{i}'><h2>{esc(t)}</h2>{b}</section>" for i, t, b in sections)
    out = R / "final_pipeline_report.html"
    out.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>RNA-seq report: {esc(s['name'])}</title>
<style>{CSS}</style></head><body><nav>{nav}</nav><main>
<h1>Bulk RNA-seq analysis report</h1><p class="muted">{esc(s['name'])} · pipeline v{__version__}</p>
{content}</main></body></html>""")
    return out
