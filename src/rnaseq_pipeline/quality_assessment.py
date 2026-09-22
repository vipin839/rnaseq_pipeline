"""Quality gate: parse FastQC data and recommend PASS / REVIEW / TRIMMING RECOMMENDED / FAIL.

Every threshold comes from config `quality_gate` (documented in docs/PIPELINE_METHODS.md).
The result is a recommendation; the user makes the final decision.
"""
import html
import json
import math
import statistics
from pathlib import Path

from . import qc_manager

ORDER = {"PASS": 0, "REVIEW": 1, "TRIMMING RECOMMENDED": 2, "FAIL": 3}


def parse_modules(data_text):
    """Return {module_name: {'status': s, 'rows': [[...]], 'header': [...], 'extra': {...}}}."""
    mods, cur = {}, None
    for line in data_text.splitlines():
        if line.startswith(">>END_MODULE"):
            cur = None
            continue
        if line.startswith(">>"):
            name, _, status = line[2:].partition("\t")
            cur = mods.setdefault(name, {"status": status.strip(), "rows": [], "header": None, "extra": {}})
            continue
        if cur is None:
            continue
        if line.startswith("#"):
            parts = line[1:].split("\t")
            if len(parts) == 2 and cur["header"] is None and not cur["rows"]:
                try:
                    cur["extra"][parts[0]] = float(parts[1])
                    continue
                except ValueError:
                    pass
            cur["header"] = parts
            continue
        cur["rows"].append(line.split("\t"))
    return mods


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def metrics_from_fastqc(zip_path):
    summary_text, data = qc_manager.read_fastqc_zip(zip_path)
    mods = parse_modules(data)
    m = {"fastqc_status": {}}
    for line in summary_text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            m["fastqc_status"][parts[1]] = parts[0]
    basic = {r[0]: r[1] for r in mods.get("Basic Statistics", {}).get("rows", []) if len(r) >= 2}
    m["total_sequences"] = int(basic.get("Total Sequences", 0))
    m["gc_percent"] = _float(basic.get("%GC"))
    m["sequence_length"] = basic.get("Sequence length", "")
    pbq = mods.get("Per base sequence quality", {}).get("rows", [])
    means = [_float(r[1]) for r in pbq if len(r) > 1]
    lq = [_float(r[3]) for r in pbq if len(r) > 3]
    if means:
        k = max(1, math.ceil(len(means) * 0.10))
        m["tail_mean_quality"] = round(statistics.fmean(means[-k:]), 2)
        m["overall_mean_quality"] = round(statistics.fmean(means), 2)
    else:
        m["tail_mean_quality"] = m["overall_mean_quality"] = math.nan
    m["min_lower_quartile"] = min((x for x in lq if not math.isnan(x)), default=math.nan)
    ad = mods.get("Adapter Content", {}).get("rows", [])
    m["adapter_max_percent"] = round(max((_float(v) for r in ad for v in r[1:]), default=0.0), 3)
    ov = mods.get("Overrepresented sequences", {}).get("rows", [])
    m["overrepresented_max_percent"] = round(max((_float(r[2]) for r in ov if len(r) > 2), default=0.0), 3)
    dup = mods.get("Sequence Duplication Levels", {}).get("extra", {})
    dd = dup.get("Total Deduplicated Percentage")
    m["duplication_percent"] = round(100 - dd, 2) if dd is not None else math.nan
    nc = mods.get("Per base N content", {}).get("rows", [])
    m["n_content_max_percent"] = round(max((_float(r[1]) for r in nc if len(r) > 1), default=0.0), 3)
    psq = mods.get("Per sequence quality scores", {}).get("rows", [])
    if psq:
        tot = sum(_float(r[1]) for r in psq)
        low = sum(_float(r[1]) for r in psq if _float(r[0]) < 20)
        m["low_quality_reads_percent"] = round(100 * low / tot, 3) if tot else 0.0
    return m


def assess(per_file, thresholds, paired):
    """per_file: {sample: {'R1': metrics, 'R2': metrics?}}. Returns {sample: {...status, reasons, metrics}}."""
    T = thresholds
    gcs = [m["gc_percent"] for s in per_file.values() for m in s.values() if not math.isnan(m["gc_percent"])]
    gc_median = statistics.median(gcs) if gcs else math.nan
    out = {}
    for sample, mates in per_file.items():
        reasons, status = [], "PASS"

        def flag(level, msg):  # called only within this iteration, so it always sees this sample's list
            nonlocal status
            reasons.append(f"[{level}] {msg}")  # noqa: B023
            if ORDER[level] > ORDER[status]:
                status = level

        if paired:
            n1 = mates.get("R1", {}).get("total_sequences")
            n2 = mates.get("R2", {}).get("total_sequences")
            if n1 != n2:
                flag("FAIL", f"R1/R2 read counts differ ({n1} vs {n2})")
        for mate, m in mates.items():
            p = f"{mate}: " if paired else ""
            a = m["adapter_max_percent"]
            if a > T["adapter_max_percent_trim"]:
                flag("TRIMMING RECOMMENDED", f"{p}adapter contamination up to {a:.1f}% "
                     f"(threshold {T['adapter_max_percent_trim']}%)")
            elif a > T["adapter_max_percent_review"]:
                flag("REVIEW", f"{p}adapter content up to {a:.1f}%")
            tq = m["tail_mean_quality"]
            if not math.isnan(tq):
                if tq < T["tail_mean_quality_trim"]:
                    flag("TRIMMING RECOMMENDED", f"{p}low-quality read tails (mean Q{tq:.1f} over last 10% "
                         f"of positions; threshold Q{T['tail_mean_quality_trim']})")
                elif tq < T["tail_mean_quality_review"]:
                    flag("REVIEW", f"{p}moderate quality drop at read ends (mean Q{tq:.1f})")
            lq = m["min_lower_quartile"]
            if not math.isnan(lq) and lq < T["per_base_lower_quartile_min"]:
                flag("REVIEW", f"{p}per-base lower quartile falls to Q{lq:.0f}")
            if m["overrepresented_max_percent"] > T["overrepresented_max_percent"]:
                flag("REVIEW", f"{p}overrepresented sequence at {m['overrepresented_max_percent']:.2f}% "
                     "(common in RNA-seq for abundant transcripts; inspect FastQC)")
            if m["n_content_max_percent"] > T["n_content_max_percent"]:
                flag("REVIEW", f"{p}N content up to {m['n_content_max_percent']:.1f}% at some positions")
            d = m["duplication_percent"]
            if not math.isnan(d) and d > T["duplication_review_percent"]:
                flag("REVIEW", f"{p}duplication {d:.0f}% (often biological in RNA-seq; not a trimming issue)")
            if not math.isnan(gc_median) and abs(m["gc_percent"] - gc_median) > T["gc_deviation_review"]:
                flag("REVIEW", f"{p}%GC {m['gc_percent']:.0f} deviates from cohort median {gc_median:.0f}")
            if m["total_sequences"] < T["min_reads_review"]:
                flag("REVIEW", f"{p}only {m['total_sequences']:,} reads (threshold {T['min_reads_review']:,})")
        out[sample] = {"status": status, "reasons": reasons, "metrics": mates}
    return out


def write_outputs(result, out_dir, thresholds, title):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    j = out_dir / "quality_assessment.json"
    j.write_text(json.dumps({"title": title, "thresholds": thresholds, "samples": result}, indent=2,
                            default=lambda x: None if isinstance(x, float) and math.isnan(x) else str(x)))
    cols = ["total_sequences", "gc_percent", "tail_mean_quality", "min_lower_quartile", "adapter_max_percent",
            "overrepresented_max_percent", "duplication_percent", "n_content_max_percent"]
    t = out_dir / "quality_assessment.tsv"
    with open(t, "w") as f:
        f.write("\t".join(["sample", "mate", "status"] + cols + ["reasons"]) + "\n")
        for s, r in result.items():
            for mate, m in r["metrics"].items():
                f.write("\t".join([s, mate, r["status"]] + [str(m.get(c, "")) for c in cols]
                                  + ["; ".join(r["reasons"])]) + "\n")
    h = out_dir / "quality_assessment.html"
    color = {"PASS": "#2e7d32", "REVIEW": "#b26a00", "TRIMMING RECOMMENDED": "#c62828", "FAIL": "#6a1b9a"}
    rows = "".join(
        f"<tr><td>{html.escape(s)}</td><td style='color:{color[r['status']]};font-weight:600'>{r['status']}</td>"
        f"<td>{'<br>'.join(html.escape(x) for x in r['reasons']) or '—'}</td></tr>"
        for s, r in result.items())
    th = "".join(f"<li><code>{html.escape(k)}</code> = {v}</li>" for k, v in thresholds.items())
    h.write_text(f"""<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2em;max-width:1100px}}table{{border-collapse:collapse;width:100%}}
td,th{{border:1px solid #ccc;padding:6px;text-align:left;vertical-align:top}}th{{background:#f3f3f3}}</style></head>
<body><h1>{html.escape(title)}</h1><p>Recommendations are advisory. FastQC WARN/FAIL flags do not by themselves
mean data are unusable (e.g. duplication and per-base content biases are expected in RNA-seq).</p>
<table><tr><th>Sample</th><th>Status</th><th>Reasons</th></tr>{rows}</table>
<h2>Thresholds used</h2><ul>{th}</ul></body></html>""")
    return [j, t, h]
