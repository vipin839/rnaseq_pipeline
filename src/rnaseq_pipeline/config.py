"""Configuration loading, merging, validation and snapshotting."""
import copy
import os
from datetime import datetime
from pathlib import Path

import yaml

from . import PACKAGE_DIR, PipelineError
from . import validators as V

DEFAULT_CONFIG = PACKAGE_DIR / "config" / "default_config.yaml"
CATALOG = PACKAGE_DIR / "config" / "reference_catalog.yaml"
USER_CONFIG = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "rnaseq-pipeline" / "config.yaml"

STRANDEDNESS = ("auto", "unstranded", "forward", "reverse")


def load_yaml(path):
    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data or {}


def save_yaml(data, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, default_flow_style=False)
    os.replace(tmp, path)


def deep_merge(base, override):
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load(project_config_path=None):
    cfg = load_yaml(DEFAULT_CONFIG)
    if project_config_path and Path(project_config_path).exists():
        cfg = deep_merge(cfg, load_yaml(project_config_path))
    return cfg


def get(cfg, dotted, default=None):
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def set_value(cfg, dotted, value):
    parts = dotted.split(".")
    cur = cfg
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = value


def resolve_threads(cfg, cores):
    t = cfg.get("threads", "auto")
    if t in (None, "auto"):
        return max(1, cores - 1)
    return V.threads(t, cores)


def resolve_memory_gb(cfg, ram_gb):
    m = cfg.get("max_memory_gb", "auto")
    if m in (None, "auto"):
        return max(1.0, round(ram_gb * 0.8, 1))
    return V.number(m, "max_memory_gb", 0.5, ram_gb)


def validate(cfg, cores=None, ram_gb=None):
    """Return a list of human-readable problems (empty if valid)."""
    errs = []

    def chk(fn, *a, **k):
        try:
            fn(*a, **k)
        except ValueError as e:
            errs.append(str(e))

    t = cfg.get("threads")
    if t != "auto":
        chk(V.threads, t, cores or os.cpu_count() or 1)
    if cfg.get("max_memory_gb") != "auto":
        chk(V.number, cfg.get("max_memory_gb"), "max_memory_gb", 0.5)
    if cfg.get("read_type") not in ("auto", "paired", "single"):
        errs.append("read_type must be auto, paired or single")
    if cfg.get("import_mode") not in ("symlink", "copy"):
        errs.append("import_mode must be symlink or copy")
    if cfg.get("strandedness") not in STRANDEDNESS:
        errs.append(f"strandedness must be one of {', '.join(STRANDEDNESS)}")
    if cfg.get("trim_adapters") not in ("ask", "always", "never"):
        errs.append("trim_adapters must be ask, always or never")
    chk(V.probability, cfg.get("alpha"), "alpha")
    chk(V.number, cfg.get("log2fc_threshold"), "log2fc_threshold", 0, 20)
    chk(V.number, cfg.get("lfc_test_threshold"), "lfc_test_threshold", 0, 20)
    chk(V.positive_int, cfg.get("top_n_genes"), "top_n_genes", 1, 100000)
    chk(V.positive_int, cfg.get("heatmap_max_genes"), "heatmap_max_genes", 2, 10000)
    mcf = cfg.get("min_count_filter") or {}
    chk(V.positive_int, mcf.get("min_count"), "min_count_filter.min_count", 0)
    if mcf.get("min_samples") != "smallest_group":
        chk(V.positive_int, mcf.get("min_samples"), "min_count_filter.min_samples", 1)
    if cfg.get("transformation") not in ("vst", "rlog"):
        errs.append("transformation must be vst or rlog")
    if cfg.get("lfc_shrinkage") not in ("apeglm", "normal", "none"):
        errs.append("lfc_shrinkage must be apeglm, normal or none")
    chk(V.design_formula, cfg.get("design_formula"))
    qg = cfg.get("quality_gate") or {}
    for k in ("adapter_max_percent_trim", "adapter_max_percent_review", "overrepresented_max_percent",
              "n_content_max_percent", "duplication_review_percent"):
        chk(V.percentage, qg.get(k), f"quality_gate.{k}")
    for k in ("tail_mean_quality_trim", "tail_mean_quality_review", "per_base_lower_quartile_min"):
        chk(V.number, qg.get(k), f"quality_gate.{k}", 0, 60)
    fp = cfg.get("fastp_parameters") or {}
    chk(V.number, fp.get("qualified_quality_phred"), "fastp.qualified_quality_phred", 0, 60)
    chk(V.percentage, fp.get("unqualified_percent_limit"), "fastp.unqualified_percent_limit")
    chk(V.positive_int, fp.get("length_required"), "fastp.length_required", 1, 10000)
    chk(V.positive_int, fp.get("cut_right_window_size"), "fastp.cut_right_window_size", 1, 1000)
    chk(V.number, fp.get("cut_right_mean_quality"), "fastp.cut_right_mean_quality", 0, 60)
    si = cfg.get("strandedness_inference") or {}
    chk(V.number, si.get("stranded_min_fraction"), "strandedness_inference.stranded_min_fraction", 0.5, 1)
    chk(V.number, si.get("unstranded_max_diff"), "strandedness_inference.unstranded_max_diff", 0, 0.5)
    fc = cfg.get("featurecounts_parameters") or {}
    chk(V.positive_int, fc.get("min_mapping_quality"), "featurecounts.min_mapping_quality", 0, 255)
    if not fc.get("feature_type") or not fc.get("attribute"):
        errs.append("featurecounts feature_type and attribute must be set")
    for key in ("fastp_parameters", "hisat2_parameters", "stringtie_parameters", "featurecounts_parameters"):
        extra = (cfg.get(key) or {}).get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, (str, int, float)) for x in extra):
            errs.append(f"{key}.extra_args must be a list of strings")
    chk(V.number, cfg.get("min_overall_alignment_rate_fail"), "min_overall_alignment_rate_fail", 0, 100)
    chk(V.number, cfg.get("min_overall_alignment_rate_warn"), "min_overall_alignment_rate_warn", 0, 100)
    n = cfg.get("ncbi") or {}
    if not isinstance(n.get("email", ""), str) or not isinstance(n.get("api_key", ""), str):
        errs.append("ncbi.email and ncbi.api_key must be strings")
    elif n.get("email") and "@" not in n["email"]:
        errs.append("ncbi.email must be an email address (or empty)")
    dl = cfg.get("download") or {}
    mr = dl.get("max_reads", 0)
    chk(V.positive_int, mr, "download.max_reads", 0)
    if isinstance(mr, int) and 0 < mr < 1000:
        errs.append("download.max_reads must be 0 (full data) or at least 1000 reads (pilot subset)")
    return errs


def require_valid(cfg, cores=None):
    errs = validate(cfg, cores)
    if errs:
        raise PipelineError("configuration is invalid:\n  - " + "\n  - ".join(errs),
                            remedy="fix the values in config/project_config.yaml or via menu 7")


def snapshot(cfg, project_config_dir, extra=None):
    """Save the exact final configuration used for this run."""
    from .secrets import redacted
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    data = redacted(cfg)
    if extra:
        data["_resolved"] = extra
    path = Path(project_config_dir) / f"config_used_{stamp}.yaml"
    save_yaml(data, path)
    return path


def load_catalog():
    return load_yaml(CATALOG)
