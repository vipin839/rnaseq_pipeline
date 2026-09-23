"""Reproducibility manifest (pipeline_manifest/)."""
import json
import platform
import sys
from datetime import datetime

from . import __version__, dependency_manager, ui
from .checkpoint import fingerprint
from .secrets import redacted
from . import config as C


# written only when the tools run from conda environments (tool versions are always in "tool_versions")
CONDA_EXPORT_FILES = ("environment.yml", "package_versions.txt")
REQUIRED = ("tool_versions", "r_version", "input_files", "reference_checksums", "parameters", "design", "commands")


def problems(project):
    """What a complete manifest must contain; used by the report stage (resume) and --validate-project."""
    f = project.path("pipeline_manifest", "manifest.json")
    if not f.exists():
        return ["manifest.json missing"]
    try:
        m = json.loads(f.read_text())
    except ValueError:
        return ["manifest.json is not valid JSON"]
    missing = [k for k in REQUIRED if not m.get(k)]
    return [f"manifest lacks {', '.join(missing)} (written by an older version or an incomplete run; "
            "re-running the report step regenerates it)"] if missing else []


def write(ctx):
    p, cfg = ctx.project, ctx.cfg
    out = p.path("pipeline_manifest")
    out.mkdir(parents=True, exist_ok=True)
    tools = dependency_manager.detect_tools()
    r_info, r_pkgs = dependency_manager.detect_r(ctx.envs.rscript())
    vers = dependency_manager.write_versions(p.path("logs", "software_versions.tsv"), tools, r_info, r_pkgs)
    env_files = ctx.envs.export(out, ctx.sysinfo)
    exe = ctx.envs.conda_bin()
    envs_found = [k for k in ("tools", "r") if ctx.envs.prefix(k)]
    conda_export = {
        "exportable": bool(exe and envs_found),
        "files": [f.name for f in env_files if f.name in CONDA_EXPORT_FILES],
        "reason": None if exe and envs_found else
        ("no conda/mamba executable found" if not exe else "the tools were not run from a conda environment"),
    }
    checksums = {}
    for stage in ctx.cp.dir.glob("*.json"):
        if ".invalid." in stage.name:
            continue
        data = ctx.cp.read(stage.stem) or {}
        for o in data.get("outputs", []):
            if "sha256" in o:
                checksums[o["path"]] = o["sha256"]
    s = p.state
    first = next(iter(p.samples.values()), {})
    inputs = {}
    for sid, rec in p.samples.items():
        files = {}
        for mate in ("r1", "r2"):
            if rec.get(mate) and p.abs(rec[mate]).exists():
                files[mate] = {"path": rec[mate], "resolved": str(p.abs(rec[mate]).resolve()),
                               **fingerprint(p.abs(rec[mate]))}
        if rec.get("provider_md5"):
            files["provider_md5_verified"] = rec["provider_md5"]
        inputs[sid] = files
    ref_manifest = p.path("reference", "reference_manifest.yaml")
    ref_info = C.load_yaml(ref_manifest) if ref_manifest.exists() else {}
    cmd_log = p.path("logs", "commands.jsonl")
    n_cmds = sum(1 for _ in open(cmd_log)) if cmd_log.exists() else 0
    manifest = {
        "project_id": s["project_id"], "project_name": s["name"], "created": s["created"],
        "manifest_date": datetime.now().isoformat(timespec="seconds"), "pipeline_version": __version__,
        "system": {"os": ctx.sysinfo.get("distribution"), "kernel": ctx.sysinfo.get("kernel"),
                   "architecture": ctx.sysinfo.get("architecture"), "cpu_cores": ctx.sysinfo.get("cpu_cores"),
                   "ram_gb": ctx.sysinfo.get("ram_total_gb"), "python": sys.version.split()[0],
                   "platform": platform.platform()},
        "r_version": r_info.get("version"),
        "tool_versions": {t["key"]: t["version"] for t in tools},
        "conda_export": conda_export,
        "tool_paths": {t["key"]: t["path"] for t in tools},
        "r_packages": {x["package"]: x["version"] for x in r_pkgs},
        "reference": s.get("reference", {}) and {k: v for k, v in s["reference"].items() if k != "package"},
        "reference_checksums": ref_info.get("sha256"),
        "reference_download": ref_info.get("files"),
        "input_files": inputs,
        "commands": {"log": "logs/commands.jsonl", "count": n_cmds,
                     "note": "every external command with arguments, exit codes and duration"},
        "data_source": s.get("data_source"), "accessions": s.get("accessions"),
        "samples": {sid: {"status": r.get("status"), "fail_reason": r.get("fail_reason"), "source": r.get("source"),
                          "accession": r.get("accession"), "bio_unit": r.get("bio_unit"), "reads": r.get("reads"),
                          "read_length": r.get("read_length"), "pilot_max_reads": r.get("pilot_max_reads"),
                          "original_files": r.get("original"), "trimmed": bool(r.get("trimmed"))}
                    for sid, r in p.samples.items()},
        "read_type": s.get("read_type"), "read_length": first.get("read_length"),
        "strandedness": {k: v for k, v in (s.get("strandedness") or {}).items() if k != "evidence"},
        "qc_decision": s.get("qc_decision"), "design": s.get("design"),
        "threads": ctx.threads, "max_memory_gb": ctx.mem_gb,
        "parameters": redacted(cfg),
        "deseq2_settings": {k: cfg.get(k) for k in ("alpha", "log2fc_threshold", "lfc_test_threshold",
                                                     "min_count_filter", "independent_filtering", "cooks_cutoff",
                                                     "lfc_shrinkage", "transformation", "design_formula")},
        "checksums_sha256": checksums,
        "checkpoints": {st.stem: (ctx.cp.read(st.stem) or {}).get("completed_at")
                        for st in ctx.cp.dir.glob("*.json") if ".invalid." not in st.name},
        "history": s.get("history", [])[-200:],
    }
    if any(r.get("pilot_max_reads") for r in p.samples.values()):
        manifest["WARNING"] = "PILOT MODE: only a subset of reads was downloaded; results are for testing only"
    j = out / "manifest.json"
    j.write_text(json.dumps(manifest, indent=2, default=str))
    y = out / "manifest.yaml"
    C.save_yaml(json.loads(json.dumps(manifest, default=str)), y)
    cfg_used = C.snapshot(cfg, p.path("config"), {"threads": ctx.threads, "max_memory_gb": ctx.mem_gb})
    ui.ok(f"manifest written: {p.rel(j)}")
    return [j, y, vers, cfg_used, *env_files]
