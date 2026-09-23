"""`rnaseq-pipeline --validate-project DIR`: non-interactive project health report (PASS / WARNING / FAIL)."""
import re

from . import environment_manager, report_manager, system_check, ui, workflow
from .project import Project

PASS, WARN, FAIL = "PASS", "WARNING", "FAIL"
_CREDENTIAL_RE = re.compile(r"api_key:\s*['\"]?[A-Za-z0-9]{12,}")


# Which health area each stage belongs to (section 47 of the stabilization plan).
AREA_OF = {"data_acquired": "DATA", "fastq_verified": "DATA", "raw_qc_completed": "QC", "quality_assessed": "QC",
           "trimming_completed": "QC", "reference_ready": "REFERENCE", "alignment_completed": "ALIGNMENT",
           "bam_qc_completed": "ALIGNMENT", "strandedness_determined": "ALIGNMENT",
           "stringtie_completed": "QUANTIFICATION", "featurecounts_completed": "QUANTIFICATION",
           "count_matrix_completed": "QUANTIFICATION", "design_confirmed": "DESIGN",
           "deseq2_completed": "DESEQ2 + RESULTS", "report_generated": "REPORT"}
AREAS = ["SYSTEM", "INSTALLATION", "CONFIGURATION", "DATA", "REFERENCE", "QC", "ALIGNMENT", "QUANTIFICATION",
         "DESIGN", "DESEQ2 + RESULTS", "REPORT", "SECURITY", "RESUME", "SCIENTIFIC VALIDATION"]
SCIENTIFIC = ["DATA", "REFERENCE", "QC", "ALIGNMENT", "QUANTIFICATION", "DESIGN", "DESEQ2 + RESULTS", "REPORT"]
_RANK = {PASS: 0, WARN: 1, FAIL: 2}


def _worst(states):
    return max(states, key=_RANK.get) if states else PASS


def check(root, deep=True):
    """Return (overall, rows, first_pending); rows are (area, item, status, detail), grouped by AREAS."""
    from . import config as C
    from . import dependency_manager
    rows = []
    p = Project.open(root)                       # also removes credentials older versions stored
    cfg = p.config()
    envs = environment_manager.Environments(cfg)
    envs.activate()
    ctx = workflow.Context(p, cfg, envs, system_check.collect(p.root))
    si = ctx.sysinfo
    warn_gb = (cfg.get("storage") or {}).get("min_free_gb_warn", 50)
    rows.append(("SYSTEM", "free space for this project", WARN if si["disk_free_gb"] < warn_gb else PASS,
                 f"{si['disk_free_gb']} GB free on {si['filesystem']} ({si['mountpoint']})"))
    # Missing or unusable tools do not make existing results invalid; they only matter for continuing.
    tools = dependency_manager.detect_tools()
    bad_tools = [f"{t['key']} ({t['status'].lower()})" for t in tools
                 if t["required"] and t["status"] not in ("AVAILABLE", "OPTIONAL")]
    rows.append(("INSTALLATION", "scientific tools", WARN if bad_tools else PASS,
                 f"needed to continue this project: {', '.join(bad_tools)}" if bad_tools else "all usable"))
    r_info, r_pkgs = dependency_manager.detect_r(envs.rscript())
    bad_r = ([] if r_info["status"] == "AVAILABLE" else ["R"]) + \
        [x["package"] for x in r_pkgs if x["required"] and x["status"] != "AVAILABLE"]
    rows.append(("INSTALLATION", "R + DESeq2 packages", WARN if bad_r else PASS,
                 f"missing: {', '.join(bad_r)}" if bad_r else f"R {r_info.get('version')}"))
    errs = C.validate(cfg, ctx.sysinfo["cpu_cores"])
    rows.append(("CONFIGURATION", "project_config.yaml", FAIL if errs else PASS, "; ".join(errs) or "valid"))
    memo, first_pending, first_invalid = {}, None, None
    for i, st in enumerate(workflow.STAGES, 1):
        status, why = workflow.stage_status(ctx, st, deep=deep, memo=memo)
        item = f"STEP {i}: {st.title}"
        if status == "VALID":
            rows.append((AREA_OF[st.key], item, PASS, "outputs re-verified"))
        elif status == "PENDING":
            first_pending = first_pending or st.title
            rows.append((AREA_OF[st.key], item, WARN, "not run yet"))
        else:
            first_invalid = first_invalid or st.title
            rows.append((AREA_OF[st.key], item, FAIL, why or "invalid"))
    bad = {s: r for s, r in p.samples.items() if r.get("status") in ("FAILED", "EXCLUDED")}
    rows.append(("DATA", f"{len(p.samples)} samples registered", WARN if bad else PASS,
                 f"{len(bad)} failed/excluded: {', '.join(list(bad)[:5])}" if bad else "all usable"))
    rep = p.path("reports", "final_pipeline_report.html")
    if rep.exists():
        probs = report_manager.validate(ctx)
        rows.append(("REPORT", rep.name, FAIL if probs else PASS, "; ".join(probs[:3]) or
                     "links resolve; key numbers (samples, input reads, alignment rate, mapped reads/%, assigned and counted reads, DE counts, thresholds, versions) match their source files"))
    if p.path("pipeline_manifest", "manifest.json").exists():
        from . import manifest
        probs = manifest.problems(p)
        rows.append(("REPORT", "reproducibility manifest", FAIL if probs else PASS, "; ".join(probs) or "complete"))
    leaked = [str(f.relative_to(p.root)) for f in p.root.rglob("*.yaml")
              if f.is_file() and _CREDENTIAL_RE.search(f.read_text(errors="replace"))]
    leaked += [str(f.relative_to(p.root)) for f in (p.root / "pipeline_manifest").glob("*.json")
               if _CREDENTIAL_RE.search(f.read_text(errors="replace").replace('"', "").replace(",", "\n"))]
    rows.append(("SECURITY", "stored credentials", FAIL if leaked else PASS,
                 f"found in {', '.join(leaked[:3])}" if leaked else "none found in project files"))
    rows.append(("RESUME", "what resuming would do", PASS,
                 f"re-run from: {first_invalid}" if first_invalid else
                 (f"continue with: {first_pending}" if first_pending else "nothing to re-run; every step is valid")))
    sci = _worst([r[2] for r in rows if r[0] in SCIENTIFIC])
    rows.append(("SCIENTIFIC VALIDATION", "independent cross-checks of all steps", sci,
                 {PASS: "every step's outputs and cross-checks re-verified",
                  WARN: "not every step has run yet (or samples were excluded)",
                  FAIL: "at least one scientific step failed re-verification (see BLOCKING ISSUES)"}[sci]))
    rows.sort(key=lambda r: AREAS.index(r[0]))
    overall = _worst([r[2] for r in rows])
    return overall, rows, first_pending


def show(root, overall, rows, first_pending):
    ui.header(f"PROJECT HEALTH — {root}")
    area = None
    for a, item, status, detail in rows:
        if a != area:
            ui.section(a)
            area = a
        ui.status({"PASS": "OK", "WARNING": "WARNING", "FAIL": "FAILED"}[status], f"{item:<48} {detail[:110]}")
    ui.section("AREA SUMMARY")
    for a in AREAS:
        st = _worst([r[2] for r in rows if r[0] == a])
        print(f"  {a:<24} {st}")
    blocking = [f"{a} — {item}: {detail}" for a, item, st, detail in rows if st == FAIL and a != "SCIENTIFIC VALIDATION"]
    nonblocking = [f"{a} — {item}: {detail}" for a, item, st, detail in rows if st == WARN and a != "SCIENTIFIC VALIDATION"]
    print("\nBLOCKING ISSUES: " + ("none" if not blocking else ""))
    for x in blocking:
        print(f"  - {x[:200]}")
    print("NON-BLOCKING ISSUES: " + ("none" if not nonblocking else ""))
    for x in nonblocking:
        print(f"  - {x[:200]}")
    print()
    ui.rule("=")
    print(f"PROJECT HEALTH: {overall}")
    if overall == FAIL:
        print("  Stages marked FAILED will be re-run when you resume the project (rnaseq-pipeline --project DIR).")
    elif first_pending:
        print(f"  The analysis is not finished; next step: {first_pending}.")
    ui.rule("=")
