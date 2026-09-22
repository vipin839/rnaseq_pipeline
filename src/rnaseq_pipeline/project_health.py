"""`rnaseq-pipeline --validate-project DIR`: non-interactive project health report (PASS / WARNING / FAIL)."""
import re

from . import environment_manager, report_manager, system_check, ui, workflow
from .project import Project

PASS, WARN, FAIL = "PASS", "WARNING", "FAIL"
_CREDENTIAL_RE = re.compile(r"api_key:\s*['\"]?[A-Za-z0-9]{12,}")


def check(root, deep=True):
    """Return (overall, rows) where rows are (area, item, status, detail)."""
    rows = []
    p = Project.open(root)                       # also removes credentials older versions stored
    cfg = p.config()
    envs = environment_manager.Environments(cfg)
    envs.activate()
    ctx = workflow.Context(p, cfg, envs, system_check.collect(p.root))
    from . import config as C
    errs = C.validate(cfg, ctx.sysinfo["cpu_cores"])
    rows.append(("Configuration", "project_config.yaml", FAIL if errs else PASS, "; ".join(errs) or "valid"))
    memo, first_pending = {}, None
    for i, st in enumerate(workflow.STAGES, 1):
        status, why = workflow.stage_status(ctx, st, deep=deep, memo=memo)
        if status == "VALID":
            rows.append(("Stages", f"{i:>2}. {st.title}", PASS, "outputs re-verified"))
        elif status == "PENDING":
            first_pending = first_pending or st.title
            rows.append(("Stages", f"{i:>2}. {st.title}", WARN, "not run yet"))
        else:
            rows.append(("Stages", f"{i:>2}. {st.title}", FAIL, why or "invalid"))
    bad = {s: r for s, r in p.samples.items() if r.get("status") in ("FAILED", "EXCLUDED")}
    rows.append(("Samples", f"{len(p.samples)} registered", WARN if bad else PASS,
                 f"{len(bad)} failed/excluded: {', '.join(list(bad)[:5])}" if bad else "all usable"))
    rep = p.path("reports", "final_pipeline_report.html")
    if rep.exists():
        probs = report_manager.validate(ctx)
        rows.append(("Report", rep.name, FAIL if probs else PASS, "; ".join(probs[:3]) or
                     "every link resolves, every number matches the result files"))
    if p.path("pipeline_manifest", "manifest.json").exists():
        from . import manifest
        probs = manifest.problems(p)
        rows.append(("Reproducibility", "manifest.json", FAIL if probs else PASS, "; ".join(probs) or "complete"))
    leaked = [str(f.relative_to(p.root)) for f in p.root.rglob("*.yaml")
              if f.is_file() and _CREDENTIAL_RE.search(f.read_text(errors="replace"))]
    leaked += [str(f.relative_to(p.root)) for f in (p.root / "pipeline_manifest").glob("*.json")
               if _CREDENTIAL_RE.search(f.read_text(errors="replace").replace('"', "").replace(",", "\n"))]
    rows.append(("Security", "stored credentials", FAIL if leaked else PASS,
                 f"found in {', '.join(leaked[:3])}" if leaked else "none found in project files"))
    states = {r[2] for r in rows}
    overall = FAIL if FAIL in states else (WARN if WARN in states else PASS)
    return overall, rows, first_pending


def show(root, overall, rows, first_pending):
    ui.header(f"PROJECT HEALTH — {root}")
    area = None
    for a, item, status, detail in rows:
        if a != area:
            ui.section(a.upper())
            area = a
        ui.status({"PASS": "OK", "WARNING": "WARNING", "FAIL": "FAILED"}[status], f"{item:<48} {detail[:110]}")
    print()
    ui.rule("=")
    print(f"PROJECT HEALTH: {overall}")
    if overall == FAIL:
        print("  Stages marked FAILED will be re-run when you resume the project (rnaseq-pipeline --project DIR).")
    elif first_pending:
        print(f"  The analysis is not finished; next step: {first_pending}.")
    ui.rule("=")
