"""Experimental design / sample metadata: interactive editing, validation, explicit confirmation.

Groups are NEVER inferred automatically. Public metadata is shown only as
reference; the user must assign (or explicitly choose a field for) every value.
"""
import csv
import hashlib
import json
import re
from datetime import datetime
from fractions import Fraction
from pathlib import Path

from . import PipelineError, ui
from . import validators as V

BASE_COLUMNS = ["sample", "condition", "replicate"]


def new_table(samples):
    return {"columns": ["condition"], "rows": {s: {"condition": ""} for s in samples}}


def with_replicates(table, var="condition"):
    counters = {}
    out = {}
    for s, row in table["rows"].items():
        lvl = row.get(var, "")
        counters[lvl] = counters.get(lvl, 0) + 1
        out[s] = counters[lvl]
    return out


def model_rank(table, variables):
    """Rank of the additive model matrix (intercept + treatment-coded factors) and its column count."""
    samples = list(table["rows"])
    cols = [[Fraction(1)] * len(samples)]
    for v in variables:
        levels = sorted({table["rows"][s][v] for s in samples})
        for lvl in levels[1:]:
            cols.append([Fraction(1 if table["rows"][s][v] == lvl else 0) for s in samples])
    rows = [list(r) for r in zip(*cols)]  # samples x coefficients
    m, n = len(rows), len(cols)
    rank, r = 0, 0
    for c in range(n):
        piv = next((i for i in range(r, m) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(m):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
        rank += 1
        if r == m:
            break
    return rank, n


def validate(table, formula, reference_level, contrasts):
    """Return list of problems (empty = valid)."""
    problems = []
    try:
        variables = V.design_formula(formula, table["columns"])
    except ValueError as e:
        return [str(e)]
    rows = table["rows"]
    if len(rows) < 2:
        return ["at least 2 samples are required"]
    for s, row in rows.items():
        for v in variables:
            val = row.get(v, "")
            if not val:
                problems.append(f"sample {s}: '{v}' is empty")
                continue
            try:
                V.factor_level(val)
            except ValueError as e:
                problems.append(f"sample {s}: {e}")
    if problems:
        return problems
    interest = variables[-1]
    levels = {}
    for row in rows.values():
        levels[row[interest]] = levels.get(row[interest], 0) + 1
    if len(levels) < 2:
        problems.append(f"'{interest}' has only one level ({list(levels)}); nothing to compare")
    single = [l for l, n in levels.items() if n < 2]
    if single:
        problems.append(f"levels without biological replicates (n<2): {single}. DESeq2 cannot estimate "
                        "within-group variability reliably; add replicates or drop these samples")
    for v in variables[:-1]:
        lv = {row[v] for row in rows.values()}
        if len(lv) < 2:
            problems.append(f"covariate '{v}' has a single level; remove it from the formula")
    if not problems:
        rank, ncoef = model_rank(table, variables)
        if rank < ncoef:
            problems.append(f"design is not full rank ({rank} < {ncoef} coefficients): "
                            f"'{interest}' is confounded with another variable (e.g. every batch contains "
                            "only one condition). The effects cannot be separated")
        elif ncoef >= len(rows):
            problems.append("no residual degrees of freedom (as many coefficients as samples)")
    if reference_level and reference_level not in levels:
        problems.append(f"reference level '{reference_level}' is not a level of '{interest}'")
    for num, den in contrasts or []:
        if num not in levels or den not in levels or num == den:
            problems.append(f"invalid contrast {num} vs {den}")
    return problems


def default_contrasts(levels, reference):
    return [(l, reference) for l in sorted(levels) if l != reference]


def show(table, info_cols=None, info=None):
    reps = with_replicates(table)
    headers = ["Sample", *[c.capitalize() for c in table["columns"]], "Replicate"]
    info_cols = info_cols or []
    rows = []
    for s, row in table["rows"].items():
        extra = [str((info or {}).get(s, {}).get(c, ""))[:30] for c in info_cols]
        rows.append([s, *[row.get(c, "") or "—" for c in table["columns"]], reps[s] if row.get("condition") else "",
                     *extra])
    ui.table(rows, headers + [f"(info) {c}" for c in info_cols], max_col=32)


def public_fields(info):
    """Metadata fields with >1 distinct value across samples (candidates the user may choose from)."""
    fields = {}
    for s, meta in info.items():
        for k, v in meta.items():
            if isinstance(v, str) and v:
                fields.setdefault(k, set()).add(v)
    return {k: v for k, v in fields.items() if 1 < len(v) < len(info) + 1 and not k.endswith("accession")}


def sanitize_level(v):
    v = re.sub(r"[^A-Za-z0-9_.]", "_", v.strip())
    if not v or not v[0].isalpha():
        v = "g_" + v
    return v[:64]


def import_tsv(table, path):
    path = V.existing_file(path)
    with open(path, newline="") as f:
        dialect = "excel-tab" if "\t" in f.readline() else "excel"
        f.seek(0)
        r = csv.DictReader(f, dialect=dialect)
        if not r.fieldnames or r.fieldnames[0].lower() not in ("sample", "sample_id", "run"):
            raise PipelineError("first column of the metadata file must be 'sample'")
        key = r.fieldnames[0]
        cols = [V.column_name(c) for c in r.fieldnames[1:] if c.lower() != "replicate"]
        seen = set()
        for rec in r:
            s = rec[key].strip()
            if s not in table["rows"]:
                raise PipelineError(f"metadata file lists unknown sample {s!r}")
            seen.add(s)
            for c in cols:
                table["rows"][s][c] = V.factor_level(rec[c]) if rec[c] else ""
        missing = set(table["rows"]) - seen
        if missing:
            raise PipelineError(f"metadata file is missing samples: {sorted(missing)}")
    for c in cols:
        if c not in table["columns"]:
            table["columns"].append(c)


def write_metadata(table, path):
    reps = with_replicates(table)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", *table["columns"], "replicate"])
        for s, row in table["rows"].items():
            w.writerow([s, *[row.get(c, "") for c in table["columns"]], reps[s]])
    return path


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_confirmed(project, table, formula, reference, contrasts):
    meta_path = project.path("counts", "sample_metadata.tsv")
    write_metadata(table, meta_path)
    variables = V.design_formula(formula, table["columns"])
    d = {"confirmed": True, "confirmed_at": datetime.now().isoformat(timespec="seconds"),
         "formula": formula, "variables": variables, "variable_of_interest": variables[-1],
         "reference_level": reference, "contrasts": [list(c) for c in contrasts],
         "metadata_file": project.rel(meta_path), "metadata_sha256": file_sha256(meta_path),
         "samples": list(table["rows"])}
    dpath = project.path("counts", "design.json")
    dpath.write_text(json.dumps(d, indent=2))
    project.state["design"] = d
    project.save()
    return meta_path, dpath


def confirmed_design(project):
    dpath = project.path("counts", "design.json")
    if not dpath.exists():
        return None
    d = json.loads(dpath.read_text())
    mp = project.abs(d["metadata_file"])
    if not d.get("confirmed") or not mp.exists() or file_sha256(mp) != d["metadata_sha256"]:
        return None
    return d


def interactive(project, samples, info, cfg):
    """Edit + confirm the design. Returns True if confirmed."""
    draft = project.state.get("design_draft")
    table = draft if draft and set(draft["rows"]) == set(samples) else new_table(samples)
    formula = (project.state.get("design") or {}).get("formula") or cfg.get("design_formula", "~ condition")
    reference = None
    contrasts = None
    fields = public_fields(info)
    info_cols = [c for c in ("sample_title", "geo_title", "geo_treatment", "geo_source_name") if c in fields][:2]

    def save_draft():
        project.state["design_draft"] = table
        project.save()

    while True:
        ui.header("EXPERIMENTAL DESIGN / SAMPLE METADATA",
                  "Groups are never guessed — please assign or confirm every value")
        show(table, info_cols, info)
        ui.kv([("Design formula", formula), ("Reference level", reference or "(not set)"),
               ("Contrasts", ", ".join(f"{a} vs {b}" for a, b in contrasts) if contrasts else "(default: each vs reference)")])
        opts = ["Assign condition to each sample (one by one)",
                "Assign condition by pattern (sample name / description contains text)",
                "Use a public metadata field as condition (you choose the field)",
                "Import metadata from a TSV/CSV file (sample, condition, batch, ...)",
                "Add or edit a covariate column (batch, sex, tissue, genotype, time, ...)",
                "Set design formula (e.g. ~ condition  or  ~ batch + condition)",
                "Set reference level and contrasts",
                "Validate and CONFIRM design",
                "Save draft and return"]
        c = ui.choose(None, opts)
        try:
            if c == 0:
                for s in table["rows"]:
                    hint = " | ".join(str(info.get(s, {}).get(k, ""))[:40] for k in info_cols)
                    cur = table["rows"][s].get("condition", "")
                    table["rows"][s]["condition"] = ui.ask(f"  {s}{'  (' + hint + ')' if hint else ''} condition",
                                                           default=cur or None, validator=V.factor_level)
            elif c == 1:
                text = ui.ask("  Text to match (case-insensitive) in sample name or description")
                level = ui.ask("  Condition to assign", validator=V.factor_level)
                hits = [s for s in table["rows"] if text.lower() in s.lower()
                        or any(text.lower() in str(v).lower() for v in info.get(s, {}).values())]
                ui.info(f"matches: {', '.join(hits) or 'none'}")
                if hits and ui.ask_yes_no(f"  Assign '{level}' to these {len(hits)} sample(s)?", True):
                    for s in hits:
                        table["rows"][s]["condition"] = level
            elif c == 2:
                if not fields:
                    ui.warn("no public metadata fields with varying values are available")
                    continue
                names = sorted(fields)
                i = ui.choose("Public metadata fields (values shown are as submitted; verify them)",
                              [f"{n}: {', '.join(sorted(fields[n]))[:70]}" for n in names])
                col = ui.ask("  Store into column", default="condition", validator=V.column_name)
                mapping = {}
                for val in sorted(fields[names[i]]):
                    mapping[val] = ui.ask(f"  Level name for '{val}'", default=sanitize_level(val),
                                          validator=V.factor_level)
                for s in table["rows"]:
                    v = info.get(s, {}).get(names[i], "")
                    table["rows"][s][col] = mapping.get(v, "")
                if col not in table["columns"]:
                    table["columns"].append(col)
            elif c == 3:
                import_tsv(table, ui.ask("  Path to metadata file"))
                ui.ok("metadata imported")
            elif c == 4:
                col = ui.ask("  Column name (e.g. batch, sex, tissue)", validator=V.column_name)
                if col not in table["columns"]:
                    table["columns"].append(col)
                for s in table["rows"]:
                    table["rows"][s][col] = ui.ask(f"  {s} {col}", default=table["rows"][s].get(col) or None,
                                                   validator=V.factor_level)
            elif c == 5:
                formula = ui.ask("  Design formula", default=formula,
                                 validator=lambda f: (V.design_formula(f, table["columns"]), f.strip())[1])
            elif c == 6:
                var = V.design_formula(formula, table["columns"])[-1]
                levels = sorted({r[var] for r in table["rows"].values() if r.get(var)})
                if len(levels) < 2:
                    ui.warn(f"assign at least two levels of '{var}' first")
                    continue
                i = ui.choose(f"Reference (baseline) level of '{var}'", levels)
                reference = levels[i]
                if ui.ask_yes_no("  Compare every other level against the reference?", True):
                    contrasts = default_contrasts(levels, reference)
                else:
                    contrasts = []
                    while True:
                        num = ui.choose("Numerator (test) level", levels)
                        den = ui.choose("Denominator (baseline) level", levels)
                        contrasts.append((levels[num], levels[den]))
                        if not ui.ask_yes_no("  Add another contrast?", False):
                            break
            elif c == 7:
                var = V.design_formula(formula, table["columns"])[-1]
                levels = sorted({r[var] for r in table["rows"].values() if r.get(var)})
                ref = reference
                cfg_ref = cfg.get("reference_level")
                if ref is None and cfg_ref not in (None, "ask") and cfg_ref in levels:
                    ref = cfg_ref
                if ref is None:
                    ui.warn("set the reference level first (option 7)")
                    continue
                cons = contrasts if contrasts is not None else default_contrasts(levels, ref)
                problems = validate(table, formula, ref, cons)
                if problems:
                    ui.error("design is not valid:")
                    for p in problems:
                        print(f"   - {p}")
                    continue
                ui.ok("design is valid")
                show(table)
                ui.kv([("Formula", formula), ("Variable of interest", var), ("Reference", ref),
                       ("Contrasts", ", ".join(f"{a} vs {b}" for a, b in cons))])
                if ui.ask_yes_no("CONFIRM this experimental design for DESeq2?", False):
                    save_confirmed(project, table, formula, ref, cons)
                    project.state.pop("design_draft", None)
                    project.save()
                    ui.ok("design confirmed and saved (counts/sample_metadata.tsv, counts/design.json)")
                    return True
            elif c == 8:
                save_draft()
                return False
            save_draft()
        except (ValueError, PipelineError) as e:
            ui.error(str(e))
