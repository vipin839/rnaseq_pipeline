"""Data sources: local FASTQ import and public accession metadata (SRA / ENA / GEO)."""
import os
import re
import shutil
from pathlib import Path

from . import PipelineError, ena_manager, geo_manager, sra_manager, ui
from . import validators as V

FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
PAIR_PATTERNS = [
    re.compile(r"^(?P<base>.+?)[_.]R(?P<mate>[12])(?P<suffix>_\d{3})?\.(fastq|fq)(\.gz)?$", re.I),
    re.compile(r"^(?P<base>.+?)[_.](?P<mate>[12])\.(fastq|fq)(\.gz)?$", re.I),
]
LANE_RE = re.compile(r"_L\d{3}$")
ILLUMINA_S_RE = re.compile(r"_S\d+$")


def sanitize_sample_name(name):
    n = re.sub(r"[^A-Za-z0-9._-]", "_", name).strip("._-")
    if not n or not n[0].isalpha():
        n = "S_" + n
    return n[:64]


def find_fastqs(paths):
    files = []
    for p in paths:
        p = Path(os.path.expanduser(p)).resolve()
        if p.is_dir():
            files += sorted(f for f in p.iterdir() if f.is_file() and FASTQ_EXT.search(f.name))
        elif p.is_file() and FASTQ_EXT.search(p.name):
            files.append(p)
        elif p.exists():
            raise PipelineError(f"not a FASTQ file (.fastq/.fq[.gz]): {p}")
        else:
            raise PipelineError(f"path not found: {p}")
    if not files:
        raise PipelineError("no FASTQ files found", remedy="files must end in .fastq, .fq, .fastq.gz or .fq.gz")
    return files


def pair_files(files):
    """Group files into samples. Returns list of dicts {id, r1, r2, bio_unit}; raises on ambiguity."""
    groups = {}
    singles = []
    for f in files:
        m = None
        for pat in PAIR_PATTERNS:
            m = pat.match(f.name)
            if m:
                break
        if m:
            key = m.group("base")
            groups.setdefault(key, {})
            mate = "r" + m.group("mate")
            if mate in groups[key]:
                raise PipelineError(f"two files map to the same sample/mate: {groups[key][mate].name}, {f.name}")
            groups[key][mate] = f
        else:
            singles.append(f)
    samples = []
    for base, g in groups.items():
        if "r1" in g and "r2" in g:
            samples.append({"raw_name": base, "r1": g["r1"], "r2": g["r2"]})
        elif "r1" in g:
            samples.append({"raw_name": base, "r1": g["r1"], "r2": None})
        else:
            raise PipelineError(f"R2 file without matching R1: {g['r2'].name}",
                                remedy="provide both mates, or rename files")
    for f in singles:
        samples.append({"raw_name": FASTQ_EXT.sub("", f.name), "r1": f, "r2": None})
    layouts = {s["r2"] is not None for s in samples}
    if len(layouts) > 1:
        raise PipelineError("mixture of paired-end and single-end files",
                            cause="Version 1 requires one library layout per project",
                            remedy="analyse paired and single-end samples in separate projects")
    seen = set()
    for s in sorted(samples, key=lambda s: s["raw_name"]):
        sid = sanitize_sample_name(s["raw_name"])
        base, i = sid, 2
        while sid in seen:
            sid = f"{base}_{i}"
            i += 1
        seen.add(sid)
        s["id"] = sid
        bio = LANE_RE.sub("", s["raw_name"])
        s["bio_unit"] = sanitize_sample_name(ILLUMINA_S_RE.sub("", bio))
    return sorted(samples, key=lambda s: s["id"])


def import_local(project, samples, mode="symlink"):
    """Link (default) or copy FASTQs into data/fastq/. Original files are never modified."""
    fq_dir = project.path("data", "fastq")
    for s in samples:
        rec = project.add_sample(s["id"], source="local", layout="PAIRED" if s["r2"] else "SINGLE",
                                 bio_unit=s.get("bio_unit", s["id"]),
                                 original={"r1": str(s["r1"]), "r2": str(s["r2"]) if s["r2"] else None})
        for mate, src in (("r1", s["r1"]), ("r2", s["r2"])):
            if src is None:
                rec[mate] = None
                continue
            ext = ".fastq.gz" if src.name.endswith(".gz") else ".fastq"
            dest = fq_dir / f"{s['id']}_{mate.upper()}{ext}"
            if dest.exists() or dest.is_symlink():
                if dest.resolve() != Path(src).resolve():
                    raise PipelineError(f"{dest} already exists and points elsewhere; refusing to overwrite")
            elif mode == "copy":
                tmp = dest.with_name(dest.name + ".part")
                shutil.copy2(src, tmp)
                os.replace(tmp, dest)
            else:
                dest.symlink_to(Path(src).resolve())
            rec[mate] = project.rel(dest) if mode == "copy" else str(Path("data/fastq") / dest.name)
    project.state["read_type"] = "paired" if samples[0]["r2"] else "single"
    project.save()


def existing_project_fastqs(project):
    """Discover FASTQ files already placed in <project>/data/fastq."""
    return pair_files(find_fastqs([project.path("data", "fastq")]))


# ------------------------------ public data ------------------------------

SC_WORDS = re.compile(r"(10x genomics|chromium (next gem )?single cell|single[- ]cell 3'|single[- ]cell 5'|"
                      r"drop-?seq|indrop|scrna-?seq|snrna-?seq|single[- ]nucle(us|i) rna|cel-?seq|"
                      r"smart-?seq|split-?seq|sci-rna)", re.I)


def single_cell_samples(project, sample_ids):
    """{sample: reason} for registered samples whose stored metadata says single-cell."""
    out = {}
    for sid in sample_ids:
        rec = project.samples[sid]
        reason = single_cell_reason({**(rec.get("metadata") or {})})
        if reason:
            out[sid] = reason
    return out


def single_cell_reason(run):
    """Why a run looks like single-cell RNA-seq (None if it looks like bulk)."""
    src = (run.get("library_source") or "").upper()
    if "SINGLE CELL" in src:
        return f"library source is '{run['library_source']}'"
    for field in ("library_construction_protocol", "experiment_title", "study_title", "library_name"):
        m = SC_WORDS.search(run.get(field) or "")
        if m:
            return f"{field.replace('_', ' ')} mentions '{m.group(0)}'"
    return None

def _runs_for_experiments(srx_list):
    """Runs for many SRA experiments: query each parent study once instead of every experiment."""
    wanted, runs, done_studies = set(srx_list), [], set()
    remaining = list(srx_list)
    while remaining:
        probe = remaining[0]
        got = ena_manager.runs_for(probe)
        study = next((r.get("study_accession") for r in got if r.get("study_accession")), None)
        if study and study not in done_studies and len(remaining) > 3:
            done_studies.add(study)
            batch = [r for r in ena_manager.runs_for(study) if r.get("experiment_accession") in wanted]
            if batch:
                runs += batch
                covered = {r["experiment_accession"] for r in batch}
                remaining = [x for x in remaining if x not in covered]
                continue
        runs += got
        remaining = remaining[1:]
    missing = wanted - {r.get("experiment_accession") for r in runs}
    for srx in sorted(missing):  # not in ENA (yet): ask NCBI directly
        runs += [_runinfo_record(i) for i in sra_manager.runinfo(srx)]
    still = wanted - {r.get("experiment_accession") for r in runs}
    if still:
        ui.warn(f"{len(still)} GEO sample(s) have no retrievable sequencing runs and were skipped: "
                f"{', '.join(sorted(still)[:5])}{'...' if len(still) > 5 else ''}")
    return runs


def _runinfo_record(i):
    """NCBI runinfo row -> the ENA-style run dict used everywhere else (no direct FASTQ URLs: SRA route)."""
    return {"run_accession": i["Run"], "experiment_accession": i.get("Experiment", ""),
            "sample_accession": i.get("BioSample", ""), "study_accession": i.get("BioProject", ""),
            "library_layout": i.get("LibraryLayout", ""), "library_strategy": i.get("LibraryStrategy", ""),
            "library_source": i.get("LibrarySource", ""), "instrument_model": i.get("Model", ""),
            "read_count": i.get("spots", ""), "scientific_name": i.get("ScientificName", ""),
            "tax_id": i.get("TaxID", ""), "sample_title": i.get("SampleName", ""),
            "fastq_bytes": "", "_files": []}


def resolve_public(accessions, source):
    """Resolve accessions to run records with metadata. Returns (runs, study_titles)."""
    runs, titles, gsm_meta = [], set(), {}
    for acc, kind in accessions:
        if kind == "geo_series":
            title, gsms, notes = geo_manager.series(acc)
            for n in notes:
                ui.info(n)
            titles.add(title)
            srx_to_gsm = {}
            for g in gsms:
                gsm_meta[g["gsm"]] = g
                for srx in g["srx"]:
                    srx_to_gsm[srx] = g["gsm"]
            if not srx_to_gsm:
                raise PipelineError(f"{acc}: none of its {len(gsms)} sample(s) link to SRA raw data",
                                    cause="the series may be microarray, or raw reads are held elsewhere (e.g. dbGaP)")
            ui.running(f"Looking up sequencing runs for {len(srx_to_gsm)} GEO sample(s)...")
            for r in _runs_for_experiments(list(srx_to_gsm)):
                r["_gsm"] = srx_to_gsm.get(r.get("experiment_accession"))
                runs.append(r)
        elif kind == "geo_sample":
            raise PipelineError("single GSM accessions are not supported; enter the GSE series or the SRR runs")
        else:
            got = ena_manager.runs_for(acc)
            if not got and source == "sra":
                got = [_runinfo_record(i) for i in sra_manager.runinfo(acc)]
            if not got:
                raise PipelineError(f"no sequencing runs found for {acc}",
                                    remedy="check that the accession is public and correct")
            runs += got
    for r in runs:
        titles.add(r.get("study_title", "")) if r.get("study_title") else None
        gsm = r.get("_gsm") or (r.get("sample_alias") if str(r.get("sample_alias", "")).startswith("GSM") else None)
        if gsm and gsm in gsm_meta:
            g = gsm_meta[gsm]
            r["_geo"] = {"gsm": gsm, "title": g.get("title", ""), "source_name": g.get("source_name", ""),
                         **g.get("characteristics", {})}
    uniq, seen = [], set()
    for r in runs:
        if r["run_accession"] not in seen:
            seen.add(r["run_accession"])
            uniq.append(r)
    return uniq, sorted(t for t in titles if t)


def show_public(runs, titles):
    ui.section("PUBLIC DATA — retrieved metadata (review carefully; not verified)")
    for t in titles:
        print(f"  Study: {t}")
    orgs = sorted({r.get("scientific_name", "?") for r in runs})
    layouts = sorted({r.get("library_layout", "?") for r in runs})
    strategies = sorted({r.get("library_strategy", "?") for r in runs})
    ui.kv([("Organism(s)", ", ".join(orgs)), ("Library layout", ", ".join(layouts)),
           ("Library strategy", ", ".join(strategies)), ("Runs", len(runs))], indent=2)
    print()
    rows = []
    for r in runs:
        desc = r.get("_geo", {}).get("title") or r.get("sample_title") or r.get("experiment_title", "")
        size = sum(f[2] or 0 for f in r.get("_files", [])) / 1e9
        rows.append((r["run_accession"], r.get("experiment_accession", ""), r.get("sample_accession", ""),
                     r.get("library_layout", ""), r.get("instrument_model", ""),
                     f"{int(r['read_count']):,}" if str(r.get("read_count", "")).isdigit() else "?",
                     f"{size:.2f}" if size else "?", desc))
    ui.table(rows, ["Run", "Experiment", "Sample", "Layout", "Instrument", "Reads", "GB", "Description"])
    warnings = []
    if len(orgs) > 1:
        warnings.append("runs come from more than one organism")
    if len(layouts) > 1:
        warnings.append("mixture of paired and single-end runs (not supported in one project)")
    if any(s not in ("RNA-Seq", "?") for s in strategies):
        warnings.append(f"library strategy is not RNA-Seq for some runs: {strategies}")
    sc = [r for r in runs if single_cell_reason(r)]
    if sc:
        warnings.append(f"{len(sc)} of {len(runs)} run(s) are SINGLE-CELL RNA-seq "
                        f"({single_cell_reason(sc[0])}); this bulk pipeline cannot analyse them")
    exps = [r.get("experiment_accession") for r in runs]
    if len(set(exps)) < len(exps):
        warnings.append("some experiments have several runs (lanes/technical replicates); they can be "
                        "summed per biological sample at the count-matrix stage")
    for w in warnings:
        ui.warn(w)
    return warnings


def register_public(project, runs, source):
    layouts = {r.get("library_layout") for r in runs}
    if len(layouts) != 1 or layouts.pop() not in ("PAIRED", "SINGLE"):
        raise PipelineError("runs must all be PAIRED or all SINGLE", remedy="select a consistent subset of runs")
    for r in runs:
        meta = {k: r.get(k, "") for k in ("experiment_accession", "sample_accession", "study_accession",
                                          "library_strategy", "library_source", "library_selection",
                                          "library_construction_protocol",
                                          "instrument_model", "scientific_name", "tax_id", "sample_title",
                                          "sample_alias", "experiment_title", "read_count")}
        meta.update({f"geo_{k}": v for k, v in r.get("_geo", {}).items()})
        project.add_sample(V.sample_name(r["run_accession"]), source=source, accession=r["run_accession"],
                           layout=r["library_layout"], bio_unit=r.get("experiment_accession") or r["run_accession"],
                           r1=None, r2=None, metadata=meta,
                           download={"files": r.get("_files", [])})
    project.state["read_type"] = "paired" if runs[0]["library_layout"] == "PAIRED" else "single"
    project.state["organism_hint"] = runs[0].get("scientific_name")
    project.save()
