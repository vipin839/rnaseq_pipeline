"""Pipeline state machine: stage definitions, per-stage screens, checkpoints, resume.

Each stage follows START -> VALIDATION (inputs) -> EXECUTION -> POST-VALIDATION -> SUCCESS/FAILURE.
A checkpoint is written only after post-validation succeeds.
"""
import csv
import json
import os
from pathlib import Path

from . import (PipelineError, Terminated, UserAbort, alignment, bam_manager, count_matrix, data_manager, dependency_manager, design,
               ena_manager, fastq_validator, featurecounts, qc_manager, quality_assessment, r_bridge,
               reference_manager, runner, sra_manager, storage, strandedness, stringtie, trimming, ui)
from . import config as C
from .checkpoint import Checkpoints


class Context:
    def __init__(self, project, cfg, envs, sysinfo, dry_run=False, auto=False):
        self.project = project
        self.cfg = cfg
        self.envs = envs
        self.sysinfo = sysinfo
        self.dry_run = dry_run
        self.auto = auto
        self.cp = Checkpoints(project)
        self.threads = C.resolve_threads(cfg, sysinfo["cpu_cores"])
        self.mem_gb = C.resolve_memory_gb(cfg, sysinfo["ram_total_gb"])

    def log(self, stage, sample=None):
        d = self.project.path("logs", stage)
        d.mkdir(parents=True, exist_ok=True)
        return d / (f"{sample}.log" if sample else f"{stage}.log")

    @property
    def samples(self):
        return self.project.active_samples()

    def ref(self):
        r = self.project.state.get("reference")
        if not r:
            raise PipelineError("no reference selected", remedy="use 'Manage Reference Databases' first")
        return r


# ============================================================================ helpers

def _fail_samples(ctx, stage, failures):
    """Report per-sample failures; ask whether to continue without them."""
    if not failures:
        return
    ui.section(f"{len(failures)} SAMPLE(S) FAILED IN {stage.upper()}")
    for sid, why in failures.items():
        print(f"  {sid}: {why}")
    remaining = [s for s in ctx.samples if s not in failures]
    for sid, why in failures.items():
        ctx.project.mark_failed(sid, why, stage)
    if len(remaining) < 2 or ctx.dry_run:
        raise PipelineError(f"{stage}: too few valid samples remain ({len(remaining)})", stage=stage)
    if not ui.ask_yes_no(f"Continue with the remaining {len(remaining)} sample(s)? (failed samples are "
                         "excluded from all downstream steps; their files are kept)", default=False):
        raise PipelineError(f"{stage}: stopped by user after sample failures", stage=stage,
                            remedy="fix the failed samples, use 'Samples' menu to re-include them, and resume")


def _validation_cache(ctx):
    p = ctx.project.path("data", "metadata", "fastq_validation_cache.json")
    try:
        return p, json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return p, {}


def _file_key(path):
    st = Path(path).stat()
    return f"{Path(path).resolve()}|{st.st_size}|{st.st_mtime_ns}"


def validate_fastqs(ctx, pairs, report_path, stage, expected=None):
    """pairs: {sid: (r1, r2)}. Streams every file; uses a cache keyed by path/size/mtime to avoid
    re-reading unchanged files after an interruption. Returns rows; raises nothing (caller decides)."""
    cache_path, cache = _validation_cache(ctx)
    fv = ctx.cfg["fastq_validation"]
    jobs, rows_by = [], {}
    for sid, (r1, r2) in pairs.items():
        key = "|".join(_file_key(p) for p in (r1, r2) if p and Path(p).exists())
        if key and key in cache and all(r["status"] != "FAILED" for r in cache[key]):
            rows_by[sid] = cache[key]
            continue
        jobs.append((sid, str(r1), str(r2) if r2 else None, fv["allowed_bases"], fv["quality_offset"]))
    if pairs and len(jobs) < len(pairs):
        ui.skipped(f"{len(pairs) - len(jobs)} sample(s) unchanged since a previous successful validation")
    workers = fv.get("workers", "auto")
    workers = min(len(jobs), max(1, ctx.threads // 2)) if workers == "auto" else int(workers)
    if jobs and not ctx.dry_run:
        ui.running(f"Validating {sum(2 if j[2] else 1 for j in jobs)} FASTQ file(s) with {max(1, workers)} worker(s) "
                   "(streaming; full read of every file)")

        def progress(sid, rows):
            st = "FAILED" if any(r["status"] == "FAILED" for r in rows) else "OK"
            reads = rows[0].get("reads", 0)
            (ui.ok if st == "OK" else ui.failed)(f"{sid}: {reads:,} reads" + ("" if st == "OK" else
                                                  f" — {rows[0]['problem'] or rows[-1]['problem']}"))

        for rows in [fastq_validator.validate_many(jobs, workers=max(1, workers), progress=progress)]:
            for r in rows:
                rows_by.setdefault(r["sample"], []).append(r)
        for sid, (r1, r2) in pairs.items():
            if sid in rows_by and all(r["status"] != "FAILED" for r in rows_by[sid]):
                key = "|".join(_file_key(p) for p in (r1, r2) if p)
                cache[key] = rows_by[sid]
        cache_path.write_text(json.dumps(cache, indent=1))
    elif jobs:
        ui.status("DRY-RUN", f"would stream-validate {len(jobs)} sample(s)")
        return []
    rows = [r for sid in pairs for r in rows_by.get(sid, [])]
    if expected:
        for r in rows:
            exp = expected.get(r["sample"])
            if r["status"] != "FAILED" and exp is not None and r.get("reads") != exp:
                r.update(status="FAILED", problem="read count mismatch",
                         reason=f"{r['reads']} reads, expected {exp} (from fastp report)",
                         recommended_action="re-run trimming")
    fastq_validator.write_report(rows, report_path)
    return rows


def _failures_from_rows(rows):
    out = {}
    for r in rows:
        if r["status"] == "FAILED" and r["sample"] not in out:
            out[r["sample"]] = f"{r['problem']} in {Path(r['file']).name} (record {r['record']}): {r['reason']}. " \
                               f"Action: {r['recommended_action']}"
    return out


def _print_validation_failure(rows):
    for r in rows:
        if r["status"] == "FAILED" and r["problem"] != "mate failed":
            print()
            ui.kv([("STATUS", "FAILED"), ("Sample", r["sample"]), ("File", r["file"]),
                   ("Problem", r["problem"]), ("Record", r["record"]), ("Reason", r["reason"]),
                   ("Recommended action", r["recommended_action"])], indent=2)


# ============================================================================ stages

class Stage:
    key = ""
    title = ""
    depends = ()
    expensive = False
    settings = ()          # shown/editable on the stage screen
    result_settings = ()   # settings whose change makes this stage's existing result stale

    def required_tools(self, ctx):
        """Executables this stage needs; checked before it starts (see run_stage)."""
        return ()

    def overview(self, ctx):
        return []

    def check_inputs(self, ctx):
        return []

    def estimate_gb(self, ctx):
        return 0.0

    def execute(self, ctx):
        raise NotImplementedError

    def revalidate(self, ctx, data):
        return []


def _set_aside_partial_pair(project, acc):
    """Before switching archive for a run, move any mate file already fetched from the other archive aside so
    both mates always come from the same source (never overwritten, never mixed)."""
    for f in project.path("data", "fastq").glob(f"{acc}*.fastq.gz"):
        target = f.with_name(f.name + ".other_archive")
        os.replace(f, target)
        ui.info(f"kept {f.name} aside as {target.name} (mates must come from one archive)")


class DataStage(Stage):
    key, title, expensive = "data_acquired", "DATA DOWNLOAD / IMPORT", True
    settings = ("download.source_preference", "download.max_reads", "import_mode")
    result_settings = ("download.max_reads",)


    def required_tools(self, ctx):
        srcs = {ctx.project.samples[x].get("source") for x in ctx.samples
                if not ctx.project.samples[x].get("r1")}
        return (("curl",) if "ena" in srcs else ()) + (
            ("prefetch", "vdb-validate", "fasterq-dump", "fastq-dump") if srcs & {"sra", "ena"} else ())
    def overview(self, ctx):
        src = ctx.project.state.get("data_source")
        pending = [s for s in ctx.samples if not ctx.project.samples[s].get("r1")]
        return [("Source", src), ("Samples", len(ctx.samples)), ("To download", len(pending)),
                ("Read type", ctx.project.state.get("read_type")),
                ("Pilot subset", f"first {ctx.cfg['download']['max_reads']:,} reads (TEST ONLY)"
                 if ctx.cfg["download"]["max_reads"] else "no (full data)")]

    def check_inputs(self, ctx):
        if not ctx.project.samples:
            return ["no samples defined — choose 'Continue pipeline' in the project menu (it asks for a data source)"]
        if not ctx.samples:
            bad = [r for r in ctx.project.samples.values() if r.get("status") in ("FAILED", "EXCLUDED")]
            return [f"all {len(bad)} sample(s) are FAILED/EXCLUDED — choose 'Continue pipeline' in the project "
                    "menu to retry them or replace the data source"]
        sc = data_manager.single_cell_samples(ctx.project, ctx.samples)
        if sc:
            first = next(iter(sc))
            return [f"{len(sc)} sample(s) are SINGLE-CELL RNA-seq ({first}: {sc[first]}); this bulk pipeline "
                    "cannot analyse them — use 'Choose or replace the data source' in the project menu"]
        want = ctx.cfg.get("read_type", "auto")
        have = ctx.project.state.get("read_type")
        if want != "auto" and have and want != have:
            return [f"config read_type is '{want}' but the data are {have}-end; fix the config or the input"]
        return []

    def estimate_gb(self, ctx):
        return storage.estimate(ctx.project, self.key, ctx.samples)

    def execute(self, ctx):
        p = ctx.project
        max_reads = int(ctx.cfg["download"].get("max_reads") or 0)
        failures = {}
        for sid in ctx.samples:
            rec = p.samples[sid]
            if rec["source"] == "local":
                continue
            if rec.get("r1") and p.abs(rec["r1"]).exists():
                continue
            try:
                files = (rec.get("download") or {}).get("files") or []
                expected = 2 if rec.get("layout") == "PAIRED" else 1
                if rec["source"] == "ena" and len(files) != expected:
                    why = "has no mirrored FASTQ" if not files else \
                        f"lists {len(files)} FASTQ file(s) for a {rec.get('layout')} run"
                    ui.info(f"{sid}: ENA {why}; using the NCBI SRA route (fasterq-dump splits the reads itself)")
                    rec["source"] = "sra"
                    p.save()
                if rec["source"] == "sra":
                    r1, r2 = sra_manager.download_run(rec["accession"], p.path("data", "raw"), p.path("data", "fastq"),
                                                      p.path("temp"), ctx.threads, rec["layout"] == "PAIRED",
                                                      max_reads=max_reads, log_file=ctx.log("data", sid))
                else:
                    run = {"run_accession": rec["accession"], "library_layout": rec["layout"],
                           "_files": [tuple(f) for f in rec["download"]["files"]]}
                    try:
                        r1, r2 = ena_manager.download_run(run, p.path("data", "fastq"), max_reads=max_reads,
                                                          log_file=ctx.log("data", sid),
                                                          retries=int(ctx.cfg["download"].get("retries", 3)))
                    except PipelineError as ena_err:
                        if runner.which("fasterq-dump") is None or runner.which("fastq-dump") is None:
                            raise
                        ui.warn(f"{sid}: ENA download failed ({ena_err}); the same run is mirrored at NCBI SRA — "
                                "trying the SRA route")
                        _set_aside_partial_pair(p, rec["accession"])
                        rec["source"] = "sra"
                        p.save()
                        r1, r2 = sra_manager.download_run(rec["accession"], p.path("data", "raw"),
                                                          p.path("data", "fastq"), p.path("temp"), ctx.threads,
                                                          rec["layout"] == "PAIRED", max_reads=max_reads,
                                                          log_file=ctx.log("data", sid))
                if not ctx.dry_run:
                    rec["r1"] = p.rel(r1)
                    rec["r2"] = p.rel(r2) if r2 else None
                    rec["pilot_max_reads"] = max_reads or None
                    if rec["source"] == "ena" and not max_reads:
                        rec["provider_md5"] = [f[1] for f in (rec.get("download") or {}).get("files", [])]
                    p.save()
                    ui.ok(f"{sid}: downloaded")
            except PipelineError as e:
                failures[sid] = str(e) + (f" ({e.cause})" if e.cause else "")
                ui.failed(f"{sid}: {e}")
        if ctx.dry_run:
            return [], {}, {}
        _fail_samples(ctx, self.key, failures)
        outs = []
        for sid in ctx.samples:
            outs += [x for x in p.fastqs(sid, trimmed=False) if x]
        missing = [str(o) for o in outs if not Path(o).exists()]
        if missing:
            raise PipelineError(f"FASTQ files missing after import: {missing[:3]}", stage=self.key)
        return outs, {"samples": ctx.samples, "max_reads": max_reads}, {"files": len(outs)}


class FastqVerifyStage(Stage):
    key, title, depends = "fastq_verified", "FASTQ VERIFICATION", ("data_acquired",)
    settings = ("fastq_validation.workers", "fastq_validation.quality_offset", "fastq_validation.allowed_bases")
    result_settings = ("fastq_validation.quality_offset", "fastq_validation.allowed_bases")

    def overview(self, ctx):
        return [("Samples", len(ctx.samples)), ("Read type", ctx.project.state.get("read_type")),
                ("Checks", "4-line records, headers, alphabet, quality, lengths, gzip, mate IDs/counts")]

    def execute(self, ctx):
        p = ctx.project
        pairs = {sid: p.fastqs(sid, trimmed=False) for sid in ctx.samples}
        report = p.path("data", "metadata", "fastq_validation_report.tsv")
        rows = validate_fastqs(ctx, pairs, report, self.key)
        if ctx.dry_run:
            return [], {}, {}
        failures = _failures_from_rows(rows)
        if failures:
            _print_validation_failure(rows)
            ui.error("Samples that fail FASTQ validation cannot be aligned.")
        _fail_samples(ctx, self.key, failures)
        warn = [r for r in rows if r["status"] == "WARNING"]
        for r in warn:
            ui.warn(f"{r['sample']} {r['mate']}: {r['problem']} — {r['recommended_action']}")
        for r in rows:
            if r["sample"] in ctx.samples:
                p.samples[r["sample"]].setdefault("reads", {})[r["mate"]] = r["reads"]
                p.samples[r["sample"]]["read_length"] = r["max_length"]
        p.save()
        ui.ok(f"FASTQ validation passed for {len(ctx.samples)} sample(s); report: {p.rel(report)}")
        return [report], {"samples": ctx.samples}, {"files": len(rows)}


class RawQCStage(Stage):
    key, title, depends = "raw_qc_completed", "RAW READ QC (FastQC + MultiQC)", ("fastq_verified",)


    def required_tools(self, ctx):
        return ("fastqc", "multiqc")
    def overview(self, ctx):
        return [("Files", sum(2 if ctx.project.fastqs(s, trimmed=False)[1] else 1 for s in ctx.samples)),
                ("Threads", ctx.threads)]

    def execute(self, ctx):
        return run_qc(ctx, self.key, trimmed=False)


def run_qc(ctx, stage, trimmed):
    p = ctx.project
    sub = "trimmed" if trimmed else "raw"
    files, expected = [], {}
    for sid in ctx.samples:
        r1, r2 = p.fastqs(sid, trimmed=trimmed)
        for mate, f in (("R1", r1), ("R2", r2)):
            if f:
                files.append(f)
                if trimmed:
                    expected[str(f)] = p.samples[sid]["trimmed"].get("expected_reads_per_file")
                else:
                    expected[str(f)] = p.samples[sid].get("reads", {}).get(mate)
    fq_dir, mq_dir = p.path("qc", f"fastqc_{sub}"), p.path("qc", f"multiqc_{sub}")
    qc_manager.run_fastqc(files, fq_dir, p.path("temp"), ctx.threads, ctx.log(stage), stage)
    if not ctx.dry_run:
        qc_manager.validate_fastqc(files, fq_dir, expected, stage)
    inputs = [fq_dir] + ([p.path("qc", "fastp")] if trimmed else [])
    qc_manager.run_multiqc(inputs, mq_dir, f"{p.name}: {sub} reads", ctx.log(stage), stage,
                           expected_sources=[fq_dir / f"{qc_manager.fastqc_basename(f)}.zip" for f in files])
    if ctx.dry_run:
        return [], {}, {}
    return [fq_dir, mq_dir], {"samples": ctx.samples}, {"files": len(files)}


def assess_quality(ctx, trimmed):
    p = ctx.project
    sub = "trimmed" if trimmed else "raw"
    per = {}
    for sid in ctx.samples:
        r1, r2 = p.fastqs(sid, trimmed=trimmed)
        per[sid] = {}
        for mate, f in (("R1", r1), ("R2", r2)):
            if f:
                z = p.path("qc", f"fastqc_{sub}", qc_manager.fastqc_basename(f) + ".zip")
                per[sid][mate] = quality_assessment.metrics_from_fastqc(z)
    res = quality_assessment.assess(per, ctx.cfg["quality_gate"], p.is_paired())
    out_dir = p.path("qc", "assessment", sub)
    files = quality_assessment.write_outputs(res, out_dir, ctx.cfg["quality_gate"],
                                             f"Quality assessment — {sub} reads")
    return res, files


def show_assessment(res):
    ui.section("QUALITY ASSESSMENT")
    ui.table([(s, r["status"]) for s, r in res.items()], ["Sample", "Status"])
    reasons = {}
    for s, r in res.items():
        for x in r["reasons"]:
            reasons.setdefault(x.split("] ", 1)[1].split(" (")[0].split(":")[-1].strip()[:70], []).append(s)
    if reasons:
        print("\nReasons:")
        for s, r in res.items():
            for x in r["reasons"]:
                print(f"  {s}: {x}")
    print("\nNote: FastQC warnings do not automatically mean data are unusable.")


# What each finding means for bulk RNA-seq: technical problems that trimming fixes, versus patterns that are normal
# for RNA-seq libraries (and are no reason to trim or drop a sample).
_REVIEW_NOTES = [
    ("adapter", "Adapter content: TECHNICAL — reads run into the adapter when the fragment is shorter than the read; "
                "adapter bases do not align. Trimming removes them."),
    ("tail", "Low-quality read ends: TECHNICAL — error-prone 3' bases lower alignment; trimming removes them."),
    ("lower quartile", "Per-base lower quartile: TECHNICAL if it drops only at the ends (trimming helps); if it is low "
                       "throughout, the run itself was poor (trimming cannot fix that)."),
    ("N content", "N content: TECHNICAL — uncalled bases at some positions (often a sequencer issue)."),
    ("duplication", "Duplication: usually normal for RNA-seq — highly expressed genes give many identical reads. Not a "
                    "trimming issue and no reason to drop a sample."),
    ("overrepresented", "Overrepresented sequences: often normal for RNA-seq (very abundant transcripts, e.g. "
                        "mitochondrial/ribosomal RNA); check them in FastQC if the share is high."),
    ("%GC", "GC deviating from the other samples: may indicate contamination or a different library type; check."),
    ("reads (threshold", "Few reads: lower statistical power; trimming cannot fix this."),
    ("R1/R2 read counts differ", "R1/R2 counts differ: the files are not proper mates (corrupt or truncated)."),
]


def show_review(ctx, res):
    """Quality gate option 'Review': the evidence per sample, what it means, and where the full reports are."""
    T = ctx.cfg["quality_gate"]
    ui.section("REVIEW — evidence per sample (thresholds from config 'quality_gate')")
    def f(v, d=1):
        return "-" if v is None or (isinstance(v, float) and v != v) else f"{v:.{d}f}"

    rows = []
    for sid, r in res.items():
        for mate, m in r["metrics"].items():
            rows.append((sid, mate, r["status"], f"{m.get('total_sequences', 0):,}",
                         f(m.get("adapter_max_percent")), f(m.get("tail_mean_quality")),
                         f(m.get("min_lower_quartile"), 0), f(m.get("duplication_percent"), 0),
                         f(m.get("gc_percent"), 0), f(m.get("n_content_max_percent")),
                         f(m.get("overrepresented_max_percent"), 2)))
    ui.table(rows, ["Sample", "Mate", "Status", "Reads", f"Adapter% (trim >{T['adapter_max_percent_trim']})",
                    f"Tail Q (trim <{T['tail_mean_quality_trim']})", "Lower-quartile Q", "Dup%", "GC%", "N%",
                    "Overrep%"])
    found = " ".join(x for r in res.values() for x in r["reasons"])
    fastqc_fails = {k for r in res.values() for m in r["metrics"].values()
                    for k, v in (m.get("fastqc_status") or {}).items() if v in ("FAIL", "WARN")}
    print("\nWhat the findings mean:")
    for key, note in _REVIEW_NOTES:
        if key in found:
            print(f"  - {note}")
    if "Per base sequence content" in fastqc_fails:
        print("  - FastQC 'Per base sequence content' warning: normal for RNA-seq (random-hexamer priming bias in "
              "the first ~12 bases); not a reason to trim.")
    if "Sequence Duplication Levels" in fastqc_fails and "duplication" not in found:
        print("  - FastQC 'Sequence Duplication Levels' warning: normal for RNA-seq (highly expressed genes).")
    p = ctx.project
    print(f"\nFull reports: {p.rel(p.path('qc', 'multiqc_raw', 'multiqc_report.html'))} · "
          f"{p.rel(p.path('qc', 'assessment', 'raw', 'quality_assessment.html'))}\n")


class QualityGateStage(Stage):
    key, title, depends = "quality_assessed", "QUALITY ASSESSMENT / QUALITY GATE", ("raw_qc_completed",)
    settings = ("trim_adapters", "quality_gate")
    result_settings = ("trim_adapters", "quality_gate")

    def execute(self, ctx):
        if ctx.dry_run:
            ui.status("DRY-RUN", "would parse FastQC results and apply thresholds from config quality_gate")
            return [], {}, {}
        res, files = assess_quality(ctx, trimmed=False)
        show_assessment(res)
        statuses = {r["status"] for r in res.values()}
        fails = {s: "; ".join(r["reasons"]) for s, r in res.items() if r["status"] == "FAIL"}
        if fails:
            _fail_samples(ctx, self.key, fails)
        recommend = "TRIMMING RECOMMENDED" in statuses
        ui.action("fastp trimming recommended" if recommend else "no trimming needed by the configured thresholds")
        mode = ctx.cfg.get("trim_adapters", "ask")
        if mode == "always":
            decision = "trim"
            ui.info("config trim_adapters=always -> trimming will run")
        elif mode == "never":
            decision = "skip"
            ui.info("config trim_adapters=never -> trimming skipped")
        else:
            ui.explain(
                "whether to trim adapters and low-quality read ends with fastp",
                "adapter sequence and low-quality 3' ends stop reads from aligning; trimming removes them, "
                "while trimming already-clean reads changes little",
                "if you trim, ALL samples are trimmed the same way (raw files are never modified); skipping "
                "despite a recommendation usually lowers alignment rates",
                "the recommendation above comes from the thresholds in config 'quality_gate'")
            options = [f"Accept recommendation ({'run fastp' if recommend else 'skip trimming'})",
                       "Run fastp", "Skip trimming", "Review the evidence per sample", "Stop pipeline"]
            while True:
                c = ui.choose("DECISION", options)
                if c != 3:
                    break
                show_review(ctx, res)
            if c == 4:
                raise UserAbort("stopped at quality gate")
            decision = {0: "trim" if recommend else "skip", 1: "trim", 2: "skip"}[c]
            if decision == "skip" and recommend:
                ui.warn("continuing without trimming despite the recommendation (your decision is recorded)")
        if "REVIEW" in statuses and decision == "skip":
            if not ui.ask_yes_no("Some samples are flagged REVIEW (see reasons). Continue with them?", default=True):
                raise UserAbort("stopped at quality gate (REVIEW samples)")
        d = {"decision": decision, "recommended": "trim" if recommend else "skip",
             "statuses": {s: r["status"] for s, r in res.items()}}
        dp = ctx.project.path("qc", "assessment", "qc_decision.json")
        dp.write_text(json.dumps(d, indent=2))
        ctx.project.state["qc_decision"] = d
        ctx.project.save()
        return files + [dp], {"samples": ctx.samples, "decision": decision}, d


class TrimmingStage(Stage):
    key, title, depends, expensive = "trimming_completed", "CONDITIONAL TRIMMING (fastp) + POST-TRIM QC", \
        ("quality_assessed",), True
    settings = ("fastp_parameters",)
    result_settings = ("fastp_parameters",)


    def required_tools(self, ctx):
        return ("fastp", "fastqc", "multiqc") if (ctx.project.state.get("qc_decision") or {}).get("decision") == "trim" else ()
    def overview(self, ctx):
        d = (ctx.project.state.get("qc_decision") or {}).get("decision")
        fp = ctx.cfg["fastp_parameters"]
        return [("Decision", d), ("Samples", len(ctx.samples)), ("Threads", ctx.threads),
                ("Parameters", f"Q{fp['qualified_quality_phred']}, cut_right={fp['cut_right']} "
                               f"(w{fp['cut_right_window_size']}/Q{fp['cut_right_mean_quality']}), "
                               f"min length {fp['length_required']}")]

    def estimate_gb(self, ctx):
        d = (ctx.project.state.get("qc_decision") or {}).get("decision")
        return storage.estimate(ctx.project, self.key, ctx.samples) if d == "trim" else 0

    def execute(self, ctx):
        p = ctx.project
        decision = (p.state.get("qc_decision") or {}).get("decision")
        if decision != "trim" and not ctx.dry_run:
            for sid in p.samples:
                p.samples[sid].pop("trimmed", None)
            p.save()
            ui.skipped("trimming not requested at the quality gate; raw reads will be aligned")
            dp = p.path("qc", "assessment", "qc_decision.json")
            return [dp], {"samples": ctx.samples, "trimmed": False}, {"trimmed": False}
        failures = {}
        for sid in ctx.samples:
            try:
                trimming.trim_sample(p, sid, ctx.cfg["fastp_parameters"], ctx.threads, ctx.log("trimming", sid))
            except PipelineError as e:
                failures[sid] = str(e)
        if ctx.dry_run:
            return [], {}, {}
        _fail_samples(ctx, self.key, failures)
        ui.running("Validating trimmed FASTQ files")
        pairs = {sid: p.fastqs(sid, trimmed=True) for sid in ctx.samples}
        expected = {sid: p.samples[sid]["trimmed"]["expected_reads_per_file"] for sid in ctx.samples}
        rows = validate_fastqs(ctx, pairs, p.path("data", "metadata", "fastq_validation_trimmed.tsv"),
                               self.key, expected)
        bad = _failures_from_rows(rows)
        if bad:
            _print_validation_failure(rows)
            raise PipelineError("trimming produced invalid FASTQ — pipeline stopped", stage=self.key,
                                remedy="inspect fastp logs; do not align these files")
        outs, _, _ = run_qc(ctx, self.key, trimmed=True)
        res, files = assess_quality(ctx, trimmed=True)
        ui.section("POST-TRIMMING ASSESSMENT")
        ui.table([(s, r["status"], f"{p.samples[s]['trimmed']['reads_after']:,} / "
                   f"{p.samples[s]['trimmed']['reads_before']:,}") for s, r in res.items()],
                 ["Sample", "Status", "Reads kept / before"])
        trimmed_files = [x for sid in ctx.samples for x in p.fastqs(sid, trimmed=True) if x]
        return trimmed_files + outs + files + [p.path("qc", "fastp")], {"samples": ctx.samples, "trimmed": True,
                                                                          "fastp": ctx.cfg["fastp_parameters"]}, \
            {"trimmed": True}


class ReferenceStage(Stage):
    key, title, expensive = "reference_ready", "REFERENCE PREPARATION", True
    settings = ("hisat2_build.use_splice_sites_in_index", "reference_store")
    result_settings = ("hisat2_build", "featurecounts_parameters.feature_type",
                       "featurecounts_parameters.attribute")


    def required_tools(self, ctx):
        return ("samtools", "hisat2-build", "hisat2-inspect", "hisat2_extract_splice_sites.py", "hisat2_extract_exons.py", "curl", "gzip")
    def overview(self, ctx):
        r = ctx.project.state.get("reference") or {}
        return [("Reference", r.get("label", "(not selected)")), ("Store", r.get("store")),
                ("Threads", ctx.threads), ("RAM", f"{ctx.sysinfo['ram_total_gb']} GB")]

    def check_inputs(self, ctx):
        return [] if ctx.project.state.get("reference") else ["no reference selected (main menu 5)"]

    def estimate_gb(self, ctx):
        r = ctx.project.state.get("reference") or {}
        if Path(r.get("store", "/nonexistent")).joinpath("index").is_dir():
            return 0
        return storage.estimate(ctx.project, self.key, [], r.get("package"))

    def execute(self, ctx):
        r = ctx.ref()
        log = ctx.log("reference")

        def on_feature_type(alt, why, counts):
            if not ui.ask_yes_no(f"Count '{alt}' features instead (and skip StringTie2, which needs exons)?",
                                 default=True):
                return False
            ctx.cfg["featurecounts_parameters"]["feature_type"] = alt
            if counts.get("exon", 0) < counts.get("gene", 0) * 0.5:
                ctx.cfg["stringtie_parameters"]["enabled"] = False
            ctx.project.save_config(ctx.cfg)
            ui.ok(f"configuration updated: feature_type={alt}"
                  + (", StringTie2 disabled" if not ctx.cfg["stringtie_parameters"].get("enabled", True) else ""))
            return True

        manifest, info = reference_manager.prepare(
            r, ctx.cfg, ctx.threads, ctx.mem_gb, ctx.sysinfo["ram_available_gb"], log,
            on_feature_type=on_feature_type,
            confirm_ram=lambda need, avail: ui.ask_yes_no("Attempt the index build anyway?", default=False),
            dry_run=ctx.dry_run)
        if ctx.dry_run:
            return [], {}, {}
        paths, fa, gtf = info["paths"], info["fa"], info["gtf"]
        proj_manifest = ctx.project.path("reference", "reference_manifest.yaml")
        reference_manager.write_manifest(proj_manifest, manifest)
        reference_manager.link_into_project(ctx.project, paths)
        summary_txt = ctx.project.path("reference", "reference_summary.txt")
        summary_txt.write_text(
            f"Reference: {r.get('label')}\nGenome: {fa['n_sequences']} sequences, {fa['total_bp']:,} bp\n"
            f"Annotation: {gtf['genes']:,} genes, {gtf['transcripts']:,} transcripts, "
            f"{gtf['exons']:,} {info['feature_type']} features\n"
            f"HISAT2 index: {paths.index_prefix} (splice sites in index: {info['use_ss']})\n")
        r.update(info["state"])
        ctx.project.save()
        outs = [proj_manifest, summary_txt, paths.gtf, paths.splice_sites, paths.bed12,
                *reference_manager.index_files(paths.index_prefix)]
        return outs, {"reference": r.get("label")}, {"genes": gtf["genes"]}

    def revalidate(self, ctx, data):
        r = ctx.project.state.get("reference") or {}
        if not r.get("index_prefix") or not reference_manager.index_files(r["index_prefix"]):
            return ["HISAT2 index files missing"]
        return []


class AlignmentStage(Stage):
    key, title, depends, expensive = "alignment_completed", "ALIGNMENT (HISAT2)", \
        ("trimming_completed", "reference_ready"), True
    result_settings = ("hisat2_parameters", "min_overall_alignment_rate_fail")
    settings = ("hisat2_parameters", "samtools_sort_memory_per_thread", "min_overall_alignment_rate_fail",
                "min_overall_alignment_rate_warn")


    def required_tools(self, ctx):
        return ("hisat2", "samtools")
    def overview(self, ctx):
        r = ctx.project.state.get("reference") or {}
        trimmed = bool((ctx.project.state.get("qc_decision") or {}).get("decision") == "trim")
        return [("Samples", len(ctx.samples)),
                ("Validated", f"{len(ctx.samples)}/{len(ctx.project.samples)}"),
                ("Input reads", "trimmed" if trimmed else "raw"),
                ("Reference", r.get("label")), ("Threads", ctx.threads),
                ("Splice sites", "in index" if r.get("index_has_splice_sites") else "supplied at alignment"),
                ("Options", "--dta" if ctx.cfg["hisat2_parameters"].get("dta") else "")]

    def estimate_gb(self, ctx):
        return storage.estimate(ctx.project, self.key, ctx.samples)

    def ram_needed_gb(self, ctx):
        """HISAT2 holds the whole index in memory; samtools sort adds its per-thread buffers."""
        files = reference_manager.index_files(ctx.ref().get("index_prefix", "")) or []
        index_gb = sum(Path(f).stat().st_size for f in files) / 1e9
        mem = alignment.sort_memory(ctx.cfg, ctx.threads, ctx.mem_gb)
        per_thread = int(mem[:-1]) / (1024 if mem.endswith("M") else 1)
        return round(index_gb * 1.15 + 0.5 + per_thread * max(1, ctx.threads - 1), 1)

    def execute(self, ctx):
        p, r = ctx.project, ctx.ref()
        paired = p.is_paired()
        if not ctx.dry_run:
            need, avail = self.ram_needed_gb(ctx), ctx.sysinfo["ram_available_gb"]
            ui.info(f"estimated memory for alignment: {need} GB (available now: {avail} GB)")
            if need > avail:
                ui.warn("alignment may run out of memory and be killed by the operating system. Close other "
                        "programs, lower 'threads', or use a machine with more RAM.")
                if not ui.ask_yes_no("Start the alignment anyway?", default=False):
                    raise UserAbort("alignment postponed: not enough free memory")
        rows, failures, outs = [], {}, []
        for sid in ctx.samples:
            final = p.path("alignment", "bam", f"{sid}.sorted.bam")
            reads = p.samples[sid]["trimmed"]["expected_reads_per_file"] if p.samples[sid].get("trimmed") \
                else p.samples[sid].get("reads", {}).get("R1")
            if not ctx.dry_run and final.exists():
                ok, probs, m = bam_manager.validate_bam(final, reads, paired, ctx.threads)
                if ok:
                    ui.skipped(f"{sid}: existing BAM is valid; not re-aligning")
                    rows.append(self._row(ctx, sid, m, reads, "PASS", []))
                    outs += [final, Path(str(final) + ".bai"), p.path("alignment", "reports", f"{sid}.hisat2.summary")]
                    continue
                bad = final.with_name(f"{sid}.invalid.bam")
                os.replace(final, bad)
                ui.warn(f"{sid}: existing BAM invalid ({'; '.join(probs)}); moved to {bad.name}; re-aligning")
            try:
                alignment.align_sample(p, sid, r, ctx.cfg, ctx.threads, ctx.mem_gb, ctx.log("alignment", sid))
            except PipelineError as e:
                failures[sid] = str(e)
                ui.failed(f"{sid}: alignment failed")
                continue
            if ctx.dry_run:
                continue
            ok, probs, m = bam_manager.validate_bam(final, reads, paired, ctx.threads)
            hs = alignment.parse_hisat2_summary(p.path("alignment", "reports", f"{sid}.hisat2.summary"))
            if hs.get("total") != reads:
                probs.append(f"HISAT2 processed {hs.get('total')} {hs.get('unit', 'reads')}, expected {reads}")
                ok = False
            rate = hs.get("overall_alignment_rate", 0)
            if rate < ctx.cfg["min_overall_alignment_rate_fail"]:
                probs.append(f"overall alignment rate {rate}% < {ctx.cfg['min_overall_alignment_rate_fail']}% "
                             "(wrong organism/reference or contaminated library?)")
                ok = False
            elif rate < ctx.cfg["min_overall_alignment_rate_warn"]:
                ui.warn(f"{sid}: overall alignment rate {rate}% is low")
            if not ok:
                failures[sid] = "; ".join(probs)
                bad = final.with_name(f"{sid}.FAILED.bam")
                os.replace(final, bad)
                if Path(str(final) + ".bai").exists():
                    os.replace(Path(str(final) + ".bai"), Path(str(bad) + ".bai"))
                ui.failed(f"{sid}: BAM failed validation — {failures[sid]}")
                rows.append(self._row(ctx, sid, m, reads, "FAILED", probs))
                continue
            ui.ok(f"{sid}: BAM validated ({rate}% aligned, {m['primary']:,} primary records)")
            rows.append(self._row(ctx, sid, m, reads, "PASS", []))
            outs += [final, Path(str(final) + ".bai"), p.path("alignment", "reports", f"{sid}.hisat2.summary")]
        if ctx.dry_run:
            return [], {}, {}
        summ = bam_manager.write_summary(rows, p.path("alignment", "reports", "alignment_summary.tsv"))
        ui.info(f"alignment summary: {p.rel(summ)}")
        _fail_samples(ctx, self.key, failures)
        outs = [o for o in outs if not any(o.name.startswith(f + ".") for f in failures)]
        return outs + [summ], {"samples": ctx.samples}, {"aligned": len(ctx.samples)}

    def _row(self, ctx, sid, m, reads, status, probs):
        hs = alignment.parse_hisat2_summary(ctx.project.path("alignment", "reports", f"{sid}.hisat2.summary")) \
            if ctx.project.path("alignment", "reports", f"{sid}.hisat2.summary").exists() else {}
        return {"sample": sid, "status": status, "input_reads": reads, "unit": hs.get("unit"),
                "overall_alignment_rate": hs.get("overall_alignment_rate"),
                "uniquely_aligned_pct": hs.get("uniquely_aligned_pct"),
                "multi_aligned_pct": hs.get("multi_aligned_pct"), **m, "problems": "; ".join(probs)}

    def revalidate(self, ctx, data):
        problems = []
        for sid in ctx.samples:
            bam = ctx.project.path("alignment", "bam", f"{sid}.sorted.bam")
            q = runner.tool_output(["samtools", "quickcheck", "-v", str(bam)])
            if q is None or q.strip():
                problems.append(f"{sid}: BAM fails quickcheck")
        return problems


class BamQCStage(Stage):
    key, title, depends = "bam_qc_completed", "BAM QC (samtools + MultiQC)", ("alignment_completed",)


    def required_tools(self, ctx):
        return ("samtools", "multiqc")
    def execute(self, ctx):
        p = ctx.project
        out = p.path("alignment", "reports", "samtools")
        files = []
        for sid in ctx.samples:
            files += bam_manager.qc_reports(p.path("alignment", "bam", f"{sid}.sorted.bam"), out, sid, ctx.threads,
                                            ctx.log("bam_qc", sid))
        mq = p.path("qc", "multiqc_alignment")
        qc_manager.run_multiqc([p.path("alignment", "reports"), p.path("qc", "fastqc_raw")], mq,
                               f"{p.name}: alignment QC", ctx.log("bam_qc"), self.key,
                               expected_sources=files + [p.path("alignment", "reports", f"{sid}.hisat2.summary")
                                                         for sid in ctx.samples])
        if ctx.dry_run:
            return [], {}, {}
        comb = p.path("alignment", "reports", "bam_qc_summary.tsv")
        keys = ["raw total sequences", "reads mapped", "reads unmapped", "reads properly paired",
                "reads duplicated", "error rate", "average length", "insert size average"]
        with open(comb, "w") as f:
            f.write("sample\t" + "\t".join(k.replace(" ", "_") for k in keys) + "\n")
            for sid in ctx.samples:
                s = bam_manager.parse_stats_summary(out / f"{sid}.samtools_stats")
                if not s:
                    raise PipelineError(f"samtools stats output empty for {sid}", stage=self.key)
                f.write(sid + "\t" + "\t".join(s.get(k, "") for k in keys) + "\n")
        ui.ok(f"BAM QC complete: {p.rel(comb)}")
        return files + [comb, mq], {"samples": ctx.samples}, {}


class StrandednessStage(Stage):
    key, title, depends = "strandedness_determined", "LIBRARY STRANDEDNESS", ("alignment_completed", "reference_ready")
    settings = ("strandedness", "strandedness_inference")
    result_settings = ("strandedness", "strandedness_inference")

    def revalidate(self, ctx, data):
        """Later stages read the decision from the project state; it must still be the recorded, checkpointed one."""
        try:
            recorded = json.loads(ctx.project.path("alignment", "reports", "strandedness.json").read_text())
        except (OSError, json.JSONDecodeError) as e:
            return [f"recorded strandedness decision unreadable: {e}"]
        live = ctx.project.state.get("strandedness") or {}
        diff = [k for k in ("value", "featurecounts_flag") if live.get(k) != recorded.get(k)]
        if diff:
            return [f"strandedness in the project state ({live.get('value')}) differs from the confirmed decision "
                    f"({recorded.get('value')}); the decision must be made again"]
        return []

    def execute(self, ctx):
        p, cfg = ctx.project, ctx.cfg
        si = cfg["strandedness_inference"]
        configured = cfg.get("strandedness", "auto")
        evidence, calls = {}, {}
        have_rseqc = runner.which("infer_experiment.py") is not None
        if have_rseqc:
            for sid in ctx.samples:
                bam = p.path("alignment", "bam", f"{sid}.sorted.bam")
                ev = strandedness.run_infer(bam, ctx.ref()["bed12"], si["sample_reads"], ctx.log("strandedness", sid), sid)
                evidence[sid] = ev
                calls[sid] = strandedness.call(ev, si["stranded_min_fraction"], si["unstranded_max_diff"],
                                               si["max_undetermined_fraction"])
        if ctx.dry_run:
            return [], {}, {}
        ui.explain(
            "the library strandedness used for counting",
            "featureCounts must know which DNA strand each read represents: a wrong choice discards about half "
            "(stranded vs unstranded) or nearly all (forward vs reverse) of the correctly assigned reads",
            "sets featureCounts -s and the StringTie strand flag; recorded with its evidence in the report",
            "unstranded / forward / reverse (dUTP and TruSeq Stranded kits are 'reverse')")
        ui.section("STRANDEDNESS EVIDENCE (RSeQC infer_experiment)")
        if calls:
            ui.table([(s, f"{e['forward']:.3f}" if e["forward"] is not None else "?",
                       f"{e['reverse']:.3f}" if e["reverse"] is not None else "?",
                       f"{e['undetermined']:.3f}" if e["undetermined"] is not None else "?",
                       calls[s][0] or "UNDETERMINED", calls[s][1]) for s, e in evidence.items()],
                     ["Sample", "Forward frac", "Reverse frac", "Undetermined", "Call", "Confidence"])
            cons, conf, note = strandedness.consensus(calls)
        else:
            ui.warn("RSeQC (infer_experiment.py) not available; cannot infer strandedness")
            cons, conf, note = None, "none", "no inference tool"
        hints = strandedness.metadata_hints(p)
        if hints:
            print("\nMetadata hints (informational):")
            for h in hints:
                print(f"  {h}")
        options = ["unstranded", "forward", "reverse"]
        if configured != "auto":
            value, source = configured, "config (user-specified)"
            if cons and cons != configured:
                ui.warn(f"configured strandedness '{configured}' disagrees with inferred '{cons}'")
                i = ui.choose("Which strandedness should be used?", [f"{configured} (configured)", f"{cons} (inferred)"])
                value, source = (configured, "config (confirmed against inference)") if i == 0 else (cons, "RSeQC (user-selected)")
            confidence = "user-specified"
        elif cons:
            ui.ok(f"inferred strandedness: {cons.upper()} ({conf} confidence; {note})")
            if ui.ask_yes_no(f"Use '{cons}' strandedness for featureCounts/StringTie2?", default=True):
                value, source, confidence = cons, "RSeQC infer_experiment", conf
            else:
                value = options[ui.choose("Select strandedness", options)]
                source, confidence = "user-selected (overrode inference)", "user-specified"
        else:
            print()
            ui.warn("Library strandedness could not be automatically established.")
            value = options[ui.choose("Select strandedness (check the library prep kit; dUTP/TruSeq Stranded = reverse)",
                                      options)]
            source, confidence = "user-selected (inference inconclusive)", "user-specified"
        rec = {"value": value, "source": source, "confidence": confidence,
               "featurecounts_flag": strandedness.FEATURECOUNTS_FLAG[value],
               "stringtie_flag": strandedness.STRINGTIE_FLAG[value],
               "evidence": evidence, "per_sample_calls": {s: list(c) for s, c in calls.items()}}
        p.state["strandedness"] = rec
        p.save()
        out = p.path("alignment", "reports", "strandedness.json")
        out.write_text(json.dumps(rec, indent=2))
        ui.ok(f"strandedness: {value} (source: {source}; confidence: {confidence})")
        return [out], {"samples": ctx.samples, "value": value}, {"value": value}


class StringTieStage(Stage):
    key, title, depends, expensive = "stringtie_completed", "TRANSCRIPT QUANTIFICATION (StringTie2)", \
        ("strandedness_determined", "alignment_completed", "reference_ready"), True
    settings = ("stringtie_parameters",)
    result_settings = ("stringtie_parameters",)


    def required_tools(self, ctx):
        return ("stringtie",) if ctx.cfg["stringtie_parameters"].get("enabled", True) else ()
    def overview(self, ctx):
        return [("Mode", "reference-guided (-e): annotated transcripts only"),
                ("Strandedness", (ctx.project.state.get("strandedness") or {}).get("value")),
                ("Samples", len(ctx.samples)), ("Note", "independent of the DESeq2 count matrix")]

    def execute(self, ctx):
        p = ctx.project
        if not ctx.cfg["stringtie_parameters"].get("enabled", True):
            note = p.path("stringtie", "SKIPPED.txt")
            note.write_text("StringTie2 was disabled in the configuration "
                            "(stringtie_parameters.enabled: false).\n"
                            "Gene counts for DESeq2 come from featureCounts and are unaffected.\n")
            ui.skipped("StringTie2 disabled in configuration")
            return [note], {"samples": ctx.samples, "enabled": False}, {"enabled": False}
        strand = (p.state.get("strandedness") or {}).get("value", "unstranded")
        outs, failures = [], {}
        for sid in ctx.samples:
            try:
                g, a = stringtie.run_sample(p, sid, p.path("alignment", "bam", f"{sid}.sorted.bam"), ctx.ref()["gtf"],
                                            strand, ctx.cfg["stringtie_parameters"], ctx.threads,
                                            ctx.log("stringtie", sid))
                outs += [g, a]
            except PipelineError as e:
                failures[sid] = str(e)
        if ctx.dry_run:
            return [], {}, {}
        if failures:
            raise PipelineError("StringTie2 failed for: " + ", ".join(failures), stage=self.key,
                                cause="; ".join(failures.values()))
        files, ng, nt = stringtie.merge_tables(p, ctx.samples)
        ui.ok(f"StringTie2 validated: {nt:,} transcripts / {ng:,} genes quantified (TPM matrices in stringtie/merged/)")
        return outs + files, {"samples": ctx.samples, "strand": strand}, {"genes": ng, "transcripts": nt}


class FeatureCountsStage(Stage):
    key, title, depends = "featurecounts_completed", "GENE COUNTING (featureCounts)", \
        ("strandedness_determined", "alignment_completed", "reference_ready")
    settings = ("featurecounts_parameters",)
    result_settings = ("featurecounts_parameters",)


    def required_tools(self, ctx):
        return ("featureCounts",)
    def overview(self, ctx):
        fc = ctx.cfg["featurecounts_parameters"]
        s = (ctx.project.state.get("strandedness") or {}).get("value")
        return [("Samples", len(ctx.samples)), ("Paired-end", ctx.project.is_paired()),
                ("Strandedness", f"{s} (-s {strandedness.FEATURECOUNTS_FLAG.get(s, '?')})"),
                ("Feature / attribute", f"{fc['feature_type']} / {fc['attribute']}"),
                ("Multimappers", "counted" if fc["count_multimapping"] else "not counted")]

    def execute(self, ctx):
        p = ctx.project
        strand = (p.state.get("strandedness") or {}).get("value")
        if strand not in strandedness.FEATURECOUNTS_FLAG:
            raise PipelineError("strandedness not established", stage=self.key)
        bams = [p.path("alignment", "bam", f"{sid}.sorted.bam") for sid in ctx.samples]
        txt, summ, cmd = featurecounts.run(p, ctx.samples, bams, ctx.ref()["gtf"], strand, p.is_paired(),
                                           ctx.cfg["featurecounts_parameters"], ctx.threads, ctx.log("featurecounts"))
        if ctx.dry_run:
            return [], {}, {}
        genes, counts = featurecounts.parse(txt, bams)
        exp = ctx.ref().get("genes")
        if exp and len(genes) != exp:
            raise PipelineError(f"featureCounts reported {len(genes)} genes but the GTF has {exp}", stage=self.key)
        stats = featurecounts.parse_summary(summ, bams, ctx.samples)
        fragments = {s: (p.samples[s]["trimmed"]["expected_reads_per_file"] if p.samples[s].get("trimmed")
                         else (p.samples[s].get("reads") or {}).get("R1")) for s in ctx.samples}
        aln_file = p.path("alignment", "reports", "alignment_summary.tsv")
        aln_rows = {r["sample"]: r for r in csv.DictReader(open(aln_file), delimiter="\t")} if aln_file.exists() else {}
        errors, warns = featurecounts.reconcile(stats, fragments, aln_rows)
        for w in warns:
            ui.warn(w)
        if errors:
            raise PipelineError("featureCounts totals are inconsistent with the alignments:\n  - "
                                + "\n  - ".join(errors), stage=self.key,
                                remedy="re-run alignment and gene counting; if it persists, report it with the logs")
        ui.ok("featureCounts totals reconcile with the BAM files (independent check)")
        ui.table([(s, f"{v['assigned']:,}", f"{v['assigned_pct']}%") for s, v in stats.items()],
                 ["Sample", "Assigned", "Assigned %"])
        low = [s for s, v in stats.items() if v["assigned_pct"] < featurecounts.LOW_ASSIGNED_WARN]
        if low:
            ui.warn(f"low assignment rate (<{featurecounts.LOW_ASSIGNED_WARN}%) for {low}: check strandedness, "
                    "annotation and rRNA contamination")
        (p.path("featurecounts", "featurecounts_stats.json")).write_text(json.dumps(stats, indent=2))
        ui.ok(f"featureCounts output validated: {len(genes):,} genes x {len(bams)} samples")
        return [txt, summ, p.path("featurecounts", "featurecounts_stats.json")], \
            {"samples": ctx.samples, "strand": strand, "command": runner.fmt(cmd)}, {"genes": len(genes)}


class CountMatrixStage(Stage):
    key, title, depends = "count_matrix_completed", "GENE COUNT MATRIX", ("featurecounts_completed",)

    def execute(self, ctx):
        if ctx.dry_run:
            ui.status("DRY-RUN", "would build counts/gene_count_matrix.tsv/.csv from featureCounts output")
            return [], {}, {}
        p = ctx.project
        bams = [p.path("alignment", "bam", f"{sid}.sorted.bam") for sid in ctx.samples]
        genes, counts = featurecounts.parse(p.path("featurecounts", "featurecounts.txt"), bams)
        b2s = {str(b): s for b, s in zip(bams, ctx.samples)}
        matrix = count_matrix.build(genes, counts, b2s, ctx.samples)
        fc_stats = json.loads(p.path("featurecounts", "featurecounts_stats.json").read_text())
        mismatch = count_matrix.check_against_assigned(matrix, {s: fc_stats[s]["assigned"] for s in ctx.samples})
        if mismatch:
            raise PipelineError("count matrix does not add up to featureCounts' assigned reads:\n  - "
                                + "\n  - ".join(mismatch), stage=self.key,
                                remedy="re-run gene counting; the featureCounts output may have been modified")
        groups = {}
        for sid in ctx.samples:
            groups.setdefault(p.samples[sid].get("bio_unit") or sid, []).append(sid)
        collapsed = False
        extra = []
        if any(len(v) > 1 for v in groups.values()):
            ui.explain(
                "whether several sequencing runs of the same library are counted as one sample",
                "runs/lanes of one library are technical replicates; treating them as separate samples inflates "
                "the replicate count and makes DESeq2 over-confident (too many false positives)",
                "summing is the standard handling; the per-run matrix is kept for inspection")
            ui.section("TECHNICAL REPLICATES / LANES DETECTED")
            for u, runs in groups.items():
                if len(runs) > 1:
                    print(f"  {u}: {', '.join(runs)}")
            print("Runs of the same library (lanes / technical replicates) are normally SUMMED before DESeq2.")
            if ui.ask_yes_no("Sum counts of runs belonging to the same biological sample?", default=True):
                count_matrix.write(genes, matrix, p.path("counts", "gene_count_matrix_per_run.tsv"),
                                   p.path("counts", "gene_count_matrix_per_run.csv"))
                matrix = count_matrix.collapse(genes, matrix, {count_matrix_safe(u): r for u, r in groups.items()})
                collapsed = True
                extra.append("Technical replicates summed: " +
                             "; ".join(f"{u} = {'+'.join(r)}" for u, r in groups.items() if len(r) > 1))
        tsv, csvp = count_matrix.write(genes, matrix, p.path("counts", "gene_count_matrix.tsv"),
                                       p.path("counts", "gene_count_matrix.csv"))
        info = count_matrix.validate_file(tsv, expected_samples=list(matrix), expected_genes=len(genes))
        expected_units = {u: sum(fc_stats[r]["assigned"] for r in runs) for u, runs in
                          ({count_matrix_safe(u): r for u, r in groups.items()} if collapsed
                           else {s: [s] for s in ctx.samples}).items()}
        written = count_matrix.check_against_assigned({s: [info["library_sizes"][s]] for s in info["samples"]},
                                                      expected_units)
        if written:
            raise PipelineError("written count matrix does not match featureCounts: " + "; ".join(written),
                                stage=self.key)
        ui.ok("count matrix column sums equal featureCounts assigned reads for every sample (independent check)")
        units = list(matrix)
        p.state["count_units"] = units
        p.state["count_unit_runs"] = {count_matrix_safe(u): r for u, r in groups.items()} if collapsed \
            else {s: [s] for s in ctx.samples}
        p.save()
        sinfo = p.path("counts", "sample_info.tsv")
        with open(sinfo, "w") as f:
            f.write("sample\truns\tsource\tdescription\n")
            for u in units:
                runs = p.state["count_unit_runs"][u]
                rec = p.samples[runs[0]]
                meta = rec.get("metadata") or {}
                desc = meta.get("geo_title") or meta.get("sample_title") or ""
                f.write(f"{u}\t{','.join(runs)}\t{rec.get('source')}\t{desc}\n")
        summ = count_matrix.write_summary(p.path("counts", "count_matrix_summary.txt"), info, extra)
        ui.ok(f"count matrix validated: {info['genes']:,} genes x {len(units)} samples "
              f"({info['all_zero_genes']:,} genes with zero counts everywhere)")
        return [tsv, csvp, summ, sinfo], {"samples": ctx.samples, "collapsed": collapsed}, \
            {"genes": info["genes"], "units": units}


def count_matrix_safe(u):
    return data_manager.sanitize_sample_name(u)


class DesignStage(Stage):
    key, title, depends = "design_confirmed", "EXPERIMENTAL DESIGN / METADATA CONFIRMATION", ("count_matrix_completed",)
    settings = ("design_formula", "reference_level")

    def check_inputs(self, ctx):
        return [] if ctx.envs.rscript() else ["Rscript not found (the R environment is needed to validate the "
                                              "design; install it via main menu 4)"]

    def execute(self, ctx):
        if ctx.dry_run:
            ui.status("DRY-RUN", "design must be confirmed interactively before DESeq2")
            return [], {}, {}
        p = ctx.project
        units = p.state["count_units"]
        info = {}
        for u in units:
            rec = p.samples[p.state["count_unit_runs"][u][0]]
            info[u] = {k: v for k, v in (rec.get("metadata") or {}).items() if isinstance(v, str)}
        if not design.interactive(p, units, info, ctx.cfg):
            raise UserAbort("design not confirmed yet")
        d = design.confirmed_design(p)
        meta = p.abs(d["metadata_file"])
        problems = count_matrix.check_metadata_match(units, d["samples"])
        if problems:
            raise PipelineError("metadata does not match count matrix: " + "; ".join(problems), stage=self.key)
        params = r_bridge.build_params(p, ctx.cfg, d)
        pp = r_bridge.write_params(p, params)
        rs = ctx.envs.rscript()
        if rs:
            r_bridge.run(rs, pp, ctx.log("deseq2"), validate_only=True)
        return [meta, p.path("counts", "design.json")], {"formula": d["formula"]}, d


class DESeq2Stage(Stage):
    key, title, depends, expensive = "deseq2_completed", "R / DESeq2 DIFFERENTIAL EXPRESSION", \
        ("design_confirmed", "count_matrix_completed"), True
    result_settings = ("alpha", "log2fc_threshold", "lfc_test_threshold", "min_count_filter", "lfc_shrinkage",
                       "transformation", "top_n_genes", "heatmap_max_genes", "independent_filtering",
                       "cooks_cutoff", "plot_formats", "plot_dpi", "annotation")
    settings = ("alpha", "log2fc_threshold", "lfc_test_threshold", "min_count_filter", "lfc_shrinkage",
                "transformation", "top_n_genes", "heatmap_max_genes", "independent_filtering", "cooks_cutoff")

    def check_inputs(self, ctx):
        probs = []
        if not ctx.envs.rscript():
            probs.append("Rscript not found (install the R environment via main menu 4)")
        if not ctx.dry_run and not design.confirmed_design(ctx.project):
            probs.append("experimental design is not confirmed (or metadata changed since confirmation)")
        return probs

    def execute(self, ctx):
        p, cfg = ctx.project, ctx.cfg
        d = design.confirmed_design(p) if not ctx.dry_run else (p.state.get("design") or {})
        s = p.state.get("strandedness") or {}
        ui.explain(
            "starting the statistical test (DESeq2) with the design and thresholds below",
            "these settings define which genes are called up- or downregulated",
            "changing the design or a threshold later re-runs only DESeq2 and the report")
        ui.section("PYTHON PIPELINE COMPLETE")
        ui.kv([("Count matrix", p.rel(p.path("counts", "gene_count_matrix.tsv"))),
               ("Metadata", d.get("metadata_file")), ("Design", d.get("formula")),
               ("Samples", len(d.get("samples", []))),
               ("Contrasts", ", ".join(f"{a} vs {b}" for a, b in d.get("contrasts", []))),
               ("Reference", (p.state.get("reference") or {}).get("label")),
               ("Strandedness", f"{s.get('value')} ({s.get('source')})"),
               ("Thresholds", f"padj < {cfg['alpha']}, |log2FC| >= {cfg['log2fc_threshold']}"),
               ("Low-count filter", f">= {cfg['min_count_filter']['min_count']} counts in >= "
                                    f"{cfg['min_count_filter']['min_samples']} samples")])
        (p.path("reports", "transition_record.txt")).write_text(
            f"PYTHON PIPELINE COMPLETE\nCount matrix: counts/gene_count_matrix.tsv\nMetadata: {d.get('metadata_file')}\n"
            f"Design: {d.get('formula')}\nSamples: {', '.join(d.get('samples', []))}\n"
            f"Reference: {(p.state.get('reference') or {}).get('label')}\nStrandedness: {s.get('value')}\n")
        if not ctx.dry_run and not ui.ask_yes_no("Proceed to DESeq2?", default=True):
            raise UserAbort("DESeq2 not started")
        if ctx.dry_run:
            runner.run([ctx.envs.rscript() or "Rscript", "--vanilla", r_bridge.R_MAIN, "deseq2_params.json"],
                       stage=self.key)
            return [], {}, {}
        orgdb = None
        r = p.state.get("reference") or {}
        if cfg.get("annotation", {}).get("enabled"):
            orgdb = C.load_catalog()["organisms"].get(r.get("organism"), {}).get("bioc_orgdb")
        params = r_bridge.build_params(p, cfg, d, orgdb)
        pp = r_bridge.write_params(p, params)
        r_bridge.run(ctx.envs.rscript(), pp, ctx.log("deseq2"))
        summary = r_bridge.validate_outputs(p, params)
        for name, c in summary["contrasts"].items():
            for kind in ("upregulated", "downregulated"):
                link = p.path("results", kind, f"{name}_{kind}.tsv")
                target = Path(c["dir"]) / f"{kind}.tsv"
                if link.is_symlink():
                    link.unlink()
                if not link.exists():
                    link.symlink_to(os.path.relpath(target, link.parent))
        ui.section("DIFFERENTIAL EXPRESSION RESULTS")
        ui.table([(n, c["genes_tested"], c["significant"], c["up"], c["down"]) for n, c in summary["contrasts"].items()],
                 ["Contrast", "Genes tested", "Significant", "Up", "Down"])
        ui.info(f"thresholds: padj < {cfg['alpha']} and |log2FC| >= {cfg['log2fc_threshold']}")
        outs = [p.path("results", "deseq2"), p.path("results", "plots"), p.path("results", "tables")]
        return outs, {"params": params}, {"contrasts": {n: {k: c[k] for k in ("significant", "up", "down")}
                                                        for n, c in summary["contrasts"].items()}}

    def revalidate(self, ctx, data):
        """Re-derive the DESeq2 results from the count matrix and design again (r_bridge.independent_checks)."""
        pf = ctx.project.path("results", "deseq2", "deseq2_params.json")
        try:
            params = json.loads(pf.read_text())
        except (OSError, json.JSONDecodeError) as e:
            return [f"DESeq2 parameter file unreadable: {e}"]
        try:
            r_bridge.validate_outputs(ctx.project, params)
        except PipelineError as e:
            return [str(e).replace("\n  - ", "; ")]
        return []


class ReportStage(Stage):
    key, title, depends = "report_generated", "FINAL REPORT + MANIFEST", ("deseq2_completed",)

    def execute(self, ctx):
        from . import manifest, report_manager
        if ctx.dry_run:
            ui.status("DRY-RUN", "would write reports/final_pipeline_report.html and pipeline_manifest/")
            return [], {}, {}
        mfiles = manifest.write(ctx)
        rep = report_manager.generate(ctx)
        problems = report_manager.validate(ctx, rep)
        if problems:
            raise PipelineError("the generated report failed validation:\n  - " + "\n  - ".join(problems),
                                stage=self.key, remedy="re-run the report step; if it persists, report it with logs")
        ui.ok("report validated: links resolve; key numbers (samples, input reads, alignment rate, mapped reads/%, assigned and counted reads, DE counts, thresholds, versions) match their source files")
        ui.ok(f"final report: {rep}")
        return [rep] + mfiles, {}, {}

    def revalidate(self, ctx, data):
        from . import manifest, report_manager
        return manifest.problems(ctx.project) + report_manager.validate(ctx)


STAGES = [DataStage(), FastqVerifyStage(), RawQCStage(), QualityGateStage(), TrimmingStage(), ReferenceStage(),
          AlignmentStage(), BamQCStage(), StrandednessStage(), StringTieStage(), FeatureCountsStage(),
          CountMatrixStage(), DesignStage(), DESeq2Stage(), ReportStage()]
BY_KEY = {s.key: s for s in STAGES}
for _i, _s in enumerate(STAGES, 1):
    # failure screens show the step the user sees; modules tag errors with the key or its first word
    ui.STAGE_TITLES[_s.key] = f"STEP {_i}: {_s.title}"
    ui.STAGE_TITLES.setdefault(_s.key.split("_")[0], ui.STAGE_TITLES[_s.key])
ui.STAGE_TITLES.update(count_matrix=ui.STAGE_TITLES["count_matrix_completed"],
                       bam_qc=ui.STAGE_TITLES["bam_qc_completed"])


# ============================================================================ status / resume

def settings_snapshot(cfg, stage):
    return {k: C.get(cfg, k) for k in stage.result_settings}


def settings_changes(ctx, stage, data):
    """Human-readable list of result-relevant settings that differ from those the checkpoint was made with."""
    recorded = (data.get("params") or {}).get("settings")
    if recorded is None:  # checkpoint written by a version that did not record settings
        return []
    now = settings_snapshot(ctx.cfg, stage)
    changes = []
    for k in sorted(set(recorded) | set(now)):
        if json.dumps(recorded.get(k), sort_keys=True, default=str) != json.dumps(now.get(k), sort_keys=True,
                                                                                  default=str):
            changes.append(f"setting '{k}' changed ({_short(recorded.get(k))} -> {_short(now.get(k))})")
    return changes


def _short(v):
    s = json.dumps(v, default=str)
    return s if len(s) <= 40 else s[:37] + "..."


def stage_status(ctx, stage, deep=False, memo=None):
    """VALID / INVALID(reason) / PENDING. A stage is VALID only if all upstream stages are VALID too."""
    memo = {} if memo is None else memo
    if stage.key in memo:
        return memo[stage.key]
    data = ctx.cp.read(stage.key)
    if data is None:
        memo[stage.key] = ("PENDING", None)
        return memo[stage.key]
    for dep in stage.depends:
        s, _ = stage_status(ctx, BY_KEY[dep], deep, memo)
        if s != "VALID":
            memo[stage.key] = ("INVALID", f"upstream stage '{dep}' is {s}")
            return memo[stage.key]
    probs = settings_changes(ctx, stage, data) + ctx.cp.verify_files(stage.key, deep=deep)
    cleaned = set(ctx.project.state.get("cleaned_files", []))
    probs = [x for x in probs if not any(x.endswith(c) for c in cleaned)]
    rec_samples = set((data.get("params") or {}).get("samples") or [])
    if rec_samples and not set(ctx.samples) <= rec_samples:
        probs.append("sample selection changed (new samples not yet processed)")
    for dep in stage.depends:
        if ctx.cp.read(dep) is None:
            probs.append(f"upstream stage {dep} has no checkpoint")
    if not probs:
        try:
            probs = stage.revalidate(ctx, data)
        except Exception as e:  # fail-safe: an unexpected error marks the stage INVALID (re-run), never VALID
            ui.log.debug("revalidation of %s raised", stage.key, exc_info=True)
            probs = [f"could not re-check outputs ({type(e).__name__}: {e}); the stage will be re-run"]
    memo[stage.key] = ("VALID", None) if not probs else ("INVALID", "; ".join(probs[:3]))
    return memo[stage.key]


def status_table(ctx, deep=False):
    rows, memo = [], {}
    for i, st in enumerate(STAGES, 1):
        s, why = stage_status(ctx, st, deep, memo)
        rows.append((i, st.title, s, why or ""))
    return rows


def show_status(ctx, deep=False):
    ui.section(f"PROJECT STATUS — {ctx.project.name}")
    rows = status_table(ctx, deep)
    ui.table(rows, ["#", "Stage", "Status", "Detail"], max_col=60)
    return rows


# ============================================================================ interactive runner

def stage_screen(ctx, idx, stage):
    """Returns 'run', 'menu' or 'cancel'."""
    while True:
        ui.header(f"STEP {idx} — {stage.title}")
        pairs = stage.overview(ctx) or [("Samples", len(ctx.samples))]
        est = stage.estimate_gb(ctx)
        if est:
            pairs.append(("Estimated storage", f"{est:.1f} GB"))
        pairs.append(("CPU / threads", f"{ctx.sysinfo['cpu_cores']} cores / {ctx.threads} threads"))
        ui.kv(pairs)
        problems = stage.check_inputs(ctx)
        for pr in problems:
            ui.error(pr)
        c = ui.choose(None, ["Start " + stage.title.split(" (")[0].lower(), "Review settings", "Change settings",
                             "Validate inputs again", "Return to the project menu", "Cancel"])
        if c == 0:
            if problems:
                ui.error("this step cannot start until the problem above is fixed — returning to the project menu")
                return "menu"
            return "run"
        if c == 1:
            show_settings(ctx, stage)
        elif c == 2:
            change_settings(ctx, stage)
        elif c == 3:
            probs = stage.check_inputs(ctx)
            for dep in stage.depends:
                s, why = stage_status(ctx, BY_KEY[dep])
                if s != "VALID":
                    probs.append(f"upstream {dep}: {s} {why or ''}")
            if probs:
                for pr in probs:
                    ui.error(pr)
            else:
                ui.ok("inputs validated")
        elif c == 4:
            return "menu"
        else:
            return "cancel"


def show_settings(ctx, stage):
    ui.section("SETTINGS")
    if not stage.settings:
        print("  (no configurable settings for this stage)")
    for key in stage.settings:
        print(f"  {key}: {json.dumps(C.get(ctx.cfg, key), default=str)}")


def change_settings(ctx, stage):
    if not stage.settings:
        ui.info("no configurable settings for this stage")
        return
    keys = [k for k in stage.settings if not isinstance(C.get(ctx.cfg, k), dict)]
    for k in stage.settings:
        v = C.get(ctx.cfg, k)
        if isinstance(v, dict):
            keys += [f"{k}.{sub}" for sub, sv in v.items() if not isinstance(sv, dict)]
    i = ui.choose("Which setting?", [f"{k} = {json.dumps(C.get(ctx.cfg, k), default=str)}" for k in keys])
    key = keys[i]
    old = C.get(ctx.cfg, key)
    raw = ui.ask(f"New value for {key}", default=json.dumps(old, default=str))
    import yaml
    try:
        new = yaml.safe_load(raw)
    except yaml.YAMLError:
        new = raw  # not YAML syntax: treat as a plain string; validation below decides if it is acceptable
    trial = json.loads(json.dumps(ctx.cfg, default=str))
    C.set_value(trial, key, new)
    errs = C.validate(trial, ctx.sysinfo["cpu_cores"])
    if errs:
        for e in errs:
            ui.error(e)
        return
    C.set_value(ctx.cfg, key, new)
    ctx.project.save_config(ctx.cfg)
    ui.ok(f"{key} set to {new!r} (saved to config/project_config.yaml)")
    if ctx.cp.exists(stage.key):
        ui.warn("this stage already has a checkpoint; re-run it for the new setting to take effect")


def run_stage(ctx, idx, stage):
    """Execute one stage with full START/VALIDATION/EXECUTION/POST-VALIDATION/checkpoint semantics."""
    ui.status("RUNNING", f"STEP {idx}: {stage.title} — START")
    for dep in stage.depends:
        s, why = stage_status(ctx, BY_KEY[dep])
        if s != "VALID" and not ctx.dry_run:
            raise PipelineError(f"upstream stage '{dep}' is {s}" + (f": {why}" if why else ""), stage=stage.key,
                                remedy="resume the project so earlier stages are completed first")
    missing = [t for t in stage.required_tools(ctx) if runner.which(t) is None]
    if missing:
        raise PipelineError(f"required software not found for this step: {', '.join(missing)}", stage=stage.key,
                            cause="the scientific tools environment is missing or incomplete",
                            remedy="run 'rnaseq-pipeline --check' for details, then main menu 4 "
                                   "(Manage Dependencies) to install them")
    unusable = dependency_manager.unusable(stage.required_tools(ctx))
    if unusable:
        raise PipelineError("installed software cannot be used for this step: " + "; ".join(unusable),
                            stage=stage.key, cause="an outdated or incompatible build is installed; results from it "
                                                   "would not match the tested pipeline",
                            remedy="run 'rnaseq-pipeline --check' for the exact fix, or main menu 4 "
                                   "(Manage Dependencies); nothing was run")
    if not ctx.dry_run:
        retry = [s for s, r in ctx.project.samples.items()
                 if r.get("status") == "FAILED" and r.get("failed_stage") == stage.key]
        if retry:
            for s in retry:
                rec = ctx.project.samples[s]
                rec.update(status="OK", previous_failure=rec.get("fail_reason"), fail_reason=None)
                rec.pop("failed_stage", None)
            ctx.project.save()
            ui.info(f"retrying {len(retry)} sample(s) that failed in this step last time")
    probs = stage.check_inputs(ctx)
    if probs:
        raise PipelineError("; ".join(probs), stage=stage.key)
    ui.ok("inputs validated")
    est = stage.estimate_gb(ctx)
    if est >= 1 and not ctx.dry_run:
        if not storage.preflight(ctx.project.root, est, ctx.cfg):
            raise UserAbort("insufficient/borderline storage")
    try:
        outputs, params, summary = stage.execute(ctx)
    except PermissionError as e:
        where = e.filename or "a project file"
        raise PipelineError(f"no permission to write {where}", stage=stage.key,
                            cause="the project folder (or part of it) is read-only for your user",
                            remedy=f"make it writable (for example: chmod -R u+w '{ctx.project.root}') or copy the "
                                   "project to a writable location, then resume; outputs of this step were not "
                                   "accepted") from e
    if ctx.dry_run:
        ui.status("DRY-RUN", f"{stage.title}: no outputs written, no checkpoint")
        return
    params = dict(params or {})
    params["settings"] = settings_snapshot(ctx.cfg, stage)
    ctx.cp.write(stage.key, outputs, params=params, depends_on=stage.depends, summary=summary)
    ctx.project.log_event(f"stage {stage.key} completed")
    ui.ok(f"STEP {idx}: {stage.title} — SUCCESS (checkpoint written)")


def run_pipeline(ctx, from_stage=None, only=None):
    """Walk the state machine; skip valid stages; stop on failure/abort. Returns True when all done."""
    start = 0 if from_stage is None else [s.key for s in STAGES].index(from_stage)
    for i, st in enumerate(STAGES[start:], start + 1):
        if only and st.key != only:
            continue
        s, why = stage_status(ctx, st)
        if s == "VALID":
            ui.ok(f"STEP {i}: {st.title} — already complete (checkpoint revalidated)")
            continue
        if s == "INVALID":
            ui.warn(f"STEP {i}: checkpoint no longer valid ({why}); stage will be re-run")
        if not ctx.dry_run and not (ctx.auto and not st.expensive):
            choice = stage_screen(ctx, i, st)
            if choice != "run":
                return False
        try:
            run_stage(ctx, i, st)
        except KeyboardInterrupt as e:
            print()
            how = f"stopped by {e.signame}" if isinstance(e, Terminated) else "interrupted"
            ui.warn(f"{how} during {st.title}; running tools were stopped and partial outputs (*.partial) will "
                    "not be used. Resume the project to continue.")
            ctx.project.log_event(f"stage {st.key} {how}")
            if isinstance(e, Terminated):
                raise
            return False
        except UserAbort as e:
            ui.info(f"stopped: {e}")
            return False
        except PipelineError as e:
            if not e.stage:
                e.stage = st.key
            ui.explain_failure(e)
            ctx.project.log_event(f"stage {st.key} failed: {e}")
            return False
    if ctx.dry_run:
        ui.ok("dry run complete — nothing was executed or written")
        return True
    return True
