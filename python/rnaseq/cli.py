"""Interactive terminal application: main menu, project wizard, resume/validate, dependencies, references."""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from . import (PIPELINE_ROOT, PipelineError, UserAbort, __version__, data_manager, dependency_manager, design,
               environment_manager, logger, reference_manager, runner, storage, system_check, ui, workflow)
from . import config as C
from . import validators as V
from .project import DEFAULT_PROJECTS_DIR, Project

BANNER = f"""
===========================================================
        BULK RNA-seq ANALYSIS PIPELINE
                 Version {__version__}
===========================================================

Language architecture:
Python -> preprocessing / QC / alignment / quantification
R      -> DESeq2 / differential expression / downstream analysis
"""


class App:
    def __init__(self, args):
        self.args = args
        self.cfg = C.load()
        if args.config:
            self.cfg = C.deep_merge(self.cfg, C.load_yaml(V.existing_file(args.config)))
            C.require_valid(self.cfg, os.cpu_count())
        self.envs = environment_manager.Environments(self.cfg)
        self.envs.activate()
        self.projects_dir = Path(args.projects_dir).expanduser() if args.projects_dir else DEFAULT_PROJECTS_DIR

    # ------------------------------------------------------------------ helpers
    def sysinfo(self, path=None):
        return system_check.collect(path or self.projects_dir, r_bin=self.envs.rscript(),
                                    conda_bin=self.envs.conda_bin())

    def context(self, project):
        cfg = project.config()
        C.require_valid(cfg, os.cpu_count())
        envs = environment_manager.Environments(cfg)
        envs.activate()
        info = self.sysinfo(project.root)
        ctx = workflow.Context(project, cfg, envs, info, dry_run=self.args.dry_run, auto=self.args.auto)
        return ctx

    def open_project(self, root):
        p = Project.open(root)
        if not self.args.dry_run:
            p.lock()
        logger.attach_project(p.path("logs"))
        ui.info(f"opened project {p.name} ({p.root})")
        return p

    # ------------------------------------------------------------------ main menu
    def main_menu(self):
        print(BANNER)
        while True:
            try:
                c = ui.choose("MAIN MENU", ["Start New RNA-seq Project", "Resume Existing Project",
                                            "Validate Existing Project", "Manage Dependencies / Environment",
                                            "Manage Reference Databases", "View Previous Reports",
                                            "Pipeline Configuration", "Exit"])
                if c == 7:
                    print("Goodbye.")
                    return 0
                [self.new_project, self.resume_project, self.validate_project, self.manage_dependencies,
                 self.manage_references, self.view_reports, self.configuration][c]()
            except UserAbort as e:
                ui.info(f"returned to main menu ({e})")
                if str(e) == "input closed":
                    return 0
            except PipelineError as e:
                ui.explain_failure(e)
            finally:
                logger.detach_project()

    # ------------------------------------------------------------------ 1. new project
    def new_project(self):
        ui.header("NEW RNA-seq PROJECT")
        info = self.sysinfo()
        system_check.display(info)
        crit, warns = system_check.assess(info, cfg=self.cfg)
        for w in warns:
            ui.warn(w)
        if crit:
            for c in crit:
                ui.error(c)
            raise PipelineError("critical infrastructure problem — cannot create a project here",
                                remedy="fix the problems above (or choose another location) and retry")
        if not self.check_dependencies(quiet=True):
            if not ui.ask_yes_no("Some required software is missing. Continue creating the project anyway? "
                                 "(tools can be installed later from main menu 4)", default=True):
                return
        name = ui.ask("Project name", validator=V.project_name)
        loc = ui.ask("Parent directory for the project", default=str(self.projects_dir),
                     validator=lambda v: str(Path(os.path.expanduser(v)).resolve()))
        linfo = self.sysinfo(loc)
        crit, warns = system_check.assess(linfo, cfg=self.cfg)
        for w in warns:
            ui.warn(w)
        if crit:
            for c in crit:
                ui.error(c)
            raise PipelineError(f"location {loc} is not suitable", remedy="choose a directory on an ext4/xfs disk")
        p = Project.create(name, loc, base_config=self.cfg)
        p.lock()
        logger.attach_project(p.path("logs"))
        (p.path("pipeline_manifest", "system_information.txt")).write_text(
            "\n".join(f"{k}: {v}" for k, v in linfo.items()) + "\n")
        ctx = self.context(p)
        ctx.cp.write("project_initialized", [p.path("pipeline_manifest", "system_information.txt")],
                     params={"name": name})
        ui.ok(f"project created: {p.root}")
        self.choose_data(p)
        self.select_reference(p)
        self.optional_databases(p)
        self.select_samples(p)
        if ui.ask_yes_no("Start the pipeline now?", default=True):
            workflow.run_pipeline(self.context(p))
        p.unlock()

    # ------------------------------------------------------------------ data input
    def choose_data(self, p):
        while True:
            c = ui.choose("DATA INPUT", ["NCBI SRA accession(s)", "NCBI GEO accession", "ENA accession(s)",
                                         "Local FASTQ files", "Existing project data"])
            try:
                if c in (0, 1, 2):
                    self._public(p, {0: "sra", 1: "geo", 2: "ena"}[c])
                elif c == 3:
                    self._local(p)
                else:
                    samples = data_manager.existing_project_fastqs(p)
                    self._confirm_local(p, samples, mode="symlink", source="existing")
                return
            except (PipelineError, ValueError) as e:
                ui.error(str(e))
                if not ui.ask_yes_no("Try again?", default=True):
                    raise UserAbort("data input cancelled")

    def _public(self, p, kind):
        label = {"sra": "SRA run/experiment/study (SRR/SRX/SRP/PRJNA...)", "geo": "GEO series (GSE...)",
                 "ena": "ENA run/experiment/study (ERR/SRR/PRJEB/PRJNA...)"}[kind]
        accs = ui.ask(f"{label}; separate multiple with spaces or commas", validator=V.accession_list)
        if kind == "geo" and any(k != "geo_series" for _, k in accs):
            raise ValueError("enter GEO series accessions (GSE...) for this option")
        if kind != "geo" and any(k.startswith("geo") for _, k in accs):
            raise ValueError("GEO accessions belong to option 2")
        route = "sra" if kind == "sra" else self.cfg["download"]["source_preference"]
        ui.running("Retrieving metadata (ENA Portal / NCBI GEO)...")
        runs, titles = data_manager.resolve_public(accs, route)
        data_manager.show_public(runs, titles)
        print()
        idx = ui.choose_many("Select runs to include", [f"{r['run_accession']}  "
                                                       f"{r.get('_geo', {}).get('title') or r.get('sample_title', '')}"
                                                       for r in runs])
        runs = [runs[i] for i in idx]
        total = sum(f[2] or 0 for r in runs for f in r.get("_files", [])) / 1e9
        ui.info(f"{len(runs)} run(s), approx. {total:.1f} GB of compressed FASTQ")
        mode = ui.choose("Download scope", ["Full data (for real analysis)",
                                            "Pilot subset: first N reads per run (TESTING ONLY; recorded in manifest)"])
        cfg = p.config()
        if mode == 1:
            n = ui.ask("Reads per run", default="500000", validator=lambda v: V.positive_int(v, "reads", 1000))
            cfg["download"]["max_reads"] = n
            ui.warn("PILOT MODE: results from a read subset are for testing only")
        else:
            cfg["download"]["max_reads"] = 0
        p.save_config(cfg)
        data_manager.register_public(p, runs, "sra" if route == "sra" else "ena")
        p.state["data_source"] = kind
        p.state["accessions"] = [a for a, _ in accs]
        p.save()
        ui.ok(f"registered {len(runs)} run(s); download happens in the DATA DOWNLOAD step")

    def _local(self, p):
        print("Enter FASTQ files or directories, one per line (blank line to finish).")
        paths = []
        while True:
            v = ui.ask("  path", allow_empty=True, default="")
            if not v:
                break
            paths.append(v)
        if not paths:
            raise ValueError("no paths entered")
        samples = data_manager.pair_files(data_manager.find_fastqs(paths))
        mode = ["symlink", "copy"][ui.choose("Import mode", ["Symlink (no duplication; recommended)",
                                                             "Copy into project (uses extra disk)"])]
        self._confirm_local(p, samples, mode, "local")

    def _confirm_local(self, p, samples, mode, source):
        ui.section("DETECTED SAMPLES")
        ui.table([(s["id"], s["bio_unit"], s["r1"].name, s["r2"].name if s["r2"] else "-") for s in samples],
                 ["Sample", "Biological unit", "R1", "R2"])
        layout = "paired-end" if samples[0]["r2"] else "single-end"
        ui.info(f"{len(samples)} sample(s), {layout}")
        if not ui.ask_yes_no("Is this sample list and read type correct?", default=True):
            if ui.ask_yes_no("Rename samples?", default=False):
                for s in samples:
                    s["id"] = ui.ask(f"  name for {s['r1'].name}", default=s["id"], validator=V.sample_name)
                if len({s['id'] for s in samples}) != len(samples):
                    raise ValueError("sample names must be unique")
            else:
                raise ValueError("sample detection rejected; rename files as <sample>_R1/_R2.fastq.gz and retry")
        if source == "existing":
            for s in samples:
                p.add_sample(s["id"], source="local", layout="PAIRED" if s["r2"] else "SINGLE", bio_unit=s["bio_unit"],
                             r1=p.rel(s["r1"]), r2=p.rel(s["r2"]) if s["r2"] else None)
            p.state["read_type"] = "paired" if samples[0]["r2"] else "single"
        else:
            data_manager.import_local(p, samples, mode)
        p.state["data_source"] = "local"
        p.save()
        ui.ok(f"{len(samples)} sample(s) registered ({mode})")

    # ------------------------------------------------------------------ references
    def select_reference(self, p):
        cat = C.load_catalog()
        ui.header("REFERENCE DATABASE SETUP")
        hint = p.state.get("organism_hint")
        if hint:
            ui.info(f"public metadata organism: {hint}")
        orgs = list(cat["organisms"])
        pcfg = p.config()
        org_default = orgs.index(pcfg["organism"]) if pcfg.get("organism") in orgs else \
            (len(orgs) if pcfg.get("organism") == "custom" else None)
        oi = ui.choose("Organism", [f"{k} ({v['scientific_name']})" for k, v in cat["organisms"].items()]
                       + ["Search NCBI for ANY organism (downloads genome + annotation)",
                          "Custom organism (local genome FASTA + GTF)"], default=org_default)
        if oi == len(orgs):
            ref = self._ncbi_reference(p)
            ref["store"] = str(reference_manager.store_dir(p.config(), ref))
            return self._finish_reference(p, ref)
        if oi < len(orgs):
            org = orgs[oi]
            sci = cat["organisms"][org]["scientific_name"]
            if hint and hint.lower() != sci.lower():
                ui.warn(f"selected organism {sci} differs from the data's organism {hint}")
                if not ui.ask_yes_no("Continue with this mismatch?", default=False):
                    return self.select_reference(p)
        else:
            org = None
        options = ["GENCODE", "Ensembl", "NCBI RefSeq", "UCSC", "Custom local genome + GTF", "Existing reference/index",
                   "Search NCBI for this organism"]
        src_key = {0: "gencode", 1: "ensembl", 2: "refseq", 3: "ucsc"}
        src_default = {"gencode": 0, "ensembl": 1, "refseq": 2, "ucsc": 3, "custom": 4, "existing": 5}.get(
            pcfg.get("reference_source"))
        si = ui.choose("Reference source", options, default=src_default) if org else 4
        if si < 4:
            pk = {k: v for k, v in reference_manager.packages_for(cat, org).items() if v["source"] == src_key[si]}
            if not pk:
                ui.warn(f"no {options[si]} package configured for {org}; add one to config/reference_catalog.yaml "
                        "or use a custom reference")
                return self.select_reference(p)
            keys = list(pk)
            ki = ui.choose("Available packages", [f"{k}: {v['assembly']} / {v['annotation_release']}" for k, v in pk.items()])
            pkg = pk[keys[ki]]
            ref = {"id": keys[ki], "kind": "catalog", "package": pkg, "organism": org, "source": pkg["source"],
                   "assembly": pkg["assembly"], "annotation_assembly": pkg["assembly"],
                   "annotation_release": pkg["annotation_release"],
                   "label": f"{pkg['assembly']} / {pkg['annotation_release']}"}
            self._show_package(pkg)
        elif si == 5:
            ref = self._existing_reference(org)
        elif si == 6:
            ref = self._ncbi_reference(p, cat["organisms"].get(org, {}).get("scientific_name", org))
        else:
            ref = self._custom_reference(org)
        ref["store"] = str(reference_manager.store_dir(p.config(), ref))
        return self._finish_reference(p, ref)

    def _finish_reference(self, p, ref):
        store_state = "prepared (will be revalidated)" if (Path(ref["store"]) / "reference_manifest.yaml").exists() \
            else "not yet prepared"
        ui.kv([("Reference", ref["label"]), ("Stored in", ref["store"]), ("Status", store_state)])
        if not ui.ask_yes_no("Use this reference for the project?", default=True):
            return self.select_reference(p)
        old = p.state.get("reference")
        if old and old.get("id") != ref["id"]:
            ui.warn("reference changed: alignment and all downstream stages will be re-run")
            workflow.Checkpoints(p).invalidate("reference_ready", "reference changed by user")
        p.state["reference"] = ref
        p.save()
        ui.ok(f"reference selected: {ref['label']}")

    def _ncbi_reference(self, p, default_query=None):
        """Search NCBI Datasets for any organism and build a downloadable reference package."""
        default_query = default_query or p.state.get("organism_hint")
        while True:
            q = ui.ask("Organism name, strain or taxid (e.g. 'Escherichia coli UTI89', 'Danio rerio', 7955)",
                       default=default_query)
            ui.running(f"Searching NCBI assemblies for {q!r}...")
            recs = reference_manager.search_ncbi(q, limit=8, reference_only=True)
            if not recs:
                recs = reference_manager.search_ncbi(q, limit=12, reference_only=False)
            recs = [r for r in recs if r["annotation"]] or recs
            recs.sort(key=lambda r: (not r["accession"].startswith("GCF"), r["level"] != "Complete Genome"))
            if not recs:
                ui.warn(f"no assemblies found for {q!r}")
                if not ui.ask_yes_no("Search again?", default=True):
                    raise UserAbort("no reference chosen")
                continue
            ui.section(f"NCBI ASSEMBLIES for {q!r}")
            ui.table([(r["accession"], r["assembly_name"], r["organism"] + (f" ({r['strain']})" if r["strain"] else ""),
                       r["level"], f"{r['genome_bp'] / 1e6:.1f}", r["annotation"] or "NO ANNOTATION")
                      for r in recs],
                     ["Accession", "Assembly", "Organism", "Level", "Mb", "Annotation"], max_col=42)
            i = ui.choose("Which assembly?", [f"{r['accession']}  {r['assembly_name']}  {r['organism']}" for r in recs]
                          + ["Search again"])
            if i == len(recs):
                default_query = q
                continue
            rec = recs[i]
            if not rec["annotation"]:
                ui.error("this assembly has no NCBI annotation (no GTF) — gene counting needs one")
                continue
            ui.running("Checking that the annotation file exists on the NCBI FTP site...")
            if not reference_manager.annotation_available(rec):
                ui.error("NCBI does not publish a GTF for this assembly; choose another one")
                continue
            pkg = reference_manager.ncbi_package(rec)
            self._show_package(pkg)
            return {"id": rec["accession"].replace(".", "_") + "_" + re.sub(r"[^A-Za-z0-9]", "", rec["assembly_name"]),
                    "kind": "catalog", "package": pkg, "organism": rec["organism"], "source": pkg["source"],
                    "assembly": rec["assembly_name"], "annotation_assembly": rec["assembly_name"],
                    "annotation_release": pkg["annotation_release"], "accession": rec["accession"],
                    "label": f"{rec['organism']} {rec['assembly_name']} ({rec['accession']})"}

    def _show_package(self, pkg):
        ui.section("REFERENCE PACKAGE")
        ui.kv([("Assembly", pkg["assembly"]), ("Annotation", pkg["annotation_release"]),
               ("Genome FASTA", f"{pkg['genome']['url']}  (~{pkg['genome']['approx_size_gb']} GB)"),
               ("Annotation GTF", f"{pkg['annotation']['url']}  (~{pkg['annotation']['approx_size_gb']} GB)"),
               ("Checksum", pkg.get("checksum", {}).get("type", "none")),
               ("HISAT2 index (built locally)", f"~{pkg.get('index_approx_size_gb', '?')} GB")])
        ui.info("Files are downloaded only when the REFERENCE PREPARATION step runs (you will confirm again).")

    def _custom_reference(self, org):
        fasta = ui.ask("Genome FASTA (.fa/.fasta[.gz])", validator=lambda v: str(V.existing_file(v)))
        gtf = ui.ask("Annotation GTF (.gtf[.gz])", validator=lambda v: str(V.existing_file(v)))
        organism = org or ui.ask("Organism name (e.g. Arabidopsis_thaliana)", validator=V.project_name)
        ga = ui.ask("Genome assembly name (e.g. TAIR10)", validator=V.project_name)
        aa = ui.ask("Assembly the annotation was made for", default=ga, validator=V.project_name)
        rel = ui.ask("Annotation release/version", default="unknown")
        override = False
        if reference_manager._norm(ga) != reference_manager._norm(aa):
            ui.error(f"genome assembly {ga} and annotation assembly {aa} differ — coordinates will not match")
            if ui.ask("Type OVERRIDE to continue anyway (NOT recommended), or press Enter to cancel",
                      default="", allow_empty=True) != "OVERRIDE":
                raise UserAbort("incompatible reference rejected")
            override = True
        ref = {"id": reference_manager.custom_id(fasta, gtf), "kind": "custom", "fasta": fasta, "gtf": gtf,
               "organism": organism, "source": "custom", "assembly": ga, "annotation_assembly": aa,
               "annotation_release": rel, "label": f"{organism} {ga} / {rel} (custom)",
               "override_compatibility": override}
        return ref

    def _existing_reference(self, org):
        store = Path(os.path.expanduser(self.cfg["reference_store"]))
        found = sorted(d for d in store.glob("*") if (d / "reference_manifest.yaml").exists()) if store.exists() else []
        opts = [f"{d.name}: {C.load_yaml(d / 'reference_manifest.yaml').get('label')}" for d in found]
        i = ui.choose("Prepared references in the shared store", opts + ["Other: give index prefix + FASTA + GTF"])
        if i < len(found):
            m = C.load_yaml(found[i] / "reference_manifest.yaml")
            ref = {"id": found[i].name, "kind": "existing", "organism": m.get("organism"), "source": m.get("source"),
                   "assembly": m.get("genome_assembly"), "annotation_assembly": m.get("annotation_assembly"),
                   "annotation_release": m.get("annotation_release"), "label": m.get("label")}
            ref["fasta"] = m["genome"]["path"]
            ref["gtf"] = m["annotation"]["path"]
            return ref
        ref = self._custom_reference(org)
        prefix = ui.ask("HISAT2 index prefix (path without .1.ht2)",
                        validator=lambda v: v if reference_manager.index_files(Path(os.path.expanduser(v)))
                        else (_ for _ in ()).throw(ValueError("no .1.ht2 ... .8.ht2 files found for this prefix")))
        ref["external_index"] = str(Path(os.path.expanduser(prefix)).resolve())
        return ref

    def optional_databases(self, p):
        ui.section("DATABASES")
        print("PRIMARY DATABASES (required)\n  Genome FASTA\n  Gene annotation GTF\n  Transcript annotation (from GTF)")
        print("\nOPTIONAL DATABASES (only if you request them)\n  Gene ID annotation (Bioconductor OrgDb)\n"
              "  GO (clusterProfiler)\n  Reactome (ReactomePA; human/mouse)")
        cfg = p.config()
        if ui.ask_yes_no("Enable optional GO/Reactome enrichment after DESeq2?", default=False):
            cfg["annotation"]["enabled"] = True
            cfg["annotation"]["go"] = ui.ask_yes_no("  GO biological process enrichment?", True)
            cfg["annotation"]["reactome"] = ui.ask_yes_no("  Reactome pathway enrichment?", False)
            ui.info("Required R packages (~0.3–0.5 GB) can be installed from main menu 4 when needed.")
        p.save_config(cfg)

    def select_samples(self, p):
        if len(p.samples) < 2:
            return
        c = ui.choose("SAMPLES TO PROCESS", ["All samples", "Selected samples (debugging/testing)"])
        if c == 0:
            p.state["selected_samples"] = None
        else:
            ids = list(p.samples)
            idx = ui.choose_many("Samples", ids)
            p.state["selected_samples"] = [ids[i] for i in idx]
        p.save()

    # ------------------------------------------------------------------ 2. resume
    def pick_project(self):
        projects = Project.list_projects(self.projects_dir)
        opts = [f"{d.name}" for d in projects] + ["Enter a project path"]
        i = ui.choose("PROJECTS", opts)
        if i == len(projects):
            return ui.ask("Project directory", validator=lambda v: str(V.existing_dir(v)))
        return projects[i]

    def resume_project(self, root=None):
        p = self.open_project(root or self.pick_project())
        try:
            self.project_menu(p)
        finally:
            p.unlock()

    def project_menu(self, p):
        while True:
            ctx = self.context(p)
            workflow.show_status(ctx)
            c = ui.choose(f"PROJECT: {p.name}", [
                "Continue pipeline (run remaining stages)", "Re-run a specific stage", "Select samples (all / subset)",
                "Samples: status / re-include failed", "Choose data source (new projects only)",
                "Change reference", "Experimental design (edit / confirm)", "Clean up intermediate files",
                "Dry run (show what would be executed)", "Return to main menu"])
            if c == 0:
                workflow.run_pipeline(ctx)
            elif c == 1:
                i = ui.choose("Stage", [s.title for s in workflow.STAGES])
                st = workflow.STAGES[i]
                if ctx.cp.exists(st.key) and ui.ask_yes_no(
                        f"Mark '{st.title}' for re-run? (outputs are kept; downstream stages will re-run)", True):
                    ctx.cp.invalidate(st.key, "re-run requested by user")
                workflow.run_pipeline(self.context(p), from_stage=st.key)
            elif c == 2:
                self.select_samples(p)
            elif c == 3:
                self.samples_menu(p)
            elif c == 4:
                if p.samples:
                    ui.warn("this project already has samples; create a new project for different data")
                else:
                    self.choose_data(p)
            elif c == 5:
                self.select_reference(p)
            elif c == 6:
                if ctx.cp.exists("design_confirmed"):
                    ctx.cp.invalidate("design_confirmed", "design edited by user")
                workflow.run_pipeline(self.context(p), from_stage="design_confirmed", only="design_confirmed")
            elif c == 7:
                self.cleanup(ctx)
            elif c == 8:
                self.dry_run(p)
            else:
                return

    def samples_menu(self, p):
        rows = [(s, r.get("status"), r.get("failed_stage") or "", (r.get("fail_reason") or "")[:60])
                for s, r in p.samples.items()]
        ui.table(rows, ["Sample", "Status", "Failed at", "Reason"])
        bad = [s for s, r in p.samples.items() if r.get("status") in ("FAILED", "EXCLUDED")]
        if bad and ui.ask_yes_no("Re-include FAILED/EXCLUDED samples (after you fixed the cause)?", False):
            idx = ui.choose_many("Samples to re-include", bad)
            for i in idx:
                rec = p.samples[bad[i]]
                rec.update(status="OK", fail_reason=None)
                rec.pop("failed_stage", None)
            p.save()
            ui.ok("samples re-included; affected stages will re-run")
        good = [s for s, r in p.samples.items() if r.get("status") == "OK"]
        if good and ui.ask_yes_no("Exclude samples from the analysis?", False):
            idx = ui.choose_many("Samples to exclude", good)
            for i in idx:
                p.samples[good[i]].update(status="EXCLUDED", fail_reason="excluded by user")
            p.save()

    def cleanup(self, ctx):
        if ctx.cfg.get("storage", {}).get("cleanup_intermediates") == "never":
            ui.info("cleanup disabled by config (storage.cleanup_intermediates: never)")
            return
        cands = storage.cleanup_candidates(ctx.project, ctx.cp)
        if not cands:
            ui.info("no removable intermediate files (raw data, references and BAMs are never offered)")
            return
        ui.section("REMOVABLE INTERMEDIATES (downstream outputs validated)")
        for i, (f, why) in enumerate(cands, 1):
            size = sum(x.stat().st_size for x in ([f] if f.is_file() else f.rglob("*")) if x.is_file()) / 1e9
            print(f"  {i}. {ctx.project.rel(f)}  ({size:.2f} GB) — {why}")
        idx = ui.choose_many("Remove which", [ctx.project.rel(f) for f, _ in cands])
        if not ui.ask_yes_no(f"Permanently delete {len(idx)} item(s)?", False):
            return
        cleaned = ctx.project.state.setdefault("cleaned_files", [])
        for i in idx:
            f = cands[i][0]
            if f.is_dir():
                shutil.rmtree(f)
            else:
                f.unlink()
            cleaned.append(ctx.project.rel(f))
        ctx.project.save()
        ui.ok("cleanup done (recorded in project.json)")

    def dry_run(self, p):
        prev = runner.DRY_RUN
        runner.DRY_RUN = True
        self.args.dry_run = True
        try:
            ui.header("DRY RUN", "commands are displayed, nothing is executed")
            ctx = self.context(p)
            system_check.display(ctx.sysinfo)
            crit, warns = system_check.assess(ctx.sysinfo, cfg=ctx.cfg)
            for w in warns:
                ui.warn(w)
            for c in crit:
                ui.error(c)
            self.check_dependencies(quiet=True)
            for i, st in enumerate(workflow.STAGES, 1):
                s, why = workflow.stage_status(ctx, st)
                if s == "VALID":
                    ui.skipped(f"STEP {i}: {st.title} (complete)")
                    continue
                ui.section(f"STEP {i}: {st.title}")
                ov = st.overview(ctx)
                if ov:
                    ui.kv(ov, indent=2)
                est = st.estimate_gb(ctx)
                if est:
                    ui.info(f"estimated new storage: {est:.1f} GB")
                for pr in st.check_inputs(ctx):
                    ui.warn(pr)
                try:
                    st.execute(ctx)
                except (PipelineError, KeyError, TypeError, FileNotFoundError) as e:
                    ui.info(f"further commands depend on outputs of earlier steps ({type(e).__name__}: {e})")
                    break
        finally:
            runner.DRY_RUN = prev
            self.args.dry_run = prev

    # ------------------------------------------------------------------ 3. validate
    def validate_project(self):
        p = self.open_project(self.pick_project())
        try:
            ctx = self.context(p)
            ui.info("deep validation: re-hashing all recorded outputs and re-running stage integrity checks")
            rows = workflow.show_status(ctx, deep=True)
            bad = [r for r in rows if r[2] == "INVALID"]
            if bad:
                ui.warn(f"{len(bad)} stage(s) invalid; 'Resume Existing Project' will re-run them")
            else:
                ui.ok("all completed stages validated")
            if ctx.cp.exists("fastq_verified") and ui.ask_yes_no(
                    "Also re-read every FASTQ file end-to-end (slow for large data)?", False):
                cache = p.path("data", "metadata", "fastq_validation_cache.json")
                if cache.exists():
                    cache.rename(cache.with_suffix(".json.bak"))
                pairs = {s: p.fastqs(s, trimmed=False) for s in ctx.samples}
                rows = workflow.validate_fastqs(ctx, pairs, p.path("data", "metadata", "fastq_revalidation.tsv"),
                                                "validate")
                fails = workflow._failures_from_rows(rows)
                (ui.ok("all FASTQ files valid") if not fails else ui.error(f"FASTQ problems: {list(fails)}"))
        finally:
            p.unlock()

    # ------------------------------------------------------------------ 4. dependencies
    def check_dependencies(self, quiet=False):
        tools = dependency_manager.detect_tools()
        r_info, r_pkgs = dependency_manager.detect_r(self.envs.rscript())
        if not quiet:
            dependency_manager.report(tools, r_info, r_pkgs)
        bad = dependency_manager.blocking(tools, r_info, r_pkgs)
        if bad:
            ui.warn(f"missing/unusable required software: {', '.join(bad)}")
        else:
            ui.ok("all required software detected")
        return not bad

    def manage_dependencies(self):
        ui.header("DEPENDENCIES / ENVIRONMENT")
        ui.kv([("Conda/mamba", self.envs.conda_bin() or "not found"),
               ("Tools env", self.envs.prefix("tools") or f"{self.envs.names['tools']} (not created)"),
               ("R env", self.envs.prefix("r") or f"{self.envs.names['r']} (not created)")])
        tools = dependency_manager.detect_tools()
        r_info, r_pkgs = dependency_manager.detect_r(self.envs.rscript())
        dependency_manager.report(tools, r_info, r_pkgs)
        c = ui.choose("OPTIONS", ["Install missing components (asks before running anything)",
                                  "Install optional annotation packages (GO/Reactome, OrgDb)",
                                  "Export environment records (environment.yml, versions, system info)",
                                  "Return"])
        log = Path.home() / ".rnaseq_pipeline" / "install.log"
        if c == 0:
            if dependency_manager.interactive_install(self.envs, tools, r_info, r_pkgs, log):
                self.envs = environment_manager.Environments(self.cfg)
                self.envs.activate()
                self.check_dependencies()
        elif c == 1:
            cat = C.load_catalog()
            orgdbs = {k: v["bioc_orgdb"] for k, v in cat["organisms"].items()}
            i = ui.choose("Organism", list(orgdbs))
            pkgs = ["bioconductor-clusterprofiler", "bioconductor-" + orgdbs[list(orgdbs)[i]].lower(),
                    "bioconductor-reactomepa"]
            cmd = self.envs.install_commands("r", pkgs)
            print("Command that would be executed:\n  " + runner.fmt(cmd))
            if ui.ask_yes_no("Proceed? (~0.5 GB)", False):
                self.envs.run_install(cmd, log)
        elif c == 2:
            out = Path.home() / ".rnaseq_pipeline" / "environment_export"
            files = self.envs.export(out, self.sysinfo())
            for f in files:
                ui.ok(str(f))

    # ------------------------------------------------------------------ 5. references
    def manage_references(self):
        ui.header("REFERENCE DATABASES")
        cat = C.load_catalog()
        ui.section("CATALOG (download on demand)")
        ui.table([(k, v["organism"], v["source"], v["assembly"], v["annotation_release"],
                   f"{v['genome']['approx_size_gb'] + v['annotation']['approx_size_gb']:.2f}")
                  for k, v in cat["packages"].items()], ["Package", "Organism", "Source", "Assembly", "Annotation", "GB"])
        store = Path(os.path.expanduser(self.cfg["reference_store"]))
        ui.section(f"PREPARED IN STORE ({store})")
        rows = []
        for d in sorted(store.glob("*")) if store.exists() else []:
            m = d / "reference_manifest.yaml"
            man = C.load_yaml(m) if m.exists() else {}
            idx = "valid" if reference_manager.index_files(d / "index" / "genome") else "missing"
            size = sum(f.stat().st_size for f in d.rglob("*") if f.is_file() and not f.is_symlink()) / 1e9
            rows.append((d.name, man.get("label", "(incomplete)"), idx, f"{size:.1f}"))
        ui.table(rows, ["ID", "Label", "HISAT2 index", "GB"])
        c = ui.choose("OPTIONS", ["Select/change the reference of a project", "Return"])
        if c == 0:
            p = self.open_project(self.pick_project())
            try:
                self.select_reference(p)
            finally:
                p.unlock()

    # ------------------------------------------------------------------ 6. reports
    def view_reports(self):
        rows = []
        for d in Project.list_projects(self.projects_dir):
            r = d / "reports" / "final_pipeline_report.html"
            rows.append((d.name, str(r) if r.exists() else "(no final report yet)"))
        ui.table(rows, ["Project", "Report"], max_col=90)
        have = [r for r in rows if r[1].endswith(".html")]
        if have and ui.ask_yes_no("Open a report in the browser?", False):
            i = ui.choose("Report", [r[0] for r in have])
            opener = shutil.which("wslview") or shutil.which("xdg-open")
            if opener:
                subprocess.Popen([opener, have[i][1]], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                ui.info(f"open this file in a browser: {have[i][1]}")

    # ------------------------------------------------------------------ 7. configuration
    def configuration(self):
        ui.header("PIPELINE CONFIGURATION")
        c = ui.choose(None, ["View default configuration", "Edit a project's configuration", "Return"])
        if c == 0:
            print(C.DEFAULT_CONFIG.read_text())
        elif c == 1:
            p = self.open_project(self.pick_project())
            try:
                cfg = p.config()
                while True:
                    key = ui.ask("Setting to change (dotted key, e.g. alpha or fastp_parameters.length_required; "
                                 "blank to finish)", allow_empty=True, default="")
                    if not key:
                        break
                    if C.get(cfg, key) is None and key.split(".")[0] not in cfg:
                        ui.error(f"unknown setting {key}")
                        continue
                    import yaml
                    raw = ui.ask(f"{key} (current {C.get(cfg, key)!r})")
                    trial = C.deep_merge(cfg, {})
                    C.set_value(trial, key, yaml.safe_load(raw))
                    errs = C.validate(trial, os.cpu_count())
                    if errs:
                        for e in errs:
                            ui.error(e)
                        continue
                    cfg = trial
                    p.save_config(cfg)
                    ui.ok(f"{key} updated")
            finally:
                p.unlock()


def parse_args(argv):
    ap = argparse.ArgumentParser(prog="rnaseq_pipeline", description="Interactive bulk RNA-seq pipeline (v1)")
    ap.add_argument("--dry-run", action="store_true", help="check tools/paths/inputs/disk and show commands only")
    ap.add_argument("--project", help="open this project directly")
    ap.add_argument("--projects-dir", help=f"where projects are listed/created (default {DEFAULT_PROJECTS_DIR})")
    ap.add_argument("--config", help="YAML file whose values override the default configuration for new projects")
    ap.add_argument("--auto", action="store_true", help="skip stage screens for inexpensive stages "
                                                        "(decisions and expensive steps still ask)")
    ap.add_argument("--check", action="store_true", help="run system + dependency checks and exit")
    ap.add_argument("--verbose", "-v", action="store_true", help="show debug output and commands")
    ap.add_argument("--version", action="version", version=f"Pipeline Version: {__version__}")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    ui.VERBOSE = args.verbose or bool(C.load().get("verbose"))
    logger.init_base_logging(args.verbose)
    runner.DRY_RUN = args.dry_run
    app = App(args)
    try:
        if args.check:
            info = app.sysinfo()
            system_check.display(info)
            crit, warns = system_check.assess(info, cfg=app.cfg)
            for w in warns:
                ui.warn(w)
            for c in crit:
                ui.error(c)
            ok = app.check_dependencies()
            return 0 if ok and not crit else 1
        if args.project:
            if args.dry_run:
                p = Project.open(args.project)
                logger.attach_project(p.path("logs"))
                app.dry_run(p)
                return 0
            app.resume_project(Path(args.project))
            return 0
        if args.dry_run:
            ui.info("DRY-RUN mode: no tool will be executed and no data downloaded")
        return app.main_menu()
    except UserAbort:
        print()
        return 0
    except KeyboardInterrupt:
        print("\n[INFO] interrupted by user")
        return 130
    except PipelineError as e:
        ui.explain_failure(e)
        return 1
