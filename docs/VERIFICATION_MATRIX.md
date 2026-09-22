# Verification Matrix

Requirement → documentation claim → implementation → tests → **observed runtime behaviour** → status.
Statuses: `VERIFIED` · `PARTIALLY VERIFIED` · `FAILS` · `MISSING` · `CONTRADICTORY` · `UNTESTED` · `RISK`.

* **BEFORE** = audit of commit `08a9c11` (127/127 tests passing), established by runtime probes on copies of a
  synthetic truth project (6 samples, 16 known DE genes) — not by reading code alone.
* **AFTER** = status after the hardening work, with the test(s) that now prove it.
* Findings are numbered `F#` and listed with severity at the end.

Baseline scientific result (synthetic truth, 40k read pairs/sample): **16/16 true DE genes recovered, 0 false
positives, 58 genes tested, 16 checkpoints**. This is the reference outcome every later change is compared against.

## 1. Workflow stages

| # | Requirement / invariant | Implementation | Existing tests | BEFORE | AFTER |
|---|---|---|---|---|---|
| 1 | Download: ENA MD5 verified, `.part` → atomic rename | `net.download`, `ena_manager.download_run` | none direct | PARTIALLY VERIFIED — MD5 mismatch kept as `.badmd5` (probe) | |
| 1 | Download: resume after interruption | `net.download` (`curl -C -`) | none | **FAILS** when server rejects ranges (curl 33) → aborts instead of restarting (F3) | |
| 1 | Download: empty/short file never accepted | `net.download` | none | **FAILS** — 0-byte file accepted when no checksum (F4) | |
| 1 | Download: actionable network errors | `runner`/`net` | none | RISK — "exit codes [22]: curl"; 404 retried (F5) | |
| 1 | SRA: prefetch → vdb-validate → fasterq-dump | `sra_manager.download_run` | none (manual evidence) | VERIFIED (manual, SRR7361181; corrupt .sra → exit 3) | |
| 1 | Single-cell data refused | `data_manager.single_cell_reason`, SRA >2 reads guard | unit (5 cases) | VERIFIED (SRR23333328 refused) | |
| 1 | One organism per project | `data_manager.show_public` | none | **RISK** — mixed organisms only warned, can be registered (F11) | |
| 2 | FASTQ structure/alphabet/lengths/gzip/mates | `fastq_validator` | 14 unit | VERIFIED | |
| 2 | Validation interruptible (Ctrl-C/SIGTERM) | `fastq_validator.validate_many` | none | RISK — ProcessPool waits for workers (F2) | |
| 3 | FastQC count == FASTQ count (layer 2) | `qc_manager.validate_fastqc` | integration | VERIFIED | |
| 4 | Advisory gate, user decides, thresholds in config | `quality_assessment`, `QualityGateStage` | integration | VERIFIED; explanation of decision minimal (F12) | |
| 5 | fastp output re-validated; counts == fastp JSON | `TrimmingStage` | integration | VERIFIED | |
| 5 | Changing fastp parameters invalidates trimming | checkpoint | none | **FAILS** — stays VALID (F1) | |
| 6 | FASTA/GTF validation, compatibility (names, lengths, assembly) | `reference_manager` | 9 unit | VERIFIED | |
| 6 | Reference change invalidates downstream | checkpoint outputs incl. store GTF/index | none | VERIFIED (probe: content change detected via hash/size) | |
| 6 | Reference checksums in manifest | `manifest` | none | MISSING — only provider checksum files in store (F9) | |
| 6 | Bacterial CDS annotation handled | `prepare`, `suggest_feature_type` | 4 unit | VERIFIED (E. coli UTI89) | |
| 7 | Pipefail + partial BAM never used | `runner.run_pipeline`, `alignment.align_sample` | failure test | VERIFIED | |
| 7 | primary records == input reads (×2 PE) (layer 2) | `bam_manager.validate_bam` | integration | VERIFIED (80,000 = 2×40,000) | |
| 7 | RAM pre-flight before alignment | — | none | MISSING (index build only) (F10) | |
| 7 | SIGTERM/SIGHUP stop children | `runner` | Ctrl-C test only | **FAILS** — orphaned hisat2/samtools after SIGTERM (F2) | |
| 8 | samtools stats/flagstat/idxstats + MultiQC | `BamQCStage` | integration | VERIFIED | |
| 9 | Strandedness inferred + confirmed, never guessed | `StrandednessStage` | unit + integration (reverse) | PARTIALLY VERIFIED — forward only on real yeast, unstranded never end-to-end (F13) | |
| 10 | StringTie `-e`, not used for DESeq2 | `stringtie` | integration | VERIFIED | |
| 11 | featureCounts gene set == annotation gene set | `FeatureCountsStage` | integration | VERIFIED | |
| 11 | featureCounts totals consistent with BAM | — | none | MISSING (F7) — measured: total ≥ input fragments, Assigned ≤ fragments | |
| 11 | Changing featureCounts parameters invalidates | checkpoint | none | **FAILS** (F1) | |
| 12 | Integer, no NA/negatives, unique IDs/samples | `count_matrix.validate_file` | 8 unit | VERIFIED | |
| 12 | Matrix column sum == featureCounts Assigned (layer 2) | — | none | MISSING (F7) — holds exactly in probe | |
| 12 | Technical-replicate collapse | `count_matrix.collapse` | none | UNTESTED | |
| 13 | Groups never invented; replicates; full rank; valid names; Python + R | `design`, R `--validate-only` | 7 unit + 4 R failure | VERIFIED | |
| 13 | Metadata edit after confirmation detected | design.json SHA-256 | none | VERIFIED (probe) | |
| 14 | Raw integer counts, formula, reference level, contrast | R `deseq2_pipeline.R` | R failure tests | VERIFIED | |
| 14 | Up/down tables satisfy thresholds (layer 2) | `r_bridge.validate_outputs` | integration | PARTIALLY VERIFIED — checks listed genes, not that *every* qualifying gene is listed (F8) | |
| 14 | Changing alpha/log2FC/filter invalidates DESeq2 | checkpoint | none | **FAILS** — results stale, report shows new values (F1) | |
| 15 | Report: links/images exist, no placeholders | `report_manager` | integration (sections only) | VERIFIED by probe (65 links, 0 broken); not enforced at runtime (F8) | |
| 15 | Report numbers == result files | — | none | MISSING runtime check (F8) | |
| 15 | Manifest from runtime (versions, config, design) | `manifest.write` | integration | PARTIALLY VERIFIED — no input fingerprints / reference checksums / command-log pointer (F9) | |

## 2. Side systems

| System | Requirement | BEFORE | AFTER |
|---|---|---|---|
| Checkpoints | output tamper, deletion, corruption → cascade | VERIFIED (9 probe scenarios) | |
| Checkpoints | configuration change → invalidation | **FAILS** (F1) | |
| Resume | Ctrl-C mid-command → no partial use, resumable | VERIFIED (failure test) | |
| Resume | SIGTERM / closed terminal | **FAILS** (F2) | |
| Security | no `shell=True`/`os.system`/`eval`; `yaml.safe_load` | VERIFIED (grep) | |
| Security | credentials never written/logged | **FAILS** — API key copied into `project_config.yaml`; would reach `config_used_*`, manifest, HTML report; appears in request URLs in error messages (F0) | |
| Errors | no broad exception hiding failures | PARTIALLY VERIFIED — 5 broad handlers, none hides a failure; should be narrowed (F6) | |
| Dependencies | missing tool detected before a stage starts | RISK — only at command start ("program not found") (F14) | |
| Doctor | `--check` functional health check | PARTIALLY VERIFIED — presence + version only (F15) | |
| Project health | non-interactive project validation | MISSING (menu only) (F15) | |
| Packaging | installable package, global command | MISSING — no pyproject, license, entry point (F16) | |
| Packaging | no dependence on checkout paths | **FAILS** — version from `VERSION` beside code; projects default inside the code dir (F16) | |
| Release | single version source, CHANGELOG, tag, CI | MISSING — no CI, no release process (F17) | |

## 3. Findings (BEFORE)

| ID | Severity | Location | Root cause | Impact |
|---|---|---|---|---|
| F0 | **HIGH (security)** | `project.Project.create`, `config.snapshot`, `manifest.write`, `report_manager`, `net.get_text` | merged config (incl. `ncbi.api_key`) persisted verbatim; key sent as URL parameter and URL echoed in errors | credential in project files, shareable report and logs |
| F1 | **HIGH (scientific)** | `checkpoint`, `workflow.stage_status` | checkpoints do not record the settings a stage depends on | resume silently keeps results computed with old parameters while the report prints new ones |
| F2 | **HIGH (operational)** | `runner`, `fastq_validator.validate_many`, `cli.main` | only SIGINT handled; children run in their own session | orphaned hisat2/samtools keep writing after `kill`/closed terminal |
| F3 | MEDIUM | `net.download` | no fallback when HTTP range requests are unsupported | interrupted downloads can never complete |
| F4 | MEDIUM | `net.download` | size/emptiness not checked when no checksum | empty/truncated files accepted |
| F5 | LOW (UX) | `net.download`, `runner` | curl exit codes not translated; 404 retried | cryptic network errors |
| F6 | LOW | 5 handlers | over-broad `except Exception` | could mask unexpected errors |
| F7 | MEDIUM (scientific) | `FeatureCountsStage`, `CountMatrixStage` | no independent count reconciliation | silent count loss would go unnoticed |
| F8 | MEDIUM (scientific) | `r_bridge.validate_outputs`, `report_manager` | completeness of up/down sets and report numbers not re-derived from tables | a filtering bug could drop genes unnoticed |
| F9 | MEDIUM (reproducibility) | `manifest` | inputs/reference not fingerprinted in manifest | analysis cannot be tied to exact input files |
| F10 | MEDIUM | `AlignmentStage` | no RAM check before alignment | OOM kill mid-run on small machines |
| F11 | MEDIUM (scientific) | `data_manager.register_public` | mixed organisms allowed | one reference cannot serve both organisms |
| F12 | LOW (UX) | decision screens | WHAT/WHY/OPTIONS/CONSEQUENCE not always shown | users decide without context |
| F13 | MEDIUM (test gap) | tests | unstranded/forward/single-end not exercised end-to-end | featureCounts `-s` mapping unproven for 2 of 3 cases |
| F14 | LOW | `workflow.run_stage` | no per-stage tool pre-flight | failure only after a stage starts |
| F15 | MEDIUM (UX) | `cli --check` | doctor checks presence, not function; no `--validate-project` | broken installs pass the check |
| F16 | **HIGH (distribution)** | repo layout | not a package; checkout-relative paths | cannot be installed or run from anywhere |
| F17 | MEDIUM | repo | no CI / release process / license | no evidence on clean machines; not redistributable |
