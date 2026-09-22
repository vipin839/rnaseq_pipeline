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
| 1 | Download: ENA MD5 verified, `.part` → atomic rename | `net.download`, `ena_manager.download_run` | none direct | PARTIALLY VERIFIED — MD5 mismatch kept as `.badmd5` (probe) | VERIFIED — `failure/test_downloads::test_verified_download`, `test_checksum_mismatch_rejected_and_kept` |
| 1 | Download: resume after interruption | `net.download` (`curl -C -`) | none | **FAILS** when server rejects ranges (curl 33) → aborts instead of restarting (F3) | VERIFIED (F3 fixed) — `test_resume_from_partial`, `test_restart_when_server_cannot_resume` |
| 1 | Download: empty/short file never accepted | `net.download` | none | **FAILS** — 0-byte file accepted when no checksum (F4) | VERIFIED (F4 fixed) — `test_empty_download_rejected`, `test_size_mismatch_rejected`, `test_existing_invalid_file_is_never_overwritten` |
| 1 | Download: actionable network errors | `runner`/`net` | none | RISK — "exit codes [22]: curl"; 404 retried (F5) | VERIFIED (F5 fixed) — `test_http_404_explained_and_not_retried`, `test_unreachable_server_explained` |
| 1 | SRA: prefetch → vdb-validate → fasterq-dump | `sra_manager.download_run` | none (manual evidence) | VERIFIED (manual, SRR7361181; corrupt .sra → exit 3) | VERIFIED (manual, SRR7361181) + used as ENA fallback on real data (SRR14208246, 1.1.0 run) |
| 1 | Pilot streaming survives transient ENA failures; mates from one archive | `ena_manager._download_head`, `DataStage._set_aside_partial_pair` | none | **FAILS** — single attempt, aborts (F18, found during 1.1.0 real-data run) | VERIFIED (F18 fixed) — `test_pilot_retries_transient_server_errors`, `test_pilot_truncated_stream_fails_with_reason`, `test_ena_failure_falls_back_to_sra_without_mixing_mates`; observed on SRR14208246 (real run: ENA → SRA fallback, all 15 stages passed) |
| 1 | Single-cell data refused | `data_manager.single_cell_reason`, SRA >2 reads guard | unit (5 cases) | VERIFIED (SRR23333328 refused) | VERIFIED — `unit::test_single_cell_detection`; SRR23333328 refused |
| 1 | One organism per project | `data_manager.show_public` | none | **RISK** — mixed organisms only warned, can be registered (F11) | VERIFIED (F11 fixed) — `unit::test_mixed_organisms_refused` |
| 2 | FASTQ structure/alphabet/lengths/gzip/mates | `fastq_validator` | 14 unit | VERIFIED | VERIFIED — unit FASTQ-validator tests |
| 2 | Validation interruptible (Ctrl-C/SIGTERM) | `fastq_validator.validate_many` | none | RISK — ProcessPool waits for workers (F2) | VERIFIED (F2 fixed) — `failure::test_fastq_validation_pool_stops_on_interrupt` (fails without the fix) |
| 3 | FastQC count == FASTQ count (layer 2) | `qc_manager.validate_fastqc` | integration | VERIFIED | VERIFIED — integration `test_trimming_and_qc` |
| 4 | Advisory gate, user decides, thresholds in config | `quality_assessment`, `QualityGateStage` | integration | VERIFIED; explanation of decision minimal (F12) | VERIFIED (F12 fixed: WHAT/WHY/OPTIONS/CONSEQUENCE) — integration |
| 5 | fastp output re-validated; counts == fastp JSON | `TrimmingStage` | integration | VERIFIED | VERIFIED — integration |
| 5 | Changing fastp parameters invalidates trimming | checkpoint | none | **FAILS** — stays VALID (F1) | VERIFIED (F1 fixed) — `test_resume_invalidation_matrix[fastp…]` |
| 6 | FASTA/GTF validation, compatibility (names, lengths, assembly) | `reference_manager` | 9 unit | VERIFIED | VERIFIED — 9 unit tests |
| 6 | Reference change invalidates downstream | checkpoint outputs incl. store GTF/index | none | VERIFIED (probe: content change detected via hash/size) | VERIFIED by probe (code unchanged); no automated test — reference file changes are not in the invalidation matrix |
| 6 | Reference checksums in manifest | `manifest` | none | MISSING — only provider checksum files in store (F9) | VERIFIED (F9 fixed) — `test_report_manifest_logs` asserts genome/annotation SHA-256 |
| 6 | Bacterial CDS annotation handled | `prepare`, `suggest_feature_type` | 4 unit | VERIFIED (E. coli UTI89) | VERIFIED — unit `test_bacterial_gtf_counts_cds`, `test_wrong_feature_type_is_rejected` |
| 7 | Pipefail + partial BAM never used | `runner.run_pipeline`, `alignment.align_sample` | failure test | VERIFIED | VERIFIED — `failure::test_pipefail_semantics`, `test_corrupted_bam_detected` |
| 7 | primary records == input reads (×2 PE) (layer 2) | `bam_manager.validate_bam` | integration | VERIFIED (80,000 = 2×40,000) | VERIFIED — `test_alignment_integrity`; real yeast pilot (see TESTING.md) |
| 7 | RAM pre-flight before alignment | — | none | MISSING (index build only) (F10) | VERIFIED (F10 fixed) — `unit::test_alignment_ram_estimate_scales_with_index` |
| 7 | SIGTERM/SIGHUP stop children | `runner` | Ctrl-C test only | **FAILS** — orphaned hisat2/samtools after SIGTERM (F2) | VERIFIED (F2 fixed) — `failure::test_termination_signal_stops_child_processes[SIGTERM/SIGHUP]` |
| 8 | samtools stats/flagstat/idxstats + MultiQC | `BamQCStage` | integration | VERIFIED | VERIFIED — integration |
| 9 | Strandedness inferred + confirmed, never guessed | `StrandednessStage` | unit + integration (reverse) | PARTIALLY VERIFIED — forward only on real yeast, unstranded never end-to-end (F13) | VERIFIED (F13 fixed) — `integration/test_library_types` (unstranded, forward, reverse single-end) + `test_wrong_strandedness_is_caught` |
| 10 | StringTie `-e`, not used for DESeq2 | `stringtie` | integration | VERIFIED | VERIFIED — integration |
| 11 | featureCounts gene set == annotation gene set | `FeatureCountsStage` | integration | VERIFIED | VERIFIED — integration |
| 11 | featureCounts totals consistent with BAM | — | none | MISSING (F7) — measured: total ≥ input fragments, Assigned ≤ fragments | VERIFIED (F7 fixed) — `unit::test_featurecounts_reconcile_invariants`; enforced at runtime |
| 11 | Changing featureCounts parameters invalidates | checkpoint | none | **FAILS** (F1) | VERIFIED (F1 fixed) — invalidation matrix |
| 12 | Integer, no NA/negatives, unique IDs/samples | `count_matrix.validate_file` | 8 unit | VERIFIED | VERIFIED — unit |
| 12 | Matrix column sum == featureCounts Assigned (layer 2) | — | none | MISSING (F7) — holds exactly in probe | VERIFIED (F7 fixed) — `unit::test_count_matrix_sums_checked_against_assigned`; enforced before and after writing |
| 12 | Technical-replicate collapse | `count_matrix.collapse` | none | UNTESTED | VERIFIED — `unit::test_technical_replicates_are_summed` |
| 13 | Groups never invented; replicates; full rank; valid names; Python + R | `design`, R `--validate-only` | 7 unit + 4 R failure | VERIFIED | VERIFIED — unit + R failure tests |
| 13 | Metadata edit after confirmation detected | design.json SHA-256 | none | VERIFIED (probe) | VERIFIED — invalidation matrix (edited metadata) |
| 14 | Raw integer counts, formula, reference level, contrast | R `deseq2_pipeline.R` | R failure tests | VERIFIED | VERIFIED — R failure tests; truth recovery 16/16, 0 FP |
| 14 | Up/down tables satisfy thresholds (layer 2) | `r_bridge.validate_outputs` | integration | PARTIALLY VERIFIED — checks listed genes, not that *every* qualifying gene is listed (F8) | VERIFIED (F8 fixed) — `unit::test_threshold_sets_*` (3 tests); enforced at runtime |
| 14 | Changing alpha/log2FC/filter invalidates DESeq2 | checkpoint | none | **FAILS** — results stale, report shows new values (F1) | VERIFIED (F1 fixed) — invalidation matrix |
| 15 | Report: links/images exist, no placeholders | `report_manager` | integration (sections only) | VERIFIED by probe (65 links, 0 broken); not enforced at runtime (F8) | VERIFIED (F8 fixed) — enforced at runtime; `test_report_tampering_detected[broken link, truncated]` |
| 15 | Report numbers == result files | — | none | MISSING runtime check (F8) | VERIFIED (F8 fixed) — `test_report_validates_against_result_files`, `test_report_tampering_detected[number]` |
| 15 | Manifest from runtime (versions, config, design) | `manifest.write` | integration | PARTIALLY VERIFIED — no input fingerprints / reference checksums / command-log pointer (F9) | VERIFIED (F9 fixed) — `test_report_manifest_logs`; `manifest.problems` used by `--validate-project` |

## 2. Side systems

| System | Requirement | BEFORE | AFTER |
|---|---|---|---|
| Checkpoints | output tamper, deletion, corruption → cascade | VERIFIED (9 probe scenarios) | VERIFIED — invalidation matrix (corrupted/deleted outputs) |
| Checkpoints | configuration change → invalidation | **FAILS** (F1) | VERIFIED (F1 fixed) — invalidation matrix (13 scenarios) |
| Resume | Ctrl-C mid-command → no partial use, resumable | VERIFIED (failure test) | VERIFIED — `failure::test_interrupted_command_is_terminated` |
| Resume | SIGTERM / closed terminal | **FAILS** (F2) | VERIFIED (F2 fixed) — signal tests |
| Security | no `shell=True`/`os.system`/`eval`; `yaml.safe_load` | VERIFIED (grep) | VERIFIED — grep + `failure::test_no_shell_expansion` |
| Security | credentials never written/logged | **FAILS** — API key copied into `project_config.yaml`; would reach `config_used_*`, manifest, HTML report; appears in request URLs in error messages (F0) | VERIFIED (F0 fixed) — `security/` (8 tests), `test_no_credentials_anywhere_in_project` |
| Errors | no broad exception hiding failures | PARTIALLY VERIFIED — 5 broad handlers, none hides a failure; should be narrowed (F6) | VERIFIED (F6 fixed) — handlers narrowed |
| Dependencies | missing tool detected before a stage starts | RISK — only at command start ("program not found") (F14) | VERIFIED (F14 fixed) — `unit::test_missing_tool_reported_before_stage_starts` |
| Doctor | `--check` functional health check | PARTIALLY VERIFIED — presence + version only (F15) | VERIFIED (F15 fixed) — `integration::test_health_check_runs_real_mini_job`; `packaging::test_check_without_tools_fails_cleanly` |
| Project health | non-interactive project validation | MISSING (menu only) (F15) | VERIFIED (F15 fixed) — `test_validate_project_passes_on_finished_project`, `test_validate_project_fails_on_edited_count_matrix` |
| Packaging | installable package, global command | MISSING — no pyproject, license, entry point (F16) | VERIFIED (F16 fixed) — `packaging/` (venv + pipx, 8 tests) |
| Packaging | no dependence on checkout paths | **FAILS** — version from `VERSION` beside code; projects default inside the code dir (F16) | VERIFIED (F16 fixed) — `packaging::test_projects_and_state_live_in_home_not_package`, `test_installed_command_from_any_directory` |
| Packaging | Bioconda recipe | MISSING | VERIFIED locally — `bioconda-utils lint`: All checks OK; `conda-build` succeeded and the recipe tests passed (`--check`: 38 passed, 1 warning [NCBI email], 0 failed). Not submitted to Bioconda (user decision) |
| Release | single version source, CHANGELOG, tag, CI | MISSING — no CI, no release process (F17) | VERIFIED (F17 fixed) — `unit/test_release.py` (6 tests); CI workflow |

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
| F18 | MEDIUM (operational) | `ena_manager._download_head`, `DataStage` | pilot streaming (curl piped to head) made one attempt; no fallback when ENA fails | real-data pilot runs aborted on ENA connection drops (curl 56, observed on SRP314352) |
| F20 | MEDIUM (installation) | `dependency_manager.detect_tools`, `doctor`, `runner` | SIGILL crash reported as "version unknown / reinstall"; CPU level not known | on CPUs without AVX2/BMI2 (user's lab machine) StringTie 3.x cannot run and the advice reinstalled the same build — fixed: CPU level shown, SIGILL explained, `stringtie=2.2.3` proposed (build disassembled: 0 AVX2/BMI2 instructions); verified: fix command applied to a clone of the tools env, `--check` 0 failed, 35/35 integration tests pass with StringTie 2.2.3; 4 unit tests (fail without the fix) |
| F19 | MEDIUM (installation) | `environment_manager.Environments.prefix/rscript` | environments under `$MAMBA_ROOT_PREFIX` not found; PATH fallback took the first `Rscript`, which can be the bare R that RSeQC pulls into the tools env | with micromamba or the tools env on PATH, DESeq2 reported missing (found by the first CI integration run) — fixed; `unit::test_rscript_on_path_prefers_the_r_with_deseq2`, `test_envs_found_under_mamba_root_prefix` (both fail without the fix) |
