# Changelog

## 1.1.0 — 2026-09-22

Hardening and packaging release. Scientific results are unchanged: on the synthetic truth dataset the count matrix is
byte-identical to 1.0.0 and the DESeq2 statistics differ by 0 (16/16 true DE genes, 0 false positives).

### Security
* **Credentials are never stored or logged.** Old: the merged configuration, including `ncbi.api_key` and
  `ncbi.email`, was written into `config/project_config.yaml`, configuration snapshots, the manifest and the HTML report,
  and request URLs containing the key could appear in error messages. New: personal settings are stripped from project
  files, redacted in snapshots/manifest/report, and scrubbed from logs, terminal output and errors. Reason: a project
  folder or report is meant to be shared. Consequence: opening a 1.0.0 project removes the stored credentials from it.

### Behaviour changes
* **Changing a scientific setting invalidates exactly the affected stages.** Old: a checkpoint did not record the
  settings it was computed with, so after changing e.g. fastp, featureCounts or DESeq2 thresholds, resume kept the old
  results while the report printed the new values. New: each stage records its result-relevant settings and becomes
  `STALE` when one changes; later stages cascade. Consequence: after a settings change, resume re-runs those stages.
* **SIGTERM and SIGHUP (kill, closed terminal) stop child processes.** Old: only Ctrl-C was handled; HISAT2/samtools
  could keep running and writing. New: all three signals stop the running process group and the FASTQ validation pool;
  no partial output is ever accepted. Exit code 143.
* **Downloads:** restart from zero when a server cannot resume (curl 33); empty or wrong-size files are rejected even
  without a checksum; HTTP 404 and other permanent errors are explained and not retried; curl errors are translated.
  Pilot (first-N-reads) streaming from ENA retries interrupted connections and, if ENA keeps failing, falls back to the
  NCBI SRA route for that run. Both mates always come from the same archive (a single ENA mate is set aside as
  `*.other_archive`).
* **Independent count reconciliation.** featureCounts totals are checked against the reads entering alignment, and every
  count-matrix column sum must equal featureCounts `Assigned`. A mismatch stops the stage.
* **DESeq2 output completeness.** The up/down/significant tables are re-derived from the full result table; missing or
  extra genes, or wrong labels, stop the stage (previously only listed genes were checked).
* **DESeq2 results are re-derived without trusting R.** Old: R's output was checked against R's own summary, so a
  swapped baseline (every fold change reversed) or a silently dropped gene passed. New: Python recomputes the
  filtered gene set, size factors, normalized counts and baseMean from the count matrix, checks p-value ranges, and
  checks every contrast's direction against the group means. Reason: an independent layer must not depend on the
  component it checks. Consequence: inconsistent results stop the DESeq2 stage, and resume or `--validate-project`
  marks them INVALID. Recorded paths are re-pointed when a project has been moved.
* **QC reports must be complete and belong to the right files.** Old: FastQC reports were checked by read count
  only (identical counts, as in balanced designs, let a report of another file pass), and MultiQC only had to exit 0
  (it silently skips unreadable files and overwrites samples with the same name). New: each FastQC report must name
  the FASTQ it analysed, and every expected input (FastQC reports; HISAT2 summaries and samtools outputs) must appear
  in MultiQC's `multiqc_sources.txt`.
* **Plot data is reconciled with the results.** Old: plots were only checked to exist. New: the data table behind
  every plot must match the result files (MA/volcano genes, values and up/down labels; heatmap gene lists and
  z-scoring; library sizes = count-matrix column sums; well-formed PCA/distance/correlation tables).
* **More report numbers are re-checked, and the claim is exact.** Old: only samples, reference, strandedness, formula,
  thresholds and per-contrast DE counts were re-checked, while the terminal said "every number matches". New:
  per-sample input reads, alignment rate, primary mapped reads and %, assigned and counted reads, and software
  versions are re-derived from the files the tools wrote; messages list what is checked. The report check also no
  longer depends on the project's original location.
* **The report is validated after it is written:** every link and image must exist, the document must be complete, and
  every number marked in it (samples, reference, strandedness, formula, thresholds, up/down/significant/tested per
  contrast) must equal the value re-derived from the result files.
* **Manifest** now records input-file fingerprints, provider MD5s, reference genome/annotation SHA-256, reference
  download details and the command log. A 1.0.0 project's report stage is re-run once on resume to add them.
* **Pre-flight checks:** required tools per stage before it starts; RAM before alignment (with an explanation and a
  choice); more than one organism in a project is refused; Rscript is required before the design stage.
* Decision screens (quality gate, trimming, strandedness, design) show WHAT / WHY / OPTIONS / CONSEQUENCE.
* **Per-area project health.** `--validate-project` (and main menu 3) now groups every check into areas — system,
  installation, configuration, data, reference, QC, alignment, quantification, design, DESeq2 + results, report,
  security, resume — prints an area summary plus a scientific-validation area that fails if any scientific area
  fails, and lists BLOCKING (FAIL) and NON-BLOCKING (WARNING) issues. Missing tools are a non-blocking warning (the
  results stay valid; the tools are needed only to continue).
* **Aligner status is shown and enforced.** A new project shows every aligner with its real status (HISAT2:
  supported and validated; STAR: not available — planned, not implemented or validated) and whether it is installed.
  The new setting `aligner` accepts only `hisat2` (anything else is refused with the reason) and is recorded with the
  alignment checkpoint, so the aligner can never change silently. A setting that did not exist when a checkpoint
  was written counts as its default (the behaviour of that run), so upgrading does not force recomputation.
* **Quality gate "Review the evidence per sample" option.** Shows each sample's metrics next to the thresholds,
  explains each finding as a technical problem (trimming helps) or a normal RNA-seq pattern (duplication,
  per-base content bias), lists the report files, and returns to the same decision. "Stop pipeline" is now option 5.
* **Every failure screen also states Impact and Retry** (what was kept or discarded; whether and when a retry is
  safe and what it repeats), and names the step as shown on screen ("STEP 6: REFERENCE PREPARATION") instead of an
  internal key. Both fields are written to `logs/pipeline_errors.log`.

### New
* Installable Python package (`src/rnaseq_pipeline`, `pyproject.toml`): commands `rnaseq-pipeline` and
  `rnaseq_pipeline`, and `python -m rnaseq_pipeline`. The version has a single source (`__version__`).
* User configuration `~/.config/rnaseq-pipeline/config.yaml`, loaded automatically; projects default to
  `~/rnaseq_projects`; state in `~/.rnaseq_pipeline`. The installed package directory is never written to.
* `--check` is a functional health check: it runs a real mini-job through HISAT2, samtools, featureCounts, StringTie,
  fastp, FastQC and DESeq2, and reports `HEALTH: PASS / WARNING / FAIL` (exit 1 only on FAIL). `--quick` skips the
  functional and network tests.
* `--validate-project DIR` re-verifies a project non-interactively: `PROJECT HEALTH: PASS / WARNING / FAIL`.
* `--runtime-env tools|r` prints the path of the shipped conda environment file.
* MIT license; GitHub Actions CI (lint, unit/security/failure tests on Ubuntu 22.04/24.04 × Python 3.10–3.13, package
  build and clean-install tests, integration tests with the real tools); Bioconda recipe prepared (not yet submitted).

### Fixed
* **No orphaned processes from quick commands.** Version probes and checks such as `samtools flagstat` ran without
  the process-group handling of pipeline commands; on Ctrl-C/SIGTERM Python stopped only the direct child, so the
  real program behind a wrapper script (hisat2, hisat2-build, fastqc) kept running. They now run in their own
  process group, which is stopped on interrupt or timeout.
* **Count matrix could not be rebuilt in a moved or copied project.** featureCounts records absolute BAM paths and
  they were compared with the project's current location ("sample columns do not match"). Columns are now matched by
  file name, in order.
* **Outdated or unusable tools are refused before a step starts** (e.g. featureCounts 1.5.0 when 2.0.0 is the
  minimum); previously only their presence was checked and they ran.
* **Clear explanations instead of tracebacks or wrong advice:** a read-only project (opening it, or writing during a
  step) names the folder and the fix and points to the read-only `--validate-project`; an R process killed by a
  signal says so (e.g. SIGKILL, often out of memory); a missing R package points to repairing the R environment
  instead of "fix the metadata/design".
* **Strandedness in the project state is checked against the confirmed decision.** Old: if the value later stages
  read (project state) differed from the checkpointed decision, featureCounts results computed with the old value
  stayed VALID, and a later re-run would silently use the new one. New: the strandedness stage is INVALID until the
  decision is made again, and everything downstream follows.
* **A correct BAM was rejected when the system clock stepped backwards.** Old: the BAM index was required to be
  newer than the BAM by file time. The WSL2 clock stepped back ~1 s between writing a BAM and its index, and a
  correct BAM was rejected ("index is older than the BAM"); conversely, a stale index with a newer time passed. New:
  the index must describe the BAM's content — the records counted by the index (`samtools idxstats`) must equal the
  records in the BAM (`samtools flagstat`, QC-passed + QC-failed).
* **Report failed when the tools were not run from conda.** Old: the report always listed `environment.yml` and
  `package_versions.txt`; without a conda environment (tools on PATH, modules, containers) they cannot exist, the
  report showed "(missing)", and its validation stopped the very last stage of a finished analysis. New: the manifest
  records whether a conda export was possible (`conda_export`), and the report states "not produced" with the
  reason and links `logs/software_versions.tsv` (tool versions are always recorded). A conda export that should
  exist and is missing still fails validation.
* **Intermittent HISAT2 index-build crash.** `hisat2-build` 2.2.3 crashes with SIGSEGV in about 0.6% of
  multithreaded builds (3 of 480 here; 0 of 395 single-threaded builds), which stopped the reference stage (and made
  one integration test fail intermittently). Now a crash by SIGSEGV/SIGBUS/SIGABRT discards the partial index and
  retries once single-threaded; the index is byte-identical to a multithreaded build and is validated as always.
  Out-of-memory kills and SIGILL are not retried.
* `tests/data/make_synthetic.py`: a variable collision made the `strand` argument change gene strands; forward,
  unstranded and single-end datasets are now generated correctly (the default dataset is byte-identical).
* Five over-broad `except Exception` handlers narrowed to the errors they are meant to handle.
* `--check` reported NCBI as unreachable when it answered with an HTTP error status.
* R was not found in micromamba installations (`$MAMBA_ROOT_PREFIX/envs`), and when several R installations are on
  PATH the first one was used even if it lacks DESeq2 (the tools environment contains a bare R, required by RSeQC).
  Now environments under `$MAMBA_ROOT_PREFIX` are found and the R that has DESeq2 is preferred.
* **Older CPUs:** every Bioconda build of StringTie 3.x is compiled for x86-64-v3 (AVX2/BMI2) and is killed with
  SIGILL on older processors; `--check` said "reinstall", which reinstalls the same build. Now `--check` shows the
  CPU instruction level, recognises a SIGILL crash as a CPU incompatibility, and names the fix
  (`stringtie=2.2.3`, whose build uses no AVX2/BMI2 instructions); main menu 4 proposes that version, and a
  pipeline step killed by SIGILL says so. StringTie only produces the transcript TPM tables; the DESeq2 input
  (featureCounts) is unaffected. The installed StringTie version is recorded in each project, as before.
* **Creating the tools environment on an older CPU** installed StringTie 3.x, which then crashed. Now
  `--runtime-env tools` and main menu 4 give such CPUs an adapted copy of the environment file (in
  `~/.rnaseq_pipeline/envs/`) that pins `stringtie=2.2.3`; the shipped file is unchanged and still used on
  modern CPUs.

### Tests
* 262 tests (127 at the start of this work, commit 08a9c11): security, download failures, signals, reconciliation invariants, report tampering, a
  13-scenario resume-invalidation matrix, unstranded / forward / single-end libraries and a deliberately wrong
  strandedness setting end to end, project-health and doctor checks, release consistency, and clean installs
  (venv and pipx).

## 1.0.0 — 2026-09-19

First release.

* Interactive terminal application (`./rnaseq_pipeline`) with main menu, new-project wizard, resume, validation,
  dependency and reference management, configuration editing, `--dry-run`, `--auto`, `--config`, `--check`.
* Data sources: NCBI SRA (prefetch / vdb-validate / fasterq-dump), NCBI GEO (GSE → SRX → runs), ENA (MD5-verified
  FASTQ), local FASTQ (symlink or copy, automatic pairing), pilot subset mode for testing.
* Streaming FASTQ validation (structure, alphabet, quality, lengths, gzip integrity, mate IDs/counts).
* FastQC + MultiQC, a configurable quality gate, conditional fastp with full post-trimming re-validation and QC.
* Reference catalog (GENCODE, Ensembl, NCBI RefSeq, UCSC for human/mouse; Ensembl yeast), custom and existing
  references, provider checksum verification, genome/annotation compatibility checks, RAM-aware HISAT2 index build.
* HISAT2 → samtools sort streaming alignment with pipefail semantics, BAM integrity validation (primary records = input
  reads), BAM QC, RSeQC-based strandedness with user confirmation.
* StringTie2 reference-guided transcript quantification (independent) and featureCounts gene counts (DESeq2 input).
* Validated count matrix, technical-replicate summing, interactive experimental design with rank/replicate checks.
* R/DESeq2: low-count filter, Wald test, BH adjustment, explicit up/down thresholds, apeglm shrinkage, PCA,
  distance/correlation heatmaps, MA, volcano, significant and top-gene heatmaps (PNG + PDF + source data),
  optional GO/Reactome enrichment.
* Checkpoints with content fingerprints and cascade invalidation; resume without recomputation.
* Logs (pipeline, errors, command history, JSONL), software versions, environment exports, R sessionInfo,
  reproducibility manifest, final HTML report.
* 111 automated tests, including an end-to-end CLI run on synthetic data with known truth.
