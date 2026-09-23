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
* **The report is validated after it is written:** every link and image must exist, the document must be complete, and
  every number marked in it (samples, reference, strandedness, formula, thresholds, up/down/significant/tested per
  contrast) must equal the value re-derived from the result files.
* **Manifest** now records input-file fingerprints, provider MD5s, reference genome/annotation SHA-256, reference
  download details and the command log. A 1.0.0 project's report stage is re-run once on resume to add them.
* **Pre-flight checks:** required tools per stage before it starts; RAM before alignment (with an explanation and a
  choice); more than one organism in a project is refused; Rscript is required before the design stage.
* Decision screens (quality gate, trimming, strandedness, design) show WHAT / WHY / OPTIONS / CONSEQUENCE.

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
* 211 tests (127 at the start of this work, commit 08a9c11): security, download failures, signals, reconciliation invariants, report tampering, a
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
