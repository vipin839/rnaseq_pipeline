# Changelog

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
