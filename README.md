# Bulk RNA-seq Pipeline (`rnaseq-pipeline`)

An interactive, terminal-based pipeline for **bulk** RNA-seq: from public accessions (SRA, GEO, ENA) or local FASTQ files
to differential-expression results, plots and an HTML report.

* **Python** does everything up to the gene count matrix: download/import, FASTQ validation, FastQC/MultiQC, the quality gate,
  fastp, reference preparation, HISAT2, BAM validation and QC, strandedness, StringTie2 and featureCounts.
* **R** runs DESeq2 and everything after it: normalisation, tests, up/down genes, and plots.
* Every stage **validates its outputs independently** (for example, FASTQ reads = HISAT2 reads = BAM primary records;
  count-matrix column sums = featureCounts assigned reads; up/down tables re-derived from the full result table).
* **Checkpoints** make interrupted runs resumable. A changed input, output or scientific setting re-runs exactly the
  affected stages.
* Scientific decisions (trimming, strandedness, experimental design) are **never guessed**; you confirm each one, with
  an explanation.

```
DATA (SRA / GEO / ENA / local, or search NCBI) -> FASTQ validation -> FastQC + MultiQC -> quality gate
  -> [fastp -> validation -> FastQC + MultiQC]  (only if you decide to trim)
  -> reference (catalog, any organism from NCBI, or your own FASTA + GTF; validated; HISAT2 index)
  -> HISAT2 | samtools sort -> BAM validation -> BAM QC -> strandedness (inferred, confirmed)
  -> StringTie2 (transcript TPM)   and   featureCounts (gene counts)      <- independent
  -> count matrix -> experimental design (you confirm) -> DESeq2 -> plots -> validated HTML report + manifest
```

## Install

```bash
pipx install git+https://github.com/vipin839/rnaseq_pipeline.git         # the command
mamba env create -f "$(rnaseq-pipeline --runtime-env tools)"               # FastQC, HISAT2, samtools, ...
mamba env create -f "$(rnaseq-pipeline --runtime-env r)"                   # R + DESeq2
rnaseq-pipeline --check                                                    # functional health check
```

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for conda setup, the Bioconda recipe (one-command install after it
is published), user settings and upgrading from 1.0.0.

## Use

```bash
rnaseq-pipeline                              # interactive menu
rnaseq-pipeline --project ~/rnaseq_projects/RNAseq_MyStudy   # resume a project
rnaseq-pipeline --validate-project DIR       # re-verify a project: PROJECT HEALTH PASS / WARNING / FAIL
rnaseq-pipeline --dry-run --project DIR      # show what would run, execute nothing
rnaseq-pipeline --check                      # installation health check
```

## What is supported

| | Status | Notes |
|---|---|---|
| Bulk RNA-seq, paired-end or single-end, Illumina | **Supported** | the layout must be the same for all samples in a project |
| Unstranded / forward / reverse libraries | **Supported** | inferred with RSeQC and confirmed; all three tested end to end |
| Human, mouse, yeast (catalog); any organism with an NCBI GTF; custom FASTA + GTF | **Supported** | genome and annotation compatibility is checked (names, lengths, assembly) |
| Bacteria (CDS-only NCBI annotation) | **Supported** | counts CDS; StringTie2 skipped; tested with E. coli UTI89 |
| Designs `~ condition`, `~ batch + condition`, several groups/contrasts | **Supported** | ≥2 biological replicates per group; confounding refused |
| Lanes / technical replicates | **Supported** | summed on request |
| Time course (LRT), interactions, continuous covariates | **Not supported** | planned |
| Mixed layouts or organisms in one project | **Refused** | use separate projects |
| Single-cell RNA-seq (10x, Drop-seq, Smart-seq, …) | **Refused** | detected from metadata and read structure; use Cell Ranger / STARsolo / alevin-fry |
| Human-scale data on a small machine | **Partially** | the index is built without splice sites below ~160 GB RAM (sites are given at alignment); ~8 GB RAM needed to align |
| GO / Reactome enrichment | **Partially** | optional; needs extra R packages; not part of the tested core |
| macOS, Windows (native) | **Not supported** | Linux and WSL2 only |
| Web interface, HPC schedulers, containers | **Not supported** | future versions |

## Documentation

| Document | Contents |
|---|---|
| [docs/INSTALLATION.md](docs/INSTALLATION.md) | install options, runtime layer, user settings, upgrading |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | step-by-step use, menus, decisions, resuming, logs |
| [docs/PIPELINE_METHODS.md](docs/PIPELINE_METHODS.md) | what each stage does, parameters, thresholds, validation |
| [docs/VERIFICATION_MATRIX.md](docs/VERIFICATION_MATRIX.md) | requirement → code → tests → evidence → status |
| [docs/TESTING.md](docs/TESTING.md) | test suites, synthetic truth data, real-data evidence |
| [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | errors and fixes |
| [docs/DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md) | architecture and extension points |
| [docs/RELEASE.md](docs/RELEASE.md) | release procedure (tag, PyPI, Bioconda) |
| [CHANGELOG.md](CHANGELOG.md) | version history |

License: MIT.
