# Bulk RNA-seq Pipeline — Version 1.0.0

An interactive, terminal-based pipeline for bulk RNA-seq, from public accessions or local FASTQ files
to differential expression results and a final HTML report.

* **Python** runs everything up to the gene count matrix: download/import, FASTQ validation, FastQC/MultiQC,
  the quality gate, fastp, reference preparation, HISAT2, BAM validation and QC, strandedness, StringTie2 and featureCounts.
* **R** runs DESeq2 and everything after it: normalisation, statistics, up/down genes, and plots.
* One launcher, `./rnaseq_pipeline`, with menus that guide you through each step.
* Every stage is validated, checkpointed and resumable, and every command is logged.

```
DATA (SRA / GEO / ENA / local) -> FASTQ validation -> FastQC + MultiQC -> quality gate
  -> [fastp -> validation -> FastQC + MultiQC]  (only if you decide to trim)
  -> reference preparation (validated genome + GTF, HISAT2 index)
  -> HISAT2 | samtools sort -> BAM validation -> BAM QC -> strandedness
  -> StringTie2 (transcript TPM)   and   featureCounts (gene counts)      <- parallel, independent
  -> count matrix -> experimental design (you confirm) -> R / DESeq2 -> plots -> HTML report
```

## Quick start

```bash
cd ~/rnaseq_pipeline
./rnaseq_pipeline            # interactive menu
./rnaseq_pipeline --check    # check the system and software, then exit
./rnaseq_pipeline --dry-run --project projects/RNAseq_MyStudy   # show what would run, execute nothing
```

Choose **1. Start New RNA-seq Project** and follow the prompts. If the program is closed or the computer
restarts, choose **2. Resume Existing Project**. Completed stages are revalidated, not recomputed.

## Documentation

| Document | Contents |
|---|---|
| [docs/INSTALLATION.md](docs/INSTALLATION.md) | Installing the software (no administrator rights needed) |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | Step-by-step use, menus, resuming, changing references/design, logs |
| [docs/PIPELINE_METHODS.md](docs/PIPELINE_METHODS.md) | What every stage does, why, parameters, thresholds and validation checks |
| [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | Common errors and how to fix them |
| [docs/DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md) | Architecture, adding stages/aligners, tests, the roadmap for Version 2 and 3 |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

## Project layout (created per analysis)

```
RNAseq_<name>/
  config/          project_config.yaml (edit this) + config_used_<time>.yaml per run
  data/            fastq/ (raw, linked or copied), trimmed/, raw/ (.sra), metadata/
  reference/       links to the shared reference store + reference_manifest.yaml
  qc/              fastqc_raw/, multiqc_raw/, fastqc_trimmed/, multiqc_trimmed/, assessment/, fastp/
  alignment/       bam/ (sorted + indexed), reports/ (HISAT2, samtools, alignment_summary.tsv)
  stringtie/       <sample>/, abundance/, merged/ (TPM matrices)
  featurecounts/   featurecounts.txt, featurecounts.summary
  counts/          gene_count_matrix.tsv/.csv, sample_metadata.tsv, design.json
  results/         deseq2/<contrast>/, upregulated/, downregulated/, plots/, tables/
  logs/            pipeline.log, pipeline_errors.log, command_history.log, commands.jsonl
  reports/         final_pipeline_report.html
  checkpoints/     one JSON per completed stage
  pipeline_manifest/  manifest.json/.yaml, environment.yml, package versions, R sessionInfo
```

## Tests

```bash
~/miniforge3/envs/rnaseq-tools/bin/python -m pytest tests -q
```

There are 111 tests: unit, failure-injection, and an end-to-end run through the real interactive CLI on a synthetic
dataset with known differential expression. The end-to-end test checks that the true DE genes are recovered
and that strandedness is inferred correctly.
