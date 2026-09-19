# User Guide

## Starting

```bash
cd ~/rnaseq_pipeline
./rnaseq_pipeline
```

| Option | Meaning |
|---|---|
| `--project DIR` | open a project directly (goes to the project menu) |
| `--dry-run` | check tools, paths, reference, inputs and disk; print the commands; execute nothing |
| `--auto` | skip the menu for quick stages; decisions and expensive steps still ask |
| `--config FILE` | YAML file whose values override the defaults for new projects (for example a lab profile) |
| `--projects-dir DIR` | where projects are listed and created (default `~/rnaseq_pipeline/projects`) |
| `--verbose` | also print every command and debug message on screen |
| `--check` | run the system and dependency checks, then exit |

Status tags: `[INFO]` information, `[OK]` validated success, `[RUNNING]` in progress, `[WARNING]` needs your attention
but not fatal, `[ERROR]`/`[FAILED]` stopped, `[SKIPPED]` not needed (already valid, or not requested), `[ACTION]` a recommendation.

## A new project, step by step

1. **System check.** OS, CPU, RAM, disk, filesystem, Python, R, conda. Critical problems (for example a FAT32 disk
   or a read-only location) stop here.
2. **Project name.** Letters, digits, `.`, `_` and `-` only. A folder `RNAseq_<name>` is created.
3. **Data input**
   * *NCBI SRA*: run, experiment, study or BioProject accessions. Downloaded with `prefetch`, checked with
     `vdb-validate`, and converted with `fasterq-dump`.
   * *NCBI GEO*: a series (GSE). Samples are linked to SRA runs, and GEO sample characteristics are shown to help you.
   * *ENA*: any ENA/SRA accession. FASTQ files are downloaded directly and their MD5 checksums verified.
   * *Local FASTQ*: files or directories. Pairs are detected from `_R1/_R2`, `_1/_2` or Illumina `_R1_001` names.
     Files are **symlinked** by default, so nothing is duplicated and the originals are never modified.
   * The metadata table (organism, layout, instrument, reads, size, description) is shown for you to review.
     You then choose **Full data**, or a **Pilot subset** (the first N reads of each run) for testing only.
     Pilot runs carry a warning in the report and the manifest.
4. **Reference.** Choose the organism and source (GENCODE, Ensembl, RefSeq, UCSC, a custom local FASTA + GTF,
   or an existing prepared reference or index). The download URLs and sizes are shown. Nothing is downloaded until
   the *Reference preparation* step, where you confirm again.
5. **Optional databases.** GO/Reactome enrichment, off by default.
6. **Samples.** Use all samples, or a subset (useful for debugging).

The pipeline then walks through the stages. Each stage shows a screen like this:

```
STEP 7 — ALIGNMENT (HISAT2)
Samples: 6   Validated: 6/6   Reference: GRCh38 / GENCODE 46   Threads: 11   Estimated storage: 40.3 GB
1. Start alignment  2. Review settings  3. Change settings  4. Validate inputs again  5. Return to main menu  6. Cancel
```

### Decisions you will be asked to make

* **Quality gate.** The assessment is PASS / REVIEW / TRIMMING RECOMMENDED / FAIL per sample, with reasons. You choose:
  accept the recommendation, run fastp, skip trimming, or stop.
* **Strandedness.** RSeQC evidence is shown per sample. You confirm the inferred value or choose one. If it cannot be
  determined you **must** choose. It is never guessed silently.
* **Technical replicates.** If several runs belong to one biological sample (lanes), you are asked whether to sum them.
* **Experimental design.** Nothing is assumed. You assign a condition to every sample, either one by one, by text pattern,
  from a TSV file, or by explicitly picking a public metadata field. You can add covariates (batch, sex, ...),
  set the formula (`~ condition` or `~ batch + condition`), the reference (baseline) level and the contrasts. The
  design is validated (replicates, confounding, valid names) and must be **confirmed** before DESeq2.
* **Proceed to DESeq2?** A transition record is shown first.

## Resuming

**Main menu → 2. Resume Existing Project** (or `./rnaseq_pipeline --project DIR`).
The status table shows each stage as VALID, INVALID (with the reason) or PENDING. "Continue pipeline" skips VALID
stages after re-checking their outputs, and re-runs anything invalid. If a stage was interrupted (Ctrl-C, crash,
reboot), its partial outputs (`*.partial.*`) are never used. The stage simply runs again. Samples that were already
aligned and still pass validation are not re-aligned.

## Validating a project

**Main menu → 3.** Re-hashes every recorded output and re-runs the integrity checks (for example BAM quickcheck).
It can optionally re-read every FASTQ file.

## Changing things later

* **Reference:** project menu → *Change reference*. Alignment and every later stage will re-run.
* **Design / contrasts:** project menu → *Experimental design*. Only DESeq2 and the report re-run.
* **Parameters:** main menu → 7, or *Change settings* on a stage screen, or edit `config/project_config.yaml`.
  Every value is validated. Re-run the affected stage (project menu → *Re-run a specific stage*); later stages follow.
* **Samples:** project menu → *Samples*. Re-include failed samples after fixing them, or exclude samples.

## Where to look

| What | File |
|---|---|
| Everything, with timestamps | `logs/pipeline.log` |
| Warnings and errors only | `logs/pipeline_errors.log` |
| Every command, exit code, duration | `logs/command_history.log` (readable), `logs/commands.jsonl` (structured) |
| Raw output of each tool | `logs/<stage>/<sample>.log` |
| Final report | `reports/final_pipeline_report.html` |
| DE results | `results/deseq2/<contrast>/` (full, significant, up, down, top tables, statistics.txt) |
| Plots (PNG + PDF) | `results/plots/` and `results/plots/<contrast>/`. The data behind each plot is in `results/tables/` and `results/deseq2/<contrast>/plot_data/` |

## Cleaning up

Project menu → *Clean up intermediate files* offers only files whose downstream outputs are validated: trimmed FASTQ
after the BAMs are validated, `.sra` caches after the FASTQs are validated, and temporary files. Raw data, references and BAMs are
never offered.
