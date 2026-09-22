# Developer Guide

## Layout

```
pyproject.toml             package metadata; console scripts rnaseq-pipeline / rnaseq_pipeline
rnaseq_pipeline, rnaseq_pipeline.py   run from a checkout without installing
src/rnaseq_pipeline/
  __init__.py              __version__ (single source), PACKAGE_DIR, PipelineError, Terminated, signal handlers
  __main__.py              python -m rnaseq_pipeline
  cli.py                   menus, wizard, argparse (--check, --validate-project, --runtime-env, ...)
  workflow.py              Stage classes, STAGES order, checkpoints/resume, settings invalidation, stage screens
  doctor.py                --check: system/config/tools + functional mini-jobs + R/DESeq2 + network
  project_health.py        --validate-project and main menu 3
  secrets.py               credential redaction (projects, snapshots, manifest, report, logs, terminal)
  ui.py / logger.py        status output, prompts, decision explanations / log files, command history
  runner.py                safe subprocess execution (argument lists, pipefail, logging, signals, dry-run)
  config.py                load/merge (default < user config < --config), validate, snapshot
  project.py               project layout, project.json state, lock
  checkpoint.py            checkpoint write/verify, content fingerprints, manifest file fingerprints
  system_check.py          OS/CPU/RAM/disk/filesystem detection
  environment_manager.py   conda env discovery, creation, exports; ENV_FILES (package data)
  dependency_manager.py    tool/R package detection, versions, install plans
  data_manager.py / ena_manager.py / sra_manager.py / geo_manager.py / entrez.py / net.py   data sources
  fastq_validator.py       streaming FASTQ validation (terminable multiprocessing pool)
  qc_manager.py / quality_assessment.py / trimming.py
  reference_manager.py     catalog, NCBI search, downloads + checksums, FASTA/GTF validation, compatibility,
                           prepare() (stage + standalone), HISAT2 index, scan_store()
  alignment.py / bam_manager.py / strandedness.py / stringtie.py / featurecounts.py / count_matrix.py
  design.py                metadata editing + design validation (rank check)
  r_bridge.py              params.json -> Rscript -> independent output validation
  storage.py / manifest.py / report_manager.py (report + its validator)
  R/                       deseq2_pipeline.R (entry) + qc_plots.R, differential_expression.R, visualization.R,
                           downstream_analysis.R
  config/                  default_config.yaml, reference_catalog.yaml
  envs/                    rnaseq-tools.yml, rnaseq-r.yml (runtime layer)
packaging/bioconda/        Bioconda recipe + sha256 helper
.github/workflows/ci.yml   lint, unit (Ubuntu 22.04/24.04 x Python 3.10-3.13), packaging, integration
tests/  unit/ security/ failure/ integration/ packaging/ data/make_synthetic.py
```

## Principles

* **Modules do not prompt the user**, except `cli.py`, `workflow.py` (stage screens and decisions) and `design.py`. The other modules are
  plain functions that take paths and parameters, so a future API or web backend (Version 3) can call them directly.
* **Every external command goes through `runner.run` / `runner.run_pipeline`.** Never use `subprocess` with `shell=True`. Empty
  arguments are rejected, so unset variables cannot silently change a command.
* **Outputs are written to `*.partial` / `*.part` / `*.building` and renamed only after validation.**
* **A stage writes its checkpoint only after post-validation.** The checkpoint fingerprint is content-based (outputs + params),
  so re-running a stage with identical results does not invalidate later stages.
* Every scientific parameter lives in `src/rnaseq_pipeline/config/default_config.yaml` and is validated in `config.validate`.

## Adding a stage

1. Subclass `workflow.Stage`: set `key` (checkpoint name), `title`, `depends`, `expensive`, `settings`,
   `result_settings` (settings whose change must invalidate the stage), `required_tools(ctx)`, and implement
   `overview`, `check_inputs`, `estimate_gb`, `execute(ctx) -> (outputs, params, summary)` and optionally `revalidate`.
   `execute` must handle `ctx.dry_run` (commands only).
2. Insert an instance into `STAGES` at the right position.
3. Add a report section in `report_manager.generate` and a test.

## Extension points for Version 2 (not implemented)

* **Other aligners (STAR) / pseudo-aligners (Salmon, kallisto):** add a module next to `alignment.py` and an `aligner`
  config key. `AlignmentStage` calls it and the BAM validation is shared. Salmon/kallisto would add a separate
  quantification stage producing a count matrix (tximport in R).
* **Execution backends (SLURM, HPC, containers):** `runner.run_pipeline` is the single execution point. Add a backend that wraps
  the argument list (for example `sbatch --wrap` or `apptainer exec`) behind a config switch.
* **Parallel sample scheduling:** per-sample loops in `AlignmentStage` / `TrimmingStage` can be moved to a worker pool.
  Outputs are already per sample and atomic.
* **Configuration profiles:** `--config FILE` already merges a YAML profile over the defaults.
* **Richer designs** (interactions, continuous covariates): extend `validators.design_formula`, `design.model_rank` and the
  R validation.

## Version 3 (web) notes

Nothing in the core modules depends on the terminal. A web layer would create a `Project`, fill the samples, reference and
design through the same functions, and run `workflow.run_stage` for each stage in a job queue, reading status via
`workflow.status_table`. The prompts in `QualityGateStage`, `StrandednessStage` and the design step would need non-interactive
decision inputs, which is the main refactor to plan for.

## Tests

```bash
pytest tests                                         # everything (~6 min; needs the tools)
pytest tests/unit tests/security tests/failure       # fast; see docs/TESTING.md
python3 tests/data/make_synthetic.py OUTDIR 150000                            # generate the synthetic dataset
```

The integration test drives the real CLI with scripted answers. If you add or reorder a prompt, update
`new_project_answers` in `tests/integration/test_end_to_end.py`.

## Versioning

Semantic versioning; the single source is `__version__` in `src/rnaseq_pipeline/__init__.py` (see docs/RELEASE.md). Record changes in `CHANGELOG.md`. Checkpoints and the manifest store the pipeline version.
