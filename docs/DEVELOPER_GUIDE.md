# Developer Guide

## Layout

```
rnaseq_pipeline            bash launcher (sets PYTHONNOUSERSITE, execs rnaseq_pipeline.py)
rnaseq_pipeline.py         entry point (Python version / PyYAML check, then rnaseq.cli.main)
python/rnaseq/
  cli.py                   menus, wizard, argparse           (the only module that talks to the user at length)
  workflow.py              Stage classes, STAGES order, checkpoint/resume logic, stage screens
  ui.py / logger.py        status output, prompts / log files, command history
  runner.py                safe subprocess execution (argument lists, pipefail, logging, Ctrl-C handling, dry-run)
  config.py                load/merge/validate/snapshot YAML config
  project.py               project layout, project.json state, lock
  checkpoint.py            checkpoint write/verify, content fingerprints
  system_check.py          OS/CPU/RAM/disk/filesystem detection
  environment_manager.py   conda env discovery, creation, exports
  dependency_manager.py    tool/R package detection, versions, install plans
  data_manager.py / ena_manager.py / sra_manager.py / geo_manager.py / net.py   data sources
  fastq_validator.py       streaming FASTQ validation (multiprocess)
  qc_manager.py / quality_assessment.py / trimming.py
  reference_manager.py     catalog, download+checksums, FASTA/GTF validation, compatibility, HISAT2 index
  alignment.py / bam_manager.py / strandedness.py / stringtie.py / featurecounts.py / count_matrix.py
  design.py                metadata editing + design validation (rank check)
  r_bridge.py              params.json -> Rscript -> output validation
  storage.py / manifest.py / report_manager.py
R/  deseq2_pipeline.R (entry) + qc_plots.R, differential_expression.R, visualization.R, downstream_analysis.R
config/ default_config.yaml, reference_catalog.yaml
envs/   rnaseq-tools.yml, rnaseq-r.yml
tests/  unit/, failure/, integration/ (+ data/make_synthetic.py)
```

## Principles

* **Modules do not prompt the user**, except `cli.py`, `workflow.py` (stage screens and decisions) and `design.py`. The other modules are
  plain functions that take paths and parameters, so a future API or web backend (Version 3) can call them directly.
* **Every external command goes through `runner.run` / `runner.run_pipeline`.** Never use `subprocess` with `shell=True`. Empty
  arguments are rejected, so unset variables cannot silently change a command.
* **Outputs are written to `*.partial` / `*.part` / `*.building` and renamed only after validation.**
* **A stage writes its checkpoint only after post-validation.** The checkpoint fingerprint is content-based (outputs + params),
  so re-running a stage with identical results does not invalidate later stages.
* Every scientific parameter lives in `config/default_config.yaml` and is validated in `config.validate`.

## Adding a stage

1. Subclass `workflow.Stage`: set `key` (checkpoint name), `title`, `depends`, `expensive`, `settings`, and implement
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
~/miniforge3/envs/rnaseq-tools/bin/python -m pytest tests -q                 # everything (~2 min)
~/miniforge3/envs/rnaseq-tools/bin/python -m pytest tests/unit tests/failure # fast
python3 tests/data/make_synthetic.py OUTDIR 150000                            # generate the synthetic dataset
```

The integration test drives the real CLI with scripted answers. If you add or reorder a prompt, update
`new_project_answers` in `tests/integration/test_end_to_end.py`.

## Versioning

Semantic versioning in `VERSION`. Record changes in `CHANGELOG.md`. Checkpoints and the manifest store the pipeline version.
