# Testing

## Running the tests

```bash
pip install -e ".[dev]"
pytest tests/unit tests/security tests/failure        # fast (~15 s); tests that need tools skip themselves
pytest tests                                          # everything (~6 min) — needs the scientific tools
RNASEQ_PACKAGING_TESTS=1 pytest tests/packaging       # clean-install tests (build wheel, venv, pipx)
ruff check src tests                                   # static analysis
```

Tools are found from `$RNASEQ_TOOLS_BIN` (the `bin/` of the tools environment), then `~/miniforge3/envs/rnaseq-tools/bin`,
then `PATH`. R is found from `$RNASEQ_RSCRIPT`, then `~/miniforge3/envs/rnaseq-r/bin/Rscript`, then `PATH`.

## What the suites prove

| Suite | Tests | What it establishes |
|---|---|---|
| `tests/unit` | 129 | validators, config rules, FASTQ validator, pairing, GEO parsing, count-matrix checks, design validation (replicates, confounding, formulas), checkpoint fingerprints, strandedness calls, reference compatibility, bacterial CDS detection, NCBI package building, store scanning, single-cell detection, featureCounts reconciliation, DESeq2 threshold-set re-derivation, manifest fingerprints, tool/RAM pre-flights, R/environment discovery, CPU-incompatible builds (SIGILL) and CPU-aware environment files, release consistency |
| `tests/security` | 8 | credentials never reach project configs, snapshots, logs, terminal or error messages; older projects are cleaned on open |
| `tests/failure` | 31 | missing tool, pipefail, empty arguments, no shell expansion, timeout, Ctrl-C, **SIGTERM/SIGHUP stop child processes**, validation pool stops at once, corrupted/truncated BAM, insufficient disk, invalid thresholds, R rejects bad designs, downloads (404, unreachable, empty, wrong size, wrong MD5, resume, restart when resume unsupported, never overwrite), pilot streaming retries/truncation, ENA→SRA fallback without mixing mates |
| `tests/integration` | 35 | the **real interactive CLI** on a synthetic dataset with known truth (16 DE genes): every stage, trimming, strandedness, truth recovery, report/manifest content, no credentials in project, resume without recomputation, dry run, corrupted FASTQ stops the run, a **13-scenario resume-invalidation matrix** (settings changes, corrupted/deleted outputs, edited metadata/results), report tampering detection, `--validate-project` PASS on a finished project and FAIL after a one-read edit, `--check` mini-job (200/200 pairs counted); plus unstranded / forward / single-end libraries and a deliberately wrong strandedness setting |
| `tests/packaging` | 8 | the built wheel installed into a fresh venv and with pipx, run with an empty HOME and minimal PATH from `~`, `/tmp` and a directory with spaces; package data present; `--check` fails cleanly without tools; the package directory is never written to |

Counts are those of release 1.1.0 (211 tests; 203 run by default, 8 packaging tests need `RNASEQ_PACKAGING_TESTS=1`).

## Distribution checks (release 1.1.0)

| Method | How it was checked | Result |
|---|---|---|
| wheel/sdist | `python -m build`, `twine check` | both PASSED, no build warnings |
| pip into a fresh venv, pipx | `tests/packaging` (8 tests) | pass |
| `pipx install git+https://github.com/vipin839/rnaseq_pipeline.git` | clean HOME and PATH | installs `rnaseq-pipeline` and `rnaseq_pipeline`; `--version` 1.1.0; package data found |
| Bioconda recipe | `bioconda-utils lint`; `conda-build` with the sdist as source; recipe tests in a fresh env | lint OK; build OK; tests pass incl. `--check` (0 failed) |
| GitHub Actions | lint, 8 unit jobs (Ubuntu 22.04/24.04 × Python 3.10–3.13), package, integration with the real tools | all 11 jobs green (commit 98e520c) |

## Synthetic truth dataset

`tests/data/make_synthetic.py OUT [pairs] [reverse|forward|unstranded] [paired|single]` writes a 3-chromosome genome, a
60-gene GTF and 6 samples (3 control, 3 treatment). GENE0001–0008 are 4× up, GENE0009–0016 4× down, the rest unchanged.
Negative-binomial noise is included, and 15% of fragments are shorter than the read length, so adapter read-through
triggers trimming. The truth table is written to `truth.tsv`.

Reference result (40,000 pairs per sample, reverse-stranded): **16/16 true DE genes recovered, 0 false positives,
58 genes tested**. Releases 1.0.0 and 1.1.0 give byte-identical count matrices and identical DESeq2 statistics
(maximum absolute difference 0).

## Real-data validation (regression evidence, not proof of general correctness)

| Dataset | What was checked | Result |
|---|---|---|
| GSE53720 (yeast, calorie restriction vs normal), pilot 500k reads | strandedness, rRNA diagnosis, DE direction | forward stranded; ~68% alignment explained by rDNA multi-mappers; glucose-repressed genes up (ADH2 +8.9, JEN1 +7.1, FBP1 +5.8 log2FC) |
| SRP314352 (yeast, low vs high glucose, wild type; 2 vs 2), pilot 200k pairs, release 1.1.0 | download robustness, reconciliation invariants on real reads, strandedness, DE direction, `--validate-project` | ENA dropped connections (curl 56): retries recovered 3 runs; SRR14208246 fell back to SRA, its lone ENA mate set aside. All 15 stages passed. Per sample: BAM primary = 2 × input pairs (e.g. 361,104 = 2 × 180,552), alignment 99.2–99.4%, featureCounts total ≥ fragments, matrix column sums = Assigned exactly. Reverse stranded (RSeQC). Glucose-repressed genes up in low glucose: HXT6 +9.3, SUC2 +6.5, CTA1 +5.4, HXT7 +5.0, JEN1 +4.9, HXT2 +4.4, ADH2 +1.9 log2FC (padj < 0.005); FBP1, PCK1 up but not significant. 430 up / 692 down. PROJECT HEALTH: PASS |
| SRR7361181 (SRA route) | prefetch → vdb-validate → fasterq-dump; a corrupted `.sra` is rejected | pass |
| E. coli UTI89 (NCBI assembly) | CDS-based annotation, plain index, `--no-spliced-alignment` | 4,954 genes |
| GSE309855 (third-party reanalysis) | GEO reanalysis links resolved | 94 runs from GSE199596 |
| SRR23333328 (10x Chromium) | single-cell refusal | refused with explanation |

One successful dataset does not prove all datasets are handled correctly. The synthetic truth tests measure
correctness; the real datasets show that behaviour on real archives and reads stays sensible.

## Writing tests

* Every bug fix gets a regression test that **fails without the fix**. Check this: the FASTQ-pool interrupt test was
  first written with files small enough to finish in time even without the fix, and was then enlarged.
* Do not weaken an assertion to make a test pass. If a test is wrong, say why in the test.
* Integration tests drive the real CLI with scripted answers. A new menu entry shifts the numbers: update
  `new_project_answers` in `tests/integration/test_end_to_end.py`.
