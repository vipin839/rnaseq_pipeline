# Release procedure

The version is defined in exactly one place, `src/rnaseq_pipeline/__init__.py` (`__version__`). `pyproject.toml` reads
it, and `tests/unit/test_release.py` fails if the CHANGELOG top entry or the Bioconda recipe disagree.

## Release readiness of 1.1.0 (P3, 25 September 2026)

| Item | Status |
|---|---|
| Core stabilization (P0, P1) and real-data acceptance (P2: Arabidopsis full data; human GRCh38) | done — see docs/TESTING.md |
| Test suite | 273 tests pass (265 by default + 8 packaging) |
| CI | 11 jobs green on Node 24 actions (checkout v7, setup-python v7, upload-artifact v7, setup-micromamba v3); no warnings |
| Clean install from GitHub (`pipx install git+…`, empty HOME, minimal PATH) | verified: commands on PATH, version 1.1.0, runtime environment files found (CPU-aware), clear errors |
| Bioconda recipe | run dependencies equal the tested environment files (enforced by `test_release.py`); dependencies solve in one environment with the tested versions; `bioconda-utils lint`: All checks OK |
| Tag `v1.1.0` and GitHub release (sdist + wheel attached) | done, 25 September 2026 (step 3 below) |
| PyPI, Bioconda submission | **not done — only on the owner's request** (steps 4–5 below) |

## 1. Prepare

1. Set `__version__` (semantic versioning: MAJOR incompatible, MINOR features, PATCH fixes).
2. Add the `## X.Y.Z — YYYY-MM-DD` entry at the top of `CHANGELOG.md`. Describe behaviour changes as
   old → new → reason → consequence.
3. Update `{% set version = "X.Y.Z" %}` in `packaging/bioconda/rnaseq-pipeline/meta.yaml` and reset `build: number` to 0.
4. Update `docs/VERIFICATION_MATRIX.md` for any requirement whose status changed.

## 2. Verify locally (all must pass)

```bash
ruff check src tests
pytest tests                                       # with the scientific tools available
RNASEQ_PACKAGING_TESTS=1 pytest tests/packaging
python -m build && twine check dist/*
rnaseq-pipeline --check                            # HEALTH: PASS or WARNING, never FAIL
```

Then run the synthetic truth dataset through the CLI and compare it with the previous release (see docs/TESTING.md):
the count matrix and DESeq2 results must be identical, unless a documented change explains the difference.

## 3. Tag and publish the source

```bash
git commit -am "Release X.Y.Z"
git tag -a vX.Y.Z -m "rnaseq-pipeline X.Y.Z"
git push origin main vX.Y.Z
```

CI (`.github/workflows/ci.yml`) runs lint, unit/security/failure tests on Ubuntu 22.04 and 24.04 × Python 3.10–3.13,
packaging tests, and the full integration suite with the real tools. Do not continue unless it is green.

Create the GitHub release from the tag and attach `dist/*` (sdist and wheel).

## 4. PyPI (optional)

```bash
twine upload dist/*            # needs a PyPI account and API token
```

After this, `pipx install rnaseq-pipeline` works. Verify it in a clean shell: `pipx install rnaseq-pipeline &&
rnaseq-pipeline --version`.

## 5. Bioconda

1. Fill the source checksum from the published tag:
   `packaging/bioconda/update_sha256.sh`
2. Fork https://github.com/bioconda/bioconda-recipes, copy `packaging/bioconda/rnaseq-pipeline/` to
   `recipes/rnaseq-pipeline/`, commit on a branch.
3. Check locally (bioconda-utils environment: `mamba create -n bioconda-utils -c conda-forge -c bioconda bioconda-utils`):
   ```bash
   bioconda-utils lint recipes config.yml --packages rnaseq-pipeline
   bioconda-utils build recipes config.yml --packages rnaseq-pipeline          # add --docker --mulled-build-and-test if Docker is available
   ```
4. Open the pull request. Bioconda CI builds the package and runs the recipe tests, including
   `rnaseq-pipeline --check`, which runs HISAT2, samtools, featureCounts, StringTie, fastp, FastQC and DESeq2 on a tiny job.
5. After it is merged: `mamba create -n rnaseq -c conda-forge -c bioconda rnaseq-pipeline`.

Later versions are usually picked up by the Bioconda autobump bot from new GitHub tags.

## 6. After publishing

Install with each published method in a clean environment and run `rnaseq-pipeline --check`. Only then list the
method as available in README.md and docs/INSTALLATION.md.
