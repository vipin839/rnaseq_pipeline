# Installation

The software has **two layers** that are installed differently:

| Layer | What | Installed with |
|---|---|---|
| Application | the `rnaseq-pipeline` command (Python, needs only PyYAML) | pip / pipx, or conda |
| Scientific runtime | FastQC, MultiQC, fastp, HISAT2, samtools, StringTie, Subread (featureCounts), SRA Toolkit, RSeQC, R + DESeq2 | conda / mamba (Bioconda) |

The scientific tools are not pip packages, so pip/pipx install **only the application**. The runtime comes from
conda. Nothing needs administrator rights, and nothing is ever installed with `sudo`.

After any method, run the health check. It runs a tiny real job through every tool and through DESeq2:

```bash
rnaseq-pipeline --check
```

## Requirements

* Linux x86_64. Tested on Ubuntu 22.04 and 24.04, including WSL2. macOS and Windows are not supported.
* Python 3.10 or newer (the test suite runs on 3.10, 3.11, 3.12 and 3.13).
* Disk: about 6 GB for the runtime, plus data. Put projects on ext4/xfs. The tool refuses FAT32 (4 GB file limit) and
  warns about Windows drives mounted in WSL (`/mnt/c`), which are slow.
* RAM: 8 GB for mammalian genomes. Small genomes (yeast, bacteria) work with less.
* CPU: any x86-64. On processors without AVX2/BMI2 (roughly pre-2013, e.g. Xeon E5 v1/v2), the Bioconda
  StringTie 3.x builds cannot run; `rnaseq-pipeline --check` detects this and gives the one-line fix
  (`mamba install -n rnaseq-tools -c conda-forge -c bioconda "stringtie=2.2.3"`).

## Option A — application with pipx, runtime with conda (recommended today)

```bash
# 1. the command (isolated environment, exposed on PATH) — verified in a clean HOME
pipx install git+https://github.com/vipin839/rnaseq_pipeline.git

# 2. conda (skip if you already have conda/mamba)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3

# 3. the scientific runtime, from the environment files shipped inside the package
~/miniforge3/bin/mamba env create -f "$(rnaseq-pipeline --runtime-env tools)"
~/miniforge3/bin/mamba env create -f "$(rnaseq-pipeline --runtime-env r)"

rnaseq-pipeline --check
```

The pipeline finds the environments named `rnaseq-tools` and `rnaseq-r` under `~/miniforge3` automatically. For another
location, set `environment.conda_root` in your user configuration (below). Step 3 can also be done interactively from
**main menu → 4. Manage Dependencies**, which shows every command and asks before running it.

## Option B — from a source checkout (developers)

```bash
git clone https://github.com/vipin839/rnaseq_pipeline.git
cd rnaseq_pipeline
python3 -m pip install -e ".[dev]"      # or: pipx install .
```

Without installing, `./rnaseq_pipeline` and `python3 rnaseq_pipeline.py` also run directly from the checkout. They
need PyYAML: `python3 -m pip install --user pyyaml`.

## Option C — everything with one conda command (Bioconda; after release)

A Bioconda recipe is prepared in `packaging/bioconda/rnaseq-pipeline/meta.yaml`. It installs the application **and**
the whole scientific runtime into one environment:

```bash
mamba create -n rnaseq -c conda-forge -c bioconda --strict-channel-priority rnaseq-pipeline
mamba activate rnaseq
rnaseq-pipeline --check
```

**Status:** the recipe passes `bioconda-utils lint` (All checks OK). The package was built locally with `conda-build`
(from the release source archive) and its recipe tests passed in a fresh environment, including
`rnaseq-pipeline --check`: `HEALTH: WARNING (38 passed, 1 warnings, 0 failed)` — the warning is the unset NCBI email;
the mini-job counted 200/200 read pairs and DESeq2 loaded. The recipe becomes installable from Bioconda only after it
has been submitted to and accepted by bioconda-recipes; see [RELEASE.md](RELEASE.md). Until then, use option A.

## Personal settings

Settings that belong to **you** rather than to a project — NCBI email and API key, where projects and references are
stored, the conda location — go in a user configuration file that is loaded automatically:

```bash
mkdir -p ~/.config/rnaseq-pipeline
cat > ~/.config/rnaseq-pipeline/config.yaml <<'EOF'
ncbi:
  email: "you@example.org"
  api_key: ""          # optional: raises NCBI's limit from 3 to 10 requests/second
projects_dir: ~/rnaseq_projects
reference_store: ~/rnaseq_references
EOF
chmod 600 ~/.config/rnaseq-pipeline/config.yaml
```

`--config FILE` adds another file on top, for example a lab profile. The API key and email are **never** copied into
projects, reports, manifests or logs.

## Where things are stored

| What | Default | Change with |
|---|---|---|
| Projects | `~/rnaseq_projects` | `projects_dir` or `--projects-dir` |
| Reference genomes and indexes (shared by all projects) | `~/rnaseq_references` | `reference_store` |
| User settings | `~/.config/rnaseq-pipeline/config.yaml` | `$XDG_CONFIG_HOME` |
| Install logs, recent state | `~/.rnaseq_pipeline` | `$RNASEQ_PIPELINE_HOME` |

The installed package directory is never written to.

## Upgrading from 1.0.0

* Projects created with 1.0.0 keep working. Opening one removes any NCBI credentials that 1.0.0 had copied into its
  `config/project_config.yaml`.
* Existing checkpoints stay valid. Checkpoints written from 1.1.0 onwards also record the settings each stage
  depended on, so changing one later re-runs exactly the affected stages.
* On resume, the report stage of a finished 1.0.0 project is re-run once, because its manifest lacks the new input
  fingerprints and reference checksums. Nothing else is recomputed.
* The source layout moved from `python/rnaseq` to `src/rnaseq_pipeline`, and new projects default to
  `~/rnaseq_projects` instead of the program folder. Existing projects anywhere can still be opened with
  `--project DIR` or "Enter a project path".

## Tested versions

FastQC 0.12.1, MultiQC 1.35, fastp 1.3.7, HISAT2 2.2.3, samtools 1.24, StringTie 3.0.3 (the current release of the StringTie2
line), Subread/featureCounts 2.1.1, SRA Toolkit 3.4.1, RSeQC 5.0.5, R 4.5.3, DESeq2 1.50.2.
