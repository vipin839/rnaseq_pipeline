# Installation

Everything installs in your home directory. **No administrator (sudo) rights are needed.**

## 1. Requirements

* Linux x86_64 (tested on Ubuntu 24.04 under WSL2)
* Python 3.10 or newer with PyYAML (`python3 -c "import yaml"`)
* About 6 GB of disk for the software, plus space for data and references
* 8 GB RAM minimum. Mammalian genomes need about 8 GB for HISAT2 alignment and index building.

Put projects on a **Linux filesystem (ext4/xfs)**, not a FAT32/exFAT USB disk or a Windows drive under `/mnt/c`.
FAT32 cannot hold files larger than 4 GB, and Windows drives mounted in WSL are slow. The pipeline
checks this and refuses FAT32.

## 2. Install Miniforge (conda + mamba)

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
```

## 3. Create the two environments

```bash
cd ~/rnaseq_pipeline
~/miniforge3/bin/mamba env create -f envs/rnaseq-tools.yml   # FastQC, MultiQC, fastp, HISAT2, samtools, StringTie, Subread, SRA Toolkit, RSeQC
~/miniforge3/bin/mamba env create -f envs/rnaseq-r.yml       # R, DESeq2, ggplot2, pheatmap, EnhancedVolcano, apeglm ...
```

You can also let the pipeline do this: **main menu → 4. Manage Dependencies / Environment →
Install missing components**. It shows exactly which commands it will run and asks before running them.

You do **not** need to "activate" the environments. The pipeline finds them in `~/miniforge3/envs/`.
To use a different location, set `environment.conda_root` in the configuration.

## 4. Check

```bash
./rnaseq_pipeline --check
```

Every component should be listed as `AVAILABLE`.

## Optional: enrichment packages

GO/Reactome enrichment is optional. Install it from main menu 4 → *Install optional annotation packages*,
or run:

```bash
~/miniforge3/bin/mamba install -n rnaseq-r -c conda-forge -c bioconda bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-reactomepa
```

## Tested versions (from the reference installation)

FastQC 0.12.1, MultiQC 1.35, fastp 1.3.7, HISAT2 2.2.3, samtools 1.24, StringTie 3.0.3 (the current release of the
StringTie2 line), Subread/featureCounts 2.1.1, SRA Toolkit 3.4.1, RSeQC 5.0.5, R 4.5.3, DESeq2 1.50.2.
