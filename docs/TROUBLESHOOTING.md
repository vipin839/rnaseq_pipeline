# Troubleshooting

Start with `logs/pipeline_errors.log`, then open the per-tool log named in the error (`logs/<stage>/<sample>.log`).
Every failure screen shows **Problem**, **Likely cause** and **Recommended action**.

| Symptom | Likely cause | What to do |
|---|---|---|
| `./rnaseq_pipeline: Permission denied` after cloning | The executable bit was lost (clone on a Windows/NTFS/exFAT drive, or a restrictive umask). `chmod +x` has no effect on `/mnt/c` style mounts | Run `bash rnaseq_pipeline` or `python3 rnaseq_pipeline.py` — both work without the bit. Better: clone into the Linux filesystem (e.g. `~/rnaseq_pipeline`), where `chmod +x rnaseq_pipeline` works |
| `bad interpreter: /usr/bin/env bash^M` | Windows line endings (CRLF) in the launcher | `sed -i 's/\r$//' rnaseq_pipeline`. The repo ships `.gitattributes` forcing LF, so a fresh clone should not have this |
| `No references found` under *Existing reference/index* | The folder has no `reference_manifest.yaml` and no FASTA+GTF pair, or it is not the configured `reference_store` | Use *Look in another folder*, or prepare it once via main menu 5 → *Download and prepare a reference now* |
| `SINGLE-CELL DATA DETECTED` / `more than two reads per spot` | The dataset is single-cell RNA-seq (e.g. 10x Genomics Chromium: an 8 bp index read, a 28 bp cell-barcode+UMI read and the cDNA read) | This pipeline is for **bulk** RNA-seq; results on single-cell reads would be wrong. Use Cell Ranger, STARsolo or alevin-fry, or choose a bulk dataset (ENA library source `TRANSCRIPTOMIC`, not `TRANSCRIPTOMIC SINGLE CELL`) |
| `all N sample(s) are FAILED/EXCLUDED` on a step screen | Every sample failed an earlier step (e.g. all downloads failed) | Project menu → *Continue pipeline*: it lists each sample's failure and offers **Retry** (the failed step runs again) or **Replace the data source**. Re-running a step always retries the samples that failed in it |
| `GEO series … has no samples with sequencing data` | Microarray series, not yet public, or raw data held elsewhere (e.g. dbGaP) | Reanalysis and SuperSeries records are followed automatically; otherwise use the accession of the study that holds the raw reads |
| `--config: file not found` / `is not valid YAML` | Wrong path or broken configuration file | Check the file named in the message (`--config FILE` or `~/.config/rnaseq-pipeline/config.yaml`): two-space indentation, values in quotes |
| `program not found: hisat2` / `required software not found for this step` | The tools env is missing or incomplete | Main menu 4 → Install missing components, or `mamba env create -f "$(rnaseq-pipeline --runtime-env tools)"`, then `rnaseq-pipeline --check`. New environments created with `--runtime-env tools` or menu 4 get the compatible version automatically |
| `rnaseq-pipeline: command not found` after `pipx install` | `~/.local/bin` is not on PATH | `pipx ensurepath`, then open a new terminal |
| `stringtie … killed by SIGILL (illegal CPU instruction)` / `exit -4` | Every StringTie 3.x build on Bioconda needs AVX2/BMI2 (CPUs from ~2013 on); this processor is older (see *CPU instruction set* in `--check`) | Reinstalling will not help. Install the compatible build: `mamba install -n rnaseq-tools -c conda-forge -c bioconda "stringtie=2.2.3"` (or main menu 4, which proposes it), then `rnaseq-pipeline --check` |
| `HEALTH: FAIL` from `--check` | A tool is installed but does not work, or R/DESeq2 is broken | The row marked FAILED names the component and the fix; usually recreate that conda environment |
| `HEALTH: WARNING` — NCBI email not set | No `ncbi.email` in your user configuration | Add it to `~/.config/rnaseq-pipeline/config.yaml` (NCBI requires an email for E-utilities) |
| Stage shows `INVALID` with `setting '…' changed (old -> new)` | You changed a setting that affects that stage's results | Expected. Continue the pipeline and the stage re-runs with the new value; or change the setting back |
| `featureCounts processed … fragments but the BAM holds …` / `column sum … != featureCounts Assigned` | Counts were lost or duplicated between steps (disk full, killed process, edited file) | The stage stopped instead of producing wrong counts. Free space and resume; if it repeats, report it with the stage log |
| `report shows … but the result files give …` | The report was edited, or results changed after the report was written | Resume; the report stage re-runs |
| `pilot download … failed after 3 attempts` then `trying the SRA route` | ENA dropped the connection repeatedly (common on some networks) | Nothing to do if the SRA route then succeeds; the single ENA mate is kept aside as `*.other_archive` so mates never mix archives |
| Terminal closed or `kill` during a run | — | Child processes are stopped and partial files are never used. Resume the project |
| `PyYAML is missing` | The system Python lacks PyYAML | `python3 -m pip install --user pyyaml`, or create the conda envs (the launcher falls back to their Python) |
| `project filesystem is vfat` | Project on a FAT32 USB disk (4 GB file limit) | Create the project on an ext4 disk (e.g. under your home directory) |
| FASTQ `truncated file` / `gzip integrity failure` | Incomplete download or copy | Delete the file and download it again. For ENA the MD5 is re-checked |
| `paired read count mismatch` / `mate ID mismatch` | R1 and R2 are not mates, or one is truncated | Check the file names and pairing, and re-download both |
| `invalid quality character` / `possible Phred+64` | Old Illumina 1.3–1.7 encoding | Convert, e.g. `seqtk seq -Q64 -V in.fq > out.fq`, then re-import |
| `genome and annotation are not compatible` | Different assemblies or providers (e.g. GRCh37 GTF with GRCh38 FASTA, or `chr1` vs `1`) | Use a catalog package, or a matching FASTA + GTF from the same release |
| `annotation looks like GFF3` | GFF3 supplied instead of GTF | `gffread annotation.gff3 -T -o annotation.gtf` |
| HISAT2 index build killed / very slow | Not enough RAM | Keep `use_splice_sites_in_index: auto` (sites are then supplied at alignment), close other programs, or use a larger machine |
| `primary alignments != input reads` | BAM truncated (disk full, killed process) | Free space and resume. The BAM was set aside as `.FAILED.bam` / `.invalid.bam` |
| `overall alignment rate … < 20%` | Wrong organism or reference, or heavy contamination | Check the organism in the metadata and the FastQC overrepresented sequences |
| Low featureCounts assignment (< 30 %) | Wrong strandedness, many multi-mappers (rRNA), or an incomplete annotation | Look at `featurecounts.summary`: large `Unassigned_NoFeatures` suggests strandedness or annotation; large `Unassigned_MultiMapping` suggests rRNA/repeats |
| `Library strandedness could not be automatically established` | Ambiguous evidence (degraded or low-depth libraries) | Check the library kit (dUTP/TruSeq Stranded = reverse; ligation-based kits are often forward) and select it |
| Design: `levels without biological replicates` | A group has only one sample | Add replicates, merge groups, or drop the sample. DESeq2 cannot estimate variance without replicates |
| Design: `not full rank` | Batch and condition are confounded (each batch has only one condition) | Remove batch from the formula. The effects cannot be separated statistically |
| `experimental design is not confirmed (or metadata changed since confirmation)` | `sample_metadata.tsv` was edited after confirmation | Project menu → Experimental design → confirm again |
| `this project is already open in another pipeline process` | Two sessions on one project | Close the other session. The lock is released automatically when it exits |
| Stage shows `INVALID` after resume | An output was modified, moved or deleted, or an upstream stage re-ran | Continue the pipeline and that stage will re-run |
| Download keeps failing | Network problem or server outage | Resume later. Downloads continue from where they stopped (`.part` files) |

## Getting more detail

```bash
rnaseq-pipeline --verbose --project DIR          # print every command
rnaseq-pipeline --validate-project DIR           # what is valid, what will re-run
less DIR/logs/command_history.log                  # all commands, exit codes and durations
```
