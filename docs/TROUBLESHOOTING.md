# Troubleshooting

Start with `logs/pipeline_errors.log`, then open the per-tool log named in the error (`logs/<stage>/<sample>.log`).
Every failure screen shows **Problem**, **Likely cause** and **Recommended action**.

| Symptom | Likely cause | What to do |
|---|---|---|
| `./rnaseq_pipeline: Permission denied` after cloning | The executable bit was lost (clone on a Windows/NTFS/exFAT drive, or a restrictive umask). `chmod +x` has no effect on `/mnt/c` style mounts | Run `bash rnaseq_pipeline` or `python3 rnaseq_pipeline.py` — both work without the bit. Better: clone into the Linux filesystem (e.g. `~/rnaseq_pipeline`), where `chmod +x rnaseq_pipeline` works |
| `bad interpreter: /usr/bin/env bash^M` | Windows line endings (CRLF) in the launcher | `sed -i 's/\r$//' rnaseq_pipeline`. The repo ships `.gitattributes` forcing LF, so a fresh clone should not have this |
| `No references found` under *Existing reference/index* | The folder has no `reference_manifest.yaml` and no FASTA+GTF pair, or it is not the configured `reference_store` | Use *Look in another folder*, or prepare it once via main menu 5 → *Download and prepare a reference now* |
| `program not found: hisat2` (or another tool) | The tools env is missing or incomplete | Main menu 4 → Install missing components, or `./rnaseq_pipeline --check` |
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
./rnaseq_pipeline --verbose --project DIR         # print every command
less DIR/logs/command_history.log                  # all commands, exit codes and durations
```
