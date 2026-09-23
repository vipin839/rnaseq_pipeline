# Pipeline Methods

For each stage: the biological purpose, what is computed, the key parameters (all in `src/rnaseq_pipeline/config/default_config.yaml`),
and the validation that must pass before the stage's checkpoint is written. A stage is never reported as successful
just because an output file exists.

## 1. Data acquisition
*Purpose:* obtain the raw sequencing reads exactly as deposited or supplied.
* ENA: `fastq_ftp` files over HTTPS, resumable (`curl -C -`), **MD5 verified** against ENA's `fastq_md5`, written to `*.part`
  and renamed only after verification.
* SRA: `prefetch` → `vdb-validate` (fails on corruption, exit ≠ 0) → `fasterq-dump --split-3 --skip-technical` in a
  temporary folder → `pigz`. Unpaired leftovers of paired runs are ignored and logged.
* GEO: GSE → GSM (title, source, characteristics) → linked SRX → runs via ENA.
* Discovery: NCBI E-utilities (`einfo` for searchable indexes, `esearch`/`esummary` for GEO series,
  `efetch` runinfo for SRA runs). Searching uses Entrez; downloading uses ENA or the SRA toolkit.
* If ENA has not mirrored a run (common for very recent studies) the SRA route is used automatically.
* Local: symlink (default) or copy. Originals are never written to.
* Pilot mode (`download.max_reads`) streams only the first N reads. It is for **testing only** and is flagged in the report and manifest.
  Interrupted streams are retried (3 attempts, increasing delay). If ENA keeps failing for a run, that run is taken
  from NCBI SRA instead (`fastq-dump -X N`); any mate already streamed from ENA is set aside as `*.other_archive`, so both
  mates always come from the same archive.
* Downloads without a provider checksum are still checked for emptiness and the size ENA reports. If a server cannot
  resume (curl 33) the download restarts from zero. Permanent errors (HTTP 404) are not retried.
* One organism per project: a selection spanning several organisms is refused (one reference cannot serve both).
* Every input file is fingerprinted in the manifest (SHA-256; for files over 50 MB, SHA-256 of size + first and last MiB).

## 2. FASTQ verification (streaming, constant memory)
Every record of every file is read: `@` header, `+` separator (and its optional ID matches), sequence alphabet
`ACGTNacgtn.`, quality characters within Phred+33 (`!`…`~`), `len(quality) == len(sequence)`, complete 4-line records
(truncation), non-empty file. gzip is decoded by an external `pigz`/`gzip` process, and its exit status detects corrupt
or truncated archives. Pairs are read **in lockstep**: identical read counts and matching mate IDs (the first token, with
`/1` `/2` removed). A Phred+64 heuristic gives a warning. Results go to `data/metadata/fastq_validation_report.tsv`.
Failed samples cannot proceed. The user can exclude them (they are kept on disk) or stop.

## 3. FastQC + MultiQC
`fastqc --threads N --noextract`. Validation: the HTML and ZIP exist, the ZIP contains `fastqc_data.txt` and `summary.txt`,
and **FastQC's "Total Sequences" equals the read count from step 2**. MultiQC aggregates the reports and its report and data
folder are checked.

## 4. Quality gate (advisory)
Metrics are parsed from `fastqc_data.txt`. Default thresholds (`quality_gate`):

| Metric | REVIEW | TRIMMING RECOMMENDED | Rationale |
|---|---|---|---|
| Adapter content (max over positions) | > 1 % | > 5 % | FastQC warns at 5 %, fails at 10 % |
| Mean quality over the last 10 % of positions | < Q28 | < Q20 | Q20 = 1 % error, the usual trimming cut-off |
| Per-base lower quartile (min) | < Q10 | – | FastQC failure criterion |
| Largest overrepresented sequence | > 1 % | – | Often biological in RNA-seq, so review only |
| N content (max) | > 5 % | – | FastQC warning level |
| Duplication | > 80 % | – | Highly expressed genes create duplicates, so it is not a trimming issue |
| %GC deviation from the cohort median | > 10 points | – | Possible contamination or rRNA |
| Reads | < 1,000,000 | – | Low depth reduces power |
| R1 ≠ R2 read counts | FAIL | | |

FastQC WARN/FAIL flags alone never make a sample unusable. **Accept recommendation** trims *all* samples if any
sample needs it, so that every sample is processed the same way.

## 5. fastp (only if chosen)
Defaults: `--detect_adapter_for_pe` (paired-end), `--qualified_quality_phred 20`, `--unqualified_percent_limit 40`,
`--cut_right --cut_right_window_size 4 --cut_right_mean_quality 20`, `--length_required 25`. Output goes to
`data/trimmed/*.partial.fastq.gz` and is renamed when complete. The trimmed FASTQs are then **fully re-validated**, and their read
count must equal fastp's JSON `after_filtering.total_reads` (divided by 2 for pairs). FastQC, MultiQC and the quality assessment are repeated.
Invalid trimmed output stops the pipeline.

## 6. Reference preparation
* Catalog references are downloaded and verified against the provider's checksums (GENCODE/UCSC/NCBI MD5,
  Ensembl BSD `sum`), and a SHA-256 is recorded. They are kept in a **shared store** (`reference_store`) and linked into projects.
* FASTA validation: header present, unique names, IUPAC alphabet, no empty sequences.
* GTF validation: 9 columns, integer coordinates with start ≤ end, strand `+ - .`, exons carry `gene_id` and
  `transcript_id`. GFF3 is detected and rejected with a conversion hint.
* **Compatibility:** chromosome names must be shared (with a chr-prefix diagnosis if not), annotation coordinates must not
  exceed chromosome lengths (this catches, for example, GRCh37 annotation on GRCh38), declared genome and annotation assemblies must agree, and the GTF
  header build must match. Any violation stops the stage unless the user typed `OVERRIDE` at selection, which is recorded.
* `samtools faidx`, BED12 from the GTF (for RSeQC), `hisat2_extract_splice_sites.py`, `hisat2_extract_exons.py`.
* **Annotations without exons** (bacteria; NCBI GTFs annotate genes as `CDS`): the feature-type counts are
  shown and you are asked whether to count `CDS` instead. StringTie2 is then skipped, a plain HISAT2 index is
  built (an empty splice-site file makes `hisat2-build` abort), alignment adds `--no-spliced-alignment`, and the
  RSeQC BED12 is derived from the counted feature.
* **HISAT2 index:** splice sites and exons are built into the index only when the estimated RAM (about 60 bytes per genome bp,
  roughly 160–200 GB for human) is available. Otherwise the plain index is built (about 8 GB for human) and the known splice
  sites are passed at alignment time with `--known-splicesite-infile`. This is standard practice with a small loss of sensitivity. The index is built in
  `index.building/` and validated with `hisat2-inspect -n`, whose sequence names must equal the FASTA names, before being moved into place. If `hisat2-build` crashes (SIGSEGV/SIGBUS/SIGABRT — an intermittent HISAT2 2.2.x multithreading crash), the partial index is discarded and the build is retried once with a single thread; the resulting index is identical. An existing
  valid index is never rebuilt.

## 7. Alignment
`hisat2 -p N --dta --new-summary --rg-id S --rg SM:S [--known-splicesite-infile] (-1 R1 -2 R2 | -U R)` piped
into `samtools sort -@ N -m MEM -T temp` → `S.partial.bam`. **Every process in the pipe is checked** (pipefail semantics).
Then `samtools quickcheck`, rename to `S.sorted.bam`, `samtools index`. No SAM file is ever written.
`--dta` is used because StringTie is downstream (it favours alignments with longer anchors).

*Validation (a BAM failing any check is renamed `S.FAILED.bam` and never enters later stages):* quickcheck (EOF block), the index must describe this BAM (records counted by `samtools idxstats` = records in the BAM; judged by content, never by file times),
header `SO:coordinate`, index present and newer than the BAM, `samtools flagstat`, **primary records = input reads
(×2 for paired-end)** (HISAT2 keeps unaligned reads, so any shortfall means a truncated BAM), HISAT2's processed total =
FASTQ read count, all primary records flagged paired (for paired-end), overall alignment rate ≥ `min_overall_alignment_rate_fail` (20 %;
warning below 70 %). Summary: `alignment/reports/alignment_summary.tsv`.

## 8. BAM QC
`samtools idxstats`, `flagstat`, `stats` per sample, a combined `bam_qc_summary.tsv`, and MultiQC over the HISAT2 + samtools reports.

## 9. Strandedness
RSeQC `infer_experiment.py` on 200,000 reads per sample against the BED12. Call: *forward* if ≥ 80 % of
informative reads are in the forward orientation (`1++,1--,2+-,2-+` / `++,--`), *reverse* if ≥ 80 % are reverse, *unstranded* if
the two fractions differ by ≤ 0.10. It is undetermined if more than 50 % of reads cannot be assigned. All samples must agree. The
user confirms, and if the result is undetermined the user must choose. Mapping: forward → featureCounts `-s 1` and StringTie `--fr`; reverse (dUTP/TruSeq
Stranded) → `-s 2` and `--rf`; unstranded → `-s 0`. Stored with its source and confidence.

## 10. StringTie2 (independent result)
`stringtie BAM -G annotation.gtf -e -o S.gtf -A S.gene_abund.tab -p N [--fr|--rf]`. `-e` = **reference-guided
quantification of annotated transcripts only** (no novel assembly), so samples are directly comparable and no merge is
needed. Validation: the GTF has transcript records with TPM, and the gene-abundance header and rows are present. Combined TPM matrices are in
`stringtie/merged/`. These values are **not** the DESeq2 input.

## 11. featureCounts (the DESeq2 input)
`featureCounts -T N -a GTF -t exon -g gene_id -s {0,1,2} -Q 0 [-p --countReadPairs -B -C]`. For paired-end data it counts
fragments, requires both ends mapped, and excludes chimeric pairs. Multi-mappers are not counted by default (`-M` off).
Validation: program header, fixed columns, sample columns equal to the BAM list in order, unique non-empty gene IDs, numeric
non-negative values, **gene count = number of genes in the GTF**. Assignment below 30 % gives a warning (check strandedness,
annotation, rRNA).

*Reconciliation with the alignment (independent layer):* for each sample, featureCounts' total must be at least the number
of fragments entering alignment (reads for single-end), `Assigned` can never exceed it, and a total above fragments +
secondary/supplementary alignments gives a warning. A violation stops the stage.

## 12. Count matrix
The annotation columns are removed and columns are renamed to sample IDs. Counts must be integers (fractional counting is refused).
If runs share a biological sample, the user may sum them, the recommended handling of lanes/technical replicates (the per-run matrix is kept).
Validation: header, unique sample and gene IDs, no annotation columns, no NA/negative/non-integer values, expected dimensions,
no zero-library samples, and **every column sum equals featureCounts `Assigned` for that sample** (checked before and after
the file is written; with summed technical replicates, the sum of their `Assigned`).

## 13. Experimental design
Stored in `counts/sample_metadata.tsv` + `counts/design.json` (with the metadata SHA-256, so DESeq2 refuses a design that was
edited after confirmation). Checks: every sample has every design variable, values are valid R names, the variable of interest has ≥ 2
levels with ≥ 2 replicates each, covariates have ≥ 2 levels, the design matrix is **full rank** (no confounding), and there are residual degrees
of freedom. The design is re-checked in R (`--validate-only`) before being confirmed.

## 14. DESeq2 (R)
1. Load the counts and metadata, then check sample identity and order, NA, negative and integer values.
2. Factors are created and the reference level set with `relevel`. The design is re-validated (replicates, `qr` rank, residual df).
3. **Low-count filter:** keep genes with ≥ `min_count` (10) reads in ≥ `min_samples` samples (default = size of the
   smallest group). This removes genes without enough information to test.
4. `DESeq()`. Normalised counts (median-of-ratios) and size factors are saved.
5. For each contrast: `results(contrast = c(var, num, den), alpha, lfcThreshold, independentFiltering, cooksCutoff)`.
   P-values are adjusted with Benjamini–Hochberg (`padj`).
6. **Significance:** up = `padj < alpha AND log2FC ≥ t`; down = `padj < alpha AND log2FC ≤ −t` (defaults alpha 0.05,
   t = 1, i.e. 2-fold). The LFC used is the unshrunken MLE. apeglm shrinkage (or `normal` when no model coefficient exists) is added
   as `log2FoldChange_shrunk` for the MA plot and ranking. The thresholds are printed in every `statistics.txt`, plot and the report.
7. Exploratory plots use a VST (or rlog) with `blind = TRUE`: PCA (top 500 variable genes), sample distance and
   correlation heatmaps, library sizes, normalised count distributions. Per contrast: MA, volcano, significant-gene heatmap
   (z-scores, capped at `heatmap_max_genes`), and a top-N heatmap. PNG (300 dpi) + PDF. The data behind every plot is saved as TSV.
8. Python then validates the outputs: every file exists, the table sizes match the summary, and the up, down and
   significant sets are **re-derived from the full result table** (padj < alpha and |log2FC| ≥ threshold, direction by sign).
   A gene missing from or wrongly added to a table, or labelled with the wrong direction, stops the stage.
9. Python also **recomputes, without trusting R**, what DESeq2 must have produced from the count matrix and the
   metadata: the tested genes must be exactly the genes passing the low-count filter; the size factors must equal
   DESeq2's median-of-ratios estimate (relative tolerance 1e-6) and the normalized counts must equal
   counts / size factor; `baseMean` must be the mean of the normalized counts; p-values and padj must lie in [0, 1]
   with padj ≥ p-value; and for every contrast at least 90% of the significant genes with |log2FC| ≥ 1 must have
   a fold-change sign that matches mean(numerator group) vs mean(baseline group) of the normalized counts (a
   swapped baseline gives ~0%; checked when there are at least 5 such genes). The same checks run again on resume
   and in `--validate-project`, with the recorded paths re-pointed if the project was moved.

## 15. Optional enrichment
If enabled and the packages are installed: clusterProfiler `enrichGO` (BP) and ReactomePA `enrichPathway` (human/mouse) on the up
and down gene sets, against the tested-gene universe. It is never required for the core results.

## Report and manifest
The HTML report marks every key number (active samples, reference, strandedness, formula, alpha, log2FC threshold, and
per contrast the tested/up/down/significant counts). After writing, the report is re-read and each marked number is
compared with the value re-derived from the result files; every link and image must exist, the document must be complete,
and no template placeholder may remain. A pilot-subset run is labelled as such at the top. The manifest records software
versions, the redacted configuration, the design, input-file fingerprints, provider MD5s, the genome and annotation
SHA-256, where the reference came from, and the command log. `--validate-project` repeats all of these checks.

## Settings and resume
Each stage records the settings its results depend on. If one changes, the stage and everything after it are marked
INVALID with the setting named and re-run on resume. Settings that cannot change results (threads, memory) invalidate
nothing.

## Storage estimates
Trimmed FASTQ ≈ 0.9× raw; sorted BAM ≈ 0.8× FASTQ.gz plus sort temp ≈ the largest sample; `fasterq-dump` temp ≈ 8× `.sra`.
You are asked to confirm when free space is below 1.5× the estimate, and the stage is refused when free space is below the estimate.
