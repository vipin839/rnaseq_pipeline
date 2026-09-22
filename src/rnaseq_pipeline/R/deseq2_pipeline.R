#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 analysis entry point (called by the Python pipeline).
#   Rscript deseq2_pipeline.R <params.json> [--validate-only]
# All thresholds/settings come from params.json (written from project config).
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: deseq2_pipeline.R <params.json> [--validate-only]")
validate_only <- "--validate-only" %in% args

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_path))
for (f in c("qc_plots.R", "differential_expression.R", "visualization.R", "downstream_analysis.R")) {
  source(file.path(script_dir, f))
}

msg <- function(tag, ...) cat(sprintf("[%s] %s\n", tag, paste0(...)))
fail <- function(...) { msg("ERROR", ...); quit(status = 2) }

main <- function() {
  p <- fromJSON(args[1], simplifyVector = TRUE)

  # ---- 1-2. load counts + metadata -----------------------------------------
  counts <- read.delim(p$counts, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  meta <- read.delim(p$metadata, check.names = FALSE, colClasses = "character", stringsAsFactors = FALSE)
  msg("INFO", sprintf("counts: %d genes x %d samples; metadata: %d rows", nrow(counts), ncol(counts), nrow(meta)))

  # ---- 3. sample matching ----------------------------------------------------
  if (!"sample" %in% colnames(meta)) fail("metadata has no 'sample' column")
  if (anyDuplicated(meta$sample)) fail("duplicate sample names in metadata")
  if (!setequal(colnames(counts), meta$sample))
    fail("samples differ between count matrix and metadata: ",
         paste(setdiff(union(colnames(counts), meta$sample), intersect(colnames(counts), meta$sample)), collapse = ", "))
  if (!identical(colnames(counts), meta$sample)) fail("sample order differs between count matrix and metadata")
  m <- as.matrix(counts)
  if (any(is.na(m))) fail("count matrix contains NA")
  if (any(m < 0)) fail("count matrix contains negative values")
  if (any(m != round(m))) fail("count matrix contains non-integer values")
  storage.mode(m) <- "integer"
  rownames(meta) <- meta$sample

  # ---- 4. design validation --------------------------------------------------
  vars <- p$variables
  miss <- setdiff(vars, colnames(meta))
  if (length(miss)) fail("design variables missing from metadata: ", paste(miss, collapse = ", "))
  for (v in vars) {
    if (any(meta[[v]] == "" | is.na(meta[[v]]))) fail("empty values in design variable '", v, "'")
    if (any(make.names(meta[[v]]) != meta[[v]])) fail("values of '", v, "' are not valid R names")
    meta[[v]] <- factor(meta[[v]])
  }
  voi <- p$variable_of_interest
  if (!p$reference_level %in% levels(meta[[voi]])) fail("reference level not present: ", p$reference_level)
  meta[[voi]] <- relevel(meta[[voi]], ref = p$reference_level)
  tab <- table(meta[[voi]])
  if (any(tab < 2)) fail("levels without replicates: ", paste(names(tab)[tab < 2], collapse = ", "))
  design <- as.formula(p$formula)
  mm <- model.matrix(design, data = meta)
  if (qr(mm)$rank < ncol(mm)) fail("design matrix is not full rank (confounded variables)")
  if (ncol(mm) >= nrow(mm)) fail("no residual degrees of freedom")
  msg("OK", "design validated: ", p$formula, " (", ncol(mm), " coefficients, ", nrow(mm), " samples)")
  if (validate_only) { msg("OK", "validate-only mode: stopping before DESeq2"); return(invisible()) }

  suppressPackageStartupMessages(library(DESeq2))
  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(p$plot_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- 5. DESeqDataSet -------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = m, colData = meta, design = design)

  # ---- 6. low-count filter (documented rule) ---------------------------------
  min_samples <- if (identical(p$min_count_filter$min_samples, "smallest_group")) min(tab) else as.integer(p$min_count_filter$min_samples)
  keep <- rowSums(counts(dds) >= p$min_count_filter$min_count) >= min_samples
  msg("INFO", sprintf("filter: keep genes with >= %d counts in >= %d samples: %d of %d genes kept",
                      p$min_count_filter$min_count, min_samples, sum(keep), length(keep)))
  if (sum(keep) < 2) fail("fewer than 2 genes pass the low-count filter")
  dds <- dds[keep, ]

  # ---- 7. DESeq --------------------------------------------------------------
  dds <- DESeq(dds, quiet = TRUE)
  msg("OK", "DESeq() model fitted")
  saveRDS(dds, file.path(p$out_dir, "dds.rds"))

  # ---- 8. normalized counts --------------------------------------------------
  nc <- counts(dds, normalized = TRUE)
  write_tsv_df(data.frame(gene_id = rownames(nc), round(nc, 4), check.names = FALSE),
               file.path(p$out_dir, "normalized_counts.tsv"))
  write_tsv_df(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
               file.path(p$out_dir, "size_factors.tsv"))

  # ---- QC / exploratory ------------------------------------------------------
  tr <- transform_counts(dds, p$transformation)
  qc <- run_qc_plots(dds, m, tr, meta, vars, voi, p)

  # ---- 9-10. contrasts -------------------------------------------------------
  summary <- list(genes_input = nrow(m), genes_tested = nrow(dds), min_samples_filter = min_samples,
                  formula = p$formula, reference_level = p$reference_level,
                  alpha = p$alpha, log2fc_threshold = p$log2fc_threshold, contrasts = list(),
                  qc_plots = qc)
  for (i in seq_len(nrow(p$contrasts))) {
    num <- p$contrasts[i, 1]; den <- p$contrasts[i, 2]
    res <- run_contrast(dds, voi, num, den, p)
    vis <- run_contrast_plots(dds, tr, res, meta, voi, p)
    res$summary$plots <- vis
    summary$contrasts[[res$name]] <- res$summary
  }

  # ---- optional downstream annotation ---------------------------------------
  summary$annotation <- run_downstream(summary, p)

  writeLines(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"),
             file.path(p$out_dir, "deseq2_summary.json"))
  writeLines(capture.output(sessionInfo()), file.path(p$out_dir, "R_sessionInfo.txt"))
  msg("OK", "DESeq2 analysis complete")
}

write_tsv_df <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

tryCatch(main(), error = function(e) {
  cat(sprintf("[ERROR] %s\n", conditionMessage(e)))
  quit(status = 1)
})
