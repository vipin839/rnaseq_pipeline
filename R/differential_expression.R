# Differential expression per contrast + up/down classification with explicit thresholds.

run_contrast <- function(dds, voi, num, den, p) {
  name <- paste0(voi, "_", num, "_vs_", den)
  cdir <- file.path(p$out_dir, name)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  res <- results(dds, contrast = c(voi, num, den), alpha = p$alpha, lfcThreshold = p$lfc_test_threshold,
                 independentFiltering = isTRUE(p$independent_filtering), cooksCutoff = isTRUE(p$cooks_cutoff))

  # LFC shrinkage (used for MA plot / ranking only; significance uses the unshrunken Wald test)
  shrink_method <- p$lfc_shrinkage
  shr <- NULL
  if (shrink_method != "none") {
    coef_name <- paste0(voi, "_", num, "_vs_", den)
    shr <- tryCatch({
      if (shrink_method == "apeglm" && coef_name %in% resultsNames(dds) && requireNamespace("apeglm", quietly = TRUE)) {
        lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
      } else {
        if (shrink_method == "apeglm") shrink_method <- "normal"  # apeglm needs a model coefficient
        lfcShrink(dds, contrast = c(voi, num, den), type = "normal", quiet = TRUE)
      }
    }, error = function(e) { msg("WARNING", "LFC shrinkage failed: ", conditionMessage(e)); shrink_method <<- "none"; NULL })
  }

  df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                   lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                   check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(shr)) df$log2FoldChange_shrunk <- shr$log2FoldChange[match(df$gene_id, rownames(shr))]
  df <- df[order(df$padj, -abs(df$log2FoldChange), na.last = TRUE), ]

  a <- p$alpha; t <- p$log2fc_threshold
  sig <- !is.na(df$padj) & df$padj < a & abs(df$log2FoldChange) >= t
  df$regulation <- ifelse(sig & df$log2FoldChange >= t, "up", ifelse(sig & df$log2FoldChange <= -t, "down", "ns"))
  up <- df[df$regulation == "up", ]
  down <- df[df$regulation == "down", ]

  write_tsv_df(df, file.path(cdir, "full_results.tsv"))
  write_tsv_df(df[sig, ], file.path(cdir, "significant_results.tsv"))
  write_tsv_df(up, file.path(cdir, "upregulated.tsv"))
  write_tsv_df(down, file.path(cdir, "downregulated.tsv"))
  n <- p$top_n_genes
  write_tsv_df(head(up, n), file.path(cdir, "top_upregulated.tsv"))
  write_tsv_df(head(down, n), file.path(cdir, "top_downregulated.tsv"))

  crit <- c(
    sprintf("Contrast: %s = %s vs %s (reference/denominator: %s)", voi, num, den, den),
    sprintf("Design: %s", p$formula),
    "Test: DESeq2 Wald test; p-values adjusted by Benjamini-Hochberg (padj)",
    sprintf("lfcThreshold in results(): %s", p$lfc_test_threshold),
    sprintf("Independent filtering: %s; Cook's distance outlier filtering: %s", isTRUE(p$independent_filtering), isTRUE(p$cooks_cutoff)),
    "",
    "SIGNIFICANCE CRITERIA",
    sprintf("  Upregulated:   padj < %s AND log2FoldChange >= %s", a, t),
    sprintf("  Downregulated: padj < %s AND log2FoldChange <= -%s", a, t),
    "  (log2FoldChange = unshrunken MLE estimate; positive = higher in the numerator level)",
    sprintf("  LFC shrinkage for plots/ranking: %s", shrink_method),
    "",
    sprintf("Genes tested (after low-count filter): %d", nrow(df)),
    sprintf("Genes with padj not NA: %d", sum(!is.na(df$padj))),
    sprintf("Significant: %d  (up %d, down %d)", sum(sig), nrow(up), nrow(down)),
    "",
    "DESeq2 summary():",
    capture.output(summary(res, alpha = a))
  )
  writeLines(crit, file.path(cdir, "statistics.txt"))
  msg("OK", sprintf("%s: %d significant (%d up, %d down)", name, sum(sig), nrow(up), nrow(down)))

  list(name = name, dir = cdir, table = df,
       summary = list(numerator = num, denominator = den, genes_tested = nrow(df),
                      padj_not_na = sum(!is.na(df$padj)), significant = sum(sig), up = nrow(up), down = nrow(down),
                      shrinkage = shrink_method, dir = cdir))
}
