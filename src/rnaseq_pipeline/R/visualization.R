# Per-contrast visualizations: MA plot, volcano plot, significant-gene heatmap, top-gene heatmap.

run_contrast_plots <- function(dds, tr, res, meta, voi, p) {
  df <- res$table
  d <- file.path(p$plot_dir, res$name); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  td <- file.path(res$dir, "plot_data"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
  out <- list()
  a <- p$alpha; t <- p$log2fc_threshold
  cols <- c(up = "#C62828", down = "#1565C0", ns = "grey70")
  subtitle <- sprintf("padj < %s and |log2FC| >= %s", a, t)

  # MA plot
  lfc_col <- if ("log2FoldChange_shrunk" %in% names(df)) "log2FoldChange_shrunk" else "log2FoldChange"
  ma <- data.frame(gene_id = df$gene_id, baseMean = df$baseMean, log2FC = df[[lfc_col]], regulation = df$regulation)
  write_tsv_df(ma, file.path(td, "ma_plot_data.tsv"))
  g <- ggplot(ma[ma$baseMean > 0, ], aes(x = baseMean, y = log2FC, colour = regulation)) +
    geom_point(size = 0.8, alpha = 0.7) + scale_x_log10() + scale_colour_manual(values = cols) +
    geom_hline(yintercept = c(-t, 0, t), linetype = c("dashed", "solid", "dashed"), colour = "grey40") +
    labs(title = paste("MA plot:", res$name), subtitle = paste0(subtitle, if (lfc_col != "log2FoldChange") "  (y: shrunken LFC)" else ""),
         x = "Mean of normalized counts", y = "log2 fold change") + theme_pub()
  out$ma <- save_gg(g, d, "ma_plot", p)

  # Volcano plot
  vd <- df[!is.na(df$padj), c("gene_id", "log2FoldChange", "pvalue", "padj", "regulation")]
  write_tsv_df(vd, file.path(td, "volcano_data.tsv"))
  if (requireNamespace("EnhancedVolcano", quietly = TRUE) && nrow(vd) > 0) {
    lab_genes <- head(vd$gene_id[vd$regulation != "ns"], 20)
    g <- EnhancedVolcano::EnhancedVolcano(vd, lab = vd$gene_id, selectLab = lab_genes, x = "log2FoldChange", y = "padj",
      pCutoff = a, FCcutoff = t, title = paste("Volcano:", res$name), subtitle = subtitle,
      ylab = bquote(~-Log[10] ~ adjusted ~ italic(P)), legendPosition = "right", labSize = 3, pointSize = 1.5,
      drawConnectors = TRUE, caption = sprintf("%d genes with padj", nrow(vd)))
  } else {
    g <- ggplot(vd, aes(x = log2FoldChange, y = -log10(padj), colour = regulation)) + geom_point(size = 0.8) +
      scale_colour_manual(values = cols) + geom_vline(xintercept = c(-t, t), linetype = "dashed") +
      geom_hline(yintercept = -log10(a), linetype = "dashed") +
      labs(title = paste("Volcano:", res$name), subtitle = subtitle, y = "-log10 adjusted p") + theme_pub()
  }
  out$volcano <- save_gg(g, d, "volcano_plot", p, width = 8, height = 6)

  x <- assay(tr)
  ann <- annotation_df(meta, p$variables)
  z <- function(m) { s <- t(scale(t(m))); s[is.na(s)] <- 0; s }

  # Heatmap of significant genes (capped)
  sig_ids <- df$gene_id[df$regulation != "ns"]
  if (length(sig_ids) >= 2) {
    ids <- head(sig_ids, p$heatmap_max_genes)
    m <- z(x[ids, , drop = FALSE])
    write_tsv_df(data.frame(gene_id = rownames(m), round(m, 4), check.names = FALSE), file.path(td, "significant_heatmap_zscores.tsv"))
    out$significant_heatmap <- save_pheatmap(list(mat = m, annotation_col = ann, annotation_colors = annotation_colors(meta, p$variables), show_rownames = length(ids) <= 60,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101), breaks = seq(-3, 3, length.out = 102),
      main = sprintf("Significant genes (n=%d%s), z-scored %s", length(ids),
                     if (length(sig_ids) > length(ids)) paste0(" of ", length(sig_ids)) else "", p$transformation)),
      d, "significant_genes_heatmap", p, height = if (length(ids) <= 60) max(5, 0.15 * length(ids) + 2) else 8)
  } else {
    msg("INFO", res$name, ": fewer than 2 significant genes; significant-gene heatmap skipped")
  }

  # Top-gene heatmap (top N by adjusted p-value, regardless of fold-change cutoff)
  ranked <- df$gene_id[!is.na(df$padj)]
  ids <- head(ranked, p$top_n_genes)
  if (length(ids) >= 2) {
    m <- z(x[ids, , drop = FALSE])
    write_tsv_df(data.frame(gene_id = rownames(m), round(m, 4), check.names = FALSE), file.path(td, "top_genes_heatmap_zscores.tsv"))
    out$top_heatmap <- save_pheatmap(list(mat = m, annotation_col = ann, annotation_colors = annotation_colors(meta, p$variables), show_rownames = length(ids) <= 60,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101), breaks = seq(-3, 3, length.out = 102),
      main = sprintf("Top %d genes by adjusted p-value (z-scored %s)", length(ids), p$transformation)),
      d, "top_genes_heatmap", p, height = max(5, 0.15 * length(ids) + 2))
  }
  out
}
