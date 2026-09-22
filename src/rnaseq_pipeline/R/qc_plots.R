# QC / exploratory plots: library sizes, count distributions, PCA, sample distances, correlation.
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          legend.position = "right")
}

save_gg <- function(plot, dir, name, p, width = 7, height = 5) {
  files <- c()
  for (fmt in p$plot_formats) {
    f <- file.path(dir, paste0(name, ".", fmt))
    if (fmt == "png") ggsave(f, plot, width = width, height = height, dpi = p$plot_dpi)
    else ggsave(f, plot, width = width, height = height)
    files <- c(files, f)
  }
  files
}

save_pheatmap <- function(args, dir, name, p, width = 7, height = 6) {
  files <- c()
  for (fmt in p$plot_formats) {
    f <- file.path(dir, paste0(name, ".", fmt))
    do.call(pheatmap, c(args, list(filename = f, width = width, height = height)))
    files <- c(files, f)
  }
  grDevices::graphics.off()
  files
}

transform_counts <- function(dds, method) {
  if (method == "rlog") return(rlog(dds, blind = TRUE))
  # vst() subsamples 1000 genes to fit the trend; fall back for small gene sets
  if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
}

annotation_df <- function(meta, vars) {
  a <- as.data.frame(meta[, vars, drop = FALSE])
  rownames(a) <- meta$sample
  a
}

# Same colours as ggplot's default discrete palette, so groups look identical in PCA and heatmaps.
annotation_colors <- function(meta, vars) {
  out <- list()
  for (v in vars) {
    lv <- levels(factor(meta[[v]]))
    out[[v]] <- setNames(scales::hue_pal()(length(lv)), lv)
  }
  out
}

run_qc_plots <- function(dds, raw, tr, meta, vars, voi, p) {
  d <- p$plot_dir; tdir <- p$table_dir
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
  out <- list()

  # library sizes
  ls <- data.frame(sample = colnames(raw), library_size = colSums(raw), group = meta[[voi]])
  write_tsv_df(ls, file.path(tdir, "library_sizes.tsv"))
  g <- ggplot(ls, aes(x = sample, y = library_size / 1e6, fill = group)) + geom_col() +
    labs(title = "Library size (reads assigned to genes)", x = NULL, y = "Million reads", fill = voi) +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  out$library_sizes <- save_gg(g, d, "library_sizes", p, width = max(6, 0.4 * ncol(raw) + 3))

  # normalized count distributions
  nc <- log2(counts(dds, normalized = TRUE) + 1)
  long <- data.frame(sample = rep(colnames(nc), each = nrow(nc)), value = as.vector(nc),
                     group = rep(meta[[voi]], each = nrow(nc)))
  qs <- t(apply(nc, 2, quantile, probs = c(0, .25, .5, .75, 1)))
  write_tsv_df(data.frame(sample = rownames(qs), qs, check.names = FALSE), file.path(tdir, "normalized_count_distribution.tsv"))
  g <- ggplot(long, aes(x = sample, y = value, fill = group)) + geom_boxplot(outlier.size = 0.3) +
    labs(title = "Normalized counts (log2 + 1)", x = NULL, y = "log2(normalized count + 1)", fill = voi) +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  out$count_distribution <- save_gg(g, d, "normalized_count_distribution", p, width = max(6, 0.4 * ncol(raw) + 3))

  # PCA (top 500 most variable genes of the transformed data)
  x <- assay(tr)
  ntop <- min(500, nrow(x))
  sel <- order(apply(x, 1, var), decreasing = TRUE)[seq_len(ntop)]
  pca <- prcomp(t(x[sel, , drop = FALSE]))
  ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pcs <- data.frame(sample = colnames(x), PC1 = pca$x[, 1], PC2 = pca$x[, min(2, ncol(pca$x))], meta[, vars, drop = FALSE],
                    check.names = FALSE)
  write_tsv_df(pcs, file.path(tdir, "pca_coordinates.tsv"))
  write_tsv_df(data.frame(PC = paste0("PC", seq_along(ve)), variance_percent = ve), file.path(tdir, "pca_variance.tsv"))
  aes_map <- aes(x = PC1, y = PC2, colour = .data[[voi]])
  g <- ggplot(pcs, aes_map) + geom_point(size = 3.5) +
    ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
    labs(title = sprintf("PCA (%s, top %d variable genes)", p$transformation, ntop),
         x = sprintf("PC1: %.1f%% variance", ve[1]), y = sprintf("PC2: %.1f%% variance", ve[min(2, length(ve))])) +
    theme_pub()
  if (length(vars) > 1) g <- g + aes(shape = .data[[vars[1]]])
  out$pca <- save_gg(g, d, "pca", p)

  ann <- annotation_df(meta, vars)
  # sample distance heatmap
  dm <- as.matrix(dist(t(x)))
  write_tsv_df(data.frame(sample = rownames(dm), dm, check.names = FALSE), file.path(tdir, "sample_distances.tsv"))
  out$sample_distance_heatmap <- save_pheatmap(list(mat = dm, clustering_distance_rows = dist(t(x)),
    clustering_distance_cols = dist(t(x)), annotation_col = ann, annotation_colors = annotation_colors(meta, p$variables),
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100), main = "Sample-to-sample Euclidean distance"),
    d, "sample_distance_heatmap", p)

  # sample correlation heatmap
  cm <- cor(x, method = "pearson")
  write_tsv_df(data.frame(sample = rownames(cm), round(cm, 5), check.names = FALSE), file.path(tdir, "sample_correlation.tsv"))
  out$sample_correlation_heatmap <- save_pheatmap(list(mat = cm, annotation_col = ann, annotation_colors = annotation_colors(meta, p$variables),
    color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100), main = sprintf("Sample correlation (Pearson, %s)", p$transformation),
    display_numbers = ncol(cm) <= 12, number_format = "%.3f"),
    d, "sample_correlation_heatmap", p)
  out
}
