# Optional biological annotation (GO / Reactome over-representation). Runs ONLY when enabled in
# config and the required packages are installed. Core DESeq2 results never depend on this module.

run_downstream <- function(summary, p) {
  ann <- p$annotation
  if (is.null(ann) || !isTRUE(ann$enabled)) {
    msg("SKIPPED", "optional functional annotation disabled (config: annotation.enabled)")
    return(list(status = "disabled"))
  }
  orgdb <- p$orgdb
  needed <- c("clusterProfiler", orgdb, if (isTRUE(ann$reactome)) "ReactomePA")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (is.null(orgdb) || length(missing)) {
    msg("WARNING", "annotation skipped; missing R packages: ", paste(c(missing, if (is.null(orgdb)) "OrgDb"), collapse = ", "))
    return(list(status = "skipped", missing = missing))
  }
  reactome_org <- c(org.Hs.eg.db = "human", org.Mm.eg.db = "mouse")[orgdb]
  out <- list(status = "done", results = list())
  for (cn in names(summary$contrasts)) {
    cdir <- summary$contrasts[[cn]]$dir
    full <- read.delim(file.path(cdir, "full_results.tsv"), stringsAsFactors = FALSE)
    ids <- sub("\\.\\d+$", "", full$gene_id)  # strip Ensembl version suffix
    keytype <- if (mean(grepl("^ENS[A-Z]*G\\d+$", ids)) > 0.5) "ENSEMBL" else "SYMBOL"
    universe <- unique(ids[!is.na(full$padj)])
    to_entrez <- function(x) {
      m <- tryCatch(clusterProfiler::bitr(x, fromType = keytype, toType = "ENTREZID", OrgDb = orgdb),
                    error = function(e) NULL)
      if (is.null(m)) character(0) else unique(m$ENTREZID)
    }
    for (dirn in c("up", "down")) {
      genes <- unique(ids[full$regulation == dirn])
      if (length(genes) < 5) { msg("INFO", cn, " ", dirn, ": <5 genes, enrichment skipped"); next }
      if (isTRUE(ann$go)) {
        ego <- tryCatch(clusterProfiler::enrichGO(genes, OrgDb = orgdb, keyType = keytype, ont = "BP",
                                                  universe = universe, pvalueCutoff = p$alpha),
                        error = function(e) { msg("WARNING", "GO failed: ", conditionMessage(e)); NULL })
        if (!is.null(ego)) {
          f <- file.path(cdir, sprintf("GO_BP_%s.tsv", dirn))
          write_tsv_df(as.data.frame(ego), f)
          out$results[[paste(cn, dirn, "GO", sep = "_")]] <- f
        }
      }
      if (isTRUE(ann$reactome)) {
        if (is.na(reactome_org)) { msg("WARNING", "Reactome supports human/mouse only; skipped"); next }
        er <- tryCatch(ReactomePA::enrichPathway(to_entrez(genes), organism = reactome_org,
                                                 universe = to_entrez(universe), pvalueCutoff = p$alpha),
                       error = function(e) { msg("WARNING", "Reactome failed: ", conditionMessage(e)); NULL })
        if (!is.null(er)) {
          f <- file.path(cdir, sprintf("Reactome_%s.tsv", dirn))
          write_tsv_df(as.data.frame(er), f)
          out$results[[paste(cn, dirn, "Reactome", sep = "_")]] <- f
        }
      }
    }
  }
  out
}
