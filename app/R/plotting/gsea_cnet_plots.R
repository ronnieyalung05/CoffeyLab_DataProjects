# R/plotting/gsea_cnet_plots.R

# builds the named descending gene vector from df + selected samples
build_gene_list <- function(df, selected_samples, gene_col) {
  numeric_cols <- base::intersect(selected_samples,
                                   base::names(df)[base::sapply(df, base::is.numeric)])
  
  if (base::length(numeric_cols) == 0) {
    base::stop("No valid numeric columns found for selected samples.")
  }
  
  if (base::length(numeric_cols) == 1) {
    genes  <- base::as.character(df[[gene_col]])
    values <- df[[numeric_cols[1]]]
  } else {
    # average expression across selected samples per gene
    all_genes <- base::list()
    for (col in numeric_cols) {
      genes  <- base::as.character(df[[gene_col]])
      values <- df[[col]]
      for (j in base::seq_along(genes)) {
        gene <- genes[j]
        if (!gene %in% base::names(all_genes)) all_genes[[gene]] <- c()
        all_genes[[gene]] <- c(all_genes[[gene]], values[j])
      }
    }
    genes  <- base::names(all_genes)
    values <- base::sapply(all_genes, mean)
  }
  
  # remove NAs, sort descending
  valid      <- !base::is.na(genes) & !base::is.na(values) & genes != ""
  genes      <- genes[valid]
  values     <- values[valid]
  gene_list  <- values
  base::names(gene_list) <- genes
  base::sort(gene_list, decreasing = TRUE)
}

run_gsea_cnet <- function(desc_gene_list, selected_ont,
                           selected_plots,
                           pvalue_cutoff  = 0.05,
                           min_gs_size    = 10,
                           max_gs_size    = 500) {
  
  gsea_res <- clusterProfiler::gseGO(
    geneList     = desc_gene_list,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = selected_ont,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    pvalueCutoff = pvalue_cutoff,
    verbose      = FALSE,
    scoreType    = "pos"
  )
  
  if (base::is.null(gsea_res) || base::nrow(gsea_res@result) == 0) {
    base::stop("No significant GSEA results found. Try adjusting pvalue cutoff or gene list.")
  }
  
  gsea_readable <- clusterProfiler::setReadable(
    gsea_res, OrgDb = org.Hs.eg.db, keyType = "SYMBOL"
  )
  
  plots <- base::list()
  
  # p1 — basic with foldChange coloring
  if ("p1" %in% selected_plots) {
    plots$p1 <- enrichplot::cnetplot(
      gsea_readable,
      foldChange = desc_gene_list
    )
  }
  
  # p2 — sized by pvalue
  if ("p2" %in% selected_plots) {
    plots$p2 <- enrichplot::cnetplot(
      gsea_readable,
      foldChange    = desc_gene_list,
      size_category = "pvalue"       # replaces categorySize in v1.28+
    )
  }
  
  # p4 — category labels only
  if ("p4" %in% selected_plots) {
    plots$p4 <- enrichplot::cnetplot(
      gsea_readable,
      node_label = "category"
    )
  }
  
  # p5 — gene labels only
  if ("p5" %in% selected_plots) {
    plots$p5 <- enrichplot::cnetplot(
      gsea_readable,
      node_label = "gene"
    )
  }
  
  # p6 — all labels
  if ("p6" %in% selected_plots) {
    plots$p6 <- enrichplot::cnetplot(
      gsea_readable,
      node_label = "all"
    )
  }
  
  # p7 — no labels, custom colors
  if ("p7" %in% selected_plots) {
    plots$p7 <- enrichplot::cnetplot(
      gsea_readable,
      node_label     = "none",
      color_category = "firebrick",
      color_item     = "steelblue"   # replaces color_gene in v1.28+
    )
  }
  
  base::return(plots)
}