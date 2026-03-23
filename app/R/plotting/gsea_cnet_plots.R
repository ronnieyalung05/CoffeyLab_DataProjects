# R/plotting/gsea_cnet_plots.R

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
  
  valid     <- !base::is.na(genes) & !base::is.na(values) & genes != ""
  genes     <- genes[valid]
  values    <- values[valid]
  gene_list <- values
  base::names(gene_list) <- genes
  base::sort(gene_list, decreasing = TRUE)
}

run_gsea_cnet <- function(desc_gene_list, selected_ont,
                           selected_plots,
                           show_category = 5,
                           pvalue_cutoff = 0.05,
                           min_gs_size   = 10,
                           max_gs_size   = 500) {
  
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
  
  n_significant <- base::nrow(gsea_res@result)
  
  # clamp show_category to actual number of significant pathways
  show_category <- base::min(show_category, n_significant)
  
  gsea_readable <- clusterProfiler::setReadable(
    gsea_res, OrgDb = org.Hs.eg.db, keyType = "SYMBOL"
  )
  
  plots <- base::list()
  
  # fold change coloring — basic
  if ("foldchange" %in% selected_plots) {
    plots$foldchange <- enrichplot::cnetplot(
      gsea_readable,
      foldChange   = desc_gene_list,
      showCategory = show_category
    )
  }
  
  # category labels only
  if ("category_labels" %in% selected_plots) {
    plots$category_labels <- enrichplot::cnetplot(
      gsea_readable,
      node_label   = "category",
      showCategory = show_category
    )
  }
  
  # gene labels only
  if ("gene_labels" %in% selected_plots) {
    plots$gene_labels <- enrichplot::cnetplot(
      gsea_readable,
      node_label   = "gene",
      showCategory = show_category
    )
  }
  
  # all labels
  if ("all_labels" %in% selected_plots) {
    plots$all_labels <- enrichplot::cnetplot(
      gsea_readable,
      node_label   = "all",
      showCategory = show_category
    )
  }
  
  # no labels, custom colors
  if ("custom_colors" %in% selected_plots) {
    plots$custom_colors <- enrichplot::cnetplot(
      gsea_readable,
      node_label     = "none",
      color_category = "firebrick",
      color_item     = "steelblue",
      showCategory   = show_category
    )
  }
  
  base::return(base::list(
    plots         = plots,
    n_significant = n_significant,
    show_category = show_category   # return actual clamped value used
  ))
}