# R/plotting/gosemsim_plots.R

run_semantic_analysis <- function(df, selected_samples, gene_col,
                                   selected_ont, top_n) {
  
  # --- build sample list from df (mirrors create_numeric_df logic) ---
  numeric_cols <- base::names(df)[base::sapply(df, base::is.numeric)]
  numeric_cols <- base::intersect(selected_samples, numeric_cols)
  
  if (base::length(numeric_cols) == 0) {
    base::stop("No valid numeric columns found for selected samples.")
  }
  
  if (base::length(numeric_cols) == 1) {
    # single sample
    sample_df  <- df[, c(gene_col, numeric_cols[1]), drop = FALSE]
    genes      <- base::as.character(sample_df[[gene_col]])
    values     <- sample_df[[numeric_cols[1]]]
    label      <- numeric_cols[1]
    
  } else {
    # multiple samples — average expression across samples per gene
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
    label  <- base::paste(numeric_cols, collapse = " + ")
  }
  
  # --- top N genes by expression ---
  top_genes <- genes[base::order(values, decreasing = TRUE)][1:top_n]
  
  # --- filter to valid org.Hs.eg.db symbols ---
  valid_symbols   <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  top_genes_valid <- base::intersect(top_genes, valid_symbols)
  omitted         <- base::setdiff(top_genes, top_genes_valid)
  
  if (base::length(omitted) > 0) {
    base::warning(paste0(base::length(omitted),
                         " gene(s) not found in org.Hs.eg.db and skipped: ",
                         base::paste(omitted, collapse = ", ")))
  }
  
  if (base::length(top_genes_valid) < 2) {
    base::stop("Too few valid genes for semantic similarity. Increase Top N.")
  }
  
  # --- compute semantic similarity ---
  hsGO <- GOSemSim::godata("org.Hs.eg.db",
                            ont     = selected_ont,
                            keytype = "SYMBOL")
  
  sim <- GOSemSim::mgeneSim(
    top_genes_valid,
    semData = hsGO,
    measure = "Wang",
    drop    = "IEA",
    combine = "BMA"
  )
  
  # --- heatmap via pheatmap, converted to ggplot ---
  heatmap_plot <- ggplotify::as.ggplot(function() {
    pheatmap::pheatmap(
      sim,
      clustering_method = "ward.D",
      color             = grDevices::colorRampPalette(c("white", "red"))(50),
      display_numbers   = TRUE,
      number_format     = "%.2f",
      main              = paste0("Semantic Similarity Heatmap: ", label,
                                 " (", selected_ont, ")")
    )
  })
  
  # --- dendrogram via ggtree ---
  hc <- stats::hclust(stats::as.dist(1 - sim), method = "ward.D")
  dendrogram_plot <- ggtree::ggtree(hc) +
    ggtree::geom_tiplab() +
    ggplot2::xlim(NA, 0.5) +
    ggplot2::labs(title = paste0("Gene Cluster Dendrogram: ", label,
                                  " (", selected_ont, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  # --- similarity matrix as dataframe for data_store ---
  sim_df <- base::as.data.frame(sim)
  
  base::return(base::list(
    heatmap_plot    = heatmap_plot,
    dendrogram_plot = dendrogram_plot,
    sim_matrix_df   = sim_df,
    label           = label
  ))
}