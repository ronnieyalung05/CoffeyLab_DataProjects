run_sample_semantic_analysis <- function(sample_list, ontologies = c("BP", "MF", "CC")) {
  # --- Step 1: Select sample(s) ---
  cat("Available samples:\n")
  sample_names <- base::names(sample_list)
  for (i in seq_along(sample_names)) cat(i, ":", sample_names[i], "\n")
  
  sample_input <- readline(prompt = "Enter the number(s) of the sample(s) you want to analyze (comma-separated for multiple, e.g., 1,2,3): ")
  sample_indices <- as.integer(strsplit(sample_input, ",")[[1]])
  
  if (any(is.na(sample_indices)) || any(sample_indices < 1) || any(sample_indices > length(sample_names))) {
    stop("Invalid selection.")
  }
  
  # --- Step 2: Combine samples if multiple selected ---
  if (length(sample_indices) == 1) {
    # Single sample: use as-is
    selected_sample <- sample_list[[sample_indices]]
    selected_sample_names <- sample_names[sample_indices]
    
    # Automatically detect the numeric column
    numeric_cols <- base::colnames(selected_sample)[base::sapply(selected_sample, base::is.numeric)]
    if (length(numeric_cols) == 0) stop("No numeric columns found for expression values.")
    expr_col <- numeric_cols[1]
    cat("Using expression column:", expr_col, "\n")
    
  } else {
    # Multiple samples: combine them
    cat("\nCombining", length(sample_indices), "samples...\n")
    
    # Collect all gene-value pairs from selected samples
    all_genes <- list()
    for (idx in sample_indices) {
      current_sample <- sample_list[[idx]]
      genes <- base::as.character(current_sample$GeneSymbol)
      
      # Automatically detect the numeric column
      numeric_cols <- base::colnames(current_sample)[base::sapply(current_sample, base::is.numeric)]
      if (length(numeric_cols) == 0) stop("No numeric columns found in sample ", sample_names[idx])
      values <- current_sample[[numeric_cols[1]]]
      
      for (j in seq_along(genes)) {
        gene <- genes[j]
        if (!gene %in% names(all_genes)) {
          all_genes[[gene]] <- c()
        }
        all_genes[[gene]] <- c(all_genes[[gene]], values[j])
      }
    }
    
    # Average values for genes appearing in multiple samples
    combined_genes <- names(all_genes)
    combined_values <- sapply(all_genes, mean)
    
    # Create combined dataframe with generic column name
    selected_sample <- data.frame(
      GeneSymbol = combined_genes,
      expr_value = combined_values,
      stringsAsFactors = FALSE
    )
    expr_col <- "expr_value"
    
    selected_sample_names <- paste(sample_names[sample_indices], collapse = " + ")
    cat("Combined", length(combined_genes), "unique genes\n\n")
  }
  
  # --- Step 2b: Use GeneSymbol as identifier ---
  gene_ids <- base::as.character(selected_sample$GeneSymbol)
  
  # --- Step 3: Select ontology ---
  cat("Available ontologies:\n")
  for (i in seq_along(ontologies)) cat(i, ":", ontologies[i], "\n")
  
  ont_index <- as.integer(readline(prompt = "Enter the number of the ontology to use (BP/MF/CC): "))
  if (is.na(ont_index) || ont_index < 1 || ont_index > length(ontologies)) stop("Invalid selection.")
  selected_ont <- ontologies[ont_index]
  
  # --- Step 4: Top N genes ---
  expr_vector <- selected_sample[[expr_col]]
  
  top_n <- as.integer(readline(prompt = "Enter the number of top expressed genes to select: "))
  top_genes <- gene_ids[base::order(expr_vector, decreasing = TRUE)][1:top_n]
  
  # --- Step 5: Filter valid IDs ---
  valid_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  top_genes_valid <- base::intersect(top_genes, valid_symbols)
  omitted <- base::setdiff(top_genes, top_genes_valid)
  if (length(omitted) > 0) {
    cat(length(omitted), "gene(s) were not found in org.Hs.eg.db and will be skipped:\n")
    print(omitted)
  }
  
  if (length(top_genes_valid) < 2) stop("Too few valid genes for semantic similarity. Increase top N.")
  
  # --- Step 6: GO semantic data ---
  hsGO <- GOSemSim::godata("org.Hs.eg.db", ont = selected_ont, keytype = "SYMBOL")
  
  # --- Step 7: Compute semantic similarity ---
  sim <- GOSemSim::mgeneSim(
    top_genes_valid, semData = hsGO, 
    measure = "Wang", drop = "IEA", combine = "BMA"
  )
  
  # --- Step 8: Heatmap ---
  pheatmap::pheatmap(
    sim, 
    clustering_method = "ward.D", 
    color = grDevices::colorRampPalette(c("white", "red"))(50),
    display_numbers = TRUE,
    number_format = "%.2f",
    main = paste("Semantic similarity heatmap:", selected_sample_names)
  )
  
  # --- Step 9: Dendrogram ---
  hc <- stats::hclust(stats::as.dist(1 - sim), method = "ward.D")
  p <- ggtree::ggtree(hc) + ggtree::geom_tiplab() + ggplot2::xlim(NA, 0.5)
  print(p)
  
  # --- Step 10: Return similarity matrix ---
  return(sim)
}