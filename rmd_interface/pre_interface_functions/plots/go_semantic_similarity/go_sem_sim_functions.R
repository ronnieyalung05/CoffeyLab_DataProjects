run_sample_semantic_analysis <- function(sample_list, ontologies = c("BP", "MF", "CC")) {
  # --- Step 1: Select sample ---
  cat("Available samples:\n")
  sample_names <- names(sample_list)
  for (i in seq_along(sample_names)) cat(i, ":", sample_names[i], "\n")
  
  sample_index <- as.integer(readline(prompt = "Enter the number of the sample you want to analyze: "))
  if (is.na(sample_index) || sample_index < 1 || sample_index > length(sample_names)) stop("Invalid selection.")
  
  selected_sample <- sample_list[[sample_index]]
  
  # --- Step 2: Use GeneSymbol as identifier ---
  gene_ids <- as.character(selected_sample$GeneSymbol)
  
  # --- Step 3: Select ontology ---
  cat("Available ontologies:\n")
  for (i in seq_along(ontologies)) cat(i, ":", ontologies[i], "\n")
  
  ont_index <- as.integer(readline(prompt = "Enter the number of the ontology to use (BP/MF/CC): "))
  if (is.na(ont_index) || ont_index < 1 || ont_index > length(ontologies)) stop("Invalid selection.")
  selected_ont <- ontologies[ont_index]
  
  # --- Step 4: Top N genes ---
  numeric_cols <- colnames(selected_sample)[sapply(selected_sample, is.numeric)]
  if (length(numeric_cols) == 0) stop("No numeric columns found for expression values.")
  
  cat("Available expression columns:\n")
  for (i in seq_along(numeric_cols)) cat(i, ":", numeric_cols[i], "\n")
  expr_index <- as.integer(readline(prompt = "Select numeric column for expression values: "))
  expr_col <- numeric_cols[expr_index]
  expr_vector <- selected_sample[[expr_col]]
  
  top_n <- as.integer(readline(prompt = "Enter the number of top expressed genes to select: "))
  top_genes <- gene_ids[order(expr_vector, decreasing = TRUE)][1:top_n]
  
  # --- Step 5: Filter valid IDs ---
  valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  top_genes_valid <- intersect(top_genes, valid_symbols)
  omitted <- setdiff(top_genes, top_genes_valid)
  if (length(omitted) > 0) {
    cat(length(omitted), "gene(s) were not found in org.Hs.eg.db and will be skipped:\n")
    print(omitted)
  }
  
  if (length(top_genes_valid) < 2) stop("Too few valid genes for semantic similarity. Increase top N.")
  
  # --- Step 6: GO semantic data ---
  hsGO <- godata("org.Hs.eg.db", ont = selected_ont, keytype = "SYMBOL")
  
  # --- Step 7: Compute semantic similarity ---
  sim <- mgeneSim(top_genes_valid, semData = hsGO, measure = "Wang", drop = "IEA", combine = "BMA")
  
  # --- Step 8: Heatmap ---
  pheatmap(sim, 
           clustering_method = "ward.D", 
           color = colorRampPalette(c("white", "red"))(50),
           display_numbers = TRUE,
           number_format = "%.2f",
           main = paste("Semantic similarity heatmap:", sample_names[sample_index]))
  
  # --- Step 9: Dendrogram ---
  hc <- hclust(as.dist(1 - sim), method = "ward.D")
  p <- ggtree(hc) + geom_tiplab() + xlim(NA, 0.5)
  print(p)
  
  # --- Step 10: Return similarity matrix ---
  return(sim)
}