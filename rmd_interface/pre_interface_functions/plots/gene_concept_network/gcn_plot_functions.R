run_gsea_cnet <- function(desc_gene_list, ontologies = c("BP", "MF", "CC")) {
  # --- Step 1: Select ontology ---
  cat("Available ontologies:\n")
  for (i in seq_along(ontologies)) cat(i, ":", ontologies[i], "\n")
  
  ont_index <- as.integer(readline(prompt = "Enter the number of the ontology to use: "))
  if (is.na(ont_index) || ont_index < 1 || ont_index > length(ontologies)) {
    stop("Invalid ontology selection.")
  }
  selected_ont <- ontologies[ont_index]
  cat("Selected ontology:", selected_ont, "\n")
  
  # --- Step 2: Run GSEA ---
  gsea_res <- gseGO(
    geneList = desc_gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = selected_ont,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    scoreType = "pos"
  )
  
  # Make gene IDs readable
  gsea_readable <- setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
  
  # --- Step 3: Generate cnetplots ---
  plots <- list()
  
  plots$p1 <- cnetplot(gsea_readable, color.params = list(foldChange = desc_gene_list))
  plots$p2 <- cnetplot(gsea_readable, categorySize = "pvalue", color.params = list(foldChange = desc_gene_list))
  plots$p3 <- cnetplot(gsea_readable, circular = TRUE, colorEdge = TRUE, color.params = list(foldChange = desc_gene_list))
  
  plots$p4 <- cnetplot(gsea_readable, node_label = "category", cex_label_category = 1.2)
  plots$p5 <- cnetplot(gsea_readable, node_label = "gene", cex_label_gene = 0.8)
  plots$p6 <- cnetplot(gsea_readable, node_label = "all")
  plots$p7 <- cnetplot(gsea_readable, node_label = "none", color_category = 'firebrick', color_gene = 'steelblue')
  
  # --- Step 4: Combine plots ---
  plots$grid1 <- cowplot::plot_grid(plots$p1, plots$p2, plots$p3, ncol = 3, labels = LETTERS[1:3], rel_widths = c(0.8, 0.8, 1.2))
  plots$grid2 <- cowplot::plot_grid(plots$p4, plots$p5, plots$p6, plots$p7, ncol = 2, labels = LETTERS[1:4])
  
  return(plots)
}