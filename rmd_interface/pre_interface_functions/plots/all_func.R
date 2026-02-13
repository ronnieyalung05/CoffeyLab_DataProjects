
# -----------------------------------------
# Menu for viewing g:Profiler (gost) plots
# -----------------------------------------

view_gost_plots_menu <- function(gost_results) {
  if (is.null(gost_results)) stop("❌ No g:Profiler results provided.")
  
  # First choice: single sample vs multi-sample
  cat("\nSelect type of g:Profiler results to view:\n")
  cat("[1] Single sample (one gost_results object)\n")
  cat("[2] Multi-sample (list of gost_results, one per sample)\n")
  
  repeat {
    choice1 <- readline(prompt = "\nEnter option number: ")
    if (!grepl("^[1-2]$", choice1)) {
      cat("Invalid input. Please enter 1 or 2.\n")
    } else break
  }
  choice1 <- as.integer(choice1)
  
  if (choice1 == 1) {
    # ----- SINGLE SAMPLE -----
    cat("\n📊 Displaying single-sample g:Profiler results...\n")
    print(gostplot(gost_results, capped = TRUE, interactive = TRUE))
    return(invisible())
  }
  
  # ----- MULTI-SAMPLE -----
  cat("\nSelect multi-sample display mode:\n")
  cat("[1] Interactive plots (no titles; separate interactive viewer windows)\n")
  cat("[2] Non-interactive plots (static, titled with sample names)\n")
  
  repeat {
    choice2 <- readline(prompt = "\nEnter option number: ")
    if (!grepl("^[1-2]$", choice2)) {
      cat("Invalid input. Please enter 1 or 2.\n")
    } else break
  }
  choice2 <- as.integer(choice2)
  
  if (choice2 == 1) {
    cat("\n📊 Displaying interactive plots for all samples...\n")
    invisible(lapply(names(gost_results), function(sample_name) {
      sample_result <- gost_results[[sample_name]]
      if (!is.null(sample_result) && !is.null(sample_result$result)) {
        cat("\n==============================\n")
        cat("📊 Sample:", sample_name, "\n")
        cat("==============================\n")
        print(gostplot(sample_result, capped = TRUE, interactive = TRUE))
      } else {
        cat("⚠️ No valid g:Profiler results for", sample_name, "\n")
      }
    }))
  } else {
    cat("\n📊 Displaying static plots with sample titles...\n")
    invisible(lapply(names(gost_results), function(sample_name) {
      sample_result <- gost_results[[sample_name]]
      if (!is.null(sample_result) && !is.null(sample_result$result)) {
        p <- gostplot(sample_result, capped = TRUE, interactive = FALSE) +
          ggtitle(paste("g:Profiler Results —", sample_name))
        print(p)
      } else {
        cat("⚠️ No valid g:Profiler results for", sample_name, "\n")
      }
    }))
  }
}

###########################
###########################
###########################
# -----------------------------------------
# 1️⃣ Single sample, multiple sources
# -----------------------------------------
plot_heatmap_single_sample <- function(gost_results, top_n = 0, per_group = FALSE) {
  if (is.null(gost_results$result) || nrow(gost_results$result) == 0) {
    stop("No valid g:Profiler results found.")
  }
  
  df <- gost_results$result %>%
    mutate(neg_log_p = -log10(p_value))
  
  if (top_n > 0) {
    if (per_group) {
      # Top N per source
      sources <- unique(df$source)
      for (s in sources) {
        sub_df <- df %>% filter(source == s) %>%
          arrange(p_value) %>%
          slice_head(n = top_n)
        
        p <- ggplot(sub_df, aes(x = source, y = reorder(term_name, neg_log_p), fill = neg_log_p)) +
          geom_tile() +
          #geom_text(aes(label = round(neg_log_p, 2)), size = 3) +  # <-- add this line
          scale_fill_gradient(low = "white", high = "red") +
          labs(title = paste0("Top ", top_n, " Terms — Source: ", s),
               x = "Source", y = "Term Name", fill = "-log10(p-value)") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8),
                axis.text.x = element_text(angle = 70, hjust = 1),
                plot.title = element_text(hjust = 0.5))
        
        print(p)
      }
      return(invisible(NULL))
    } else {
      # Top N overall
      df <- df %>% arrange(p_value) %>% slice_head(n = top_n)
    }
  }
  
  # Single heatmap (all sources together)
  p <- ggplot(df, aes(x = source, y = reorder(term_name, neg_log_p), fill = neg_log_p)) +
    geom_tile() +
    #geom_text(aes(label = round(neg_log_p, 2)), size = 3) +  # <-- add this line
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = "Enriched Terms - Single Sample",
         x = "Source", y = "Term Name", fill = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 70, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  print(p)
  invisible(NULL)
}

# -----------------------------------------
# 2️⃣ Multi-sample, ONE source
# -----------------------------------------
plot_heatmap_multi_sample_one_source <- function(gost_results_list, top_n = 0, per_group = FALSE) {
  # Find first valid source name
  first_valid <- gost_results_list[[which(!sapply(gost_results_list, is.null))[1]]]
  source_name <- first_valid$result$source[1]
  
  # Combine results into one df
  combined_df <- do.call(rbind, lapply(names(gost_results_list), function(sample_name) {
    res <- gost_results_list[[sample_name]]$result
    if (!is.null(res)) {
      res$sample <- sample_name
      return(res)
    }
    NULL
  }))
  
  combined_df <- combined_df %>% mutate(neg_log_p = -log10(p_value))
  
  if (top_n > 0) {
    if (per_group) {
      # Top N per sample
      samples <- unique(combined_df$sample)
      for (s in samples) {
        sub_df <- combined_df %>% filter(sample == s) %>%
          arrange(p_value) %>%
          slice_head(n = top_n)
        
        p <- ggplot(sub_df, aes(x = sample, y = reorder(term_name, neg_log_p), fill = neg_log_p)) +
          geom_tile() +
          #geom_text(aes(label = round(neg_log_p, 2)), size = 3) +  # <-- add this line
          scale_fill_gradient(low = "white", high = "red") +
          labs(title = paste0("Top ", top_n, " Terms — Sample: ", s, " (", source_name, ")"),
               x = "Sample", y = "Term Name", fill = "-log10(p-value)") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 70, hjust = 1),
                axis.text.y = element_text(size = 7),
                plot.title = element_text(hjust = 0.5, face = "bold"))
        
        print(p)
      }
      return(invisible(NULL))
    } else {
      # Top N overall
      combined_df <- combined_df %>% arrange(p_value) %>% slice_head(n = top_n)
    }
  }
  
  # Single heatmap for all samples
  p <- ggplot(combined_df, aes(x = sample, y = reorder(term_name, neg_log_p), fill = neg_log_p)) +
    geom_tile() +
    #geom_text(aes(label = round(neg_log_p, 2)), size = 3) +  # <-- add this line
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = paste0("Top Enriched Terms — ONE source (", source_name, ")"),
         x = "Sample", y = "Term Name", fill = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = 9),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  invisible(NULL)
}

# -----------------------------------------
# 3️⃣ Multi-sample, ONE source, Top N per sample COMBINED
# -----------------------------------------
plot_heatmap_multi_sample_one_source_combined_topn <- function(gost_results_list, top_n = 0) {
  # Combine results into one df
  combined_df <- do.call(rbind, lapply(names(gost_results_list), function(sample_name) {
    res <- gost_results_list[[sample_name]]$result
    if (!is.null(res)) {
      res$sample <- sample_name
      return(res)
    }
    NULL
  }))
  
  if (nrow(combined_df) == 0) stop("No valid results found in any sample.")
  
  combined_df <- combined_df %>% mutate(neg_log_p = -log10(p_value))
  
  # Top N per sample, then union terms
  top_terms <- do.call(c, lapply(unique(combined_df$sample), function(s) {
    combined_df %>% filter(sample == s) %>%
      arrange(p_value) %>% slice_head(n = top_n) %>% pull(term_name)
  }))
  top_terms <- unique(top_terms)
  
  # Filter combined_df to include only these top terms
  df_plot <- combined_df %>% filter(term_name %in% top_terms)
  
  # Plot single heatmap
  p <- ggplot(df_plot, aes(x = sample, y = reorder(term_name, neg_log_p), fill = neg_log_p)) +
    geom_tile() +
    #geom_text(aes(label = round(neg_log_p, 2)), size = 3) +  # <-- add this line
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = paste0("Top ", top_n, " Terms per Sample — Combined Heatmap"),
         x = "Sample", y = "Term Name", fill = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = 9),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  invisible(NULL)
}

# -----------------------------------------
# Menu function
# -----------------------------------------
plot_heatmap_gost_menu <- function(gost_results) {
  options <- c(
    "Single sample, multiple sources",
    "All samples, ONE source, results separate/combined"
  )
  
  cat("\nSelect a heatmap to view by entering the number:\n")
  for (i in seq_along(options)) cat(sprintf("[%d] %s\n", i, options[i]))
  
  # --- Choose main option ---
  repeat {
    choice <- readline(prompt = "\nEnter option number: ")
    if (!grepl("^[0-9]+$", choice)) { cat("Invalid input. Please enter a number.\n"); next }
    choice <- as.integer(choice)
    if (choice < 1 || choice > length(options)) { cat("Invalid option. Try again.\n"); next }
    break
  }
  
  # --- Ask Top N ---
  top_n <- as.integer(readline(prompt = "Enter Top N terms to display (0 for all): "))
  if (is.na(top_n) || top_n < 0) top_n <- 0
  
  # --- Ask per-group ---
  per_group <- FALSE
  display_type <- "separate"
  if (choice == 2) {
    per_group_input <- readline(prompt = "Filter Top N per group? (y/n): ")
    per_group <- tolower(per_group_input) == "y"
    
    # If per-group, ask combined vs separate
    if (per_group) {
      cat("\nDisplay options for Top N per sample:\n[1] Separate plots per sample\n[2] Combined heatmap\n")
      repeat {
        disp_choice <- readline(prompt = "Enter option number: ")
        if (!grepl("^[0-9]+$", disp_choice)) { cat("Invalid input. Please enter a number.\n"); next }
        disp_choice <- as.integer(disp_choice)
        if (!(disp_choice %in% 1:2)) { cat("Invalid option. Try again.\n"); next }
        display_type <- ifelse(disp_choice == 1, "separate", "combined")
        break
      }
    }
  }
  
  # --- Map choice to function and execute ---
  tryCatch({
    if (choice == 1) {
      # Single sample, multiple sources
      plot_heatmap_single_sample(gost_results, top_n = top_n, per_group = per_group)
      
    } else if (choice == 2) {
      # Multi-sample, ONE source
      if (per_group && display_type == "combined") {
        # Combined Top N per sample
        plot_heatmap_multi_sample_one_source_combined_topn(gost_results, top_n = top_n)
      } else {
        # Separate plots or Top N overall
        plot_heatmap_multi_sample_one_source(gost_results, top_n = top_n, per_group = per_group)
      }
    }
  }, error = function(e) {
    cat("❌ Error generating heatmap:\n", e$message, "\n")
    return(NULL)
  })
}

###########################
###########################
###########################

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

###########################
###########################
###########################
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


###########################
###########################
###########################