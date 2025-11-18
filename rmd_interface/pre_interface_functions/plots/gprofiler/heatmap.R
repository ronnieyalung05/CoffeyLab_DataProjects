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