# -----------------------------------------
# Single Sample, Multiple Sources with Top N
# -----------------------------------------
run_gost_single_sample <- function(g_protein_list, top_n = 0) {
  cat("Available samples:\n")
  print(names(g_protein_list))
  
  protein_input <- readline(prompt = "\nEnter ONE sample name from the list above: ")
  if (!(protein_input %in% names(g_protein_list))) {
    stop(paste0("Error: '", protein_input, "' is not a valid sample name.\n",
                "Valid options: ", paste(names(g_protein_list), collapse = ", ")))
  }
  
  # Sources
  cat("\nAvailable sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, CORUM, TF, MIRNA, HPA, HP\n")
  source_input <- readline(prompt = "Enter one or more sources (comma-separated): ")
  sources <- trimws(strsplit(source_input, ",")[[1]])
  
  valid_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", 
                     "CORUM", "TF", "MIRNA", "HPA", "HP")
  if (!all(sources %in% valid_sources)) {
    invalid <- sources[!sources %in% valid_sources]
    stop(paste0("Invalid source(s): ", paste(invalid, collapse = ", "),
                "\nValid options: ", paste(valid_sources, collapse = ", ")))
  }
  
  # Take top N genes based on sample column
  gpr_filter_df <- g_protein_list[[protein_input]]
  if (top_n > 0 && nrow(gpr_filter_df) > top_n) {
    gpr_filter_df <- gpr_filter_df %>% arrange(desc(.data[[protein_input]])) %>% slice_head(n = top_n)
  }
  
  genes_to_run <- gpr_filter_df$GeneSymbol
  
  safe_gost <- function(query, sources, max_attempts = 10) {
    for (i in 1:max_attempts) {
      result <- try(
        gost(query = query, organism = "hsapiens", ordered_query = TRUE,
             significant = TRUE, sources = sources),
        silent = TRUE
      )
      if (!inherits(result, "try-error") && !is.null(result)) return(result)
      cat("Attempt", i, "failed — retrying in 2 seconds...\n")
      Sys.sleep(2)
    }
    stop("g:Profiler request failed after multiple attempts")
  }
  
  result <- safe_gost(genes_to_run, sources)
  cat("\n✅ g:Profiler analysis complete for sample", protein_input, "\n")
  return(result)
}

# -----------------------------------------
# Multi-Sample, One Source with Top N per sample
# -----------------------------------------
run_gost_multi_sample_one_source <- function(g_protein_list, top_n = 0) {
  cat("\nAvailable sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, CORUM, TF, MIRNA, HPA, HP\n")
  source_input <- readline(prompt = "Enter ONE source from the list above: ")
  source_input <- trimws(source_input)
  
  valid_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP",
                     "CORUM", "TF", "MIRNA", "HPA", "HP")
  if (!(source_input %in% valid_sources)) {
    stop(paste0("Invalid source: ", source_input,
                "\nValid options: ", paste(valid_sources, collapse = ", ")))
  }
  
  results_list <- list()
  
  for (sample_name in names(g_protein_list)) {
    cat("Processing sample:", sample_name, "...\n")
    df <- g_protein_list[[sample_name]]
    
    # Top N based on that sample's column
    if (top_n > 0 && nrow(df) > top_n) {
      df <- df %>% arrange(desc(.data[[sample_name]])) %>% slice_head(n = top_n)
    }
    
    genes_to_run <- df$GeneSymbol
    
    safe_gost <- function(query, source, max_attempts = 10) {
      for (i in 1:max_attempts) {
        result <- try(
          gost(query = query, organism = "hsapiens", ordered_query = TRUE,
               significant = TRUE, sources = source),
          silent = TRUE
        )
        if (!inherits(result, "try-error") && !is.null(result)) return(result)
        cat("Attempt", i, "failed — retrying in 2 seconds...\n")
        Sys.sleep(2)
      }
      warning(paste0("g:Profiler request failed for sample ", sample_name))
      return(NULL)
    }
    
    results_list[[sample_name]] <- safe_gost(genes_to_run, source_input)
  }
  
  cat("\n✅ Multi-sample g:Profiler analysis complete.\n")
  return(results_list)
}

# -----------------------------------------
# Menu for running g:Profiler with Top N option
# -----------------------------------------
run_gost_menu_with_top_n <- function(g_protein_list) {
  options <- c(
    "Single sample, multiple sources",
    "All samples, ONE source, results separate"
  )
  
  cat("\nSelect a g:Profiler analysis option:\n")
  for (i in seq_along(options)) cat(sprintf("[%d] %s\n", i, options[i]))
  
  repeat {
    choice <- readline(prompt = "\nEnter option number: ")
    if (!grepl("^[0-9]+$", choice)) { cat("Invalid input. Please enter a number.\n"); next }
    choice <- as.integer(choice)
    if (choice < 1 || choice > length(options)) { 
      cat("Invalid option. Select a number between 1 and", length(options), "\n"); next
    }
    break
  }
  
  # Prompt for Top N
  repeat {
    top_n_input <- readline(prompt = "Enter Top N Genes (by expression values) to include (0 for all): ")
    if (!grepl("^[0-9]+$", top_n_input)) { cat("Invalid input. Enter a number.\n"); next }
    top_n <- as.integer(top_n_input)
    break
  }
  
  result <- NULL
  
  if (choice == 1) {
    # Single sample, multiple sources
    result <- tryCatch({
      run_gost_single_sample(g_protein_list, top_n = top_n)
    }, error = function(e) {
      cat("\n❌ Error in single-sample g:Profiler:\n", e$message, "\n")
      return(NULL)
    })
  } else if (choice == 2) {
    # Multi-sample, one source
    result <- tryCatch({
      run_gost_multi_sample_one_source(g_protein_list, top_n = top_n)
    }, error = function(e) {
      cat("\n❌ Error in multi-sample g:Profiler:\n", e$message, "\n")
      return(NULL)
    })
  }
  
  return(result)
}