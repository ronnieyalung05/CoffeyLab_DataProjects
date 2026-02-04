# -----------------------------------------
# Single Sample, Multiple Sources with Top N
# -----------------------------------------
run_gost_single_sample <- function(g_protein_list, top_n = 0) {
  base::cat("Available samples:\n")
  base::print(base::names(g_protein_list))
  
  protein_input <- base::readline(prompt = "\nEnter ONE sample name from the list above: ")
  if (!(protein_input %in% base::names(g_protein_list))) {
    base::stop(base::paste0("Error: '", protein_input, "' is not a valid sample name.\n",
                            "Valid options: ", base::paste(base::names(g_protein_list), collapse = ", ")))
  }
  
  # Sources
  base::cat("\nAvailable sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, CORUM, TF, MIRNA, HPA, HP\n")
  source_input <- base::readline(prompt = "Enter one or more sources (comma-separated): ")
  sources <- base::trimws(base::strsplit(source_input, ",")[[1]])
  
  valid_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP", 
                     "CORUM", "TF", "MIRNA", "HPA", "HP")
  if (!all(sources %in% valid_sources)) {
    invalid <- sources[!sources %in% valid_sources]
    base::stop(base::paste0("Invalid source(s): ", base::paste(invalid, collapse = ", "),
                            "\nValid options: ", base::paste(valid_sources, collapse = ", ")))
  }
  
  # Take top N genes based on sample column
  gpr_filter_df <- g_protein_list[[protein_input]]
  if (top_n > 0 && base::nrow(gpr_filter_df) > top_n) {
    gpr_filter_df <- gpr_filter_df %>%
      dplyr::arrange(dplyr::desc(.data[[protein_input]])) %>%
      dplyr::slice_head(n = top_n)
  }
  
  genes_to_run <- gpr_filter_df$GeneSymbol
  genes_to_run <- genes_to_run[!is.na(genes_to_run) & genes_to_run != ""]
  
  if (length(genes_to_run) == 0) {
    base::stop("No valid Gene Symbols found for this sample after filtering.")
  }
  
  safe_gost <- function(query, sources, max_attempts = 5) {
    for (i in 1:max_attempts) {
      result <- base::try(
        gost(query = query, organism = "hsapiens", ordered_query = TRUE,
             significant = TRUE, sources = sources),
        silent = FALSE
      )
      if (!base::inherits(result, "try-error") && !base::is.null(result)) return(result)
      base::cat("Attempt", i, "failed — retrying in 2 seconds...\n")
      base::Sys.sleep(3)
    }
    base::stop("g:Profiler request failed after multiple attempts")
  }
  
  result <- safe_gost(genes_to_run, sources)
  base::cat("\n✅ g:Profiler analysis complete for sample", protein_input, "\n")
  return(result)
}

# -----------------------------------------
# Multi-Sample, One Source with Top N per sample
# -----------------------------------------
run_gost_multi_sample_one_source <- function(g_protein_list, top_n = 0) {
  base::cat("\nAvailable sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, CORUM, TF, MIRNA, HPA, HP\n")
  source_input <- base::readline(prompt = "Enter ONE source from the list above: ")
  source_input <- base::trimws(source_input)
  
  valid_sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP",
                     "CORUM", "TF", "MIRNA", "HPA", "HP")
  if (!(source_input %in% valid_sources)) {
    base::stop(base::paste0("Invalid source: ", source_input,
                            "\nValid options: ", base::paste(valid_sources, collapse = ", ")))
  }
  
  results_list <- base::list()
  
  for (sample_name in base::names(g_protein_list)) {
    base::cat("Processing sample:", sample_name, "...\n")
    df <- g_protein_list[[sample_name]]
    
    # Top N based on that sample's column
    if (top_n > 0 && base::nrow(df) > top_n) {
      df <- df %>%
        dplyr::arrange(dplyr::desc(.data[[sample_name]])) %>%
        dplyr::slice_head(n = top_n)
    }
    
    genes_to_run <- df$GeneSymbol
    
    safe_gost <- function(query, source, max_attempts = 5) {
      for (i in 1:max_attempts) {
        result <- base::try(
          gost(query = query, organism = "hsapiens", ordered_query = TRUE,
               significant = TRUE, sources = source),
          silent = FALSE
        )
        if (!base::inherits(result, "try-error") && !base::is.null(result)) return(result)
        base::cat("Attempt", i, "failed — retrying in 2 seconds...\n")
        base::Sys.sleep(3)
      }
      base::warning(base::paste0("g:Profiler request failed for sample ", sample_name))
      return(NULL)
    }
    
    results_list[[sample_name]] <- safe_gost(genes_to_run, source_input)
  }
  
  base::cat("\n✅ Multi-sample g:Profiler analysis complete.\n")
  return(results_list)
}

# -----------------------------------------
# Menu for running g:Profiler with Top N option
# -----------------------------------------
run_gost_menu <- function(g_protein_list) {
  options <- c(
    "Single sample, multiple sources",
    "All samples, ONE source, results separate"
  )
  
  base::cat("\nSelect a g:Profiler analysis option:\n")
  for (i in base::seq_along(options)) base::cat(base::sprintf("[%d] %s\n", i, options[i]))
  
  repeat {
    choice <- base::readline(prompt = "\nEnter option number: ")
    if (!base::grepl("^[0-9]+$", choice)) { base::cat("Invalid input. Please enter a number.\n"); next }
    choice <- base::as.integer(choice)
    if (choice < 1 || choice > base::length(options)) { 
      base::cat("Invalid option. Select a number between 1 and", base::length(options), "\n"); next
    }
    break
  }
  
  # Prompt for Top N
  repeat {
    top_n_input <- base::readline(prompt = "Enter Top N Genes (by expression values) to include (0 for all): ")
    if (!base::grepl("^[0-9]+$", top_n_input)) { base::cat("Invalid input. Enter a number.\n"); next }
    top_n <- base::as.integer(top_n_input)
    break
  }
  
  result <- NULL
  
  if (choice == 1) {
    # Single sample, multiple sources
    result <- base::tryCatch({
      run_gost_single_sample(g_protein_list, top_n = top_n)
    }, error = function(e) {
      base::cat("\n❌ Error in single-sample g:Profiler:\n", e$message, "\n")
      return(NULL)
    })
  } else if (choice == 2) {
    # Multi-sample, one source
    result <- base::tryCatch({
      run_gost_multi_sample_one_source(g_protein_list, top_n = top_n)
    }, error = function(e) {
      base::cat("\n❌ Error in multi-sample g:Profiler:\n", e$message, "\n")
      return(NULL)
    })
  }
  
  return(result)
}
