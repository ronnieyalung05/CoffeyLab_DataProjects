# R/gost_analysis.R

VALID_GOST_SOURCES <- c("GO:BP", "GO:MF", "GO:CC")

create_numeric_df <- function(df, keep_cols = c("GeneSymbol")) {
  numeric_cols <- base::names(df)[base::sapply(df, base::is.numeric)]
  
  gpr_list <- base::lapply(numeric_cols, function(colname) {
    df %>%
      dplyr::select(dplyr::all_of(c(colname, keep_cols))) %>%
      dplyr::filter(.data[[colname]] != 0) %>%
      dplyr::arrange(dplyr::desc(.data[[colname]]))
  })
  
  base::names(gpr_list) <- numeric_cols
  base::return(gpr_list)
}

safe_gost <- function(query, sources, significant, max_attempts) {
  for (i in seq_len(max_attempts)) {
    result <- base::tryCatch(
      gprofiler2::gost(
        query          = query,
        organism       = "hsapiens",
        ordered_query  = TRUE,
        significant    = significant,
        sources        = sources
      ),
      error = function(e) {
        # connection-level errors
        if (grepl("Could not resolve host|Failed to connect|timeout|Couldn't connect",
                  e$message, ignore.case = TRUE)) {
          base::stop("CONNECTION_ERROR: Cannot reach g:Profiler API. Check your internet connection.")
        }
        NULL
      }
    )
    if (!base::is.null(result) && !base::inherits(result, "try-error")) return(result)
    if (i < max_attempts) base::Sys.sleep(3)
  }
  base::stop(paste0("g:Profiler request failed after ", max_attempts, " attempts."))
}

run_gost_single <- function(numeric_df_list, sample_name, sources,
                             top_n = 0, significant = TRUE, max_attempts = 5) {
  if (!sample_name %in% base::names(numeric_df_list)) {
    base::stop(paste0("'", sample_name, "' not found in numeric df list."))
  }
  
  df <- numeric_df_list[[sample_name]]
  
  if (top_n > 0 && base::nrow(df) > top_n) {
    df <- df %>%
      dplyr::arrange(dplyr::desc(.data[[sample_name]])) %>%
      dplyr::slice_head(n = top_n)
  }
  
  genes <- df$GeneSymbol
  genes <- genes[!base::is.na(genes) & genes != ""]
  
  if (base::length(genes) == 0) base::stop("No valid Gene Symbols found after filtering.")
  
  result <- safe_gost(genes, sources, significant, max_attempts)
  base::return(result)
}

run_gost_multi <- function(numeric_df_list, source, top_n = 0,
                            significant = TRUE, max_attempts = 5) {
  results_list <- base::list()
  
  for (sample_name in base::names(numeric_df_list)) {
    df <- numeric_df_list[[sample_name]]
    
    if (top_n > 0 && base::nrow(df) > top_n) {
      df <- df %>%
        dplyr::arrange(dplyr::desc(.data[[sample_name]])) %>%
        dplyr::slice_head(n = top_n)
    }
    
    genes <- df$GeneSymbol
    genes <- genes[!base::is.na(genes) & genes != ""]
    
    result <- base::tryCatch(
      safe_gost(genes, source, significant, max_attempts),
      error = function(e) {
        base::warning(paste0("Failed for sample '", sample_name, "': ", e$message))
        NULL
      }
    )
    
    results_list[[sample_name]] <- result
  }
  
  base::return(results_list)
}

# flattens a gost result into a saveable dataframe
gost_to_df <- function(result, sample_name = "result") {
  if (base::is.null(result)) {
    return(base::data.frame(sample = sample_name, note = "No result returned"))
  }
  
  df <- result$result
  if (base::is.null(df) || base::nrow(df) == 0) {
    return(base::data.frame(sample = sample_name, note = "No significant terms found"))
  }
  
  # parents column is a list-col — flatten to semicolon string for storage
  if ("parents" %in% base::names(df)) {
    df$parents <- base::sapply(df$parents, function(x) base::paste(x, collapse = ";"))
  }
  
  df$sample <- sample_name
  base::return(df)
}