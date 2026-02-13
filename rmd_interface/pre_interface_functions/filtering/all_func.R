filter_proteins_by_gost <- function(sample_data_with_go, 
                                    gost_results, 
                                    sample_col = NULL,
                                    cols_to_keep = c("Identified_Proteins", "Accession_Number", "UniProtID", 
                                                     "EntryName", "GeneSymbol", "goIds", "goNames", 
                                                     "ontology", "matched_goIds")) {
  
  sample_data_with_go <- base::as.data.frame(sample_data_with_go)
  
  # 1. Progress Bar
  pb <- progress::progress_bar$new(
    total = base::nrow(sample_data_with_go),
    format = "[:bar] :percent | :current/:total proteins | eta: :eta"
  )
  
  # 2. Extract GO IDs from gost_results
  gost_go_ids <- base::unique(gost_results$result$term_id)
  
  # 3. Setup matched_goIds
  sample_data_with_go$matched_goIds <- base::vector("list", base::nrow(sample_data_with_go))
  lost_proteins <- base::data.frame()
  
  # 4. The Loop
  for (i in base::seq_len(base::nrow(sample_data_with_go))) {
    pb$tick()
    protein_goids <- sample_data_with_go$goIds[[i]]
    matched <- base::intersect(protein_goids, gost_go_ids)
    
    if (base::length(matched) > 0) {
      sample_data_with_go$matched_goIds[[i]] <- matched
    } else {
      lost_proteins <- base::rbind(lost_proteins, sample_data_with_go[i, ])
    }
  }
  
  # 5. Filter the data
  filtered_data <- base::subset(sample_data_with_go, base::lengths(matched_goIds) > 0)
  
  # 6. DYNAMIC COLUMN SELECTION
  # Add the sample column to the list of desired columns if provided
  final_cols <- if (!base::is.null(sample_col)) c(cols_to_keep, sample_col) else cols_to_keep
  
  # Safety Check: Only keep columns that actually exist in the dataframe
  # This prevents the "undefined columns selected" error
  existing_cols <- base::intersect(final_cols, base::colnames(sample_data_with_go))
  
  filtered_data <- filtered_data[, existing_cols, drop = FALSE]
  
  if (base::nrow(lost_proteins) > 0) {
    lost_proteins <- lost_proteins[, existing_cols, drop = FALSE]
  }
  
  return(base::list(
    filtered_data = filtered_data,
    lost_proteins = lost_proteins,
    kept_count = base::nrow(filtered_data),
    lost_count = base::nrow(lost_proteins)
  ))
}

filter_all_sample_by_go <- function(sample_data_with_go, gost_results_list, cores = parallel::detectCores() - 2) {
  # library(parallel)  # assume already loaded
  
  # Helper to run filter_proteins_by_gost and pass sample column name
  run_filter <- function(sample_name, gost_item) {
    filter_proteins_by_gost(
      sample_data_with_go,
      gost_item,
      sample_col = sample_name
    )
  }
  
  # If the list is named, use those names as sample_col
  if (!base::is.null(base::names(gost_results_list))) {
    sample_names <- base::names(gost_results_list)
    
    results <- parallel::mcmapply(
      run_filter,
      sample_name = sample_names,
      gost_item = gost_results_list,
      SIMPLIFY = FALSE,
      mc.cores = cores
    )
    
    base::names(results) <- sample_names
  } else {
    # unnamed list: no sample_col
    results <- parallel::mclapply(
      gost_results_list,
      function(gost_item) filter_proteins_by_gost(sample_data_with_go, gost_item),
      mc.cores = cores
    )
  }
  
  return(results)
}
