create_gcn_data <- function(gpr_list) {
  gcn_data <- lapply(names(gpr_list), function(sample_name) {
    sample_df <- gpr_list[[sample_name]]
    
    # Keep only rows with valid GeneSymbol and non-NA values
    df <- sample_df %>%
      dplyr::filter(!is.na(GeneSymbol) & !is.na(.data[[sample_name]])) %>%
      dplyr::distinct(GeneSymbol, .keep_all = TRUE)
    
    # Extract values
    vals <- as.numeric(df[[sample_name]])
    names(vals) <- df$GeneSymbol
    
    # Identify duplicates and perturb slightly
    duplicated_vals <- vals[duplicated(vals)]
    if (length(duplicated_vals) > 0) {
      for (v in duplicated_vals) {
        idx <- which(vals == v)
        vals[idx] <- vals[idx] + seq_along(idx) * 1e-6
      }
    }
    
    # Sort decreasing
    vals <- sort(vals, decreasing = TRUE)
    return(vals)
  })
  
  # Name the list elements
  names(gcn_data) <- names(gpr_list)
  
  return(gcn_data)
}