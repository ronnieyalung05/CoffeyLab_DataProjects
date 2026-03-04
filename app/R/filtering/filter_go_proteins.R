# R/filter_proteins.R

filter_by_go_crossref <- function(reference_df, gost_df, cols_to_keep) {
  
  # extract all GO IDs from gost results
  gost_go_ids <- base::unique(gost_df$term_id)
  
  if (base::length(gost_go_ids) == 0) {
    base::stop("No GO term IDs found in gost results. Make sure the dataset has a 'term_id' column.")
  }
  
  if (!"goIds" %in% base::colnames(reference_df)) {
    base::stop("Reference dataset does not have a 'goIds' column.")
  }
  
  # for each protein, find which of its GO IDs match the gost results
  reference_df$matched_goIds <- base::lapply(reference_df$goIds, function(protein_goids) {
    base::intersect(protein_goids, gost_go_ids)
  })
  
  # split into matched and unmatched
  has_match     <- base::lengths(reference_df$matched_goIds) > 0
  filtered_data <- reference_df[has_match,  cols_to_keep, drop = FALSE]
  lost_data     <- reference_df[!has_match, cols_to_keep, drop = FALSE]
  
  base::return(base::list(
    filtered_data = filtered_data,
    lost_data     = lost_data,
    kept_count    = base::nrow(filtered_data),
    lost_count    = base::nrow(lost_data)
  ))
}