# R/convert_ids.R

# All valid bitr identifier types for humans
VALID_ID_TYPES <- c("UNIPROT", "SYMBOL", "ENTREZID")

convert_ids <- function(df, from_col, from_type, to_types) {
  
  if (!from_col %in% base::colnames(df)) {
    base::stop(paste("Column '", from_col, "' not found in the data frame."))
  }
  
  ids <- df[[from_col]]
  ids <- ids[!base::is.na(ids)]
  
  if (base::length(ids) == 0) {
    base::stop("No valid IDs found in column '", from_col, "'.")
  }
  
  mapping <- clusterProfiler::bitr(
    ids,
    fromType = from_type,
    toType   = to_types,
    OrgDb    = org.Hs.eg.db
  ) %>%
    dplyr::group_by(.data[[from_type]]) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  df_out <- df %>%
    dplyr::left_join(mapping, by = stats::setNames(from_type, from_col))
  
  base::return(df_out)
}

split_unmapped_ids <- function(df, check_cols) {
  # check_cols: character vector of columns that must be non-NA to be "mapped"
  
  is_unmapped <- base::Reduce(`|`, base::lapply(check_cols, function(col) base::is.na(df[[col]])))
  
  base::return(base::list(
    valid    = df[!is_unmapped, ],
    unmapped = df[is_unmapped, ]
  ))
}