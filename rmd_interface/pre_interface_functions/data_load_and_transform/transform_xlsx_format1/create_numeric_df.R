create_numeric_df <- function(df, keep_cols = c("GeneSymbol")) {
  # Select numeric columns
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  
  # Create list of dataframes, one per numeric column
  gpr_list <- lapply(numeric_cols, function(colname) {
    df %>%
      dplyr::select(all_of(c(colname, keep_cols))) %>%
      filter(.data[[colname]] != 0) %>%
      arrange(dplyr::desc(.data[[colname]]))
  })
  
  # Name each element in the list by the numeric column
  names(gpr_list) <- numeric_cols
  return(gpr_list)
}
