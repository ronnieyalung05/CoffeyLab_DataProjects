create_numeric_df <- function(df, keep_cols = c("GeneSymbol")) {
  # Select numeric columns
  numeric_cols <- base::names(df)[base::sapply(df, base::is.numeric)]
  
  # Create list of dataframes, one per numeric column
  gpr_list <- base::lapply(numeric_cols, function(colname) {
    df %>%
      dplyr::select(dplyr::all_of(c(colname, keep_cols))) %>%
      dplyr::filter(.data[[colname]] != 0) %>%
      dplyr::arrange(dplyr::desc(.data[[colname]]))
  })
  
  # Name each element in the list by the numeric column
  base::names(gpr_list) <- numeric_cols
  base::return(gpr_list)
}
