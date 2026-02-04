#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(dplyr)
convert_uniprot_ids <- function(df, uniprot_col = "UniProtID") {
  # Check that the column exists
  if (!uniprot_col %in% base::colnames(df)) {
    base::stop(base::paste("Column", uniprot_col, "not found in the data frame."))
  }
  
  # Extract UniProt IDs, remove NAs
  uni_ids <- df[[uniprot_col]]
  uni_ids <- uni_ids[!base::is.na(uni_ids)]
  
  # Stop if no valid IDs
  if (base::length(uni_ids) == 0) {
    base::stop("No valid UniProt IDs found in column '", uniprot_col, "'. Execution halted.")
  }
  
  if (base::length(uni_ids) == 0) {
    base::warning("No valid UniProt IDs found.")
    df$GeneSymbol <- NA
    df$EntrezID   <- NA
    return(df)
  }
  
  # Map UniProt -> Gene Symbol + EntrezID
  mapping <- clusterProfiler::bitr(
    uni_ids,
    fromType = "UNIPROT",
    toType   = c("SYMBOL", "ENTREZID"),
    OrgDb    = org.Hs.eg.db
  ) %>%
    dplyr::group_by(UNIPROT) %>%
    dplyr::slice_head(n = 1) %>%  # Take first mapping only
    dplyr::ungroup()
  
  # Merge back into original df
  df_out <- df %>%
    dplyr::left_join(mapping, by = c("UniProtID" = "UNIPROT")) %>%
    dplyr::rename(
      GeneSymbol = SYMBOL,
      EntrezID   = ENTREZID
    )
  
  base::return(df_out)
}

split_unmapped_ids <- function(data_list, cleaned_name = "cleaned",
                               gene_col = "GeneSymbol", entrez_col = "EntrezID") {
  # Extract the cleaned data frame
  df <- data_list[[cleaned_name]]
  
  # Identify unmapped rows
  unmapped <- df %>%
    dplyr::filter(dplyr::if_else(base::is.na(.data[[gene_col]]) | base::is.na(.data[[entrez_col]]), TRUE, FALSE))
  
  # Keep only valid rows
  valid <- df %>%
    dplyr::filter(dplyr::if_else(!base::is.na(.data[[gene_col]]) & !base::is.na(.data[[entrez_col]]), TRUE, FALSE))
  
  # Update the data list in place
  data_list[[cleaned_name]] <- valid
  data_list$unmapped <- unmapped
  
  base::return(data_list)
}
