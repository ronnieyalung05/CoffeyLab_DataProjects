#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(dplyr)
convert_uniprot_ids <- function(df) {
  # Extract UniProt IDs, remove NA
  uni_ids <- df$UniProtID
  uni_ids <- uni_ids[!is.na(uni_ids)]
  
  if (length(uni_ids) == 0) {
    warning("No valid UniProt IDs found.")
    df$GeneSymbol <- NA
    df$EntrezID   <- NA
    return(df)
  }
  
  # Map UniProt -> Gene Symbol + EntrezID
  mapping <- bitr(
    uni_ids,
    fromType = "UNIPROT",
    toType   = c("SYMBOL", "ENTREZID"),
    OrgDb    = org.Hs.eg.db
  ) %>%
    group_by(UNIPROT) %>%
    slice(1) %>%  # Take first mapping only
    ungroup()
  
  # Merge back into original df
  df_out <- df %>%
    left_join(mapping, by = c("UniProtID" = "UNIPROT")) %>%
    rename(
      GeneSymbol = SYMBOL,
      EntrezID   = ENTREZID
    )
  
  return(df_out)
}

split_unmapped_ids <- function(data_list, cleaned_name = "cleaned",
                                   gene_col = "GeneSymbol", entrez_col = "EntrezID") {
  # Extract the cleaned data frame
  df <- data_list[[cleaned_name]]
  
  # Identify unmapped rows
  unmapped <- df %>%
    filter(is.na(.data[[gene_col]]) | is.na(.data[[entrez_col]]))
  
  # Keep only valid rows
  valid <- df %>%
    filter(!is.na(.data[[gene_col]]) & !is.na(.data[[entrez_col]]))
  
  # Update the data list in place
  data_list[[cleaned_name]] <- valid
  data_list$unmapped <- unmapped
  
  return(data_list)
}