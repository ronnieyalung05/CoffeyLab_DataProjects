#library(dplyr)
#library(stringr)
clean_data <- function(df) {
  
  # 1️⃣ Remove rows where "#" is NA (not tracking these yet)
  if ("#" %in% base::colnames(df)) {
    df <- df[!base::is.na(df$`#`), ]
  }
  
  # 2️⃣ Remove unwanted columns. Doesn't break if these don't exist
  cols_to_remove <- c("#", "Visible?", "Starred?", "Alternate ID",
                      "Protein Grouping Ambiguity", "Taxonomy", "Molecular Weight")
  df <- df %>% dplyr::select(-dplyr::any_of(cols_to_remove))
  
  # 3️⃣ Rename first two columns
  base::colnames(df)[1:2] <- c("Identified_Proteins", "Accession_Number_unclean")
  
  # 6️⃣ Clean AccessionNumber, extract UniProtID & EntryName
  df <- df %>%
    dplyr::mutate(
      Accession_Number = stringr::str_remove(Accession_Number_unclean, "\\s*\\(\\+.*\\)"),
      UniProtID            = stringr::str_extract(Accession_Number, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
      EntryName            = stringr::str_extract(Accession_Number, "[A-Z0-9_]+(?=\\|?$)")
    )
  
  # 7️⃣ Track omitted rows **after renaming and removing columns**
  # Keep a copy before filtering
  pre_filter_df <- df
  
  # Filter: human proteins & rows with UniProtID
  keep_rows <- stringr::str_detect(df$EntryName, "_HUMAN$") & !base::is.na(df$UniProtID)
  
  # Rows removed go to omitted_rows
  omitted_rows <- pre_filter_df[!keep_rows, ]
  
  # Keep only the filtered rows in df
  df <- df[keep_rows, ]
  
  # Remove AccessionNumber and leave just the Cleaned version
  df <- df %>% dplyr::select(-Accession_Number_unclean)
  
  # 8️⃣ Return cleaned + omitted
  base::return(base::list(
    cleaned_data = df,
    omitted_data = omitted_rows
  ))
}

normalize_log <- function(df, numeric_cols, base = 2) {
  df <- df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(numeric_cols), 
        ~ dplyr::if_else(. == 0, 0, base::log(., base = base))
      )
    )
  
  base::return(df)
}

