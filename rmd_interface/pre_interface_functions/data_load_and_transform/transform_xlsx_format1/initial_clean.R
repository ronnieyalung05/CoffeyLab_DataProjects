#library(dplyr)
#library(stringr)
clean_and_normalize_sample_data <- function(df) {
  
  # 1️⃣ Remove rows where "#" is NA (not tracking these yet)
  df <- df[!is.na(df$`#`), ]
  
  # 2️⃣ Remove unwanted columns
  cols_to_remove <- c("#", "Visible?", "Starred?", "Alternate ID",
                      "Protein Grouping Ambiguity", "Taxonomy", "Molecular Weight")
  df <- df %>% select(-any_of(cols_to_remove))
  
  # 3️⃣ Rename first two columns
  colnames(df)[1:2] <- c("IdentifiedProteins", "AccessionNumber")
  
  # 4️⃣ Identify numeric columns: all columns after AccessionNumber
  numeric_cols <- 3:ncol(df)
  
  # 5️⃣ Normalize numeric columns (log2, 0 stays 0)
  df <- df %>%
    mutate(across(all_of(numeric_cols), ~ if_else(. == 0, 0, log2(.))))
  
  # 6️⃣ Clean AccessionNumber, extract UniProtID & EntryName
  df <- df %>%
    mutate(
      AccessionNumberClean = str_remove(AccessionNumber, "\\s*\\(\\+.*\\)"),
      UniProtID            = str_extract(AccessionNumberClean, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
      EntryName            = str_extract(AccessionNumberClean, "[A-Z0-9_]+(?=\\|?$)")
    )
  
  # 7️⃣ Track omitted rows **after renaming and removing columns**
  # Keep a copy before filtering
  pre_filter_df <- df
  
  # Filter: human proteins & rows with UniProtID
  keep_rows <- str_detect(df$EntryName, "_HUMAN$") & !is.na(df$UniProtID)
  
  # Rows removed go to omitted_rows
  omitted_rows <- pre_filter_df[!keep_rows, ]
  
  # Keep only the filtered rows in df
  df <- df[keep_rows, ]
  
  # Remove AccessionNumber and leave just the Cleaned version
  df <- df %>% select(-AccessionNumber)
  
  # 8️⃣ Return cleaned + omitted
  return(list(
    cleaned = df,
    omitted = omitted_rows
  ))
}