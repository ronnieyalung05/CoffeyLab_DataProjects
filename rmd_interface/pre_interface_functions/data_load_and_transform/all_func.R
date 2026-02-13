#library(readxl)
#library(tools)
load_dataset <- function(data_list = list(), name = "") {
  file_path <- base::file.choose()
  
  ext <- base::tolower(tools::file_ext(file_path))
  
  df <- base::switch(ext,
                     "xlsx" = readxl::read_excel(file_path),
                     "csv"  = utils::read.csv(file_path),
                     base::stop("Unsupported file type; requires .xlsx or .csv")
  )
  
  if (name == "") {
    name <- tools::file_path_sans_ext(base::basename(file_path))
  }
  
  data_list[[name]] <- df
  base::return(data_list)
}


##########
##########
##########

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


##########
##########
##########

#library(dplyr)
#library(httr)
#library(jsonlite)
#library(furrr)
#library(purrr)
#library(tibble)
# ---- Helper functions ----
`%||%` <- function(a, b) if (!base::is.null(a)) a else b

fetch_go_for_protein <- function(uniprot_id) {
  if (base::missing(uniprot_id) || !base::nzchar(uniprot_id)) base::stop("Provide a valid UniProt ID.")
  
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
  all_results <- base::list()
  page <- 1L
  limit <- 200L
  
  repeat {
    params <- base::list(
      geneProductId = uniprot_id,
      limit = limit,
      page = page,
      includeFields = base::paste(c("goName","name"), collapse = ",")
    )
    
    resp <- base::tryCatch(
      httr::GET(
        base_url,
        query = params,
        httr::add_headers(Accept = "application/json"),
        httr::timeout(60)
      ),
      error = function(e) base::stop("Request failed: ", e$message)
    )
    
    if (httr::http_error(resp)) {
      base::stop(
        base::sprintf(
          "QuickGO request failed: %s %s",
          httr::status_code(resp),
          httr::content(resp, "text", encoding = "UTF-8")
        )
      )
    }
    
    parsed <- httr::content(resp, as = "parsed", type = "application/json")
    results <- parsed$results
    results <- results[purrr::map_lgl(results, ~ !base::is.null(.x) && base::length(.x) > 0)]
    if (base::length(results) > 0) all_results <- base::append(all_results, results)
    
    total <- parsed$pageInfo$total %||% 0
    per_page <- parsed$pageInfo$resultsPerPage %||% limit
    if (page * per_page >= total) break
    page <- page + 1L
  }
  
  if (base::length(all_results) == 0) {
    return(tibble::tibble(
      UniProtID = uniprot_id,
      goIds = base::list(base::character(0)),
      goNames = base::list(NA_character_),
      ontology = base::list(NA_character_)
    ))
  }
  
  df <- jsonlite::fromJSON(
    jsonlite::toJSON(all_results, auto_unbox = TRUE),
    flatten = TRUE
  )
  df <- df %>%
    dplyr::mutate(
      goId = base::as.character(goId %||% NA_character_),
      goName = base::as.character(goName %||% NA_character_),
      goAspect = base::as.character(goAspect %||% NA_character_)
    ) %>%
    dplyr::select(goId, goName, goAspect) %>%
    dplyr::distinct(goId, .keep_all = TRUE)
  
  tibble::tibble(
    UniProtID = uniprot_id,
    goIds = base::list(df$goId),
    goNames = base::list(df$goName),
    ontology = base::list(df$goAspect)
  )
}

safe_fetch_go <- function(uid) {
  base::tryCatch(
    fetch_go_for_protein(uid),
    error = function(e) {
      base::warning(base::sprintf("Failed for UniProtID %s: %s", uid, e$message))
      tibble::tibble(
        UniProtID = uid,
        goIds = base::list(base::character(0)),
        goNames = base::list(NA_character_),
        ontology = base::list(NA_character_)
      )
    }
  )
}

# ---- Main GO mapping function ----
# Requires UniProtIDs
add_go_annotations <- function(df, uniprot_col = "UniProtID") {
  # Set up parallel processing
  future::plan(future::multisession, workers = parallel::detectCores() - 1)
  
  uids <- df[[uniprot_col]]
  
  # Fetch GO annotations in parallel with progress
  results_list <- furrr::future_map(uids, safe_fetch_go, .progress = TRUE)
  go_df <- dplyr::bind_rows(results_list)
  
  # Merge GO annotations back
  df_out <- df %>%
    dplyr::left_join(go_df, by = stats::setNames("UniProtID", uniprot_col))
  
  base::return(df_out)
}


##########
##########
##########

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


##########
##########
##########


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
