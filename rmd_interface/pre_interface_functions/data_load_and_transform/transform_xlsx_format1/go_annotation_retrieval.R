#library(dplyr)
#library(httr)
#library(jsonlite)
#library(furrr)
#library(purrr)
#library(tibble)
# ---- Helper functions ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

fetch_go_for_protein <- function(uniprot_id) {
  if (missing(uniprot_id) || !nzchar(uniprot_id)) stop("Provide a valid UniProt ID.")
  
  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
  all_results <- list()
  page <- 1L
  limit <- 200L
  
  repeat {
    params <- list(
      geneProductId = uniprot_id,
      limit = limit,
      page = page,
      includeFields = paste(c("goName","name"), collapse = ",")
    )
    
    resp <- tryCatch(
      GET(base_url, query = params, add_headers(Accept = "application/json"), timeout(60)),
      error = function(e) stop("Request failed: ", e$message)
    )
    
    if (http_error(resp)) {
      stop(sprintf("QuickGO request failed: %s %s", status_code(resp),
                   content(resp, "text", encoding = "UTF-8")))
    }
    
    parsed <- httr::content(resp, as = "parsed", type = "application/json")
    results <- parsed$results
    results <- results[map_lgl(results, ~ !is.null(.x) && length(.x) > 0)]
    if (length(results) > 0) all_results <- append(all_results, results)
    
    total <- parsed$pageInfo$total %||% 0
    per_page <- parsed$pageInfo$resultsPerPage %||% limit
    if (page * per_page >= total) break
    page <- page + 1L
  }
  
  if (length(all_results) == 0) {
    return(tibble(
      UniProtID = uniprot_id,
      goIds = list(character(0)),
      goNames = list(NA_character_),
      ontology = list(NA_character_)
    ))
  }
  
  df <- fromJSON(toJSON(all_results, auto_unbox = TRUE), flatten = TRUE)
  df <- df %>%
    mutate(
      goId = as.character(goId %||% NA_character_),
      goName = as.character(goName %||% NA_character_),
      goAspect = as.character(goAspect %||% NA_character_)
    ) %>%
    select(goId, goName, goAspect) %>%
    distinct(goId, .keep_all = TRUE)
  
  tibble(
    UniProtID = uniprot_id,
    goIds = list(df$goId),
    goNames = list(df$goName),
    ontology = list(df$goAspect)
  )
}

safe_fetch_go <- function(uid) {
  tryCatch(
    fetch_go_for_protein(uid),
    error = function(e) {
      warning(sprintf("Failed for UniProtID %s: %s", uid, e$message))
      tibble(
        UniProtID = uid,
        goIds = list(character(0)),
        goNames = list(NA_character_),
        ontology = list(NA_character_)
      )
    }
  )
}

# ---- Main GO mapping function ----
add_go_annotations <- function(df, uniprot_col = "UniProtID") {
  # Set up parallel processing
  plan(multisession, workers = parallel::detectCores() - 1)
  
  uids <- df[[uniprot_col]]
  
  # Fetch GO annotations in parallel with progress
  results_list <- future_map(uids, safe_fetch_go, .progress = TRUE)
  go_df <- bind_rows(results_list)
  
  # Merge GO annotations back
  df_out <- df %>%
    left_join(go_df, by = setNames("UniProtID", uniprot_col))
  
  return(df_out)
}