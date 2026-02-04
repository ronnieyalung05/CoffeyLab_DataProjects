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
