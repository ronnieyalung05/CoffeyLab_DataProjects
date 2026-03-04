# R/plot_store.R

PLOT_CACHE_DIR <- "session_cache/plots"

save_plot_store <- function(plot_list) {
  tryCatch({
    if (!dir.exists(PLOT_CACHE_DIR)) dir.create(PLOT_CACHE_DIR, recursive = TRUE)
    saveRDS(plot_list, file.path(PLOT_CACHE_DIR, "plot_store.rds"))
  }, error = function(e) {
    warning(paste("Could not save plot store:", e$message))
  })
}

load_plot_store <- function() {
  tryCatch({
    path <- file.path(PLOT_CACHE_DIR, "plot_store.rds")
    if (file.exists(path)) readRDS(path) else list()
  }, error = function(e) {
    warning(paste("Could not load plot store:", e$message))
    list()
  })
}