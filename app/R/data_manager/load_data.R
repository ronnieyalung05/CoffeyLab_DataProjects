load_dataset <- function(data_list = list(), name = "", file_path = NULL, original_name = NULL) {
  
  if (is.null(file_path)) stop("No file path provided")
  
  name_source <- if (!is.null(original_name)) original_name else file_path
  ext <- base::tolower(tools::file_ext(name_source))
  
  df <- base::switch(ext,
    "xlsx" = readxl::read_excel(file_path),
    "csv"  = utils::read.csv(file_path),
    "rds"  = {
      obj <- base::readRDS(file_path)
      if (is.data.frame(obj)) {
        obj
      } else if (is.list(obj)) {
        return(c(data_list, obj))
      } else {
        base::stop("RDS file does not contain a dataframe or workspace")
      }
    },
    base::stop("Unsupported file type; requires .xlsx, .csv, or .rds")
  )
  
  if (name == "") {
    name <- tools::file_path_sans_ext(base::basename(name_source))
  }
  
  data_list[[name]] <- df
  return(data_list)
}

CACHE_DIR <- "session_cache"

save_data_store <- function(data_list) {
  if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR)
  saveRDS(data_list, file.path(CACHE_DIR, "data_store.rds"))
}

load_data_store <- function() {
  path <- file.path(CACHE_DIR, "data_store.rds")
  if (file.exists(path)) readRDS(path) else list()
}