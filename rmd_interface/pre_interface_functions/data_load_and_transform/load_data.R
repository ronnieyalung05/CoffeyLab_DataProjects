#library(readxl)
#library(tools)
load_dataset <- function(data_list = list(), name = "") {
  file_path <- file.choose()
  
  ext <- tolower(tools::file_ext(file_path))
  
  df <- switch(ext,
               "xlsx" = readxl::read_excel(file_path),
               stop("Unsupported file type")
  )
  
  if (name == "") {
    name <- tools::file_path_sans_ext(basename(file_path))
  }
  
  data_list[[name]] <- df
  return(data_list)
}