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
