dataLoadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "Choose CSV, XLSX, or RDS",
              accept = c(".csv", ".xlsx", ".rds")),
    textInput(ns("name"), "Dataset name (optional)", value = ""),
    hr(),
    uiOutput(ns("status")),      # success/error message goes here
    tableOutput(ns("preview"))
  )
}

dataLoadServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    status_message <- reactiveVal(NULL)   # tracks the last status
    
    observeEvent(input$file, {
      req(input$file)
      
      tryCatch({
        updated_list <- load_dataset(
          data_list     = data_store(),
          name          = input$name,
          file_path     = input$file$datapath,
          original_name = input$file$name
        )
        
        data_store(updated_list)
        save_data_store(updated_list)
        
        loaded_name <- tail(names(updated_list), 1)
        df <- updated_list[[loaded_name]]
        
        status_message(list(
          type = "success",
          text = paste0("✓ Loaded '", loaded_name, "' — ", nrow(df), " rows × ", ncol(df), " cols")
        ))
        
      }, error = function(e) {
        status_message(list(
          type = "error",
          text = paste0("✗ Failed to load: ", e$message)
        ))
      })
    })
    
    output$status <- renderUI({
      req(status_message())
      msg <- status_message()
      
      if (msg$type == "success") {
        div(style = "color: green; font-weight: bold; padding: 8px;", msg$text)
      } else {
        div(style = "color: red; font-weight: bold; padding: 8px;", msg$text)
      }
    })
    
  })
}