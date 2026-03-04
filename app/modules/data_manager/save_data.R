# modules/save_data.R

dataSaveUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("Save a single dataset"),
    selectInput(ns("save_single_choice"), "Choose dataset", choices = NULL),
    selectInput(ns("single_format"), "File format",
                choices = c("RDS (.rds)" = "rds", "CSV (.csv)" = "csv")),
    
    # warning only shown for csv
    conditionalPanel(
      condition = paste0("input['", ns("single_format"), "'] == 'csv'"),
      div(style = "color: #e67e22; font-size: 0.85em; padding: 4px 8px;",
          "⚠ CSV cannot store list-type columns (e.g. GO term columns). Those will be dropped on export.")
    ),
    
    downloadButton(ns("download_single"), "Download dataset"),
    
    hr(),
    
    h4("Save entire workspace"),
    p("Saves all currently loaded datasets into one file you can reload later."),
    p(style = "color: #888; font-size: 0.85em;",
      "Note: workspace export is RDS only — CSV cannot store multiple datasets in one file."),
    downloadButton(ns("download_all"), "Download full workspace as .rds")
  )
}

dataSaveServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    observe({
      updateSelectInput(session, "save_single_choice", choices = names(data_store()))
    })
    
    output$download_single <- downloadHandler(
      
      filename = function() {
        ext <- if (input$single_format == "csv") ".csv" else ".rds"
        paste0(input$save_single_choice, ext)
      },
      
      content = function(file) {
        df <- data_store()[[input$save_single_choice]]
        
        if (input$single_format == "csv") {
          # drop list-type columns — they can't be written to csv
          list_cols <- sapply(df, is.list)
          if (any(list_cols)) {
            warning(paste("Dropping list columns:", paste(names(df)[list_cols], collapse = ", ")))
            df <- df[, !list_cols, drop = FALSE]
          }
          utils::write.csv(df, file, row.names = FALSE)
        } else {
          saveRDS(df, file)
        }
      }
    )
    
    output$download_all <- downloadHandler(
      filename = function() paste0("workspace_", Sys.Date(), ".rds"),
      content  = function(file) saveRDS(data_store(), file)
    )
    
  })
}