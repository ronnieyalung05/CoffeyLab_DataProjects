# modules/rename_data.R

renameDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Rename Dataset"),
    selectInput(ns("rename_source"), "Select dataset to rename", choices = NULL),
    textInput(ns("rename_to"), "New name", placeholder = "e.g. my_data_final"),
    actionButton(ns("run_rename"), "Rename",
                 style = "color: white; background-color: #8e44ad;"),
    uiOutput(ns("rename_status"))
  )
}

renameDataServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    observe({
      updateSelectInput(session, "rename_source", choices = names(data_store()))
    })
    
    observeEvent(input$run_rename, {
      req(input$rename_source, input$rename_to)
      tryCatch({
        new_name <- trimws(input$rename_to)
        
        if (new_name == "") stop("New name cannot be blank")
        if (new_name %in% names(data_store())) stop("A dataset with that name already exists")
        
        updated <- data_store()
        updated[[new_name]]            <- updated[[input$rename_source]]
        updated[[input$rename_source]] <- NULL
        data_store(updated)
        save_data_store(updated)
        
        output$rename_status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Renamed '", input$rename_source, "' → '", new_name, "'"))
        })
      }, error = function(e) {
        output$rename_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })
    
  })
}