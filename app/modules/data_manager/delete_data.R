dataDeleteUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Remove a loaded dataset (refresh browser after running)"),
    selectInput(ns("delete_choice"), "Select dataset to remove", choices = NULL),
    actionButton(ns("delete_single"), "Remove Selected", 
                 style = "color: white; background-color: #e74c3c;"),
    hr(),
    h4("Clear everything (refresh browser after running)"),
    p("Removes all loaded datasets from the app and cache."),
    actionButton(ns("delete_all"), "Clear All Datasets",
                 style = "color: white; background-color: #c0392b;"),
    hr(),
    uiOutput(ns("delete_status"))
  )
}

dataDeleteServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    # keep dropdown in sync with what's loaded
    observe({
      updateSelectInput(session, "delete_choice", choices = names(data_store()))
    })
    
    # remove single dataset
    observeEvent(input$delete_single, {
      req(input$delete_choice)
      
      updated_list <- data_store()
      updated_list[[input$delete_choice]] <- NULL
      
      data_store(updated_list)
      save_data_store(updated_list)
      
      output$delete_status <- renderUI({
        div(style = "color: green; font-weight: bold; padding: 8px;",
            paste0("✓ Removed '", input$delete_choice, "'"))
      })
    })
    
    # clear everything
    observeEvent(input$delete_all, {
      data_store(list())
      save_data_store(list())
      
      output$delete_status <- renderUI({
        div(style = "color: green; font-weight: bold; padding: 8px;",
            "✓ All datasets cleared")
      })
    })
  })
}