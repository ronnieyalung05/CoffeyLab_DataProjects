# modules/clean_data.R

cleanDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("Clean Dataset"),
    selectInput(ns("clean_source"), "Select dataset to clean (only for data formatted as shown in example dataset)", choices = NULL),
    actionButton(ns("run_clean"), "Run Cleaning",
                 style = "color: white; background-color: #2980b9;"),
    uiOutput(ns("clean_status")),
    
    hr(),
    
    h4("Log Normalize"),
    selectInput(ns("norm_source"), "Select dataset to normalize", choices = NULL),
    uiOutput(ns("col_selector")),
    selectInput(ns("log_base"), "Log base",
                choices = c("2" = 2, "10" = 10, "e" = exp(1)), selected = 2),
    textInput(ns("norm_output_name"), "Save normalized data as (optional)",
              placeholder = "leave blank for auto name (__[original-name]__log[base])"),
    actionButton(ns("run_norm"), "Run Normalization",
                 style = "color: white; background-color: #2980b9;"),
    uiOutput(ns("norm_status"))
  )
}

cleanDataServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    observe({
      nms <- names(data_store())
      updateSelectInput(session, "clean_source", choices = nms)
      updateSelectInput(session, "norm_source",  choices = nms)
    })
    
    output$col_selector <- renderUI({
      req(input$norm_source, data_store()[[input$norm_source]])
      df       <- data_store()[[input$norm_source]]
      num_cols <- names(df)[sapply(df, is.numeric)]
      checkboxGroupInput(session$ns("norm_cols"), "Columns to normalize",
                         choices = num_cols, selected = num_cols)
    })
    
    observeEvent(input$run_clean, {
      req(input$clean_source)
      tryCatch({
        df       <- data_store()[[input$clean_source]]
        result   <- clean_data(df)
        base_name  <- input$clean_source
        clean_name <- paste0("__", base_name, "__cleaned")
        omit_name  <- paste0("__", base_name, "__omitted")
        
        updated <- data_store()
        updated[[clean_name]] <- result$cleaned_data
        updated[[omit_name]]  <- result$omitted_data
        data_store(updated)
        save_data_store(updated)
        
        output$clean_status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ '", base_name, "' cleaned (saved as __[name]__cleaned and __[name]__omitted→ ",
                     nrow(result$cleaned_data), " rows kept, ",
                     nrow(result$omitted_data), " rows omitted"))
        })
      }, error = function(e) {
        output$clean_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })
    
    observeEvent(input$run_norm, {
      req(input$norm_source, input$norm_cols)
      tryCatch({
        df     <- data_store()[[input$norm_source]]
        result <- normalize_log(df, input$norm_cols, base = as.numeric(input$log_base))
        
        norm_name <- if (input$norm_output_name == "") {
          paste0("__", input$norm_source, "__log", input$log_base)
        } else {
          input$norm_output_name
        }
        
        updated <- data_store()
        updated[[norm_name]] <- result
        data_store(updated)
        save_data_store(updated)
        
        output$norm_status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Normalized '", input$norm_source, "' → '", norm_name, "'"))
        })
      }, error = function(e) {
        output$norm_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })
    
  })
}