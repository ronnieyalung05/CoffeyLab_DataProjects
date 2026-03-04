# modules/data_load.R (dataViewer section)

dataViewerUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("selected_df"), "Select a loaded dataset", choices = NULL),
    verbatimTextOutput(ns("df_info")),
    
    # save filtered view
    fluidRow(
      column(4,
        textInput(ns("filtered_name"), "Save current view as",
                  placeholder = "name for filtered dataset")
      ),
      column(2,
        br(),
        actionButton(ns("save_filtered"), "Save Current View",
                     style = "color: white; background-color: #27ae60;")
      ),
      column(6,
        br(),
        uiOutput(ns("save_status"))
      )
    ),
    
    hr(),
    DT::dataTableOutput(ns("df_table"))
  )
}

dataViewerServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    observe({
      updateSelectInput(session, "selected_df", choices = names(data_store()))
    })
    
    output$df_info <- renderPrint({
      req(input$selected_df)
      df <- data_store()[[input$selected_df]]
      cat("Rows:", nrow(df), "\n")
      cat("Cols:", ncol(df), "\n")
      cat("Column names:", paste(names(df), collapse = ", "), "\n")
    })
    
    output$df_table <- DT::renderDataTable({
      req(input$selected_df)
      data_store()[[input$selected_df]]
    }, options = list(pageLength = 15, scrollX = TRUE))
    
    observeEvent(input$save_filtered, {
      req(input$selected_df)
      
      tryCatch({
        # df_table_rows_all = all rows passing current filter, in display order
        # df_table_rows_current = only the current page
        # we want all filtered rows, not just the current page
        filtered_rows <- input$df_table_rows_all
        
        if (base::is.null(filtered_rows) || base::length(filtered_rows) == 0) {
          stop("No rows in current view — adjust your filter and try again.")
        }
        
        full_df       <- data_store()[[input$selected_df]]
        filtered_df   <- full_df[filtered_rows, ]
        
        save_name <- if (trimws(input$filtered_name) == "") {
          paste0("__", input$selected_df, "__filtered")
        } else {
          trimws(input$filtered_name)
        }
        
        if (save_name %in% names(data_store())) {
          stop(paste0("'", save_name, "' already exists — choose a different name."))
        }
        
        updated <- data_store()
        updated[[save_name]] <- filtered_df
        data_store(updated)
        save_data_store(updated)
        
        output$save_status <- renderUI({
          div(style = "color: green; font-weight: bold;",
              paste0("✓ Saved '", save_name, "' (",
                     nrow(filtered_df), " of ", nrow(full_df), " rows)"))
        })
        
      }, error = function(e) {
        output$save_status <- renderUI({
          div(style = "color: red; font-weight: bold;", paste0("✗ ", e$message))
        })
      })
    })
    
  })
}