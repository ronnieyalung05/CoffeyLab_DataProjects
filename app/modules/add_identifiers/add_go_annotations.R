# modules/go_annotations.R

goAnnotationsUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("Fetch GO Annotations"),
    p("Queries the EBI QuickGO API for each protein in parallel. One request per row."),
    
    selectInput(ns("source_df"),   "Select dataset",    choices = NULL),
    selectInput(ns("uniprot_col"), "UniProt ID column", choices = NULL),
    
    hr(),
    
    h4("Parallelism"),
    numericInput(ns("n_workers"), "Number of parallel workers",
                 value = max(1, parallel::detectCores() - 1),
                 min = 1, max = parallel::detectCores(), step = 1),
    p(style = "color: #888; font-size: 0.85em;",
      paste0("Your machine has ", parallel::detectCores(), " cores. ",
             "Recommended: leave 1 free for the Shiny session.")),
    
    hr(),
    
    h4("Output name"),
    textInput(ns("output_name"), "Save annotated dataset as",
              placeholder = "defaults to original name"),
    
    actionButton(ns("run"), "Fetch GO Annotations",
                 style = "color: white; background-color: #2980b9;"),
    
    br(), br(),
    uiOutput(ns("row_count_warning")),
    
    # spinner wraps the status output — spins while result is NULL
    shinycssloaders::withSpinner(
      uiOutput(ns("status")),
      type    = 8,        # circle spinner style
      color   = "#2980b9",
      size    = 0.7,
      hide.ui = FALSE     # keep layout stable while spinning
    )
  )
}

goAnnotationsServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    # NULL = idle, TRUE = running, character = done/error message
    running <- reactiveVal(FALSE)
    
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })
    
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      cols        <- names(data_store()[[input$source_df]])
      default_col <- if ("UniProtID" %in% cols) "UniProtID" else cols[1]
      updateSelectInput(session, "uniprot_col", choices = cols, selected = default_col)
    })
    
    output$row_count_warning <- renderUI({
      req(input$source_df, data_store()[[input$source_df]])
      n        <- nrow(data_store()[[input$source_df]])
      workers  <- max(1, input$n_workers)
      est_low  <- round((n * 0.5) / workers / 60, 1)
      est_high <- round((n * 2)   / workers / 60, 1)
      
      if (n > 200) {
        div(style = "color: #e67e22; font-weight: bold; padding: 8px;",
            paste0("⚠ ", n, " rows across ", workers, " workers — estimated ",
                   est_low, "–", est_high, " minutes."))
      } else {
        div(style = "color: #888; padding: 8px;",
            paste0(n, " rows across ", workers, " workers — estimated ",
                   round(est_low * 60), "–", round(est_high * 60), " seconds."))
      }
    })
    
    # status only renders when NOT running — this is what drives the spinner
    # when running() is TRUE, status returns NULL, which triggers the spinner
    output$status <- renderUI({
      if (running()) return(NULL)   # spinner shows while this is NULL
      
      div(style = "color: #888; padding: 8px;",
          "Ready — select a dataset and click Fetch GO Annotations.")
    })
    
    observeEvent(input$run, {
      req(input$source_df, input$uniprot_col)
      
      df      <- data_store()[[input$source_df]]
      n       <- nrow(df)
      workers <- max(1, input$n_workers)
      
      running(TRUE)   # triggers spinner
      
      tryCatch({
        result <- add_go_annotations(
          df          = df,
          uniprot_col = input$uniprot_col,
          workers     = workers
        )
        
        out_name <- if (trimws(input$output_name) == "") {
          input$source_df
        } else {
          trimws(input$output_name)
        }
        
        updated <- data_store()
        updated[[out_name]] <- result
        data_store(updated)
        save_data_store(updated)
        
        n_with_go <- sum(sapply(result$goIds, function(x) length(x) > 0))
        
        running(FALSE)  # stops spinner
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Done — '", out_name, "': ",
                     n_with_go, " of ", n, " proteins had GO terms"))
        })
        
      }, error = function(e) {
        running(FALSE)  # always stop spinner on error too
        
        is_connection_error <- grepl("CONNECTION_ERROR", e$message)
        clean_msg <- gsub("CONNECTION_ERROR: ", "", e$message)
        
        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              if (is_connection_error) {
                paste0("🌐 Connection error: ", clean_msg)
              } else {
                paste0("✗ Error: ", clean_msg)
              }
          )
        })
      })
    })
    
  })
}