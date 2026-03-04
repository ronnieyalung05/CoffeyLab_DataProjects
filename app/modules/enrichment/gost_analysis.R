# modules/gost_analysis.R

gostAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("g:Profiler Enrichment Analysis"),
    
    # --- Step 1: dataset + gene symbol col ---
    selectInput(ns("source_df"), "Dataset", choices = NULL),
    selectInput(ns("gene_col"),  "Gene Symbol column", choices = NULL),
    
    hr(),
    
    # --- Step 2: mode ---
    radioButtons(ns("mode"), "Analysis mode",
                 choices = c(
                   "Single sample, multiple sources" = "single",
                   "Multiple samples, one source"    = "multi"
                 ),
                 selected = "single",
                 inline   = TRUE),
    
    hr(),
    
    # --- Step 3a: single mode ---
    conditionalPanel(
      condition = paste0("input['", ns("mode"), "'] == 'single'"),
      selectInput(ns("single_sample"), "Sample (numeric column)", choices = NULL),
      checkboxGroupInput(ns("single_sources"), "Sources",
                         choices  = VALID_GOST_SOURCES,
                         selected = c("GO:BP", "GO:MF", "GO:CC"),
                         inline   = TRUE)
    ),
    
    # --- Step 3b: multi mode ---
    conditionalPanel(
      condition = paste0("input['", ns("mode"), "'] == 'multi'"),
      uiOutput(ns("multi_sample_selector")),
      selectInput(ns("multi_source"), "Source",
                  choices  = VALID_GOST_SOURCES,
                  selected = "GO:BP")
    ),
    
    hr(),
    
    # --- Shared options ---
    fluidRow(
      column(3,
        numericInput(ns("top_n"), "Top N genes by expression values (0 = all)",
                     value = 0, min = 0, step = 10)
      ),
      column(3,
        checkboxInput(ns("significant"), "Significant terms only (p < 0.05)", value = TRUE)
      ),
      column(3,
        numericInput(ns("max_attempts"), "Max retry attempts (adding more increases chances of success, but also increases the chance of not catching connectivity errors)",
                     value = 5, min = 1, max = 20, step = 1)
      )
    ),
    
    hr(),
    
    # --- Output name + run ---
    textInput(ns("output_name"), "Save results as",
              placeholder = "leave blank for auto name"),
    
    actionButton(ns("run_gost"), "Run g:Profiler",
                 style = "color: white; background-color: #2980b9; width: 100%;"),
    
    br(), br(),
    
    shinycssloaders::withSpinner(
      uiOutput(ns("gost_status")),
      type    = 8,
      color   = "#2980b9",
      size    = 0.7,
      hide.ui = FALSE
    )
  )
}

gostAnalysisServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    running <- reactiveVal(FALSE)
    
    # sync dataset dropdown
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })
    
    # when dataset changes, update column dropdowns
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      df       <- data_store()[[input$source_df]]
      all_cols <- names(df)
      num_cols <- all_cols[sapply(df, is.numeric)]
      
      # gene symbol col — prefer GeneSymbol if it exists
      gene_default <- if ("GeneSymbol" %in% all_cols) "GeneSymbol" else all_cols[1]
      updateSelectInput(session, "gene_col", choices = all_cols, selected = gene_default)
      
      # single sample dropdown — only numeric cols
      updateSelectInput(session, "single_sample", choices = num_cols)
    })
    
    # multi mode: checkbox list of numeric cols so user picks which to include
    output$multi_sample_selector <- renderUI({
      req(input$source_df, data_store()[[input$source_df]])
      df       <- data_store()[[input$source_df]]
      num_cols <- names(df)[sapply(df, is.numeric)]
      
      checkboxGroupInput(session$ns("multi_samples"),
                         "Samples to include (numeric columns)",
                         choices  = num_cols,
                         selected = num_cols,
                         inline   = TRUE)
    })
    
    # idle message
    output$gost_status <- renderUI({
      if (running()) return(NULL)
      div(style = "color: #888; padding: 8px;",
          "Ready — configure options above and click Run g:Profiler.")
    })
    
    observeEvent(input$run_gost, {
      req(input$source_df, input$gene_col)
      
      df       <- data_store()[[input$source_df]]
      num_cols <- names(df)[sapply(df, is.numeric)]
      
      # build numeric df internally — user never sees this step
      numeric_df <- create_numeric_df(df, keep_cols = input$gene_col)
      
      running(TRUE)
      
      tryCatch({
        
        if (input$mode == "single") {
          req(input$single_sample, input$single_sources)
          
          # rename GeneSymbol col inside the list to match what run_gost_single expects
          numeric_df_renamed <- lapply(numeric_df, function(d) {
            names(d)[names(d) == input$gene_col] <- "GeneSymbol"
            d
          })
          
          result <- run_gost_single(
            numeric_df_list = numeric_df_renamed,
            sample_name     = input$single_sample,
            sources         = input$single_sources,
            top_n           = input$top_n,
            significant     = input$significant,
            max_attempts    = input$max_attempts
          )
          
          out_name <- if (trimws(input$output_name) == "") {
            paste0("__gost__", input$source_df, "__", input$single_sample)
          } else {
            trimws(input$output_name)
          }
          
          result_df       <- gost_to_df(result, sample_name = input$single_sample)
          updated         <- data_store()
          updated[[out_name]] <- result_df
          updated[[paste0("__raw__", out_name)]] <- result
          data_store(updated)
          save_data_store(updated)
          
          n_terms <- if (!is.null(result$result)) nrow(result$result) else 0
          running(FALSE)
          
          output$gost_status <- renderUI({
            div(style = "color: green; font-weight: bold; padding: 8px;",
                paste0("✓ Done — '", out_name, "': ", n_terms, " significant terms"))
          })
          
        } else {
          req(input$multi_samples, input$multi_source)
          
          # filter numeric_df to only selected samples
          numeric_df_filtered <- numeric_df[input$multi_samples]
          
          numeric_df_renamed <- lapply(numeric_df_filtered, function(d) {
            names(d)[names(d) == input$gene_col] <- "GeneSymbol"
            d
          })
          
          results <- run_gost_multi(
            numeric_df_list = numeric_df_renamed,
            source          = input$multi_source,
            top_n           = input$top_n,
            significant     = input$significant,
            max_attempts    = input$max_attempts
          )
          
          out_name <- if (trimws(input$output_name) == "") {
            paste0("__gost__multi__", input$source_df, "__", input$multi_source)
          } else {
            trimws(input$output_name)
          }
          
          result_df <- dplyr::bind_rows(
            base::mapply(gost_to_df, results, names(results), SIMPLIFY = FALSE)
          )
          
          updated         <- data_store()
          updated[[out_name]] <- result_df
          updated[[paste0("__raw__", out_name)]] <- results
          data_store(updated)
          save_data_store(updated)
          
          n_terms   <- nrow(result_df)
          n_samples <- length(results)
          running(FALSE)
          
          output$gost_status <- renderUI({
            div(style = "color: green; font-weight: bold; padding: 8px;",
                paste0("✓ Done — '", out_name, "': ",
                       n_terms, " total terms across ", n_samples, " samples"))
          })
        }
        
      }, error = function(e) {
        running(FALSE)
        is_connection <- grepl("CONNECTION_ERROR", e$message)
        clean_msg     <- gsub("CONNECTION_ERROR: ", "", e$message)
        
        output$gost_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              if (is_connection) paste0("🌐 Connection error: ", clean_msg)
              else paste0("✗ Error: ", clean_msg))
        })
      })
    })
    
  })
}