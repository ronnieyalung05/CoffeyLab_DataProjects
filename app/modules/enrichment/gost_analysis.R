# modules/enrichment/gost_analysis.R

gostAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("g:Profiler Enrichment Analysis"),
    tool_description(
      "Runs functional enrichment analysis using the g:Profiler web service. ",
      "Tests a ranked list of genes against Gene Ontology databases (BP, MF, CC), ",
      "KEGG pathways, Reactome, and WikiPathways. ",
      "Supports single-sample (multiple sources) or multi-sample (one source) modes."
    ),
    citation_link(
      "Kolberg L, et al. g:Profiler — interoperable web service for functional enrichment analysis. Nucleic Acids Research, 2023",
      "https://doi.org/10.1093/nar/gkad347"
    ),

    hr(),

    selectInput(ns("source_df"), "Dataset", choices = NULL),
    selectInput(ns("gene_col"),  "Gene Symbol column", choices = NULL),

    hr(),

    radioButtons(ns("mode"), "Analysis mode",
                 choices = c(
                   "Single sample, multiple sources" = "single",
                   "Multiple samples, one source"    = "multi"
                 ),
                 selected = "single",
                 inline   = TRUE),

    hr(),

    conditionalPanel(
      condition = paste0("input['", ns("mode"), "'] == 'single'"),
      selectInput(ns("single_sample"), "Sample (numeric column)", choices = NULL),
      checkboxGroupInput(ns("single_sources"), "Sources",
                         choices  = VALID_GOST_SOURCES,
                         selected = c("GO:BP", "GO:MF", "GO:CC"),
                         inline   = TRUE)
    ),

    conditionalPanel(
      condition = paste0("input['", ns("mode"), "'] == 'multi'"),
      uiOutput(ns("multi_sample_selector")),
      selectInput(ns("multi_source"), "Source",
                  choices  = VALID_GOST_SOURCES,
                  selected = "GO:BP")
    ),

    hr(),

    fluidRow(
      column(3,
        numericInput(ns("top_n"), "Top N genes by expression (0 = all)",
                     value = 0, min = 0, step = 10)
      ),
      column(3,
        checkboxInput(ns("significant"), "Significant terms only (p < 0.05)", value = TRUE)
      ),
      column(3,
        numericInput(ns("max_attempts"), "Max retry attempts",
                     value = 5, min = 1, max = 20, step = 1)
      )
    ),

    hr(),

    textInput(ns("output_name"), "Save results as",
              placeholder = "leave blank for auto name"),

    action_row(
      tool_button(ns("run_gost"), "Run g:Profiler")
    ),

    br(),
    uiOutput(ns("gost_status"))
  )
}

gostAnalysisServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {

    gost_msg <- reactiveVal(list(type = "idle", text = "Ready — configure options above and click Run g:Profiler."))

    output$gost_status <- renderUI({
      msg <- gost_msg()
      style <- switch(msg$type,
        idle    = "color: #888; padding: 8px;",
        success = "color: green; font-weight: bold; padding: 8px;",
        error   = "color: red; font-weight: bold; padding: 8px;"
      )
      div(style = style, msg$text)
    })

    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })

    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      df       <- data_store()[[input$source_df]]
      all_cols <- names(df)
      num_cols <- all_cols[sapply(df, is.numeric)]

      gene_default <- if ("GeneSymbol" %in% all_cols) "GeneSymbol" else all_cols[1]
      updateSelectInput(session, "gene_col", choices = all_cols, selected = gene_default)
      updateSelectInput(session, "single_sample", choices = num_cols)
    })

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

    observeEvent(input$run_gost, {
      req(input$source_df, input$gene_col)

      df         <- data_store()[[input$source_df]]
      numeric_df <- create_numeric_df(df, keep_cols = input$gene_col)

      button_loading(session, "run_gost")

      withProgress(message = "Running g:Profiler…", value = NULL, {

        tryCatch({

          if (input$mode == "single") {
            req(input$single_sample, input$single_sources)

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

            result_df           <- gost_to_df(result, sample_name = input$single_sample)
            updated             <- data_store()
            updated[[out_name]] <- result_df
            updated[[paste0("__raw__", out_name)]] <- result
            data_store(updated)
            save_data_store(updated)

            n_terms <- if (!is.null(result$result)) nrow(result$result) else 0
            gost_msg(list(type = "success",
                          text = paste0("✓ Done — '", out_name, "': ", n_terms, " significant terms")))

          } else {
            req(input$multi_samples, input$multi_source)

            numeric_df_filtered <- numeric_df[input$multi_samples]
            numeric_df_renamed  <- lapply(numeric_df_filtered, function(d) {
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

            updated             <- data_store()
            updated[[out_name]] <- result_df
            updated[[paste0("__raw__", out_name)]] <- results
            data_store(updated)
            save_data_store(updated)

            n_terms   <- nrow(result_df)
            n_samples <- length(results)
            gost_msg(list(type = "success",
                          text = paste0("✓ Done — '", out_name, "': ",
                                        n_terms, " total terms across ", n_samples, " samples")))
          }

        }, error = function(e) {
          is_connection <- grepl("CONNECTION_ERROR", e$message)
          clean_msg     <- gsub("CONNECTION_ERROR: ", "", e$message)
          gost_msg(list(type = "error",
                        text = if (is_connection) paste0("🌐 Connection error: ", clean_msg)
                               else paste0("✗ Error: ", clean_msg)))
        })

      }) # end withProgress

      button_reset(session, "run_gost")
    })

  })
}
