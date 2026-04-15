# modules/add_identifiers/add_go_annotations.R

goAnnotationsUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("Fetch GO Annotations"),
    tool_description(
      "Queries the EBI QuickGO REST API for each protein in the selected dataset, ",
      "fetching associated Gene Ontology term IDs, names, and aspects (BP/MF/CC). ",
      "Requests are parallelised across workers to reduce wall-clock time. ",
      "One HTTP request is made per protein row."
    ),
    citation_link(
      "Binns D, et al. QuickGO: a web-based tool for Gene Ontology searching. Bioinformatics, 2009",
      "https://doi.org/10.1093/bioinformatics/btp536"
    ),

    hr(),

    selectInput(ns("source_df"),   "Select dataset",    choices = NULL),
    selectInput(ns("uniprot_col"), "UniProt ID column", choices = NULL),

    hr(),

    h5("Parallelism"),
    numericInput(ns("n_workers"), "Number of parallel workers",
                 value = max(1, parallel::detectCores() - 1),
                 min = 1, max = parallel::detectCores(), step = 1),
    p(style = "color: #888; font-size: 0.85em;",
      paste0("Your machine has ", parallel::detectCores(), " cores. ",
             "Recommended: leave 1 free for the Shiny session.")),

    hr(),

    h5("Output name"),
    textInput(ns("output_name"), "Save annotated dataset as",
              placeholder = "defaults to original name"),

    uiOutput(ns("row_count_warning")),

    action_row(
      tool_button(ns("run"), "Fetch GO Annotations")
    ),

    br(),
    shinycssloaders::withSpinner(
      uiOutput(ns("status")),
      type    = 6,
      color   = "#2980b9",
      size    = 0.6,
      hide.ui = FALSE
    )
  )
}

goAnnotationsServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {

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

    output$status <- renderUI({
      if (running()) return(NULL)
      div(style = "color: #888; padding: 8px;",
          "Ready — select a dataset and click Fetch GO Annotations.")
    })

    observeEvent(input$run, {
      req(input$source_df, input$uniprot_col)

      df      <- data_store()[[input$source_df]]
      workers <- max(1, input$n_workers)

      running(TRUE)
      button_loading(session, "run")

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
        n         <- nrow(df)

        running(FALSE)
        button_reset(session, "run")

        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Done — '", out_name, "': ",
                     n_with_go, " of ", n, " proteins had GO terms"))
        })

      }, error = function(e) {
        running(FALSE)
        button_reset(session, "run")

        is_connection_error <- grepl("CONNECTION_ERROR", e$message)
        clean_msg <- gsub("CONNECTION_ERROR: ", "", e$message)

        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              if (is_connection_error) {
                paste0("🌐 Connection error: ", clean_msg)
              } else {
                paste0("✗ Error: ", clean_msg)
              })
        })
      })
    })

  })
}
