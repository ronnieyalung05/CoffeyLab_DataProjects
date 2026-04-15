# modules/filtering/filter_go_proteins.R

filterProteinsUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("GO Cross-reference"),
    tool_description(
      "Filters a protein dataset by intersecting each protein's GO term IDs ",
      "(fetched via QuickGO) with the enriched GO terms returned by g:Profiler. ",
      "Produces two output datasets: proteins whose GO annotations overlap with ",
      "at least one enriched term (matched), and those that do not (unmatched)."
    ),

    hr(),

    selectInput(ns("reference_df"),    "Reference dataset (contains goIds column)", choices = NULL),
    selectInput(ns("gost_results_df"), "g:Profiler results dataset (contains term_id column)", choices = NULL),

    hr(),

    h5("Columns to keep in output"),
    uiOutput(ns("cols_selector")),

    hr(),

    fluidRow(
      column(6,
        textInput(ns("filtered_name"), "Save matched proteins as",
                  placeholder = "leave blank for auto name")
      ),
      column(6,
        textInput(ns("lost_name"), "Save unmatched proteins as",
                  placeholder = "leave blank for auto name")
      )
    ),

    action_row(
      tool_button(ns("run"), "Run Cross-reference")
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

filterProteinsServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {

    running <- reactiveVal(FALSE)

    observe({
      nms <- base::names(data_store())
      updateSelectInput(session, "reference_df",    choices = nms)
      updateSelectInput(session, "gost_results_df", choices = nms)
    })

    output$cols_selector <- renderUI({
      req(input$reference_df, data_store()[[input$reference_df]])
      cols <- c(base::names(data_store()[[input$reference_df]]), "matched_goIds")

      checkboxGroupInput(session$ns("cols_to_keep"),
                         label    = NULL,
                         choices  = cols,
                         selected = cols,
                         inline   = TRUE)
    })

    output$status <- renderUI({
      if (running()) return(NULL)
      div(style = "color: #888; padding: 8px;",
          "Ready — select datasets and click Run Cross-reference.")
    })

    observeEvent(input$run, {
      req(input$reference_df, input$gost_results_df, input$cols_to_keep)

      running(TRUE)
      button_loading(session, "run")

      tryCatch({
        ref_df  <- data_store()[[input$reference_df]]
        gost_df <- data_store()[[input$gost_results_df]]

        if (!"term_id" %in% base::names(gost_df)) {
          base::stop("g:Profiler results dataset does not have a 'term_id' column.")
        }

        result <- filter_by_go_crossref(
          reference_df = ref_df,
          gost_df      = gost_df,
          cols_to_keep = input$cols_to_keep
        )

        filtered_name <- if (trimws(input$filtered_name) == "") {
          paste0("__crossref_matched__", input$reference_df)
        } else {
          trimws(input$filtered_name)
        }

        lost_name <- if (trimws(input$lost_name) == "") {
          paste0("__crossref_unmatched__", input$reference_df)
        } else {
          trimws(input$lost_name)
        }

        updated                  <- data_store()
        updated[[filtered_name]] <- result$filtered_data
        updated[[lost_name]]     <- result$lost_data
        data_store(updated)
        save_data_store(updated)

        running(FALSE)
        button_reset(session, "run")

        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              base::paste0("✓ Done — ", result$kept_count, " proteins matched, ",
                           result$lost_count, " unmatched"))
        })

      }, error = function(e) {
        running(FALSE)
        button_reset(session, "run")

        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              base::paste0("✗ Error: ", e$message))
        })
      })
    })

  })
}
