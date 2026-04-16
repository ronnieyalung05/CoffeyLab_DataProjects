# modules/plotting/expression_heatmap.R

expressionHeatmapUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("Expression Heatmap"),
    tool_description(
      "Visualize protein or gene abundance across samples as a clustered heatmap. ",
      "Select any dataset, choose which samples to include, and optionally filter ",
      "to the most variable rows. Rows = proteins/genes, columns = samples."
    ),

    hr(),

    fluidRow(
      column(6,
        selectInput(ns("source_df"),  "Dataset",           choices = NULL),
        selectInput(ns("label_col"),  "Row label column",  choices = NULL)
      ),
      column(6,
        numericInput(ns("top_n_var"), "Top N rows by variance (0 = all)",
                     value = 50, min = 0, step = 10),
        selectInput(ns("scale"), "Scale",
                    choices = c("None" = "none", "Row" = "row", "Column" = "column"),
                    selected = "row")
      )
    ),

    h5("Sample columns to include"),
    p(style = "color: #888; font-size: 0.85em;",
      "Un-check samples to exclude them from the heatmap."),
    uiOutput(ns("sample_picker")),

    hr(),

    textInput(ns("plot_name"), "Save plot as",
              placeholder = "leave blank for auto name"),
    action_row(
      tool_button(ns("run"), "Generate Heatmap", color = "#27ae60")
    ),
    br(),
    uiOutput(ns("status"))
  )
}

expressionHeatmapServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {

    status_msg <- reactiveVal(list(
      type = "idle",
      text = "Select a dataset and click Generate Heatmap."
    ))

    output$status <- renderUI({
      msg   <- status_msg()
      style <- switch(msg$type,
        idle    = "color: #888; padding: 8px;",
        success = "color: green; font-weight: bold; padding: 8px;",
        error   = "color: red; font-weight: bold; padding: 8px;"
      )
      div(style = style, msg$text)
    })

    # ── Dataset picker ────────────────────────────────────────────────────────
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })

    # ── Label column picker (character/factor cols) ───────────────────────────
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      df        <- data_store()[[input$source_df]]
      char_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
      updateSelectInput(session, "label_col",
                        choices  = char_cols,
                        selected = if (length(char_cols) > 0) char_cols[1] else NULL)
    })

    # ── Sample checkbox group (numeric cols) ──────────────────────────────────
    num_cols <- reactive({
      req(input$source_df, data_store()[[input$source_df]])
      df <- data_store()[[input$source_df]]
      names(df)[sapply(df, is.numeric)]
    })

    output$sample_picker <- renderUI({
      cols <- num_cols()
      if (length(cols) == 0)
        return(div(style = "color: red; font-size: 0.85em;", "No numeric columns found."))
      checkboxGroupInput(session$ns("sample_cols"), label = NULL,
                         choices = cols, selected = cols, inline = TRUE)
    })

    # ── Generate ──────────────────────────────────────────────────────────────
    observeEvent(input$run, {
      req(input$source_df, input$label_col, input$sample_cols)
      button_loading(session, "run")

      withProgress(message = "Generating expression heatmap…", value = NULL, {
        tryCatch({
          df <- data_store()[[input$source_df]]

          p <- make_expression_heatmap(
            df          = df,
            label_col   = input$label_col,
            sample_cols = input$sample_cols,
            top_n_var   = as.integer(input$top_n_var),
            scale       = input$scale
          )

          plot_name <- if (trimws(input$plot_name) == "") {
            paste0("expr_heatmap_", input$source_df)
          } else {
            trimws(input$plot_name)
          }

          updated              <- plot_store()
          updated[[plot_name]] <- p
          plot_store(updated)
          save_plot_store(updated)

          status_msg(list(
            type = "success",
            text = paste0("✓ Heatmap saved as '", plot_name, "' — view it in View Plots.")
          ))
        }, error = function(e) {
          status_msg(list(type = "error", text = paste0("✗ Error: ", e$message)))
        })
      })

      button_reset(session, "run")
    })

  })
}
