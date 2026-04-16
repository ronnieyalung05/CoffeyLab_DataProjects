# modules/plotting/corshift_heatmap.R

corshiftHeatmapUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("CorShift Heatmap"),
    tool_description(
      "Visualize significant protein pair co-abundance patterns from a CorShift result. ",
      "Rows = protein pairs, columns = individual samples (Group A then Group B). ",
      "Color = average z-scored co-abundance of the two proteins per sample. ",
      "Blue = co-low, Red = co-high. Group annotation bar shown at top."
    ),

    hr(),

    fluidRow(
      column(6,
        selectInput(ns("result_df"),  "CorShift results dataset",  choices = NULL),
        selectInput(ns("source_df"),  "Abundance matrix dataset",  choices = NULL)
      ),
      column(6,
        selectInput(ns("uniprot_col"), "UniProt ID column", choices = NULL),
        numericInput(ns("top_n"), "Top N pairs by correlation shift (0 = all)",
                     value = 0, min = 0, step = 10),
        radioButtons(ns("label_type"), "Row label",
                     choices  = c("Gene symbol" = "symbol", "UniProt ID" = "uniprot"),
                     selected = "symbol", inline = TRUE)
      )
    ),

    hr(),

    h5("Assign samples to groups"),
    p(style = "color: #888; font-size: 0.85em;",
      "Select which numeric columns belong to each group. Must match the groups used in CorShift."),

    fluidRow(
      column(6,
        tags$label("Group A", style = "font-weight: 600; color: #2980b9;"),
        uiOutput(ns("group_a_picker"))
      ),
      column(6,
        tags$label("Group B", style = "font-weight: 600; color: #c0392b;"),
        uiOutput(ns("group_b_picker"))
      )
    ),

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

corshiftHeatmapServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {

    status_msg <- reactiveVal(list(
      type = "idle",
      text = "Select a CorShift results dataset and the original abundance matrix, assign groups, then click Generate Heatmap."
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

    # ── Sync dataset dropdowns ────────────────────────────────────────────────
    observe({
      nms <- names(data_store())
      updateSelectInput(session, "result_df", choices = nms)
      updateSelectInput(session, "source_df", choices = nms)
    })

    # ── UniProt column picker ─────────────────────────────────────────────────
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      cols        <- names(data_store()[[input$source_df]])
      default_col <- if ("UniProtID" %in% cols) "UniProtID" else cols[1]
      updateSelectInput(session, "uniprot_col", choices = cols, selected = default_col)
    })

    # ── Numeric columns for group pickers ─────────────────────────────────────
    num_cols <- reactive({
      req(input$source_df, data_store()[[input$source_df]])
      df <- data_store()[[input$source_df]]
      names(df)[sapply(df, is.numeric)]
    })

    output$group_a_picker <- renderUI({
      cols <- num_cols()
      if (length(cols) == 0)
        return(div(style = "color: red; font-size: 0.85em;", "No numeric columns found."))
      checkboxGroupInput(session$ns("group_a"), label = NULL,
                         choices = cols, selected = NULL, inline = FALSE)
    })

    output$group_b_picker <- renderUI({
      cols <- num_cols()
      if (length(cols) == 0)
        return(div(style = "color: red; font-size: 0.85em;", "No numeric columns found."))
      checkboxGroupInput(session$ns("group_b"), label = NULL,
                         choices = cols, selected = NULL, inline = FALSE)
    })

    # ── Generate Heatmap ──────────────────────────────────────────────────────
    observeEvent(input$run, {
      req(input$result_df, input$source_df, input$uniprot_col,
          input$group_a, input$group_b)
      button_loading(session, "run")

      withProgress(message = "Generating CorShift heatmap…", value = NULL, {
        tryCatch({
          result_df <- data_store()[[input$result_df]]
          df        <- data_store()[[input$source_df]]

          # Build numeric matrix from the original abundance dataset
          all_sample_cols <- c(input$group_a, input$group_b)
          mat             <- as.matrix(df[, all_sample_cols, drop = FALSE])
          rownames(mat)   <- df[[input$uniprot_col]]
          mat             <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

          p <- make_corshift_heatmap(
            result_df    = result_df,
            mat          = mat,
            group_a_cols = input$group_a,
            group_b_cols = input$group_b,
            top_n        = as.integer(input$top_n),
            label_type   = input$label_type
          )
          if (is.null(p)) stop("Heatmap could not be generated — no valid pairs matched.")

          plot_name <- if (trimws(input$plot_name) == "") {
            paste0("corshift_heatmap_", input$result_df)
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
