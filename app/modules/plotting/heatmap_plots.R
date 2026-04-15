# modules/plotting/heatmap_plots.R

heatmapPlotUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("GO Enrichment Heatmap"),
    tool_description(
      "Generates a clustered heatmap from flat g:Profiler results, ",
      "using -log10(p-value) as the color intensity. ",
      "Rows and columns are hierarchically clustered, producing dendrograms ",
      "that reveal which terms and samples are most similar. ",
      "Requires a flat (non-raw) g:Profiler results dataset."
    ),
    citation_link(
      "Kolberg L, et al. g:Profiler — interoperable web service for functional enrichment analysis. Nucleic Acids Research, 2023",
      "https://doi.org/10.1093/nar/gkad347"
    ),

    hr(),

    selectInput(ns("gost_df"), "g:Profiler results dataset", choices = NULL),

    uiOutput(ns("sample_picker")),

    hr(),

    radioButtons(ns("view_type"), "Output view",
                 choices = c(
                   "Combined — all selected samples in one heatmap" = "combined",
                   "Per sample — one heatmap per selected sample"   = "per_sample"
                 ),
                 selected = "combined"),

    hr(),

    fluidRow(
      column(5,
        numericInput(ns("top_n"), "Top N terms (0 = all)",
                     value = 20, min = 0, step = 5)
      ),
      column(7,
        conditionalPanel(
          condition = paste0("input['", ns("view_type"), "'] == 'combined'"),
          radioButtons(ns("per_sample_topn"), "How to apply Top N",
                       choices = c(
                         "Globally — top N terms across all selected samples" = "global",
                         "Per sample — each sample contributes its own top N (union shown)" = "per_sample"
                       ),
                       selected = "global")
        )
      )
    ),

    hr(),

    textInput(ns("plot_name"), "Save plot(s) as",
              placeholder = "leave blank for auto name"),

    action_row(
      tool_button(ns("run"), "Generate & Save Heatmap")
    ),

    br(),
    uiOutput(ns("status"))
  )
}

heatmapPlotServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {

    status_msg <- reactiveVal(list(type = "idle", text = "Ready — select a dataset and click Generate."))

    output$status <- renderUI({
      msg <- status_msg()
      style <- switch(msg$type,
        idle    = "color: #888; padding: 8px;",
        success = "color: green; font-weight: bold; padding: 8px;",
        error   = "color: red; font-weight: bold; padding: 8px;"
      )
      div(style = style, msg$text)
    })

    observe({
      all_names  <- base::names(data_store())
      flat_names <- all_names[!startsWith(all_names, "__raw__")]
      updateSelectInput(session, "gost_df", choices = flat_names)
    })

    output$sample_picker <- renderUI({
      req(input$gost_df, data_store()[[input$gost_df]])
      df <- data_store()[[input$gost_df]]

      if ("sample" %in% base::names(df)) {
        samples <- base::unique(df$sample)
        checkboxGroupInput(session$ns("selected_samples"),
                           "Select samples to include",
                           choices  = samples,
                           selected = samples,
                           inline   = FALSE)
      } else {
        div(style = "color: #888; font-size: 0.85em; padding: 4px 8px;",
            "Single sample result — no sample column detected.")
      }
    })

    observeEvent(input$run, {
      req(input$gost_df, data_store()[[input$gost_df]])

      button_loading(session, "run")

      withProgress(message = "Generating heatmap…", value = NULL, {

        tryCatch({
          df             <- data_store()[[input$gost_df]]
          has_sample_col <- "sample" %in% base::names(df)

          if (has_sample_col) {
            req(input$selected_samples)
            df <- df %>% dplyr::filter(sample %in% input$selected_samples)
            if (base::nrow(df) == 0) base::stop("No data found for selected samples.")
          } else {
            df$sample <- base::gsub("__gost__", "", input$gost_df)
          }

          per_sample_topn <- !is.null(input$per_sample_topn) &&
                             input$per_sample_topn == "per_sample"

          plots <- if (input$view_type == "combined") {
            make_heatmap_combined(df, top_n = input$top_n, per_sample = per_sample_topn)
          } else {
            make_heatmap_per_sample(df, top_n = input$top_n)
          }

          updated <- plot_store()
          for (plot_name in base::names(plots)) {
            out_name <- if (trimws(input$plot_name) == "") {
              paste0("__heatmap__", input$gost_df, "__", plot_name)
            } else if (base::length(plots) > 1) {
              paste0(trimws(input$plot_name), "__", plot_name)
            } else {
              trimws(input$plot_name)
            }
            updated[[out_name]] <- plots[[plot_name]]
          }

          plot_store(updated)
          save_plot_store(updated)

          status_msg(list(type = "success",
                          text = paste0("✓ Saved ", base::length(plots), " heatmap(s) to plot store")))

        }, error = function(e) {
          status_msg(list(type = "error", text = paste0("✗ Error: ", e$message)))
        })

      }) # end withProgress

      button_reset(session, "run")
    })
  })
}
