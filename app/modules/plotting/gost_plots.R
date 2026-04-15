# modules/plotting/gost_plots.R

gostPlotUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("g:Profiler Plot Generator"),
    tool_description(
      "Generates a Manhattan-style enrichment plot from raw g:Profiler results, ",
      "showing -log10(adjusted p-value) across GO terms and other sources. ",
      "Requires a raw results object (datasets prefixed with ", tags$code("__raw__"), ")."
    ),
    citation_link(
      "Kolberg L, et al. g:Profiler — interoperable web service for functional enrichment analysis. Nucleic Acids Research, 2023",
      "https://doi.org/10.1093/nar/gkad347"
    ),

    hr(),

    selectInput(ns("gost_df"), "g:Profiler raw results dataset", choices = NULL),

    uiOutput(ns("sample_picker")),

    radioButtons(ns("interactive"), "Plot type",
                 choices  = c("Interactive" = "interactive", "Static" = "static"),
                 selected = "interactive", inline = TRUE),

    hr(),

    textInput(ns("plot_name"), "Save plot(s) as",
              placeholder = "leave blank for auto name"),

    action_row(
      tool_button(ns("run"), "Generate & Save Plot")
    ),

    br(),
    uiOutput(ns("status"))
  )
}

gostPlotServer <- function(id, data_store, plot_store) {
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
      all_names <- base::names(data_store())
      raw_names <- all_names[startsWith(all_names, "__raw__")]
      updateSelectInput(session, "gost_df", choices = raw_names)
    })

    is_multi <- reactive({
      req(input$gost_df, data_store()[[input$gost_df]])
      raw_obj <- data_store()[[input$gost_df]]
      !is.null(raw_obj[[1]]$result)
    })

    output$sample_picker <- renderUI({
      req(input$gost_df, data_store()[[input$gost_df]])

      if (is_multi()) {
        raw_obj      <- data_store()[[input$gost_df]]
        sample_names <- base::names(raw_obj)
        checkboxGroupInput(session$ns("selected_samples"),
                           "Select samples to plot",
                           choices  = sample_names,
                           selected = sample_names,
                           inline   = TRUE)
      } else {
        div(style = "color: #888; font-size: 0.85em; padding: 4px 8px;",
            "Single sample result detected.")
      }
    })

    observeEvent(input$run, {
      req(input$gost_df, data_store()[[input$gost_df]])

      button_loading(session, "run")

      withProgress(message = "Generating g:Profiler plot…", value = NULL, {

        tryCatch({
          raw_obj     <- data_store()[[input$gost_df]]
          interactive <- input$interactive == "interactive"

          plots <- if (!is_multi()) {
            p <- make_gostplot_single(raw_obj, interactive = interactive)
            base::list(combined = p)
          } else {
            req(input$selected_samples)
            selected <- input$selected_samples
            result   <- base::lapply(selected, function(s) {
              make_gostplot_single(raw_obj[[s]], interactive = interactive)
            })
            base::names(result) <- selected
            result
          }

          updated <- plot_store()
          for (plot_name in base::names(plots)) {
            out_name <- if (trimws(input$plot_name) == "") {
              paste0("__gostplot__", input$gost_df, "__", plot_name)
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
                          text = paste0("✓ Saved ", base::length(plots), " plot(s) to plot store")))

        }, error = function(e) {
          status_msg(list(type = "error", text = paste0("✗ Error: ", e$message)))
        })

      }) # end withProgress

      button_reset(session, "run")
    })
  })
}
