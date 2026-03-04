# modules/plotting/heatmap_plots.R

heatmapPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("GO Enrichment Heatmap Generator"),
    p(style = "color: #888; font-size: 0.85em;",
      "Requires a flat g:Profiler/gost results dataset (no __raw__ prefix). Uses pvalues for heatmap generation"),
    
    selectInput(ns("gost_df"), "g:Profiler results dataset", choices = NULL),
    
    # dynamic sample checkboxes
    uiOutput(ns("sample_picker")),
    
    hr(),
    
    # combined vs per-sample view
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
        # only relevant for combined view
        conditionalPanel(
          condition = paste0("input['", ns("view_type"), "'] == 'combined'"),
          radioButtons(ns("per_sample_topn"), "How to apply Top N",
                       choices = c(
                         "Globally — top N terms across all selected samples. If all top N belong to one sample, only that sample appears." = "global",
                         "Per sample — each sample contributes its own top N terms. Union of all those terms shown together." = "per_sample"
                       ),
                       selected = "global")
        )
      )
    ),
    
    hr(),
    
    textInput(ns("plot_name"), "Save plot(s) as",
              placeholder = "leave blank for auto name"),
    
    actionButton(ns("run"), "Generate & Save Heatmap",
                 style = "color: white; background-color: #2980b9; width: 100%;"),
    
    br(), br(),
    
    shinycssloaders::withSpinner(
      uiOutput(ns("status")),
      type    = 8,
      color   = "#2980b9",
      size    = 0.7,
      hide.ui = FALSE
    )
  )
}

heatmapPlotServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {
    
    running <- reactiveVal(FALSE)
    
    # only show non-raw flat datasets
    observe({
      all_names  <- base::names(data_store())
      flat_names <- all_names[!startsWith(all_names, "__raw__")]
      updateSelectInput(session, "gost_df", choices = flat_names)
    })
    
    # populate sample checkboxes from the selected dataset
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
        # single sample result — no sample column
        div(style = "color: #888; font-size: 0.85em; padding: 4px 8px;",
            "Single sample result — no sample column detected.")
      }
    })
    
    output$status <- renderUI({
      if (running()) return(NULL)
      div(style = "color: #888; padding: 8px;",
          "Ready — select a dataset and click Generate.")
    })
    
    observeEvent(input$run, {
      req(input$gost_df, data_store()[[input$gost_df]])
      
      running(TRUE)
      
      tryCatch({
        df      <- data_store()[[input$gost_df]]
        has_sample_col <- "sample" %in% base::names(df)
        
        # filter to selected samples if multi-sample
        if (has_sample_col) {
          req(input$selected_samples)
          df <- df %>% dplyr::filter(sample %in% input$selected_samples)
          
          if (base::nrow(df) == 0) {
            base::stop("No data found for selected samples.")
          }
        } else {
          # single sample — add a dummy sample column so functions work uniformly
          df$sample <- base::gsub("__gost__", "", input$gost_df)
        }
        
        per_sample_topn <- !is.null(input$per_sample_topn) &&
                           input$per_sample_topn == "per_sample"
        
        plots <- if (input$view_type == "combined") {
          make_heatmap_combined(df,
                                top_n      = input$top_n,
                                per_sample = per_sample_topn)
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
        
        running(FALSE)
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Saved ", base::length(plots), " heatmap(s) to plot store"))
        })
        
      }, error = function(e) {
        running(FALSE)
        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })
  })
}