# modules/plotting/gost_plots.R

gostPlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("g:Profiler Plot Generator"),
    p(style = "color: #888; font-size: 0.85em;",
      "Requires a raw g:Profiler results object (datasets prefixed with __raw__)."),
    
    selectInput(ns("gost_df"), "g:Profiler raw results dataset", choices = NULL),
    
    # sample picker — populated after dataset selected
    uiOutput(ns("sample_picker")),
    
    radioButtons(ns("interactive"), "Plot type",
                 choices  = c("Interactive" = "interactive", "Static" = "static"),
                 selected = "interactive", inline = TRUE),
    
    hr(),
    
    textInput(ns("plot_name"), "Save plot(s) as",
              placeholder = "leave blank for auto name"),
    
    actionButton(ns("run"), "Generate & Save Plot",
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

gostPlotServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {
    
    running <- reactiveVal(FALSE)
    
    # only show __raw__ datasets
    observe({
      all_names  <- base::names(data_store())
      raw_names  <- all_names[startsWith(all_names, "__raw__")]
      updateSelectInput(session, "gost_df", choices = raw_names)
    })
    
    # detect if raw object is single or multi-sample
    # single = list with $result and $meta
    # multi  = named list of such objects
    is_multi <- reactive({
      req(input$gost_df, data_store()[[input$gost_df]])
      raw_obj <- data_store()[[input$gost_df]]
      # if first element has $result it's a multi list, otherwise it's single
      !is.null(raw_obj[[1]]$result)
    })
    
    # show sample checkboxes for multi, nothing for single
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
        # single sample — no picker needed
        div(style = "color: #888; font-size: 0.85em; padding: 4px 8px;",
            "Single sample result detected.")
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
        raw_obj     <- data_store()[[input$gost_df]]
        interactive <- input$interactive == "interactive"
        
        plots <- if (!is_multi()) {
          # single sample
          p <- make_gostplot_single(raw_obj, interactive = interactive)
          base::list(combined = p)
          
        } else {
          # multi — only plot selected samples
          req(input$selected_samples)
          selected <- input$selected_samples
          
          result <- base::lapply(selected, function(s) {
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
        
        running(FALSE)
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Saved ", base::length(plots), " plot(s) to plot store"))
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