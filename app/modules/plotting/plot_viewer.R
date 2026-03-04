# modules/plotting/plot_viewer.R

plotViewerUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("View Saved Plots"),
    
    # wider dropdown using CSS to handle long names
    tags$style(HTML(paste0(
      "#", ns("selected_plot"), " { width: 100%; max-width: 800px; }"
    ))),
    selectInput(ns("selected_plot"), "Select a plot",
                choices = NULL, width = "100%"),
    
    hr(),
    
    # size controls
    fluidRow(
      column(4,
        sliderInput(ns("plot_width"), "Plot width (px)",
                    min = 400, max = 2000, value = 900, step = 50)
      ),
      column(4,
        sliderInput(ns("plot_height"), "Plot height (px)",
                    min = 300, max = 1600, value = 600, step = 50)
      )
    ),
    
    # render area
    uiOutput(ns("plot_area")),
    
    hr(),
    
    h4("Rename Plot"),
    fluidRow(
      column(5,
        textInput(ns("rename_to"), "New name",
                  placeholder = "e.g. my_heatmap_final")
      ),
      column(3, br(),
        actionButton(ns("run_rename"), "Rename",
                     style = "color: white; background-color: #8e44ad;")
      )
    ),
    uiOutput(ns("rename_status")),
    
    hr(),
    
    h4("Delete Plot"),
    actionButton(ns("delete_plot"), "Delete Selected Plot",
                 style = "color: white; background-color: #e74c3c;"),
    uiOutput(ns("delete_status"))
  )
}

plotViewerServer <- function(id, plot_store) {
  moduleServer(id, function(input, output, session) {
    
    observe({
      updateSelectInput(session, "selected_plot", choices = names(plot_store()))
    })
    
    # detect type and render with user-controlled dimensions
    output$plot_area <- renderUI({
      req(input$selected_plot, plot_store()[[input$selected_plot]])
      ns <- session$ns
      p  <- plot_store()[[input$selected_plot]]
      w  <- paste0(input$plot_width,  "px")
      h  <- paste0(input$plot_height, "px")
      
      if (inherits(p, "plotly")) {
        plotly::plotlyOutput(ns("plotly_render"), width = w, height = h)
      } else {
        plotOutput(ns("ggplot_render"), width = w, height = h)
      }
    })
    
    output$plotly_render <- plotly::renderPlotly({
      req(input$selected_plot)
      p <- plot_store()[[input$selected_plot]]
      req(inherits(p, "plotly"))
      p
    })
    
    output$ggplot_render <- renderPlot({
      req(input$selected_plot)
      p <- plot_store()[[input$selected_plot]]
      req(inherits(p, "gg") || inherits(p, "ggplot"))
      print(p)
    })
    
    # --- rename ---
    observeEvent(input$run_rename, {
      req(input$selected_plot, input$rename_to)
      tryCatch({
        new_name <- trimws(input$rename_to)
        if (new_name == "")                        stop("New name cannot be blank.")
        if (new_name %in% names(plot_store()))     stop("A plot with that name already exists.")
        
        updated                        <- plot_store()
        updated[[new_name]]            <- updated[[input$selected_plot]]
        updated[[input$selected_plot]] <- NULL
        plot_store(updated)
        save_plot_store(updated)
        
        output$rename_status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Renamed '", input$selected_plot, "' → '", new_name, "'"))
        })
      }, error = function(e) {
        output$rename_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ ", e$message))
        })
      })
    })
    
    # --- delete ---
    observeEvent(input$delete_plot, {
      req(input$selected_plot)
      tryCatch({
        deleted            <- input$selected_plot
        updated            <- plot_store()
        updated[[deleted]] <- NULL
        plot_store(updated)
        save_plot_store(updated)
        
        output$delete_status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Deleted '", deleted, "'"))
        })
      }, error = function(e) {
        output$delete_status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ ", e$message))
        })
      })
    })
    
  })
}