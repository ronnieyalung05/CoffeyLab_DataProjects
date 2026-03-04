# modules/plotting/gosemsim_plots.R

goSemSimUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("GO Semantic Similarity Analysis"),
    p(style = "color: #888; font-size: 0.85em;",
      "Computes gene-level semantic similarity using GOSemSim. ",
      "Requires a dataset with a GeneSymbol column and numeric sample columns."),
    
    selectInput(ns("source_df"),  "Select dataset", choices = NULL),
    selectInput(ns("gene_col"),   "GeneSymbol column", choices = NULL),
    
    hr(),
    
    uiOutput(ns("sample_picker")),
    
    hr(),
    
    fluidRow(
      column(4,
        selectInput(ns("ontology"), "Ontology",
                    choices  = c("BP" = "BP", "MF" = "MF", "CC" = "CC"),
                    selected = "BP")
      ),
      column(4,
        numericInput(ns("top_n"), "Top N genes by expression",
                     value = 50, min = 2, step = 10)
      )
    ),
    
    hr(),
    
    h5("Output names"),
    fluidRow(
      column(4,
        textInput(ns("heatmap_name"), "Similarity heatmap name",
                  placeholder = "leave blank for auto name")
      ),
      column(4,
        textInput(ns("dendro_name"), "Dendrogram name",
                  placeholder = "leave blank for auto name")
      ),
      column(4,
        textInput(ns("matrix_name"), "Similarity matrix name",
                  placeholder = "leave blank for auto name")
      )
    ),
    
    actionButton(ns("run"), "Run Semantic Analysis",
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

goSemSimServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {
    
    running <- reactiveVal(FALSE)
    
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })
    
    # update gene col dropdown when dataset changes
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      cols        <- base::names(data_store()[[input$source_df]])
      gene_default <- if ("GeneSymbol" %in% cols) "GeneSymbol" else cols[1]
      updateSelectInput(session, "gene_col", choices = cols, selected = gene_default)
    })
    
    # populate sample checkboxes from numeric columns
    output$sample_picker <- renderUI({
      req(input$source_df, data_store()[[input$source_df]])
      df       <- data_store()[[input$source_df]]
      num_cols <- base::names(df)[base::sapply(df, base::is.numeric)]
      
      if (base::length(num_cols) == 0) {
        div(style = "color: red; padding: 4px 8px;",
            "No numeric columns found in this dataset.")
      } else {
        checkboxGroupInput(session$ns("selected_samples"),
                           "Select samples to include",
                           choices  = num_cols,
                           selected = num_cols,
                           inline   = FALSE)
      }
    })
    
    output$status <- renderUI({
      if (running()) return(NULL)
      div(style = "color: #888; padding: 8px;",
          "Ready — configure options above and click Run.")
    })
    
    observeEvent(input$run, {
      req(input$source_df, input$gene_col, input$selected_samples)
      
      running(TRUE)
      
      tryCatch({
        df <- data_store()[[input$source_df]]
        
        result <- run_semantic_analysis(
          df               = df,
          selected_samples = input$selected_samples,
          gene_col         = input$gene_col,
          selected_ont     = input$ontology,
          top_n            = input$top_n
        )
        
        base_label <- result$label
        
        # --- auto names ---
        heatmap_name <- if (trimws(input$heatmap_name) == "") {
          paste0("__semsim_heatmap__", input$source_df, "__", input$ontology)
        } else trimws(input$heatmap_name)
        
        dendro_name <- if (trimws(input$dendro_name) == "") {
          paste0("__semsim_dendro__", input$source_df, "__", input$ontology)
        } else trimws(input$dendro_name)
        
        matrix_name <- if (trimws(input$matrix_name) == "") {
          paste0("__semsim_matrix__", input$source_df, "__", input$ontology)
        } else trimws(input$matrix_name)
        
        # --- save plots to plot_store ---
        updated_plots                  <- plot_store()
        updated_plots[[heatmap_name]]  <- result$heatmap_plot
        updated_plots[[dendro_name]]   <- result$dendrogram_plot
        plot_store(updated_plots)
        save_plot_store(updated_plots)
        
        # --- save matrix to data_store ---
        updated_data                  <- data_store()
        updated_data[[matrix_name]]   <- result$sim_matrix_df
        data_store(updated_data)
        save_data_store(updated_data)
        
        running(FALSE)
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Done — heatmap saved as '", heatmap_name,
                     "', dendrogram as '", dendro_name,
                     "', matrix as '", matrix_name, "'"))
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