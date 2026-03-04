# modules/plotting/gsea_cnet_plots.R

CNET_PLOT_OPTIONS <- c(
  "p1 — basic cnetplot with foldChange coloring"       = "p1",
  "p2 — node size scaled by p-value"                   = "p2",
  "p4 — pathway/category labels only"                  = "p4",
  "p5 — gene labels only"                              = "p5",
  "p6 — all labels shown"                              = "p6",
  "p7 — no labels, firebrick/steelblue color scheme"   = "p7"
)

gseaCnetUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("GSEA Concept Network (cnetplot)"),
    p(style = "color: #888; font-size: 0.85em;",
      "Runs Gene Set Enrichment Analysis and generates concept network plots. ",
      "Requires a dataset with a GeneSymbol column and numeric sample columns."),
    
    selectInput(ns("source_df"), "Select dataset", choices = NULL),
    selectInput(ns("gene_col"),  "GeneSymbol column", choices = NULL),
    
    hr(),
    
    uiOutput(ns("sample_picker")),
    
    hr(),
    
    fluidRow(
      column(3,
        selectInput(ns("ontology"), "Ontology",
                    choices  = c("BP" = "BP", "MF" = "MF", "CC" = "CC"),
                    selected = "BP")
      ),
      column(3,
        numericInput(ns("pvalue_cutoff"), "p-value cutoff",
                     value = 0.05, min = 0.001, max = 1, step = 0.01)
      ),
      column(3,
        numericInput(ns("min_gs"), "Min gene set size",
                     value = 10, min = 1, step = 5)
      ),
      column(3,
        numericInput(ns("max_gs"), "Max gene set size",
                     value = 500, min = 10, step = 50)
      )
    ),
    
    hr(),
    
    h5("Select plots to generate"),
    checkboxGroupInput(ns("selected_plots"), label = NULL,
                       choices  = CNET_PLOT_OPTIONS,
                       selected = c("p1", "p4", "p6"),
                       inline   = FALSE),
    
    hr(),
    
    textInput(ns("plot_name_prefix"), "Plot name prefix",
              placeholder = "leave blank for auto name"),
    
    actionButton(ns("run"), "Run GSEA & Generate Plots",
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

gseaCnetServer <- function(id, data_store, plot_store) {
  moduleServer(id, function(input, output, session) {
    
    running <- reactiveVal(FALSE)
    
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })
    
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      cols         <- base::names(data_store()[[input$source_df]])
      gene_default <- if ("GeneSymbol" %in% cols) "GeneSymbol" else cols[1]
      updateSelectInput(session, "gene_col", choices = cols, selected = gene_default)
    })
    
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
          "Ready — configure options and click Run.")
    })
    
    observeEvent(input$run, {
      req(input$source_df, input$gene_col,
          input$selected_samples, input$selected_plots)
      
      running(TRUE)
      
      tryCatch({
        df <- data_store()[[input$source_df]]
        
        # build sorted descending gene list
        desc_gene_list <- build_gene_list(
          df               = df,
          selected_samples = input$selected_samples,
          gene_col         = input$gene_col
        )
        
        plots <- run_gsea_cnet(
          desc_gene_list = desc_gene_list,
          selected_ont   = input$ontology,
          selected_plots = input$selected_plots,
          pvalue_cutoff  = input$pvalue_cutoff,
          min_gs_size    = input$min_gs,
          max_gs_size    = input$max_gs
        )
        
        # save each selected plot to plot_store
        updated <- plot_store()
        
        for (plot_key in base::names(plots)) {
          out_name <- if (trimws(input$plot_name_prefix) == "") {
            paste0("__cnet__", input$source_df, "__",
                   input$ontology, "__", plot_key)
          } else {
            paste0(trimws(input$plot_name_prefix), "__", plot_key)
          }
          updated[[out_name]] <- plots[[plot_key]]
        }
        
        plot_store(updated)
        save_plot_store(updated)
        
        running(FALSE)
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;",
              paste0("✓ Saved ", base::length(plots), " plot(s) to plot store: ",
                     base::paste(base::names(plots), collapse = ", ")))
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