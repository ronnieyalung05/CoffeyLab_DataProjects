# modules/plotting/gsea_cnet_plots.R

CNET_PLOT_OPTIONS <- c(
  "FoldChange coloring — node color reflects expression level"         = "foldchange",
  "Category labels only — only pathway names shown"                    = "category_labels",
  "Gene labels only — only gene names shown"                           = "gene_labels",
  "All labels — both pathway and gene names shown"                     = "all_labels",
  "Custom colors — no labels, firebrick pathways / steelblue genes"    = "custom_colors"
)

gseaCnetUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("GSEA Concept Network (cnetplot)"),
    tool_description(
      "Runs Gene Set Enrichment Analysis (GSEA) on ranked gene lists and visualises ",
      "results as concept network plots, linking enriched GO pathways to the genes ",
      "that drive them. Multiple plot styles are available. ",
      "Requires a dataset with a GeneSymbol column and numeric expression columns."
    ),
    citation_link(
      "Yu G, et al. clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. OMICS, 2012",
      "https://doi.org/10.1089/omi.2011.0118"
    ),

    hr(),

    selectInput(ns("source_df"), "Select dataset",    choices = NULL),
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

    fluidRow(
      column(6,
        numericInput(ns("show_category"),
                     "Top N pathways to show across all plots",
                     value = 5, min = 1, max = 30, step = 1)
      )
    ),

    hr(),

    h5("Select plots to generate"),
    checkboxGroupInput(ns("selected_plots"), label = NULL,
                       choices  = CNET_PLOT_OPTIONS,
                       selected = c("foldchange", "all_labels"),
                       inline   = FALSE),

    hr(),

    textInput(ns("plot_name_prefix"), "Plot name prefix",
              placeholder = "leave blank for auto name"),

    action_row(
      tool_button(ns("run"), "Run GSEA & Generate Plots")
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
      button_loading(session, "run")

      tryCatch({
        df <- data_store()[[input$source_df]]

        desc_gene_list <- build_gene_list(
          df               = df,
          selected_samples = input$selected_samples,
          gene_col         = input$gene_col
        )

        result <- run_gsea_cnet(
          desc_gene_list = desc_gene_list,
          selected_ont   = input$ontology,
          selected_plots = input$selected_plots,
          show_category  = if (!is.null(input$show_category)) input$show_category else 5,
          pvalue_cutoff  = input$pvalue_cutoff,
          min_gs_size    = input$min_gs,
          max_gs_size    = input$max_gs
        )

        plots         <- result$plots
        n_sig         <- result$n_significant
        used_category <- result$show_category

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
        button_reset(session, "run")

        clamped <- used_category < input$show_category

        output$status <- renderUI({
          tagList(
            div(style = "color: #2980b9; font-weight: bold; padding: 4px 8px;",
                paste0("ℹ ", n_sig, " significant pathways found at p < ",
                       input$pvalue_cutoff,
                       if (clamped) {
                         paste0(" — requested ", input$show_category,
                                " pathways but only ", n_sig,
                                " found, showing all ", n_sig)
                       } else {
                         paste0(" — showing top ", used_category)
                       }
                )),
            div(style = "color: green; font-weight: bold; padding: 4px 8px;",
                paste0("✓ Saved ", base::length(plots), " plot(s): ",
                       base::paste(base::names(plots), collapse = ", ")))
          )
        })

      }, error = function(e) {
        running(FALSE)
        button_reset(session, "run")
        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })

  })
}
