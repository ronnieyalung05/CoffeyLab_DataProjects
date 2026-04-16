# modules/corshift/corshift_analysis.R

corshiftUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("CorShift — Protein Complex Correlation Analysis"),
    tool_description(
      "Detects protein pairs whose co-expression pattern shifts significantly between ",
      "two sample groups. A strong correlation in one group but not the other indicates ",
      "complex assembly or disassembly that would be invisible to standard differential ",
      "abundance analysis. Requires a log-normalized abundance matrix ",
      "(rows = proteins, columns = samples) and a protein interaction pairs table ",
      "(e.g. imported from CORUM via the Add Identifiers tab)."
    ),
    citation_link(
      "Buljan M, et al. A computational framework for the inference of protein complex remodeling from whole-proteome measurements. Nature Methods, 2023",
      "https://doi.org/10.1038/s41592-023-02011-w"
    ),

    hr(),

    # ── Dataset inputs ────────────────────────────────────────────────────────
    fluidRow(
      column(6,
        selectInput(ns("source_df"),   "Abundance matrix dataset", choices = NULL),
        selectInput(ns("uniprot_col"), "UniProt ID column",        choices = NULL)
      ),
      column(6,
        selectInput(ns("pairs_df"), "Protein interaction pairs dataset", choices = NULL),
        p(style = "color: #888; font-size: 0.85em; margin-top: -6px;",
          "2-column UniProt pair table — import from CORUM or upload manually.")
      )
    ),

    hr(),

    # ── Sample group assignment ───────────────────────────────────────────────
    h5("Assign samples to groups"),
    p(style = "color: #888; font-size: 0.85em;",
      "Select which numeric columns belong to each group. ",
      "A column must not appear in both groups. Minimum 3 per group."),

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

    # ── Thresholds ────────────────────────────────────────────────────────────
    h5("Thresholds"),
    fluidRow(
      column(4,
        numericInput(ns("shift_threshold"),
                     "Min correlation shift |rA − rB|",
                     value = 0.6, min = 0, max = 1, step = 0.05)
      ),
      column(4,
        numericInput(ns("min_cor"),
                     "Min within-group correlation",
                     value = 0.6, min = 0, max = 1, step = 0.05)
      ),
      column(4,
        numericInput(ns("cor_signif"),
                     "Max p-value for correlation",
                     value = 0.01, min = 0.001, max = 0.1, step = 0.005)
      )
    ),

    hr(),

    textInput(ns("output_name"), "Save results as",
              placeholder = "leave blank for auto name"),

    action_row(
      tool_button(ns("run"), "Run CorShift")
    ),

    br(),
    uiOutput(ns("status")),

    br(),
    DT::dataTableOutput(ns("results_table"))
  )
}

corshiftServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {

    status_msg <- reactiveVal(list(
      type = "idle",
      text = "Ready — select datasets, assign sample groups, and click Run CorShift."
    ))

    results <- reactiveVal(NULL)

    # ── Status display ────────────────────────────────────────────────────────
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
      updateSelectInput(session, "source_df", choices = nms)
      updateSelectInput(session, "pairs_df",  choices = nms)
    })

    # ── Update UniProt column picker ──────────────────────────────────────────
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

    # ── Results table ─────────────────────────────────────────────────────────
    output$results_table <- DT::renderDataTable({
      req(results())
      res <- results()

      if ("note" %in% names(res))
        return(DT::datatable(res, options = list(dom = "t"), rownames = FALSE))

      round_cols <- base::intersect(
        c("cor_group_A", "pval_group_A", "cor_group_B", "pval_group_B", "correlation_shift"),
        names(res)
      )

      dt <- DT::datatable(
        res,
        options  = list(scrollX = TRUE, pageLength = 15,
                        order = list(list(ncol(res) - 1L, "desc"))),
        rownames = FALSE
      )
      if (length(round_cols) > 0) dt <- DT::formatRound(dt, columns = round_cols, digits = 4)
      dt
    })

    # ── Run CorShift ──────────────────────────────────────────────────────────
    observeEvent(input$run, {
      req(input$source_df, input$uniprot_col,
          input$pairs_df,
          input$group_a, input$group_b)

      button_loading(session, "run")

      withProgress(message = "Running CorShift…", value = NULL, {
        tryCatch({

          df        <- data_store()[[input$source_df]]
          int_pairs <- data_store()[[input$pairs_df]]

          result <- run_corshift(
            df                 = df,
            uniprot_col        = input$uniprot_col,
            group_a_cols       = input$group_a,
            group_b_cols       = input$group_b,
            int_pairs          = int_pairs,
            shift_threshold    = input$shift_threshold,
            min_cor_in_samples = input$min_cor,
            cor_signif         = input$cor_signif
          )

          # ── Add gene symbol + entrez ID columns ──────────────────────────
          if (!"note" %in% names(result)) {
            prot_a  <- sub("-.*$",    "", result$protein_pair)
            prot_b  <- sub("^[^-]+-", "", result$protein_pair)
            all_ids <- unique(c(prot_a, prot_b))

            mapping <- tryCatch(
              clusterProfiler::bitr(all_ids,
                                    fromType = "UNIPROT",
                                    toType   = c("SYMBOL", "ENTREZID"),
                                    OrgDb    = org.Hs.eg.db),
              error = function(e)
                data.frame(UNIPROT  = character(),
                           SYMBOL   = character(),
                           ENTREZID = character(),
                           stringsAsFactors = FALSE)
            )

            sym_map    <- if (nrow(mapping) > 0) setNames(mapping$SYMBOL,   mapping$UNIPROT) else character()
            entrez_map <- if (nrow(mapping) > 0) setNames(mapping$ENTREZID, mapping$UNIPROT) else character()

            sym_a    <- ifelse(prot_a %in% names(sym_map),    sym_map[prot_a],    "[unmapped]")
            sym_b    <- ifelse(prot_b %in% names(sym_map),    sym_map[prot_b],    "[unmapped]")
            entrez_a <- ifelse(prot_a %in% names(entrez_map), entrez_map[prot_a], "[unmapped]")
            entrez_b <- ifelse(prot_b %in% names(entrez_map), entrez_map[prot_b], "[unmapped]")

            result$symbol_pair <- paste(sym_a, sym_b, sep = "-")
            result$entrez_pair <- paste(entrez_a, entrez_b, sep = "-")

            # Reorder: protein_pair, symbol_pair, entrez_pair, then remaining columns
            rest_cols <- base::setdiff(names(result), c("protein_pair", "symbol_pair", "entrez_pair"))
            result    <- result[, c("protein_pair", "symbol_pair", "entrez_pair", rest_cols)]
          }

          out_name <- if (trimws(input$output_name) == "") {
            paste0("__corshift__", input$source_df)
          } else {
            trimws(input$output_name)
          }

          updated             <- data_store()
          updated[[out_name]] <- result
          data_store(updated)
          save_data_store(updated)

          results(result)

          if ("note" %in% names(result)) {
            status_msg(list(type = "idle", text = result$note[1]))
          } else {
            status_msg(list(
              type = "success",
              text = paste0("✓ Done — ", nrow(result),
                            " significant pair(s) found. Saved as '", out_name, "'.")
            ))
          }

        }, error = function(e) {
          status_msg(list(type = "error", text = paste0("✗ Error: ", e$message)))
        })
      })

      button_reset(session, "run")
    })

  })
}
