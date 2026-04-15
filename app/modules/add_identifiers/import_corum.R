# modules/add_identifiers/import_corum.R

importCorumUI <- function(id) {
  ns <- NS(id)
  tagList(

    h4("Import CORUM Protein Complex Pairs"),
    tool_description(
      "Parses a CORUM ", tags$code("allComplexes.txt"), " file and expands each human protein complex ",
      "into all pairwise UniProt ID combinations. ",
      "The resulting table is the interaction pairs input required by CorShift analysis."
    ),
    citation_link(
      "Giurgiu M, et al. CORUM: the comprehensive resource of mammalian protein complexes. Nucleic Acids Research, 2019",
      "https://doi.org/10.1093/nar/gky973"
    ),

    hr(),

    # Download instructions
    div(
      style = "background: #eaf4fb; border-left: 3px solid #2980b9; border-radius: 3px; padding: 10px 14px; margin-bottom: 12px;",
      tags$b("Step 1 — Download the CORUM file"),
      tags$ol(
        style = "margin: 6px 0 0 0; padding-left: 18px; font-size: 0.9em;",
        tags$li(tags$a("Open the CORUM download page",
                       href   = "https://mips.helmholtz-muenchen.de/corum/download",
                       target = "_blank",
                       rel    = "noopener noreferrer")),
        tags$li("Download ", tags$code("allComplexes.txt"), " (plain text, not ZIP)"),
        tags$li("Upload it below")
      )
    ),

    tags$b("Step 2 — Upload the file", style = "font-size: 0.9em;"),
    fileInput(ns("corum_file"), label = NULL,
              accept = c(".txt", "text/plain"),
              placeholder = "allComplexes.txt"),

    h5("Output name"),
    textInput(ns("output_name"), label = NULL,
              placeholder = "defaults to 'corum_human_pairs'"),

    uiOutput(ns("status"))
  )
}

importCorumServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {

    status_msg <- reactiveVal(list(
      type = "idle",
      text = "Upload the allComplexes.txt file above to import CORUM pairs."
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

    # Auto-process as soon as a file is uploaded (no button needed)
    observeEvent(input$corum_file, {
      req(input$corum_file)

      withProgress(message = "Parsing CORUM file…", value = NULL, {
        tryCatch({
          pairs <- parse_corum_file(input$corum_file$datapath, organism = "Human")

          out_name <- if (trimws(input$output_name) == "") {
            "corum_human_pairs"
          } else {
            trimws(input$output_name)
          }

          updated             <- data_store()
          updated[[out_name]] <- pairs
          data_store(updated)
          save_data_store(updated)

          status_msg(list(
            type = "success",
            text = paste0("✓ Done — '", out_name, "': ",
                          nrow(pairs), " unique protein pairs from human CORUM complexes")
          ))

        }, error = function(e) {
          status_msg(list(type = "error", text = paste0("✗ Error: ", e$message)))
        })
      })
    })

  })
}
