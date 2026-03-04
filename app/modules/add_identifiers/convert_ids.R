# modules/convert_ids.R

convertIdsUI <- function(id) {
  ns <- NS(id)
  tagList(
    
    h4("ID Conversion"),
    
    selectInput(ns("source_df"), "Select dataset", choices = NULL),
    
    fluidRow(
      column(4,
        selectInput(ns("from_type"), "Convert from type",
                    choices  = VALID_ID_TYPES,
                    selected = "UNIPROT")
      ),
      column(4,
        selectInput(ns("from_col"), "Column containing those IDs", choices = NULL)
      ),
      column(4,
        selectInput(ns("to_types"), "Convert to type(s)",
                    choices  = VALID_ID_TYPES,
                    selected = c("SYMBOL", "ENTREZID"),
                    multiple = TRUE)
      )
    ),
    
    hr(),
    
    checkboxInput(ns("do_split"), "Split unmapped rows into separate dataset", value = TRUE),
    
    hr(),
    
    h4("Output names"),
    textInput(ns("output_name"), "Save converted dataset as",
              placeholder = "defaults to original name"),
    conditionalPanel(
      condition = paste0("input['", ns("do_split"), "'] == true"),
      textInput(ns("unmapped_name"), "Save unmapped dataset as",
                placeholder = "defaults to __originalname__unmapped")
    ),
    
    actionButton(ns("run"), "Run Conversion",
                 style = "color: white; background-color: #2980b9;"),
    uiOutput(ns("status"))
  )
}

convertIdsServer <- function(id, data_store) {
  moduleServer(id, function(input, output, session) {
    
    # sync dataset dropdown
    observe({
      updateSelectInput(session, "source_df", choices = names(data_store()))
    })
    
    # when dataset changes, update column dropdown
    observeEvent(input$source_df, {
      req(input$source_df, data_store()[[input$source_df]])
      cols <- names(data_store()[[input$source_df]])
      
      # guess a sensible default column based on from_type
      default_col <- switch(input$from_type,
        "UNIPROT"  = if ("UniProtID"  %in% cols) "UniProtID"  else cols[1],
        "SYMBOL"   = if ("GeneSymbol" %in% cols) "GeneSymbol" else cols[1],
        "ENTREZID" = if ("EntrezID"   %in% cols) "EntrezID"   else cols[1],
        cols[1]
      )
      
      updateSelectInput(session, "from_col", choices = cols, selected = default_col)
    })
    
    # also update column dropdown when from_type changes
    observeEvent(input$from_type, {
      req(input$source_df, data_store()[[input$source_df]])
      cols <- names(data_store()[[input$source_df]])
      
      default_col <- switch(input$from_type,
        "UNIPROT"  = if ("UniProtID"  %in% cols) "UniProtID"  else cols[1],
        "SYMBOL"   = if ("GeneSymbol" %in% cols) "GeneSymbol" else cols[1],
        "ENTREZID" = if ("EntrezID"   %in% cols) "EntrezID"   else cols[1],
        cols[1]
      )
      
      updateSelectInput(session, "from_col", choices = cols, selected = default_col)
    })
    
    observeEvent(input$run, {
      req(input$source_df, input$from_col, input$from_type, input$to_types)
      
      tryCatch({
        df        <- data_store()[[input$source_df]]
        converted <- convert_ids(
          df        = df,
          from_col  = input$from_col,
          from_type = input$from_type,
          to_types  = input$to_types
        )
        
        out_name <- if (trimws(input$output_name) == "") {
          input$source_df
        } else {
          trimws(input$output_name)
        }
        
        updated <- data_store()
        
        if (input$do_split) {
          # split on whichever columns were added by the conversion
          split_result  <- split_unmapped_ids(converted, check_cols = input$to_types)
          
          unmapped_name <- if (trimws(input$unmapped_name) == "") {
            paste0("__", input$source_df, "__unmapped")
          } else {
            trimws(input$unmapped_name)
          }
          
          updated[[out_name]]      <- split_result$valid
          updated[[unmapped_name]] <- split_result$unmapped
          
          msg <- paste0(
            "✓ Converted '", input$source_df, "' → '", out_name, "' (",
            nrow(split_result$valid), " mapped, ",
            nrow(split_result$unmapped), " unmapped → '", unmapped_name, "')"
          )
          
        } else {
          updated[[out_name]] <- converted
          msg <- paste0("✓ Converted '", input$source_df, "' → '", out_name, "'")
        }
        
        data_store(updated)
        save_data_store(updated)
        
        output$status <- renderUI({
          div(style = "color: green; font-weight: bold; padding: 8px;", msg)
        })
        
      }, error = function(e) {
        output$status <- renderUI({
          div(style = "color: red; font-weight: bold; padding: 8px;",
              paste0("✗ Error: ", e$message))
        })
      })
    })
    
  })
}