#install.packages("shiny")

library(shiny)
library(readxl)

ui <- fluidPage(
  tabsetPanel(
    tabPanel("Upload",
             fileInput("file", "Upload Excel"),
             textInput("dataset_name", "Dataset name (optional)", "dataset1"),
             actionButton("load_btn", "Load Dataset")
    ),
    
    tabPanel("Data",
             selectInput("selected_dataset", "Choose dataset:",
                         choices = NULL),
             tableOutput("data_preview")
    )
  )
)

server <- function(input, output, session) {
  
  # Persistent store for all datasets
  rv <- reactiveValues(data_list = list())
  
  # Load dataset when button clicked
  observeEvent(input$load_btn, {
    req(input$file)
    
    df <- read_excel(input$file$datapath)
    
    # Store the dataset
    rv$data_list[[input$dataset_name]] <- df
    
    # Update dropdown choices
    updateSelectInput(session, "selected_dataset",
                      choices = names(rv$data_list))
  })
  
  # Show preview of selected dataset
  output$data_preview <- renderTable({
    req(input$selected_dataset)
    head(rv$data_list[[input$selected_dataset]], 10)
  })
}

shinyApp(ui, server)