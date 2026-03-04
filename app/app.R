source("sources.R")

ui <- page_navbar(
    title = "Functional Analysis of Proteomic Data",

    tags$head(tags$script(HTML(
        "window.onbeforeunload = function() { return 'Data may be lost on refresh. Are you sure?'; };"
    ))),

    nav_panel(
        "Data Manager",
        # Inner tabs
        navset_tab(
            nav_panel("Load Data", dataLoadUI("dataLoad")),
            nav_panel("View/Search Data", dataViewerUI("dataViewer")),
            nav_panel("Rename Dataset", renameDataUI("renameData")),
            nav_panel("Clean & Normalize", cleanDataUI("cleanData")),
            nav_panel("Delete Dataset", dataDeleteUI("dataDelete")),
            nav_panel("Save Data", dataSaveUI("dataSave")),
        )
    ),

    nav_panel(
        "Add Identifiers",
        # Inner tabs
        navset_tab(
            nav_panel("Add External IDs", convertIdsUI("convertIds")),
            nav_panel("Add GO Annotations", goAnnotationsUI("goAnnotations"))
        )
    ),

    nav_panel("Enrichment Analysis", gostAnalysisUI("gostAnalysis")),

    nav_panel("GO Cross-reference", filterProteinsUI("filterProteins")),

    nav_panel("Plots, Diagrams, Visualizations", 
        # Inner tabs 
        navset_tab(
            nav_panel("View Plots", plotViewerUI("plotViewer")),
            nav_panel("Create Gost Plots", gostPlotUI("gostPlot")),
            nav_panel("Create Heatmaps", heatmapPlotUI("heatmapPlot")),
            nav_panel("GO Semantic Similarity", goSemSimUI("goSemSim")),
            nav_panel("GSEA Cnetplot", gseaCnetUI("gseaCnet"))
        )
    )

)

server <- function(input, output, session) {
    # THIS is your persistent dataframe + plot store — a named list of dataframes and plots
    # .rds files are found in session_cache; it is best not to mess with those files directly
    # modules receives this and can read or update it
    data_store <- reactiveVal(load_data_store())
    plot_store <- reactiveVal(load_plot_store())
  
    # data manager
    dataLoadServer("dataLoad", data_store)
    dataViewerServer("dataViewer", data_store)
    renameDataServer("renameData", data_store)
    cleanDataServer("cleanData",   data_store)
    dataDeleteServer("dataDelete", data_store)
    dataSaveServer("dataSave", data_store) 

    # add identifiers
    convertIdsServer("convertIds", data_store) 
    goAnnotationsServer("goAnnotations", data_store)

    # enrichment
    gostAnalysisServer("gostAnalysis", data_store)

    # filtering
    filterProteinsServer("filterProteins", data_store)

    # plotting
    gostPlotServer("gostPlot", data_store, plot_store)
    heatmapPlotServer("heatmapPlot", data_store, plot_store) 
    plotViewerServer("plotViewer", plot_store) 
    goSemSimServer("goSemSim", data_store, plot_store)
    gseaCnetServer("gseaCnet", data_store, plot_store)
}

shinyApp(ui, server)