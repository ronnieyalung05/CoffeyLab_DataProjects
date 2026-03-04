source("libraries.R")

### r/
# data manager
source("R/data_manager/load_data.R")
source("R/data_manager/clean_data.R")

# add identifiers
source("R/add_identifiers/convert_ids.R")
source("R/add_identifiers/add_go_annotations.R")

# enrichment
source("R/enrichment/gost_analysis.R")

# filtering
source("R/filtering/filter_go_proteins.R")

# plotting
source("R/plotting/gost_plots.R")
source("R/plotting/heatmap_plots.R")
source("R/plotting/plot_store.R")
source("R/plotting/gosemsim_plots.R")
source("R/plotting/gsea_cnet_plots.R")


#---------------------------------


### modules/
# data manager
source("modules/data_manager/load_data.R")
source("modules/data_manager/view_data.R")
source("modules/data_manager/save_data.R")
source("modules/data_manager/delete_data.R")
source("modules/data_manager/rename_data.R")
source("modules/data_manager/clean_data.R")

# add identifiers
source("modules/add_identifiers/convert_ids.R")
source("modules/add_identifiers/add_go_annotations.R")

# enrichment
source("modules/enrichment/gost_analysis.R")

# filtering
source("modules/filtering/filter_go_proteins.R")

# plotting
source("modules/plotting/gost_plots.R")
source("modules/plotting/heatmap_plots.R")
source("modules/plotting/plot_viewer.R")
source("modules/plotting/gosemsim_plots.R")
source("modules/plotting/gsea_cnet_plots.R")