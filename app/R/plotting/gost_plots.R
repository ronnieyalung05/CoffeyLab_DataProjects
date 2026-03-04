# R/plots_gost.R

reconstruct_gost_result <- function(df) {
  base::list(result = df)
}

reconstruct_gost_list <- function(df) {
  if (!"sample" %in% base::names(df)) {
    base::stop("Dataset does not have a 'sample' column — cannot reconstruct multi-sample results.")
  }
  sample_names <- base::unique(df$sample)
  result <- base::lapply(sample_names, function(s) {
    base::list(result = df[df$sample == s, ])
  })
  base::names(result) <- sample_names
  result
}

make_gostplot_single <- function(gost_result, interactive = TRUE) {
  gprofiler2::gostplot(gost_result, capped = TRUE, interactive = interactive)
}

make_gostplot_multi_static <- function(gost_results_list) {
  result <- base::lapply(base::names(gost_results_list), function(sample_name) {
    res <- gost_results_list[[sample_name]]
    if (!base::is.null(res) && !base::is.null(res$result)) {
      gprofiler2::gostplot(res, capped = TRUE, interactive = FALSE) +
        ggplot2::ggtitle(paste("g:Profiler Results —", sample_name))
    } else {
      NULL
    }
  })
  base::names(result) <- base::names(gost_results_list)
  result
}