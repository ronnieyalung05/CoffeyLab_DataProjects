# R/plotting/heatmap_plots.R

# builds one combined heatmap for all provided samples
# per_sample = FALSE: top N globally across all samples
# per_sample = TRUE:  top N per sample, union of terms shown together
make_heatmap_combined <- function(gost_result_df, top_n = 0, per_sample = FALSE) {
  if (base::nrow(gost_result_df) == 0) base::stop("No valid results found.")
  
  df <- gost_result_df %>%
    dplyr::mutate(neg_log_p = -log10(p_value))
  
  if (top_n > 0) {
    if (per_sample) {
      # top N per sample — union of all those terms shown in one heatmap
      # each sample contributes its own top N, result may have more than top_n rows
      top_terms <- base::unique(base::unlist(
        base::lapply(base::unique(df$sample), function(s) {
          df %>%
            dplyr::filter(sample == s) %>%
            dplyr::arrange(p_value) %>%
            dplyr::slice_head(n = top_n) %>%
            dplyr::pull(term_name)
        })
      ))
      df <- df %>% dplyr::filter(term_name %in% top_terms)
    } else {
      # top N globally — if top 8 are all from one sample, only that sample shows
      df <- df %>%
        dplyr::arrange(p_value) %>%
        dplyr::slice_head(n = top_n)
    }
  }
  
  source_label <- if ("source" %in% base::names(df)) df$source[1] else ""
  
  p <- ggplot2::ggplot(df, ggplot2::aes(
      x    = sample,
      y    = reorder(term_name, neg_log_p),
      fill = neg_log_p
    )) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::labs(
      title = if (per_sample) {
        paste0("Top ", top_n, " Terms per Sample — Combined (", source_label, ")")
      } else {
        paste0("Top ", top_n, " Terms Globally (", source_label, ")")
      },
      x    = "Sample",
      y    = "Term Name",
      fill = "-log10(p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = ggplot2::element_text(size = 7),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.border = ggplot2::element_rect(color = "grey80", fill = NA)
    )
  
  base::list(combined = p)
}

# builds one heatmap per sample, returned as named list
# per_sample = FALSE: top N globally within that sample
# per_sample = TRUE:  same as FALSE for individual sample (per_sample has no extra meaning here)
make_heatmap_per_sample <- function(gost_result_df, top_n = 0) {
  if (base::nrow(gost_result_df) == 0) base::stop("No valid results found.")
  
  df      <- gost_result_df %>% dplyr::mutate(neg_log_p = -log10(p_value))
  samples <- base::unique(df$sample)
  
  plots <- base::lapply(samples, function(s) {
    sub_df <- df %>% dplyr::filter(sample == s)
    
    if (top_n > 0) {
      sub_df <- sub_df %>%
        dplyr::arrange(p_value) %>%
        dplyr::slice_head(n = top_n)
    }
    
    source_label <- if ("source" %in% base::names(sub_df)) sub_df$source[1] else ""
    
    ggplot2::ggplot(sub_df, ggplot2::aes(
        x    = source,
        y    = reorder(term_name, neg_log_p),
        fill = neg_log_p
      )) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient(low = "white", high = "red") +
      ggplot2::labs(
        title = paste0("Top ", top_n, " Terms — Sample: ", s),
        x     = "Source",
        y     = "Term Name",
        fill  = "-log10(p-value)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y  = ggplot2::element_text(size = 7),
        plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.border = ggplot2::element_rect(color = "grey80", fill = NA)
      )
  })
  
  base::names(plots) <- samples
  plots
}