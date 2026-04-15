# R/plotting/heatmap_plots.R
# Uses pheatmap() for clustering and dendrograms.

# ── Shared helper ──────────────────────────────────────────────────────────────

# Build a -log10(p_value) matrix (terms × samples) from a long-format gost df.
# Missing combinations are filled with 0.
.build_heatmap_matrix <- function(df, row_col, col_col) {
  df %>%
    dplyr::select(dplyr::all_of(c(row_col, col_col, "neg_log_p"))) %>%
    # keep highest value if a term appears more than once for a given sample
    dplyr::group_by(dplyr::across(dplyr::all_of(c(row_col, col_col)))) %>%
    dplyr::summarise(neg_log_p = max(neg_log_p), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = dplyr::all_of(col_col),
      values_from = neg_log_p,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames(row_col) %>%
    as.matrix()
}

# ── Combined heatmap ───────────────────────────────────────────────────────────

# Builds one combined heatmap (terms × samples) with clustering and dendrogram.
# per_sample = FALSE : top N globally across all samples
# per_sample = TRUE  : top N per sample; union of those terms shown together
make_heatmap_combined <- function(gost_result_df, top_n = 0, per_sample = FALSE) {
  if (base::nrow(gost_result_df) == 0) base::stop("No valid results found.")

  df <- gost_result_df %>%
    dplyr::mutate(neg_log_p = -log10(p_value))

  if (top_n > 0) {
    if (per_sample) {
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
      top_terms <- df %>%
        dplyr::arrange(p_value) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(term_name) %>%
        base::unique()
      df <- df %>% dplyr::filter(term_name %in% top_terms)
    }
  }

  source_label <- if ("source" %in% base::names(df)) df$source[1] else ""

  title <- if (per_sample) {
    paste0("Top ", top_n, " Terms per Sample — Combined (", source_label, ")")
  } else if (top_n > 0) {
    paste0("Top ", top_n, " Terms Globally (", source_label, ")")
  } else {
    paste0("All Terms — Combined (", source_label, ")")
  }

  mat <- .build_heatmap_matrix(df, row_col = "term_name", col_col = "sample")

  # cluster rows only when there are enough rows; cluster cols only when > 1
  cluster_r <- nrow(mat) > 1
  cluster_c <- ncol(mat) > 1

  ph <- pheatmap::pheatmap(
    mat,
    color         = grDevices::colorRampPalette(c("#f9f9f9", "#e74c3c"))(60),
    cluster_rows  = cluster_r,
    cluster_cols  = cluster_c,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row  = 8,
    fontsize_col  = 9,
    main          = title,
    border_color  = "white",
    silent        = TRUE
  )

  base::list(combined = ggplotify::as.ggplot(ph))
}

# ── Per-sample heatmap ─────────────────────────────────────────────────────────

# Builds one pheatmap per sample; each shows terms × GO sources.
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

    source_label <- if ("source" %in% base::names(sub_df)) sub_df$source[1] else s

    # If multiple GO sources exist, show terms × sources; else a single-column matrix
    if ("source" %in% base::names(sub_df) && length(unique(sub_df$source)) > 1) {
      mat <- .build_heatmap_matrix(sub_df, row_col = "term_name", col_col = "source")
    } else {
      mat <- matrix(
        sub_df$neg_log_p,
        ncol     = 1,
        dimnames = list(sub_df$term_name, source_label)
      )
    }

    cluster_r <- nrow(mat) > 1
    cluster_c <- ncol(mat) > 1

    ph <- pheatmap::pheatmap(
      mat,
      color         = grDevices::colorRampPalette(c("#f9f9f9", "#e74c3c"))(60),
      cluster_rows  = cluster_r,
      cluster_cols  = cluster_c,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row  = 8,
      fontsize_col  = 9,
      main          = paste0(if (top_n > 0) paste0("Top ", top_n, " Terms") else "All Terms",
                             " — Sample: ", s),
      border_color  = "white",
      silent        = TRUE
    )

    ggplotify::as.ggplot(ph)
  })

  base::names(plots) <- samples
  plots
}
