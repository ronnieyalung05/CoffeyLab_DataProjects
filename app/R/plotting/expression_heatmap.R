# R/plotting/expression_heatmap.R
#
# General-purpose expression heatmap using pheatmap.
# Rows = proteins/genes, Columns = user-selected samples.

make_expression_heatmap <- function(df, label_col, sample_cols,
                                    top_n_var = 0, scale = "none") {

  mat            <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat)  <- df[[label_col]]

  # Drop all-NA rows
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (nrow(mat) == 0) stop("No rows with data after removing all-NA rows.")

  # Filter to top N most variable rows
  if (top_n_var > 0 && nrow(mat) > top_n_var) {
    row_vars <- apply(mat, 1, var, na.rm = TRUE)
    mat      <- mat[order(row_vars, decreasing = TRUE)[seq_len(top_n_var)], , drop = FALSE]
  }

  # Diverging palette for scaled data, sequential for raw values
  colors <- if (scale %in% c("row", "column")) {
    colorRampPalette(c("#2166ac", "white", "#d6604d"))(60)
  } else {
    colorRampPalette(c("white", "#d6604d"))(60)
  }

  ph <- pheatmap::pheatmap(
    mat,
    scale        = scale,
    cluster_rows = nrow(mat) > 1,
    cluster_cols = ncol(mat) > 1,
    color        = colors,
    silent       = TRUE
  )

  ggplotify::as.ggplot(ph)
}
