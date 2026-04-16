# R/corshift/corshift_analysis.R
#
# Wrapper around AlteredPQR::CorShift().
# Handles data preparation (matrix conversion, index mapping, pair filtering)
# so that the Shiny module only passes column names and data frames.

run_corshift <- function(df,
                         uniprot_col,
                         group_a_cols,
                         group_b_cols,
                         int_pairs,
                         shift_threshold    = 0.6,
                         min_cor_in_samples = 0.6,
                         cor_signif         = 0.01) {

  # ── Validate inputs ──────────────────────────────────────────────────────────
  if (length(base::intersect(group_a_cols, group_b_cols)) > 0)
    stop("The same column cannot appear in both Group A and Group B.")

  if (length(group_a_cols) < 3 || length(group_b_cols) < 3)
    stop("Each group needs at least 3 samples for meaningful correlation estimates.")

  # ── Build numeric matrix (proteins × samples) ─────────────────────────────
  all_sample_cols <- c(group_a_cols, group_b_cols)
  mat             <- as.matrix(df[, all_sample_cols, drop = FALSE])
  rownames(mat)   <- df[[uniprot_col]]

  # Remove rows that are entirely NA
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

  # ── Column index vectors (1-based, into mat) ──────────────────────────────
  n_a   <- length(group_a_cols)
  n_b   <- length(group_b_cols)
  idx_a <- seq_len(n_a)
  idx_b <- seq_len(n_b) + n_a

  # ── Filter int_pairs to proteins present in mat ───────────────────────────
  present <- rownames(mat)
  ip_filt <- int_pairs[int_pairs[[1]] %in% present & int_pairs[[2]] %in% present, ]

  n_matched <- nrow(ip_filt)
  if (n_matched == 0)
    stop("None of the protein pairs in the interaction table matched proteins in the dataset. ",
         "Check that the UniProt ID column uses the same identifiers as the pairs table.")

  message("CorShift: ", n_matched, " pairs matched between dataset and interaction table.")

  # CorShift requires columns named exactly ProtA and ProtB
  ip_filt <- ip_filt[, 1:2, drop = FALSE]
  names(ip_filt) <- c("ProtA", "ProtB")

  # ── Run CorShift ──────────────────────────────────────────────────────────
  result <- tryCatch(
    AlteredPQR::CorShift(
      samplesA             = idx_a,
      samplesB             = idx_b,
      quant_data_all_local = mat,
      int_pairs_local      = ip_filt,
      shift_threshold      = shift_threshold,
      min_cor_in_samples   = min_cor_in_samples,
      cor_signif           = cor_signif,
      writeTable           = FALSE
    ),
    error = function(e) NULL   # empty result / no pairs passed thresholds
  )

  # ── Format output ─────────────────────────────────────────────────────────
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame(
      note = paste0("No significant correlation shifts found across ", n_matched,
                    " matched pairs. Try relaxing the thresholds: ",
                    "lower 'Min correlation shift' (e.g. 0.4) or raise 'Max p-value' (e.g. 0.05). ",
                    "Each group also needs >3 non-NA values per protein pair.")
    ))
  }

  result <- as.data.frame(result)
  colnames(result) <- make.unique(colnames(result))  # CorShift returns duplicate "cor_p_value" col names
  result <- tibble::rownames_to_column(result, var = "protein_pair")

  # Friendly column names — only rename if the expected 8 columns are present
  expected_names <- c(
    "protein_pair",
    "cor_group_A", "pval_group_A",
    "cor_group_B", "pval_group_B",
    "n_samples_A", "n_samples_B",
    "correlation_shift"
  )
  if (ncol(result) == length(expected_names)) {
    names(result) <- expected_names
  }

  result
}

# ── CorShift heatmap ──────────────────────────────────────────────────────────
#
# Rows = significant protein pairs (labeled by symbol_pair if available),
# Cols = individual samples (Group A then Group B), annotated with group bar.
# Cell value = average z-scored co-abundance: (z_ProtA + z_ProtB) / 2.

make_corshift_heatmap <- function(result_df, mat, group_a_cols, group_b_cols,
                                   top_n = 0, label_type = "symbol") {
  if (is.null(result_df) || "note" %in% names(result_df) || nrow(result_df) == 0)
    return(NULL)

  # Row labels: use gene symbols or UniProt IDs based on label_type
  labels <- if (label_type == "symbol" && "symbol_pair" %in% names(result_df))
               result_df$symbol_pair
             else
               result_df$protein_pair

  # Split protein_pair (format "ProtA-ProtB") at first dash
  prot_a <- sub("-.*$",    "", result_df$protein_pair)
  prot_b <- sub("^[^-]+-", "", result_df$protein_pair)

  # Z-score each protein across all samples (scale() works on columns → transpose)
  mat_z <- t(scale(t(mat)))

  # Only include sample columns that actually exist in mat
  all_cols <- base::intersect(c(group_a_cols, group_b_cols), colnames(mat_z))

  score_mat <- matrix(
    NA_real_, nrow = nrow(result_df), ncol = length(all_cols),
    dimnames = list(labels, all_cols)
  )

  for (i in seq_len(nrow(result_df))) {
    a_id <- prot_a[i]; b_id <- prot_b[i]
    if (!(a_id %in% rownames(mat_z)) || !(b_id %in% rownames(mat_z))) next
    score_mat[i, ] <- (mat_z[a_id, all_cols] + mat_z[b_id, all_cols]) / 2
  }

  # Sort rows by correlation_shift descending
  if ("correlation_shift" %in% names(result_df)) {
    ord       <- order(result_df$correlation_shift, decreasing = TRUE)
    score_mat <- score_mat[ord, , drop = FALSE]
  }

  # Filter to top N pairs if requested
  if (top_n > 0 && nrow(score_mat) > top_n)
    score_mat <- score_mat[seq_len(top_n), , drop = FALSE]

  # Column annotation bar (Group A = blue, Group B = red)
  n_a_actual <- length(base::intersect(group_a_cols, all_cols))
  n_b_actual <- length(base::intersect(group_b_cols, all_cols))
  ann_col <- data.frame(
    Group = factor(c(rep("Group A", n_a_actual), rep("Group B", n_b_actual))),
    row.names = all_cols
  )
  ann_colors <- list(Group = c("Group A" = "#2980b9", "Group B" = "#c0392b"))

  lim <- max(abs(score_mat), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  ph <- pheatmap::pheatmap(
    score_mat,
    cluster_rows      = nrow(score_mat) > 1,
    cluster_cols      = FALSE,
    annotation_col    = ann_col,
    annotation_colors = ann_colors,
    color  = colorRampPalette(c("#2166ac", "white", "#d6604d"))(60),
    breaks = seq(-lim, lim, length.out = 61),
    silent = TRUE
  )

  ggplotify::as.ggplot(ph)
}
