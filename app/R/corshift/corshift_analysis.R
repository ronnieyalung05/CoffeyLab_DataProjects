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
