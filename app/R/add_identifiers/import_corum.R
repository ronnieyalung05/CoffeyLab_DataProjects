# R/add_identifiers/import_corum.R
#
# Parses an uploaded CORUM allComplexes.txt file and expands each complex
# into all pairwise protein combinations, producing the 2-column int_pairs
# table that CorShift() expects.

parse_corum_file <- function(file_path, organism = "Human") {

  # ── Read ────────────────────────────────────────────────────────────────────
  df <- utils::read.delim(file_path, sep = "\t", stringsAsFactors = FALSE,
                           encoding = "UTF-8", check.names = TRUE)

  # ── Filter to organism ──────────────────────────────────────────────────────
  org_col <- grep("^organism$", names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(org_col)) stop("Could not find Organism column. Is this a CORUM allComplexes.txt file?")

  df <- df[df[[org_col]] == organism, ]
  if (nrow(df) == 0) stop("No complexes found for organism '", organism, "' in this file.")

  # ── Find UniProt subunit column ─────────────────────────────────────────────
  uniprot_col <- grep("subunits_uniprot", names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(uniprot_col))  # fallback for older CORUM versions
    uniprot_col <- grep("uniprot", names(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(uniprot_col)) stop("Could not find a UniProt subunits column in this file.")

  # ── Expand each complex → all pairwise combinations ────────────────────────
  pairs_list <- base::lapply(seq_len(nrow(df)), function(i) {
    ids <- base::trimws(base::strsplit(df[[uniprot_col]][i], ";")[[1]])
    ids <- ids[nchar(ids) > 0]
    if (length(ids) < 2) return(NULL)
    combs <- t(utils::combn(ids, 2))
    data.frame(ProtA = combs[, 1], ProtB = combs[, 2],
               stringsAsFactors = FALSE)
  })

  pairs <- do.call(rbind, Filter(Negate(is.null), pairs_list))
  unique(pairs)
}
