# R/ui_helpers.R — shared UI helpers, constants, and button utilities
# Loaded first (see sources.R) so all modules can depend on it.

# ── Constants ──────────────────────────────────────────────────────────────────

VALID_ID_TYPES <- c("UNIPROT", "SYMBOL", "ENTREZID")

VALID_GOST_SOURCES <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP")

# ── UI component helpers ───────────────────────────────────────────────────────

# Styled action button — NOT full width; sits naturally at the bottom of a form.
# color   : hex background color
# width   : e.g. "200px", or NULL for natural (content) width
# compact : TRUE for smaller delete/rename buttons; FALSE (default) for primary actions
tool_button <- function(ns_id, label, color = "#2980b9", width = NULL, compact = FALSE) {
  padding   <- if (compact) "5px 12px" else "6px 18px"
  font_size <- if (compact) "font-size: 0.85em; " else ""
  style <- paste0(
    "color: white; background-color: ", color, "; ",
    "padding: ", padding, "; ",
    font_size,
    "border-radius: 4px; font-weight: 500; margin-top: 6px; white-space: nowrap;",
    if (!is.null(width)) paste0(" width: ", width, ";") else ""
  )
  actionButton(ns_id, label, style = style)
}

# Small info-icon link (opens a modal or scrolls to a citation)
info_button <- function(ns, modal_id) {
  actionLink(
    ns(modal_id),
    label = shiny::icon("circle-info"),
    style = "color: #2980b9; margin-left: 6px; font-size: 1.1em; vertical-align: middle;"
  )
}

# Renders a citation as a linked footnote below a tool heading.
# text : citation string shown to the user
# url  : DOI or paper URL
citation_link <- function(text, url) {
  tags$p(
    tags$a(
      href   = url,
      target = "_blank",
      rel    = "noopener noreferrer",
      class  = "citation-link",
      shiny::icon("book-open"), " ", text
    )
  )
}

# Wraps content in a subtle description block.
tool_description <- function(...) {
  div(class = "tool-description", ...)
}

# Wraps the run button + status spinner in a consistent bottom action area.
action_row <- function(...) {
  div(class = "action-row", ...)
}

# ── Button loading state helpers (require shinyjs) ─────────────────────────────

# Grey out a button while an operation is running.
button_loading <- function(session, input_id) {
  shinyjs::disable(input_id)
}

# Re-enable a button once the operation completes or errors.
button_reset <- function(session, input_id) {
  shinyjs::enable(input_id)
}
