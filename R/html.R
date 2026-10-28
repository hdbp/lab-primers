# Writes the interactive primer page: the annotated table as JSON, dropped
# into web/template.html. Base R only.

json_string <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  x <- gsub("</", "<\\/", x, fixed = TRUE)   # can't close the <script> block
  paste0("\"", x, "\"")
}

json_value <- function(v) {
  if (is.numeric(v)) ifelse(is.na(v), "null", format(v, trim = TRUE, scientific = FALSE, drop0trailing = TRUE))
  else json_string(v)
}

# One JSON object per row, one row per line so diffs of the page stay readable.
primers_json <- function(df) {
  vals <- lapply(df, json_value)
  keys <- json_string(names(df))
  rows <- vapply(seq_len(nrow(df)), function(i) {
    paste0("{", paste0(keys, ":", vapply(vals, `[`, "", i), collapse = ","), "}")
  }, "")
  if (!length(rows)) return("[]")
  paste0("[\n", paste(rows, collapse = ",\n"), "\n]")
}

write_html <- function(df, template, path) {
  page <- paste(readLines(template, warn = FALSE), collapse = "\n")
  marker <- "/*__PRIMERS_JSON__*/[]"
  if (!grepl(marker, page, fixed = TRUE)) stop("template has no ", marker, call. = FALSE)
  page <- sub(marker, primers_json(df), page, fixed = TRUE)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(page, path, useBytes = TRUE)
}
