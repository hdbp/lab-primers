# Selecting primers and writing IDT order files. Base R only; the .xlsx
# writer builds the workbook XML directly and zips it with utils::zip.

IDT_SCALES <- c("25nm", "100nm", "250nm", "1um", "5um", "10um")
IDT_PURIFICATIONS <- c("STD", "PAGE", "HPLC", "IEHPLC", "RNASE", "DUALHPLC", "PAGEHPLC")

# "HP0001-HP0004,HP0019" -> HP0001 HP0002 HP0003 HP0004 HP0019.
# A range needs both ends to share a prefix and zero-padding width.
expand_ids <- function(spec) {
  parts <- trimws(unlist(strsplit(spec, ",")))
  parts <- parts[nzchar(parts)]
  unlist(lapply(parts, function(p) {
    ends <- strsplit(p, "-", fixed = TRUE)[[1]]
    if (length(ends) != 2) return(p)
    m <- regmatches(ends, regexec("^(.*?)([0-9]+)$", ends))
    if (any(lengths(m) != 3) || m[[1]][2] != m[[2]][2] ||
        nchar(m[[1]][3]) != nchar(m[[2]][3])) {
      stop("can't read ID range '", p, "' (use e.g. HP0001-HP0006)", call. = FALSE)
    }
    nums <- seq(as.integer(m[[1]][3]), as.integer(m[[2]][3]))
    sprintf("%s%0*d", m[[1]][2], nchar(m[[1]][3]), nums)
  }), use.names = FALSE)
}

# Rows matching every filter given (NULL filters are ignored). Patterns are
# case-insensitive regular expressions. With partners = TRUE, the partner of
# each selected primer is added so pairs are never ordered half-complete.
select_primers <- function(df, ids = NULL, name = NULL, target = NULL,
                           notes = NULL, partners = TRUE) {
  keep <- rep(TRUE, nrow(df))
  if (!is.null(ids)) {
    want <- expand_ids(ids)
    unknown <- setdiff(want, df$primer_id)
    if (length(unknown)) stop("unknown primer_id: ", paste(unknown, collapse = ", "), call. = FALSE)
    keep <- keep & df$primer_id %in% want
  }
  if (!is.null(name)) keep <- keep & grepl(name, df$name, ignore.case = TRUE)
  if (!is.null(target)) keep <- keep & grepl(target, df$target, ignore.case = TRUE)
  if (!is.null(notes)) keep <- keep & grepl(notes, df$notes, ignore.case = TRUE)
  if (partners) keep <- keep | df$primer_id %in% df$partner_id[keep]
  df[keep, , drop = FALSE]
}

idt_names <- function(df) paste(df$primer_id, df$name, sep = "_")

idt_tubes <- function(df, scale = "25nm", purification = "STD") {
  if (!scale %in% IDT_SCALES) {
    stop("scale must be one of: ", paste(IDT_SCALES, collapse = ", "), call. = FALSE)
  }
  if (!purification %in% IDT_PURIFICATIONS) {
    stop("purification must be one of: ", paste(IDT_PURIFICATIONS, collapse = ", "), call. = FALSE)
  }
  data.frame(Name = idt_names(df), Sequence = df$sequence, Scale = scale,
             Purification = purification, check.names = FALSE)
}

# Assign 96-well positions (A1, B1, ... when filling by column). Returns one
# data frame per plate.
idt_plates <- function(df, fill = c("cols", "rows")) {
  fill <- match.arg(fill)
  wells <- if (fill == "cols") {
    paste0(rep(LETTERS[1:8], times = 12), rep(1:12, each = 8))
  } else {
    paste0(rep(LETTERS[1:8], each = 12), rep(1:12, times = 8))
  }
  plate <- (seq_len(nrow(df)) - 1) %/% 96
  lapply(split(seq_len(nrow(df)), plate), function(i) {
    data.frame(`Well Position` = wells[seq_along(i)], Name = idt_names(df)[i],
               Sequence = df$sequence[i], check.names = FALSE)
  })
}

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

sheet_xml <- function(tab) {
  rows <- rbind(names(tab), as.matrix(tab))
  body <- vapply(seq_len(nrow(rows)), function(r) {
    cells <- sprintf('<c r="%s%d" t="inlineStr"><is><t>%s</t></is></c>',
                     LETTERS[seq_len(ncol(rows))], r, xml_escape(rows[r, ]))
    sprintf('<row r="%d">%s</row>', r, paste(cells, collapse = ""))
  }, "")
  paste0('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
         '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
         '<sheetData>', paste(body, collapse = ""), '</sheetData></worksheet>')
}

# Minimal .xlsx: one worksheet per element of `sheets` (a named list of data
# frames), text cells only.
write_xlsx <- function(sheets, path) {
  tmp <- tempfile("xlsx")
  dir.create(file.path(tmp, "_rels"), recursive = TRUE)
  dir.create(file.path(tmp, "xl", "_rels"), recursive = TRUE)
  dir.create(file.path(tmp, "xl", "worksheets"))
  n <- length(sheets)
  ns <- "http://schemas.openxmlformats.org"
  writeLines(paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="', ns, '/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    paste(sprintf('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>', seq_len(n)), collapse = ""),
    '</Types>'), file.path(tmp, "[Content_Types].xml"))
  writeLines(paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="', ns, '/package/2006/relationships">',
    '<Relationship Id="rId1" Type="', ns, '/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
    '</Relationships>'), file.path(tmp, "_rels", ".rels"))
  writeLines(paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<workbook xmlns="', ns, '/spreadsheetml/2006/main" xmlns:r="', ns, '/officeDocument/2006/relationships"><sheets>',
    paste(sprintf('<sheet name="%s" sheetId="%d" r:id="rId%d"/>', xml_escape(names(sheets)), seq_len(n), seq_len(n)), collapse = ""),
    '</sheets></workbook>'), file.path(tmp, "xl", "workbook.xml"))
  writeLines(paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="', ns, '/package/2006/relationships">',
    paste(sprintf('<Relationship Id="rId%d" Type="%s/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>', seq_len(n), ns, seq_len(n)), collapse = ""),
    '</Relationships>'), file.path(tmp, "xl", "_rels", "workbook.xml.rels"))
  for (i in seq_len(n)) {
    writeLines(sheet_xml(sheets[[i]]), file.path(tmp, "xl", "worksheets", sprintf("sheet%d.xml", i)))
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (file.exists(path)) file.remove(path)
  old <- setwd(tmp)
  on.exit({ setwd(old); unlink(tmp, recursive = TRUE) })
  files <- c("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml",
             "xl/_rels/workbook.xml.rels", sprintf("xl/worksheets/sheet%d.xml", seq_len(n)))
  status <- utils::zip(path, files, flags = "-X -q")
  if (status != 0) stop("zip failed writing ", path, call. = FALSE)
  invisible(path)
}
