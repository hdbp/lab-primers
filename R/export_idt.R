#!/usr/bin/env Rscript
# Select primers from primers.tsv and write an IDT order file.
#
# Usage: Rscript R/export_idt.R [options]
#   --ids=HP0001-HP0006,HP0019   primer IDs; ranges allowed
#   --name=REGEX                 match on name     (case-insensitive)
#   --target=REGEX               match on target
#   --notes=REGEX                match on notes
#   --no-partners                don't add the partner of each selected primer
#   --format=tubes|tubes-xlsx|plate
#                                tubes:      tab-separated text for IDT bulk input (default)
#                                tubes-xlsx: the same as an Excel file
#                                plate:      Excel file for IDT 96-well plate upload
#   --scale=25nm                 tubes: 25nm 100nm 250nm 1um 5um 10um
#   --purification=STD           tubes: STD PAGE HPLC IEHPLC RNASE DUALHPLC PAGEHPLC
#   --fill=cols|rows             plate: fill A1,B1,... (cols, default) or A1,A2,... (rows)
#   --out=FILE                   default: orders/idt_<date>.<txt|xlsx>
#   --table=primers.tsv          primer table to read
#   --list                       show the selection; write nothing
# Filters combine with AND; with no filter, every primer is selected.

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "primer_utils.R"))
source(file.path(script_dir, "idt.R"))

args <- commandArgs(TRUE)
opt <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[length(hit)]) else default
}
flag <- function(key) paste0("--", key) %in% args

known <- c("ids", "name", "target", "notes", "format", "scale", "purification",
           "fill", "out", "table", "no-partners", "list")
bad <- args[!sub("=.*$", "", sub("^--", "", args)) %in% known | !startsWith(args, "--")]
if (length(bad)) {
  message("unknown option: ", paste(bad, collapse = " "), " (see the top of R/export_idt.R)")
  quit(status = 1)
}

df <- read_primers(opt("table", "primers.tsv"))
res <- validate_primers(df)
for (e in res$errors) message("ERROR: ", e)
if (length(res$errors)) {
  message("Not exporting: fix the errors above first.")
  quit(status = 1)
}

sel <- tryCatch(
  select_primers(df, ids = opt("ids"), name = opt("name"), target = opt("target"),
                 notes = opt("notes"), partners = !flag("no-partners")),
  error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })
if (!nrow(sel)) {
  message("No primers match.")
  quit(status = 1)
}
sel <- sel[order(sel$primer_id), , drop = FALSE]

if (flag("list")) {
  print(sel[, c("primer_id", "name", "sequence", "direction", "partner_id")], row.names = FALSE)
  message(sprintf("%d primers selected.", nrow(sel)))
  quit(status = 0)
}

fmt <- opt("format", "tubes")
if (!fmt %in% c("tubes", "tubes-xlsx", "plate")) {
  message("--format must be tubes, tubes-xlsx or plate")
  quit(status = 1)
}
ext <- if (fmt == "tubes") "txt" else "xlsx"
out <- opt("out", file.path("orders", sprintf("idt_%s_%s.%s", fmt, format(Sys.Date()), ext)))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

if (fmt == "plate") {
  plates <- idt_plates(sel, fill = opt("fill", "cols"))
  names(plates) <- sprintf("Plate%d", seq_along(plates))
  write_xlsx(plates, out)
  message(sprintf("Wrote %d primers on %d plate(s) to %s", nrow(sel), length(plates), out))
  if (nrow(sel) %% 96 && nrow(sel) %% 96 < 24) {
    message("Note: the last plate has fewer than 24 oligos; IDT plate orders may need more.")
  }
} else {
  tubes <- tryCatch(idt_tubes(sel, scale = opt("scale", "25nm"),
                              purification = opt("purification", "STD")),
                    error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })
  if (fmt == "tubes") {
    utils::write.table(tubes, out, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    write_xlsx(list(Tubes = tubes), out)
  }
  message(sprintf("Wrote %d primers (%s, %s) to %s", nrow(sel), tubes$Scale[1],
                  tubes$Purification[1], out))
}
