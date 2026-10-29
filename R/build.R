#!/usr/bin/env Rscript
# Validate primers.tsv, then build the site in site/ (not committed; the Pages
# workflow runs this on every push to main):
#   site/build/primers_annotated.tsv  - table + length, GC% and Tm of each primer, Tm difference
#   site/build/primers.fasta          - every primer (both of each pair) as FASTA
#   site/index.html + web/*.js        - the web app (search, edit, IDT export, Claude assistant)
# Usage: Rscript R/build.R [primers.tsv] [out_dir] [page]

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "primer_utils.R"))
source(file.path(script_dir, "html.R"))

args <- commandArgs(TRUE)
path <- if (length(args) >= 1) args[1] else "primers.tsv"
out_dir <- if (length(args) >= 2) args[2] else file.path("site", "build")
page <- if (length(args) >= 3) args[3] else file.path("site", "index.html")

df <- read_primers(path)
res <- validate_primers(df)
for (w in res$warnings) message("WARNING: ", w)
for (e in res$errors) message("ERROR: ", e)
if (length(res$errors)) {
  message("Not building: fix the errors above first.")
  quit(status = 1)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ann <- annotate_primers(df)
ann <- ann[order(ann$pair_id), , drop = FALSE]
utils::write.table(ann, file.path(out_dir, "primers_annotated.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
write_fasta(ann, file.path(out_dir, "primers.fasta"))
web <- file.path(script_dir, "..", "web")
write_html(ann, file.path(web, "template.html"), page)
invisible(file.copy(list.files(web, pattern = "\\.js$", full.names = TRUE), dirname(page), overwrite = TRUE))
message(sprintf("Wrote %d primer pairs to %s/ and %s", nrow(ann), out_dir, page))
