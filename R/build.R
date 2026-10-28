#!/usr/bin/env Rscript
# Validate primers.tsv, then regenerate the derived files in build/:
#   build/primers_annotated.tsv  - table + length, GC%, Tm, pair Tm difference
#   build/primers.fasta          - all primers as FASTA
# Usage: Rscript R/build.R [primers.tsv] [out_dir]

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "primer_utils.R"))

args <- commandArgs(TRUE)
path <- if (length(args) >= 1) args[1] else "primers.tsv"
out_dir <- if (length(args) >= 2) args[2] else "build"

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
ann <- ann[order(ann$primer_id), , drop = FALSE]
utils::write.table(ann, file.path(out_dir, "primers_annotated.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
write_fasta(ann, file.path(out_dir, "primers.fasta"))
message(sprintf("Wrote %d primers to %s/", nrow(ann), out_dir))
