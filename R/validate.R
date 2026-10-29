#!/usr/bin/env Rscript
# Validate the primer table. Exits non-zero if any errors are found.
# Usage: Rscript R/validate.R [primers.tsv]

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "primer_utils.R"))

args <- commandArgs(TRUE)
path <- if (length(args)) args[1] else "primers.tsv"

df <- read_primers(path)
res <- validate_primers(df)

for (w in res$warnings) message("WARNING: ", w)
for (e in res$errors) message("ERROR: ", e)
message(sprintf("%s: %d primer pairs, %d errors, %d warnings",
                path, nrow(df), length(res$errors), length(res$warnings)))
if (length(res$errors)) quit(status = 1)
