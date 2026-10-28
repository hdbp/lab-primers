#!/usr/bin/env Rscript
# Run from the repo root: Rscript tests/run_tests.R
source("R/primer_utils.R")

check <- function(desc, cond) {
  if (!isTRUE(cond)) stop("FAILED: ", desc, call. = FALSE)
  message("ok - ", desc)
}

# Reference Tm values from Biopython 1.88:
#   MeltingTemp.Tm_NN(seq, nn_table=DNA_NN3, saltcorr=5, Na=50, dnac1=125, dnac2=125)
#   i.e. 250 nM total oligo, effective concentration CT/4
ref <- read.delim("tests/tm_reference.tsv", colClasses = c("character", "numeric"))
got <- tm_nn(ref$sequence)
check("tm_nn matches Biopython within 0.1 C", all(abs(got - ref$tm_c) <= 0.1))
check("tm_nn is NA for ambiguous bases", is.na(tm_nn("ACGTNACGTACGTACGT")))

check("gc_percent basic", identical(gc_percent(c("GGCC", "AATT", "ACGT")), c(100, 0, 50)))
check("gc_percent counts S and W", identical(gc_percent("SSWW"), 50))
check("gc_percent NA for N", is.na(gc_percent("ACGN")))

ok <- read_primers("tests/fixture_valid.tsv")
res <- validate_primers(ok)
check("valid fixture has no errors", length(res$errors) == 0)
check("valid fixture warns about ambiguity", any(grepl("ambiguity", res$warnings)))
ann <- annotate_primers(ok)
check("annotation adds length", identical(ann$length, c(20L, 20L, 20L)))
check("pair Tm diff only for paired primers",
      !is.na(ann$pair_tm_diff[1]) && is.na(ann$pair_tm_diff[3]))

bad <- validate_primers(read_primers("tests/fixture_invalid.tsv"))$errors
expect <- c("sequence must be uppercase", "duplicate primer_id 'B002'",
            "direction must be 'F' or 'R'", "partner_id 'B999' does not exist",
            "does not list 'B004' back", "has the same direction")
for (e in expect) check(paste("invalid fixture reports:", e), any(grepl(e, bad, fixed = TRUE)))

hdr <- read_primers("primers.tsv")
check("primers.tsv has the expected columns", identical(names(hdr), PRIMER_COLUMNS))

message("All tests passed.")
