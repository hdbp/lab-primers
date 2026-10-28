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

source("R/idt.R")
check("expand_ids expands a range", identical(expand_ids("HP0008-HP0011,T3"),
      c("HP0008", "HP0009", "HP0010", "HP0011", "T3")))
check("expand_ids rejects mixed prefixes",
      inherits(try(expand_ids("A01-B03"), silent = TRUE), "try-error"))
check("select_primers adds the partner", identical(select_primers(ok, ids = "T001")$primer_id, c("T001", "T002")))
check("select_primers without partners", identical(select_primers(ok, ids = "T001", partners = FALSE)$primer_id, "T001"))
check("select_primers ANDs filters", identical(select_primers(ok, name = "_f$", target = "gene2")$primer_id, "T003"))
check("select_primers rejects unknown IDs",
      inherits(try(select_primers(ok, ids = "T999"), silent = TRUE), "try-error"))
tubes <- idt_tubes(ok[1:2, ], scale = "100nm", purification = "PAGE")
check("idt_tubes has IDT bulk-input columns", identical(names(tubes), c("Name", "Sequence", "Scale", "Purification")))
check("idt_tubes names carry the primer_id", identical(tubes$Name, c("T001_test_F", "T002_test_R")))
check("idt_tubes rejects unknown scale", inherits(try(idt_tubes(ok, scale = "50nm"), silent = TRUE), "try-error"))
big <- ok[rep(1, 100), ]
plates <- idt_plates(big)
check("idt_plates splits at 96 wells", identical(vapply(plates, nrow, 1L), c(`0` = 96L, `1` = 4L)))
check("idt_plates fills by column", identical(plates[[1]]$`Well Position`[c(1, 2, 9, 96)], c("A1", "B1", "A2", "H12")))
check("idt_plates fills by row", identical(idt_plates(big, "rows")[[1]]$`Well Position`[c(2, 13)], c("A2", "B1")))
xl <- tempfile(fileext = ".xlsx")
write_xlsx(list(Plate1 = plates[[1]], Plate2 = plates[[2]]), xl)
parts <- utils::unzip(xl, list = TRUE)$Name
check("write_xlsx writes one worksheet per plate",
      all(c("xl/worksheets/sheet1.xml", "xl/worksheets/sheet2.xml", "[Content_Types].xml") %in% parts))
sheet <- readLines(unz(xl, "xl/worksheets/sheet1.xml"), warn = FALSE)
check("write_xlsx sheet has headers and data",
      grepl("<t>Well Position</t>", sheet) && grepl("<t>T001_test_F</t>", sheet) && grepl('r="C97"', sheet))

hdr <- read_primers("primers.tsv")
check("primers.tsv has the expected columns", identical(names(hdr), PRIMER_COLUMNS))

message("All tests passed.")
