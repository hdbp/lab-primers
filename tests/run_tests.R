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
check("annotation adds length of each primer",
      identical(ann$forward_length, c(20L, 20L)) && identical(ann$reverse_length, c(20L, 20L)))
check("Tm diff is NA when a primer has no Tm",
      !is.na(ann$tm_diff[1]) && is.na(ann$tm_diff[2]) && !is.na(ann$reverse_tm_c[2]))
check("pair_oligos lists forward then reverse",
      identical(pair_oligos(ok[1, ])$oligo, c("T001_test_F", "T001_test_R")))
dup <- validate_primers(transform(ok[1, ], reverse = forward))$warnings
check("warns when forward and reverse are identical", any(grepl("T001 forward, T001 reverse", dup, fixed = TRUE)))

bad <- validate_primers(read_primers("tests/fixture_invalid.tsv"))$errors
expect <- c("forward must be uppercase", "duplicate pair_id 'B002'", "(B004): empty reverse")
for (e in expect) check(paste("invalid fixture reports:", e), any(grepl(e, bad, fixed = TRUE)))

source("R/idt.R")
check("expand_ids expands a range", identical(expand_ids("HP008-HP011,T3"),
      c("HP008", "HP009", "HP010", "HP011", "T3")))
check("expand_ids rejects mixed prefixes",
      inherits(try(expand_ids("A01-B03"), silent = TRUE), "try-error"))
check("select_pairs by ID", identical(select_pairs(ok, ids = "T001")$pair_id, "T001"))
check("select_pairs ANDs filters", identical(select_pairs(ok, name = "^deg", target = "gene2")$pair_id, "T002"))
check("select_pairs ANDs filters to nothing", nrow(select_pairs(ok, name = "^deg", target = "gene1")) == 0)
check("select_pairs rejects unknown IDs",
      inherits(try(select_pairs(ok, ids = "T999"), silent = TRUE), "try-error"))
tubes <- idt_tubes(ok[1, ], scale = "100nm", purification = "PAGE")
check("idt_tubes has IDT bulk-input columns", identical(names(tubes), c("Name", "Sequence", "Scale", "Purification")))
check("idt_tubes orders both primers of a pair",
      identical(tubes$Name, c("T001_test_F", "T001_test_R")) &&
      identical(tubes$Sequence, c(ok$forward[1], ok$reverse[1])))
check("idt_tubes rejects unknown scale", inherits(try(idt_tubes(ok, scale = "50nm"), silent = TRUE), "try-error"))
big <- ok[rep(1, 50), ]   # 100 oligos
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

source("R/html.R")
check("json_string escapes quotes, backslashes, tabs and </script>",
      identical(json_string('a"b\\c\td</script>'), '"a\\"b\\\\c\\td<\\/script>"'))
check("primers_json writes null for missing numbers",
      identical(primers_json(data.frame(id = "X1", tm = NA_real_)), '[\n{"id":"X1","tm":null}\n]'))
check("primers_json of no rows is []", identical(primers_json(ok[0, ]), "[]"))
page <- tempfile(fileext = ".html")
write_html(annotate_primers(ok), "web/template.html", page)
page_txt <- paste(readLines(page, warn = FALSE), collapse = "\n")
check("page embeds the primer table",
      grepl('"pair_id":"T001"', page_txt, fixed = TRUE) && !grepl("__PRIMERS_JSON__", page_txt, fixed = TRUE))

hdr <- read_primers("primers.tsv")
check("primers.tsv has the expected columns", identical(names(hdr), PRIMER_COLUMNS))

message("All tests passed.")
