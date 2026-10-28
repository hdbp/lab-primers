# lab-primers

The lab's primer registry. One curated table, `primers.tsv`, is the source of truth; everything else is derived from it.

## Layout

| Path | What it is |
|---|---|
| `primers.tsv` | The primer table. **Edit this.** |
| `build/primers_annotated.tsv` | Generated: the table plus length, GC%, Tm, and the Tm difference to the partner primer. |
| `build/primers.fasta` | Generated: all primers as FASTA (`>primer_id name\|target\|species\|direction`). |
| `R/primer_utils.R` | Validation and property functions (base R, no packages needed). |
| `R/validate.R` | Checks `primers.tsv`; exits non-zero on errors. |
| `R/build.R` | Validates, then regenerates `build/`. |
| `R/export_idt.R`, `R/idt.R` | Select primers and write IDT order files (see [Ordering from IDT](#ordering-from-idt)). |
| `tests/` | Unit tests and fixtures. |

## Columns in `primers.tsv`

Tab-separated, one primer per row, all columns present (leave optional ones empty).

| Column | Required | Rules |
|---|---|---|
| `primer_id` | yes | Unique. Letters, digits, `_ . -` only (e.g. `HP0001`). Never reuse an ID, even after a primer is retired. |
| `name` | yes | Human-readable name, e.g. `H1-0_qPCR_F`. |
| `sequence` | yes | 5'→3', uppercase IUPAC DNA (`ACGTRYSWKMBDHVN`), no spaces. Put modifications in `notes`. |
| `target` | yes | Gene symbol or locus, e.g. `H1-0`, `GAPDH`, `chr7:5,530,000`. |
| `species` | yes | Binomial name, e.g. `Homo sapiens`, `Mus musculus`. |
| `direction` | yes | `F` or `R`. |
| `partner_id` | no | `primer_id` of the pair partner. Must be reciprocal and opposite direction. |
| `notes` | no | Free text. |

`length`, `gc_pct`, and `tm_c` are **not** typed in by hand. `R/build.R` computes them, so they can't go stale.

## Adding primers

```bash
# 1. add rows to primers.tsv (Excel/LibreOffice: save as tab-delimited text)
# 2. validate and regenerate derived files
Rscript R/validate.R
Rscript R/build.R
# 3. commit both primers.tsv and build/
git add primers.tsv build/
git commit -m "Add H1-0 qPCR primers"
```

CI runs the tests, validates the table, and fails if `build/` was not regenerated.

### What validation checks

Errors (block the build): missing or extra columns, empty required fields, duplicate IDs, invalid characters in a sequence, `direction` not `F`/`R`, a `partner_id` that doesn't exist or isn't reciprocal, partners with the same direction.

Warnings: sequences shorter than 10 nt, ambiguity codes (Tm is left blank), identical sequences under different IDs, partners with a different target or species.

## Ordering from IDT

`R/export_idt.R` selects primers from `primers.tsv` and writes a file to upload to IDT.

```bash
# preview a selection (writes nothing)
Rscript R/export_idt.R --name=esrrb --list
# tubes: tab-separated text to paste or upload into IDT's bulk input
Rscript R/export_idt.R --ids=HP0001-HP0006,HP0019
# the same as an Excel file, at 100 nmol with PAGE purification
Rscript R/export_idt.R --target="lefty|xist" --format=tubes-xlsx --scale=100nm --purification=PAGE
# a 96-well plate upload (Excel, one worksheet per plate)
Rscript R/export_idt.R --notes="ChIP target: Oct4" --format=plate
```

Selection: `--ids` takes IDs and ranges; `--name`, `--target` and `--notes` take case-insensitive regular expressions. Filters combine with AND, and with none every primer is selected. The partner of each selected primer is added automatically so pairs are never ordered half-complete; turn that off with `--no-partners`.

| Format | File | Columns |
|---|---|---|
| `tubes` (default) | `.txt`, tab-separated | Name, Sequence, Scale, Purification |
| `tubes-xlsx` | `.xlsx` | Name, Sequence, Scale, Purification |
| `plate` | `.xlsx`, one sheet per 96-well plate | Well Position, Name, Sequence |

- Names are `<primer_id>_<name>`, so tubes arriving from IDT can be matched back to the registry.
- Scale: `25nm` (default), `100nm`, `250nm`, `1um`, `5um`, `10um`. Purification: `STD` (default), `PAGE`, `HPLC`, `IEHPLC`, `RNASE`, `DUALHPLC`, `PAGEHPLC`. Scale and purification apply to tube orders; plates take those options in IDT's plate form.
- Plates fill down columns (A1, B1, … H1, A2) by default; use `--fill=rows` for A1, A2, … A12, B1. IDT may require a minimum number of oligos per plate (24 is commonly quoted), so the script warns about small plates.
- Output goes to `orders/idt_<format>_<date>.<ext>` unless you pass `--out`. `orders/` is git-ignored.
- The table is validated first; nothing is exported while it has errors.

These layouts follow IDT's bulk-input and plate-upload templates. If IDT changes its forms, check the current template on the order page.

## Tm calculation

Nearest-neighbour model (SantaLucia 1998 unified parameters) with SantaLucia's entropy salt correction: 50 mM Na⁺, 250 nM total primer, effective strand concentration CT/4 (Primer3's convention). The tests check it against Biopython's `Tm_NN` to within 0.1 °C. There's no Mg²⁺/dNTP correction, so values come out a few degrees below Primer3's PCR-buffer defaults. Use them to compare primers and to check that a pair is matched (`pair_tm_diff`), not as an annealing temperature. Sequences with ambiguity codes get no Tm.

## Using it from R

```r
source("R/primer_utils.R")
p <- read.delim("build/primers_annotated.tsv", colClasses = c(primer_id = "character"))
subset(p, target == "GAPDH" & species == "Homo sapiens")
```
