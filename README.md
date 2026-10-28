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

## Tm calculation

Nearest-neighbour model (SantaLucia 1998 unified parameters) with SantaLucia's entropy salt correction: 50 mM Na⁺, 250 nM total primer, effective strand concentration CT/4 (Primer3's convention). The tests check it against Biopython's `Tm_NN` to within 0.1 °C. There's no Mg²⁺/dNTP correction, so values come out a few degrees below Primer3's PCR-buffer defaults. Use them to compare primers and to check that a pair is matched (`pair_tm_diff`), not as an annealing temperature. Sequences with ambiguity codes get no Tm.

## Using it from R

```r
source("R/primer_utils.R")
p <- read.delim("build/primers_annotated.tsv", colClasses = c(primer_id = "character"))
subset(p, target == "GAPDH" & species == "Homo sapiens")
```
