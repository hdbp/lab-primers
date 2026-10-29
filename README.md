# lab-primers

The lab's primer registry. One curated table, `primers.tsv`, is the source of truth; everything else is derived from it. Each row is a primer pair: the forward and reverse primer side by side.

## Layout

| Path | What it is |
|---|---|
| `primers.tsv` | The primer table. **Edit this.** |
| `R/primer_utils.R` | Validation and property functions (base R, no packages needed). |
| `R/validate.R` | Checks `primers.tsv`; exits non-zero on errors. |
| `R/build.R` | Validates, then builds the site into `site/` (not committed): the web app plus `site/build/primers_annotated.tsv` (the table with length, GC% and Tm of each primer and their Tm difference) and `site/build/primers.fasta` (every primer, two records per pair: `><pair_id>_<name>_F target\|species\|F`). |
| `web/template.html` | The web app; `R/build.R` fills in a snapshot of the primers. |
| `web/primer_core.js` | The app's validation, GC% and Tm, mirroring `R/primer_utils.R` (tests check they agree). |
| `web/github.js`, `web/assistant.js` | The app's GitHub sync and Claude assistant. |
| `R/html.R` | Writes the page. |
| `.github/workflows/` | `validate` runs the tests on every push and PR; `pages` builds and publishes the site on every push to `main`. |
| `R/export_idt.R`, `R/idt.R` | Select primers and write IDT order files (see [Ordering from IDT](#ordering-from-idt-command-line)). |
| `tests/` | Unit tests and fixtures. |

## Columns in `primers.tsv`

Tab-separated, one primer pair per row, all columns present (leave `notes` empty if there's nothing to say).

| Column | Required | Rules |
|---|---|---|
| `pair_id` | yes | Unique. Letters, digits, `_ . -` only (e.g. `HP001`). Never reuse an ID, even after a pair is retired. |
| `name` | yes | Human-readable name for the pair, e.g. `H1-0_qPCR`. Ordered primers are named `<pair_id>_<name>_F` and `_R`. |
| `target` | yes | Gene symbol or locus, e.g. `H1-0`, `GAPDH`, `chr7:5,530,000`. |
| `species` | yes | Binomial name, e.g. `Homo sapiens`, `Mus musculus`. |
| `forward` | yes | Forward primer, 5'→3', uppercase IUPAC DNA (`ACGTRYSWKMBDHVN`), no spaces. Put modifications in `notes`. |
| `reverse` | yes | Reverse primer, same rules. |
| `notes` | no | Free text. |

Length, GC% and Tm are **not** typed in by hand. `R/build.R` computes them, so they can't go stale.

## Adding primers

In the web app: **Add pair** (or a card from the Claude assistant), then **Review & open pull request**. See [The web app](#the-web-app).

With git:

```bash
# 1. add rows to primers.tsv (Excel/LibreOffice: save as tab-delimited text)
# 2. validate
Rscript R/validate.R
# 3. commit and open a pull request
git add primers.tsv
git commit -m "Add H1-0 qPCR primer pair"
```

Either way, CI runs the tests and validates the table. When the change is merged, the `pages` workflow rebuilds the site; nothing generated is committed.

### What validation checks

Errors (block the build): missing or extra columns, empty required fields (including a missing forward or reverse primer), duplicate IDs, invalid characters in a sequence.

Warnings: primers shorter than 10 nt, ambiguity codes (Tm is left blank), the same sequence used twice (in two pairs, or as both primers of one pair).

## The web app

The site is published with GitHub Pages at `https://<owner>.github.io/lab-primers/` (one-time setup: Settings → Pages → Source: **GitHub Actions**). It loads `primers.tsv` from `main` every time it opens, so it always shows the current table; if GitHub can't be reached it falls back to the copy built into the page. To try it locally, run `Rscript R/build.R` and serve `site/` (for example `python3 -m http.server --directory site`).

**Browse and order**

- **Search:** matches ID, name, target and both sequences. Several words must all match; use `"quoted phrases"`. Tick *Include notes* to search notes too.
- **Filter and sort:** filter by species; click a column header to sort.
- **Notes:** each pair's notes button shows its source paper, coordinates and validation status.
- **Select:** tick pairs. Both primers of each ticked pair are ordered. *Select shown* ticks everything the current search shows.
- **Order:** in the bottom bar, pick the IDT format, scale and purification, then **Download**. Formats and column layouts are the same as `R/export_idt.R` below. For tube text, **Copy** puts the lines on the clipboard to paste straight into IDT's bulk input.
- **Warnings:** pairs whose primers' Tm differs by more than 5 °C are highlighted; the bar warns if a plate has fewer than 24 oligos.

**Edit and sync to GitHub**

- **Add pair**, or **Edit** on a row, opens a form that shows length, GC% and Tm as you type and runs the same checks as `R/validate.R` (including duplicate IDs and sequences already in the table). Saved edits go into a draft, highlighted in the table, kept in your browser until you send or discard them.
- **Review & open pull request** shows every change, then creates a branch (`primer-app/…`) with the new `primers.tsv` and opens a pull request against `main`. If someone changed the table meanwhile, your draft is re-applied to the new version; a pair edited both here and there stops the save so nothing is overwritten.
- This needs a GitHub token, entered under **Settings** and stored only in that browser: a [fine-grained personal access token](https://github.com/settings/personal-access-tokens/new) for this repository only, with **Contents** and **Pull requests** set to read and write. Reading needs no token.

**Ask Claude**

The assistant searches papers and primer databases for published, validated pairs, or designs new ones, and shows each candidate as a card with Tm and GC% computed by the app and any problems flagged (duplicates of pairs already in the table, short primers, mismatched Tm). It is told never to call a pair published unless it read it in a source it can cite, and designed pairs are labelled unvalidated. **Review & add to draft** opens the pair in the edit form; nothing reaches GitHub until you open the pull request.

- It needs an Anthropic API key (Settings; stored only in that browser). Requests go directly from the browser to the Anthropic API and are billed to that key; set a spending limit on the key. A Claude.ai subscription can't be used from a web page.
- It uses Claude Opus 5.5 with web search and web fetch. It can't run in-silico PCR or BLAST, so check new pairs against the genome (UCSC In-Silico PCR, Primer-BLAST) before ordering; the cards' notes say what is still unchecked.

## Ordering from IDT (command line)

`R/export_idt.R` selects primer pairs from `primers.tsv` and writes a file to upload to IDT, with both primers of each pair.

```bash
# preview a selection (writes nothing)
Rscript R/export_idt.R --name=esrrb --list
# tubes: tab-separated text to paste or upload into IDT's bulk input
Rscript R/export_idt.R --ids=HP001-HP003,HP010
# the same as an Excel file, at 100 nmol with PAGE purification
Rscript R/export_idt.R --target="lefty|xist" --format=tubes-xlsx --scale=100nm --purification=PAGE
# a 96-well plate upload (Excel, one worksheet per plate)
Rscript R/export_idt.R --notes="ChIP target: Oct4" --format=plate
```

Selection: `--ids` takes pair IDs and ranges; `--name`, `--target` and `--notes` take case-insensitive regular expressions. Filters combine with AND, and with none every pair is selected.

| Format | File | Columns |
|---|---|---|
| `tubes` (default) | `.txt`, tab-separated | Name, Sequence, Scale, Purification |
| `tubes-xlsx` | `.xlsx` | Name, Sequence, Scale, Purification |
| `plate` | `.xlsx`, one sheet per 96-well plate | Well Position, Name, Sequence |

- Names are `<pair_id>_<name>_F` and `_R`, forward first, so tubes arriving from IDT can be matched back to the registry.
- Scale: `25nm` (default), `100nm`, `250nm`, `1um`, `5um`, `10um`. Purification: `STD` (default), `PAGE`, `HPLC`, `IEHPLC`, `RNASE`, `DUALHPLC`, `PAGEHPLC`. Scale and purification apply to tube orders; plates take those options in IDT's plate form.
- Plates fill down columns (A1, B1, … H1, A2) by default; use `--fill=rows` for A1, A2, … A12, B1. IDT may require a minimum number of oligos per plate (24 is commonly quoted), so the script warns about small plates.
- Output goes to `orders/idt_<format>_<date>.<ext>` unless you pass `--out`. `orders/` is git-ignored.
- The table is validated first; nothing is exported while it has errors.

These layouts follow IDT's bulk-input and plate-upload templates. If IDT changes its forms, check the current template on the order page.

## Tm calculation

Nearest-neighbour model (SantaLucia 1998 unified parameters) with SantaLucia's entropy salt correction: 50 mM Na⁺, 250 nM total primer, effective strand concentration CT/4 (Primer3's convention). The tests check it against Biopython's `Tm_NN` to within 0.1 °C. There's no Mg²⁺/dNTP correction, so values come out a few degrees below Primer3's PCR-buffer defaults. Use them to compare primers and to check that a pair is matched (`tm_diff`), not as an annealing temperature. Sequences with ambiguity codes get no Tm.

## Using it from R

```r
source("R/primer_utils.R")
p <- annotate_primers(read_primers("primers.tsv"))
subset(p, target == "GAPDH" & species == "Homo sapiens")
```
