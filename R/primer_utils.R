# Core primer utilities: reading, validation, and derived properties
# (length, GC%, Tm). Base R only, so the repo runs anywhere R is installed.

PRIMER_COLUMNS <- c("primer_id", "name", "sequence", "target", "species",
                    "direction", "partner_id", "notes")

IUPAC_DNA <- "ACGTRYSWKMBDHVN"

# SantaLucia (1998) unified nearest-neighbour parameters.
# dH in kcal/mol, dS in cal/(K*mol). Keys are 5'->3' dinucleotides on the
# top strand; reverse-complement keys share the same values.
NN_PARAMS <- local({
  p <- rbind(
    AA = c(-7.9, -22.2), AT = c(-7.2, -20.4), TA = c(-7.2, -21.3),
    CA = c(-8.5, -22.7), GT = c(-8.4, -22.4), CT = c(-7.8, -21.0),
    GA = c(-8.2, -22.2), CG = c(-10.6, -27.2), GC = c(-9.8, -24.4),
    GG = c(-8.0, -19.9)
  )
  rc <- c(TT = "AA", AC = "GT", TG = "CA", AG = "CT", TC = "GA", CC = "GG")
  p <- rbind(p, p[rc, , drop = FALSE])
  rownames(p)[(nrow(p) - length(rc) + 1):nrow(p)] <- names(rc)
  colnames(p) <- c("dH", "dS")
  p
})

read_primers <- function(path) {
  df <- utils::read.delim(path, colClasses = "character", na.strings = character(0),
                          quote = "", comment.char = "", check.names = FALSE,
                          strip.white = TRUE)
  df
}

gc_percent <- function(seq) {
  vapply(seq, function(s) {
    b <- strsplit(s, "")[[1]]
    if (any(!b %in% c("A", "C", "G", "T"))) {
      # S (G/C) counts as GC; W (A/T) as not; other ambiguity codes -> NA
      if (any(!b %in% c("A", "C", "G", "T", "S", "W"))) return(NA_real_)
    }
    round(100 * sum(b %in% c("G", "C", "S")) / length(b), 1)
  }, numeric(1), USE.NAMES = FALSE)
}

# Nearest-neighbour Tm (SantaLucia 1998) with SantaLucia's entropy salt
# correction. primer_nM is the total oligo concentration CT and the
# effective strand concentration is CT/4 (Primer3 convention; equals
# Biopython Tm_NN(nn_table=DNA_NN3, saltcorr=5, dnac1 = dnac2 = CT/2)).
# No Mg2+/dNTP correction, so values run a few degrees below Primer3's
# PCR-buffer defaults; use them to compare primers, not as absolute Ta.
# Returns NA for sequences containing ambiguity codes.
tm_nn <- function(seq, na_mM = 50, primer_nM = 250) {
  vapply(seq, function(s) {
    if (!grepl("^[ACGT]+$", s) || nchar(s) < 2) return(NA_real_)
    b <- strsplit(s, "")[[1]]
    n <- length(b)
    pairs <- paste0(b[-n], b[-1])
    dH <- sum(NN_PARAMS[pairs, "dH"])
    dS <- sum(NN_PARAMS[pairs, "dS"])
    for (end in c(b[1], b[n])) {
      if (end %in% c("G", "C")) { dH <- dH + 0.1; dS <- dS - 2.8 }
      else                      { dH <- dH + 2.3; dS <- dS + 4.1 }
    }
    dS <- dS + 0.368 * (n - 1) * log(na_mM / 1000)
    ct <- primer_nM * 1e-9
    tm <- (1000 * dH) / (dS + 1.987 * log(ct / 4)) - 273.15
    round(tm, 1)
  }, numeric(1), USE.NAMES = FALSE)
}

# Returns a list(errors = character, warnings = character).
validate_primers <- function(df) {
  errors <- character(0)
  warnings <- character(0)
  err <- function(...) errors <<- c(errors, paste0(...))
  warn <- function(...) warnings <<- c(warnings, paste0(...))

  missing_cols <- setdiff(PRIMER_COLUMNS, names(df))
  extra_cols <- setdiff(names(df), PRIMER_COLUMNS)
  if (length(missing_cols)) err("missing columns: ", paste(missing_cols, collapse = ", "))
  if (length(extra_cols)) err("unexpected columns: ", paste(extra_cols, collapse = ", "))
  if (length(missing_cols) || length(extra_cols)) return(list(errors = errors, warnings = warnings))
  if (nrow(df) == 0) return(list(errors = errors, warnings = warnings))

  row <- seq_len(nrow(df)) + 1L  # file line number (header is line 1)
  where <- function(i) sprintf("line %d (%s)", row[i], df$primer_id[i])

  for (col in c("primer_id", "name", "sequence", "target", "species", "direction")) {
    for (i in which(df[[col]] == "")) err(where(i), ": empty ", col)
  }

  bad_id <- which(df$primer_id != "" & !grepl("^[A-Za-z0-9_.-]+$", df$primer_id))
  for (i in bad_id) err(where(i), ": primer_id may only contain letters, digits, '_', '.', '-'")
  for (id in unique(df$primer_id[duplicated(df$primer_id) & df$primer_id != ""])) {
    err("duplicate primer_id '", id, "' on lines ",
        paste(row[df$primer_id == id], collapse = ", "))
  }

  bad_seq <- which(df$sequence != "" & !grepl(paste0("^[", IUPAC_DNA, "]+$"), df$sequence))
  for (i in bad_seq) {
    err(where(i), ": sequence must be uppercase IUPAC DNA (", IUPAC_DNA,
        ") with no spaces or modifications")
  }
  short <- which(df$sequence != "" & nchar(df$sequence) < 10)
  for (i in short) warn(where(i), ": sequence is only ", nchar(df$sequence[i]), " nt")
  ambig <- which(grepl("[^ACGT]", df$sequence) & !seq_len(nrow(df)) %in% bad_seq &
                   df$sequence != "")
  for (i in ambig) warn(where(i), ": sequence contains ambiguity codes; Tm not computed")

  for (i in which(df$direction != "" & !df$direction %in% c("F", "R"))) {
    err(where(i), ": direction must be 'F' or 'R', got '", df$direction[i], "'")
  }

  dup_seq <- unique(df$sequence[duplicated(df$sequence) & df$sequence != ""])
  for (s in dup_seq) {
    warn("identical sequence shared by: ",
         paste(df$primer_id[df$sequence == s], collapse = ", "))
  }

  ids <- df$primer_id
  for (i in which(df$partner_id != "")) {
    p <- df$partner_id[i]
    j <- match(p, ids)
    if (is.na(j)) { err(where(i), ": partner_id '", p, "' does not exist"); next }
    if (j == i) { err(where(i), ": primer lists itself as partner"); next }
    if (df$partner_id[j] != ids[i]) {
      err(where(i), ": partner '", p, "' does not list '", ids[i], "' back as its partner")
    }
    if (df$direction[i] == df$direction[j] && df$direction[i] %in% c("F", "R")) {
      err(where(i), ": partner '", p, "' has the same direction (", df$direction[i], ")")
    }
    if (df$target[i] != df$target[j]) {
      warn(where(i), ": partner '", p, "' has a different target ('", df$target[j], "')")
    }
    if (df$species[i] != df$species[j]) {
      warn(where(i), ": partner '", p, "' has a different species ('", df$species[j], "')")
    }
  }

  list(errors = errors, warnings = warnings)
}

annotate_primers <- function(df, na_mM = 50, primer_nM = 250) {
  df$length <- nchar(df$sequence)
  df$gc_pct <- gc_percent(df$sequence)
  df$tm_c <- tm_nn(df$sequence, na_mM = na_mM, primer_nM = primer_nM)
  partner_tm <- df$tm_c[match(df$partner_id, df$primer_id)]
  df$pair_tm_diff <- round(abs(df$tm_c - partner_tm), 1)
  df
}

write_fasta <- function(df, path) {
  hdr <- sprintf(">%s %s|%s|%s|%s", df$primer_id, df$name, df$target,
                 df$species, df$direction)
  writeLines(c(rbind(hdr, df$sequence)), path)
}
