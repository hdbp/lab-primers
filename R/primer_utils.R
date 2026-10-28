# Core primer utilities: reading, validation, and derived properties
# (length, GC%, Tm). Base R only, so the repo runs anywhere R is installed.

# One row per primer pair.
PRIMER_COLUMNS <- c("pair_id", "name", "target", "species", "forward", "reverse", "notes")
SEQ_COLUMNS <- c("forward", "reverse")

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
  where <- function(i) sprintf("line %d (%s)", row[i], df$pair_id[i])

  for (col in setdiff(PRIMER_COLUMNS, "notes")) {
    for (i in which(df[[col]] == "")) err(where(i), ": empty ", col)
  }

  bad_id <- which(df$pair_id != "" & !grepl("^[A-Za-z0-9_.-]+$", df$pair_id))
  for (i in bad_id) err(where(i), ": pair_id may only contain letters, digits, '_', '.', '-'")
  for (id in unique(df$pair_id[duplicated(df$pair_id) & df$pair_id != ""])) {
    err("duplicate pair_id '", id, "' on lines ",
        paste(row[df$pair_id == id], collapse = ", "))
  }

  for (col in SEQ_COLUMNS) {
    s <- df[[col]]
    bad_seq <- which(s != "" & !grepl(paste0("^[", IUPAC_DNA, "]+$"), s))
    for (i in bad_seq) {
      err(where(i), ": ", col, " must be uppercase IUPAC DNA (", IUPAC_DNA,
          ") with no spaces or modifications")
    }
    for (i in which(s != "" & nchar(s) < 10)) warn(where(i), ": ", col, " is only ", nchar(s[i]), " nt")
    ambig <- which(s != "" & grepl("[^ACGT]", s) & !seq_len(nrow(df)) %in% bad_seq)
    for (i in ambig) warn(where(i), ": ", col, " contains ambiguity codes; Tm not computed")
  }

  # The same oligo under two pairs (or as both primers of one pair).
  oligo <- c(df$forward, df$reverse)
  label <- c(paste0(df$pair_id, " forward"), paste0(df$pair_id, " reverse"))
  for (s in unique(oligo[duplicated(oligo) & oligo != ""])) {
    warn("identical sequence shared by: ", paste(label[oligo == s], collapse = ", "))
  }

  list(errors = errors, warnings = warnings)
}

annotate_primers <- function(df, na_mM = 50, primer_nM = 250) {
  for (col in SEQ_COLUMNS) {
    df[[paste0(col, "_length")]] <- nchar(df[[col]])
    df[[paste0(col, "_gc_pct")]] <- gc_percent(df[[col]])
    df[[paste0(col, "_tm_c")]] <- tm_nn(df[[col]], na_mM = na_mM, primer_nM = primer_nM)
  }
  df$tm_diff <- round(abs(df$forward_tm_c - df$reverse_tm_c), 1)
  df
}

# One oligo per row (forward then reverse for each pair), named
# <pair_id>_<name>_F / _R. This is what gets ordered and written to FASTA.
pair_oligos <- function(df) {
  n <- nrow(df)
  i <- rep(seq_len(n), each = 2)
  dir <- rep(c("F", "R"), times = n)
  data.frame(oligo = paste(df$pair_id[i], df$name[i], dir, sep = "_"),
             sequence = c(rbind(df$forward, df$reverse)),
             pair_id = df$pair_id[i], direction = dir,
             target = df$target[i], species = df$species[i])
}

write_fasta <- function(df, path) {
  o <- pair_oligos(df)
  hdr <- sprintf(">%s %s|%s|%s", o$oligo, o$target, o$species, o$direction)
  writeLines(c(rbind(hdr, o$sequence)), path)
}
