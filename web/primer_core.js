// Primer table logic shared by the web app and tests/test_core.js:
// TSV read/write, validation, length, GC% and Tm. Mirrors R/primer_utils.R;
// tests/test_core.js checks the two give the same numbers.
(function (root) {
  "use strict";

  // One row per primer pair.
  const PRIMER_COLUMNS = ["pair_id", "name", "target", "species", "forward", "reverse", "notes"];
  const SEQ_COLUMNS = ["forward", "reverse"];
  const IUPAC_DNA = "ACGTRYSWKMBDHVN";

  // SantaLucia (1998) unified nearest-neighbour parameters: [dH kcal/mol, dS cal/(K*mol)].
  const NN = { AA: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3], CA: [-8.5, -22.7],
    GT: [-8.4, -22.4], CT: [-7.8, -21.0], GA: [-8.2, -22.2], CG: [-10.6, -27.2],
    GC: [-9.8, -24.4], GG: [-8.0, -19.9] };
  for (const [k, v] of Object.entries({ TT: "AA", AC: "GT", TG: "CA", AG: "CT", TC: "GA", CC: "GG" })) NN[k] = NN[v];

  // R's round(x, 1) (fround.c, R >= 4.0.0): take whichever of the two
  // neighbouring one-decimal values is nearer in double arithmetic; ties go to
  // the even digit.
  function round1(x) {
    const x10 = x * 10, lo = Math.floor(x10), hi = Math.ceil(x10);
    const xd = lo / 10, xu = hi / 10, dd = x - xd, du = xu - x;
    return du < dd || (du === dd && hi % 2 === 0) ? xu : xd;
  }

  function gcPercent(s) {
    if (!s || /[^ACGTSW]/.test(s)) return null;
    return round1(100 * (s.match(/[GCS]/g) || []).length / s.length);
  }

  // Same model and defaults as tm_nn() in R: 50 mM Na+, 250 nM total oligo (CT/4).
  function tmNN(s, naMM = 50, primerNM = 250) {
    if (!/^[ACGT]+$/.test(s || "") || s.length < 2) return null;
    let dH = 0, dS = 0;
    for (let i = 0; i < s.length - 1; i++) { const p = NN[s.slice(i, i + 2)]; dH += p[0]; dS += p[1]; }
    for (const end of [s[0], s[s.length - 1]]) {
      if (end === "G" || end === "C") { dH += 0.1; dS -= 2.8; } else { dH += 2.3; dS += 4.1; }
    }
    dS += 0.368 * (s.length - 1) * Math.log(naMM / 1000);
    const ct = primerNM * 1e-9;
    return round1((1000 * dH) / (dS + 1.987 * Math.log(ct / 4)) - 273.15);
  }

  function annotate(p) {
    const out = { ...p };
    for (const c of SEQ_COLUMNS) {
      out[c + "_length"] = (p[c] || "").length;
      out[c + "_gc_pct"] = gcPercent(p[c]);
      out[c + "_tm_c"] = tmNN(p[c]);
    }
    out.tm_diff = out.forward_tm_c != null && out.reverse_tm_c != null
      ? round1(Math.abs(out.forward_tm_c - out.reverse_tm_c)) : null;
    return out;
  }

  // Same rules and messages as validate_primers() in R. Rows are plain objects.
  function validate(rows) {
    const errors = [], warnings = [];
    const where = i => `line ${i + 2} (${rows[i].pair_id})`;
    rows.forEach((r, i) => {
      for (const c of PRIMER_COLUMNS) {
        if (c !== "notes" && !r[c]) errors.push(`${where(i)}: empty ${c}`);
        if (/[\t\r\n]/.test(r[c] || "")) errors.push(`${where(i)}: ${c} contains a tab or line break`);
      }
      if (r.pair_id && !/^[A-Za-z0-9_.-]+$/.test(r.pair_id))
        errors.push(`${where(i)}: pair_id may only contain letters, digits, '_', '.', '-'`);
    });
    const seen = new Map();
    rows.forEach((r, i) => { if (r.pair_id) seen.set(r.pair_id, [...(seen.get(r.pair_id) || []), i + 2]); });
    for (const [id, lines] of seen) if (lines.length > 1) errors.push(`duplicate pair_id '${id}' on lines ${lines.join(", ")}`);
    const iupac = new RegExp(`^[${IUPAC_DNA}]+$`);
    for (const c of SEQ_COLUMNS) {
      rows.forEach((r, i) => {
        const s = r[c] || "";
        if (!s) return;
        if (!iupac.test(s)) { errors.push(`${where(i)}: ${c} must be uppercase IUPAC DNA (${IUPAC_DNA}) with no spaces or modifications`); return; }
        if (s.length < 10) warnings.push(`${where(i)}: ${c} is only ${s.length} nt`);
        if (/[^ACGT]/.test(s)) warnings.push(`${where(i)}: ${c} contains ambiguity codes; Tm not computed`);
      });
    }
    const labels = new Map();
    for (const c of SEQ_COLUMNS) for (const r of rows) {
      if (r[c]) labels.set(r[c], [...(labels.get(r[c]) || []), `${r.pair_id} ${c}`]);
    }
    for (const l of labels.values()) if (l.length > 1) warnings.push(`identical sequence shared by: ${l.join(", ")}`);
    return { errors, warnings };
  }

  // primers.tsv has no quoting; fields are trimmed like read_primers() in R.
  function parseTSV(text) {
    const lines = text.replace(/\r\n?/g, "\n").split("\n").filter(l => l.length);
    const head = lines.shift().split("\t").map(s => s.trim());
    const missing = PRIMER_COLUMNS.filter(c => !head.includes(c));
    if (missing.length) throw new Error("primers.tsv is missing columns: " + missing.join(", "));
    return lines.map(l => {
      const f = l.split("\t");
      return Object.fromEntries(head.map((h, i) => [h, (f[i] ?? "").trim()]));
    });
  }

  function toTSV(rows) {
    return [PRIMER_COLUMNS, ...rows.map(r => PRIMER_COLUMNS.map(c => r[c] ?? ""))]
      .map(f => f.join("\t")).join("\n") + "\n";
  }

  // Next ID in the dominant prefix: HP014 -> HP015 (keeps the zero padding).
  function nextPairId(rows, taken = []) {
    const m = [...rows.map(r => r.pair_id), ...taken].map(id => /^([A-Za-z]+)(\d+)$/.exec(id || "")).filter(Boolean);
    const prefix = m.some(x => x[1] === "HP") || !m.length ? "HP" : m[0][1];
    const mine = m.filter(x => x[1] === prefix);
    const width = Math.max(3, ...mine.map(x => x[2].length));
    return prefix + String(Math.max(0, ...mine.map(x => +x[2])) + 1).padStart(width, "0");
  }

  function cleanSeq(s) { return String(s || "").toUpperCase().replace(/[\s\-]|5'|3'|5′|3′/g, ""); }

  const api = { PRIMER_COLUMNS, SEQ_COLUMNS, IUPAC_DNA, round1, gcPercent, tmNN, annotate, validate,
    parseTSV, toTSV, nextPairId, cleanSeq };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.PrimerCore = api;
})(typeof self !== "undefined" ? self : this);
