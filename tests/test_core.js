// Run from the repo root: node tests/test_core.js
// Checks web/primer_core.js against the Biopython reference and against
// R's annotate_primers() output (written by tests/run_tests.R to tests/r_annotated.tsv
// when present, otherwise computed here by calling Rscript).
const fs = require("fs");
const { execFileSync } = require("child_process");
const core = require("../web/primer_core.js");

function check(desc, cond) {
  if (!cond) { console.error("FAILED: " + desc); process.exit(1); }
  console.log("ok - " + desc);
}

const ref = core.parseTSV.call(null, "pair_id\tname\ttarget\tspecies\tforward\treverse\tnotes\n"); // header-only parses
check("parseTSV of header only is empty", ref.length === 0);

const tmRef = fs.readFileSync("tests/tm_reference.tsv", "utf8").trim().split("\n").slice(1).map(l => l.split("\t"));
check("tmNN matches Biopython within 0.1 C", tmRef.every(([s, t]) => Math.abs(core.tmNN(s) - +t) <= 0.1));
check("tmNN is null for ambiguous bases", core.tmNN("ACGTNACGTACGTACGT") === null);
check("gcPercent basic", core.gcPercent("GGCC") === 100 && core.gcPercent("AATT") === 0 && core.gcPercent("ACGT") === 50);
check("gcPercent counts S and W", core.gcPercent("SSWW") === 50);
check("gcPercent null for N", core.gcPercent("ACGN") === null);
check("round1 rounds half to even like R", core.round1(0.25) === 0.2 && core.round1(0.35) === 0.3 && core.round1(31.85) === 31.9 && core.round1(2.36) === 2.4);

// Same numbers as R for every pair in primers.tsv.
const text = fs.readFileSync("primers.tsv", "utf8");
const rows = core.parseTSV(text);
check("toTSV round-trips primers.tsv byte for byte", core.toTSV(rows) === text);
const rOut = execFileSync("Rscript", ["-e",
  'source("R/primer_utils.R"); a <- annotate_primers(read_primers("primers.tsv")); ' +
  'write.table(a[, c("pair_id", grep("_(length|gc_pct|tm_c)$|^tm_diff$", names(a), value = TRUE))], ' +
  'stdout(), sep = "\\t", quote = FALSE, row.names = FALSE, na = "")'], { encoding: "utf8" });
const [hdr, ...lines] = rOut.trim().split("\n").map(l => l.split("\t"));
const rAnn = Object.fromEntries(lines.map(f => [f[0], Object.fromEntries(hdr.map((h, i) => [h, f[i]]))]));
const mismatch = [];
for (const r of rows) {
  const a = core.annotate(r);
  for (const k of hdr.slice(1)) {
    const want = rAnn[r.pair_id][k] === "" ? null : +rAnn[r.pair_id][k];
    if (a[k] !== want) mismatch.push(`${r.pair_id} ${k}: js ${a[k]} vs R ${want}`);
  }
}
check("annotate matches R for every pair" + (mismatch.length ? ": " + mismatch.join("; ") : ""), mismatch.length === 0);

// Validation mirrors validate_primers(): same fixtures, same messages.
const ok = core.parseTSV(fs.readFileSync("tests/fixture_valid.tsv", "utf8"));
const okRes = core.validate(ok);
check("valid fixture has no errors", okRes.errors.length === 0);
check("valid fixture warns about ambiguity", okRes.warnings.some(w => w.includes("ambiguity")));
const bad = core.validate(core.parseTSV(fs.readFileSync("tests/fixture_invalid.tsv", "utf8"))).errors;
for (const e of ["forward must be uppercase", "duplicate pair_id 'B002'", "(B004): empty reverse"])
  check("invalid fixture reports: " + e, bad.some(m => m.includes(e)));
check("rejects tabs inside a field", core.validate([{ ...ok[0], notes: "a\tb" }]).errors.length === 1);
check("warns when forward and reverse are identical",
  core.validate([{ ...ok[0], reverse: ok[0].forward }]).warnings.some(w => w.includes("T001 forward, T001 reverse")));

check("nextPairId continues the HP series", core.nextPairId(rows) === "HP015");
check("nextPairId skips IDs already taken", core.nextPairId(rows, ["HP015"]) === "HP016");
check("nextPairId on an empty table", core.nextPairId([]) === "HP001");
check("cleanSeq strips spaces, dashes and 5'/3' labels", core.cleanSeq("5'-acg tTG-3'") === "ACGTTG");

console.log("All JS tests passed.");
