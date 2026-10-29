// Run from the repo root: node tests/test_github.js
// Checks web/github.js against a fake GitHub API: reading primers.tsv, and the
// branch -> commit -> pull request sequence (and cleanup when the commit fails).
global.self = global;
require("../web/github.js");
const { GitHub } = global;

function check(desc, cond) {
  if (!cond) { console.error("FAILED: " + desc); process.exit(1); }
  console.log("ok - " + desc);
}

let calls = [], fail = {};
const tsv = "pair_id\tname\nHP001\tµ-test ✓\n";
global.fetch = async (url, opt) => {
  const path = url.replace("https://api.github.com/repos/o/r", ""), method = opt.method;
  calls.push({ method, path, body: opt.body && JSON.parse(opt.body), auth: opt.headers.Authorization });
  const key = `${method} ${path.split("?")[0]}`;
  const reply = (status, data) => ({ ok: status < 300, status, statusText: "", json: async () => data });
  if (fail[key]) return reply(fail[key], { message: "nope" });
  if (key === "GET ") return reply(200, { default_branch: "main", permissions: { push: true } });
  if (key === "GET /contents/primers.tsv") return reply(200, { sha: "blob1", content: Buffer.from(tsv).toString("base64").replace(/(.{60})/g, "$1\n") });
  if (key === "GET /git/ref/heads/main") return reply(200, { object: { sha: "commit1" } });
  if (key === "POST /git/refs") return reply(201, {});
  if (key === "PUT /contents/primers.tsv") return reply(200, {});
  if (key === "POST /pulls") return reply(201, { number: 7, html_url: "https://github.com/o/r/pull/7" });
  if (method === "DELETE") return reply(204, null);
  return reply(404, { message: "Not Found" });
};

(async () => {
  const f = await new GitHub("o/r").readPrimers();
  check("readPrimers decodes UTF-8 base64 with line breaks", f.text === tsv && f.sha === "blob1" && f.branch === "main");
  check("reads without a token send no Authorization header", calls.every(c => !c.auth));

  calls = [];
  const g = new GitHub("o/r", "tok");
  const pr = await g.proposeChange({ text: tsv, baseSha: "blob1", base: "main", title: "T", body: "B", message: "M" });
  const seq = calls.map(c => `${c.method} ${c.path}`);
  check("proposeChange: ref, branch, commit, pull request in order",
    seq.join() === "GET /git/ref/heads/main,POST /git/refs,PUT /contents/primers.tsv,POST /pulls");
  const branch = calls[1].body.ref.replace("refs/heads/", "");
  check("branch is cut from the base head under primer-app/", branch.startsWith("primer-app/") && calls[1].body.sha === "commit1");
  check("commit sends the base blob sha and the encoded table",
    calls[2].body.sha === "blob1" && calls[2].body.branch === branch && Buffer.from(calls[2].body.content, "base64").toString() === tsv);
  check("pull request goes from the new branch to base", calls[3].body.head === branch && calls[3].body.base === "main" && pr.number === 7);
  check("token is sent as a bearer token", calls.every(c => c.auth === "Bearer tok"));

  calls = []; fail = { "PUT /contents/primers.tsv": 409 };
  let err;
  try { await g.proposeChange({ text: tsv, baseSha: "old", base: "main", title: "T", body: "B", message: "M" }); } catch (e) { err = e; }
  check("a stale base sha reports that primers.tsv changed", err && /changed on GitHub/.test(err.message));
  check("the half-made branch is deleted", calls.at(-1).method === "DELETE" && calls.at(-1).path.startsWith("/git/refs/heads/primer-app/"));

  let noTok;
  try { await new GitHub("o/r").proposeChange({ text: tsv, baseSha: "x", base: "main" }); } catch (e) { noTok = e; }
  check("opening a pull request without a token fails clearly", noTok && /token/.test(noTok.message));
  console.log("All GitHub tests passed.");
})();
