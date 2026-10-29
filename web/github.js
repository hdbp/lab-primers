// GitHub REST calls used by the app: read primers.tsv from the default branch,
// and save an edited copy as a new branch + pull request. Uses a fine-grained
// personal access token (Contents and Pull requests: read and write) when one
// is set; reading a public repo works without one.
(function (root) {
  "use strict";
  const API = "https://api.github.com";
  const FILE = "primers.tsv";

  function b64encode(text) {
    const bytes = new TextEncoder().encode(text);
    let bin = "";
    for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
    return btoa(bin);
  }
  function b64decode(b64) {
    const bin = atob(b64.replace(/\s/g, ""));
    return new TextDecoder().decode(Uint8Array.from(bin, c => c.charCodeAt(0)));
  }

  class GitHub {
    constructor(repo, token) { this.repo = repo; this.token = token || ""; }

    async call(method, path, body) {
      const headers = { Accept: "application/vnd.github+json", "X-GitHub-Api-Version": "2022-11-28" };
      if (this.token) headers.Authorization = "Bearer " + this.token;
      if (body) headers["Content-Type"] = "application/json";
      const res = await fetch(`${API}/repos/${this.repo}${path}`,
        { method, headers, body: body ? JSON.stringify(body) : undefined, cache: "no-store" });
      const data = res.status === 204 ? null : await res.json().catch(() => null);
      if (!res.ok) {
        const err = new Error(`GitHub ${res.status}: ${data && data.message || res.statusText}`);
        err.status = res.status;
        throw err;
      }
      return data;
    }

    async info() { return this.call("GET", ""); }

    // -> { text, sha, branch }
    async readPrimers(branch) {
      branch = branch || (await this.info()).default_branch;
      const f = await this.call("GET", `/contents/${FILE}?ref=${encodeURIComponent(branch)}`);
      return { text: b64decode(f.content), sha: f.sha, branch };
    }

    async openAppPulls() {
      const prs = await this.call("GET", "/pulls?state=open&per_page=50");
      return prs.filter(p => p.head.ref.startsWith("primer-app/"));
    }

    // Commits `text` as primers.tsv on a new branch cut from `base` and opens a PR.
    // `baseSha` is the blob sha the edit was made against; if base has moved on,
    // GitHub rejects the commit and the caller should reload and re-apply.
    async proposeChange({ text, baseSha, base, title, body, message }) {
      if (!this.token) throw new Error("Add a GitHub token in Settings to open pull requests.");
      const head = (await this.call("GET", `/git/ref/heads/${encodeURIComponent(base)}`)).object.sha;
      const stamp = new Date().toISOString().replace(/[-:]/g, "").replace(/\..*/, "").replace("T", "-");
      const branch = `primer-app/${stamp}-${Math.random().toString(36).slice(2, 6)}`;
      await this.call("POST", "/git/refs", { ref: `refs/heads/${branch}`, sha: head });
      try {
        await this.call("PUT", `/contents/${FILE}`, { message, content: b64encode(text), sha: baseSha, branch });
        return await this.call("POST", "/pulls", { title, head: branch, base, body });
      } catch (e) {
        await this.call("DELETE", `/git/refs/heads/${branch}`).catch(() => {});
        if (e.status === 409) e.message = "primers.tsv changed on GitHub since you loaded it. Reload, then open the pull request again.";
        throw e;
      }
    }
  }

  root.GitHub = GitHub;
})(self);
