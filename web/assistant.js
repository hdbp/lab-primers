// Claude assistant for the primer app: finds published primer pairs (web search
// and fetch) or designs new ones, and returns each candidate through the
// propose_primer_pairs tool so the page can check it and offer "Add to draft".
// Runs in the browser with the user's own Anthropic API key.
import Anthropic from "https://cdn.jsdelivr.net/npm/@anthropic-ai/sdk@0.133.0/+esm";

const MODEL = "claude-opus-5-5";

const SYSTEM = `You help a molecular biology lab maintain its primer database (one row per primer pair: name, target, species, forward, reverse, notes). The user asks you to find published, validated primer pairs or to design new ones.

How to work:
- Prefer published pairs that were used successfully in the same assay type and organism (qPCR, ChIP-qPCR, RT-qPCR, genotyping...). Search the literature and primer databases (PrimerBank, RTPrimerDB, papers' methods and supplementary tables). Fetch the source and copy sequences exactly as printed, 5'->3'. Never present a pair as published unless you read it in a source you can cite.
- If nothing suitable is published, or the user asks you to design, propose new pairs and say clearly they are designed and unvalidated. Aim for 18-25 nt, 40-60% GC, primer Tm within ~2 °C of each other, no 3' runs of G/C, amplicon size suited to the assay (qPCR 70-200 bp).
- Check the current database below and do not propose a pair that is already in it; mention it instead.
- You cannot run in-silico PCR or BLAST here. For every pair, tell the user what still needs checking (e.g. UCSC In-Silico PCR or Primer-BLAST against the right assembly) unless the source already reports genomic coordinates and specificity.

Returning pairs:
- Call propose_primer_pairs with every candidate pair. The page computes length, GC% and Tm (SantaLucia nearest-neighbour, 50 mM Na+, 250 nM oligo; runs a few °C below Primer3) and returns them with any validation problems or duplicates. If a pair has problems, fix or drop it and call the tool again.
- name: short and in the database's style, e.g. ChIP_H3K27me3_Hoxb3-CGI or qPCR_Gapdh. target: gene/locus and region, with coordinates and assembly when known. species: binomial, e.g. Mus musculus.
- notes: write the full notes field in the same style as existing rows: assay; ChIP target if any; amplicon coordinates/size and assembly if known; evidence (published and how it was validated, or "Designed, unvalidated"); source as PMID/PMC/doi; what still needs checking.
- After the tool call, give a short summary for the user. Keep it brief; the cards on the page show the details.`;

const PROPOSE_TOOL = {
  name: "propose_primer_pairs",
  description: "Show candidate primer pairs to the user as cards they can add to their draft database. The result returns computed length, GC%, Tm, Tm difference, validation problems and duplicates for each pair.",
  strict: true,
  eager_input_streaming: true,
  input_schema: {
    type: "object",
    additionalProperties: false,
    required: ["pairs"],
    properties: {
      pairs: {
        type: "array",
        items: {
          type: "object",
          additionalProperties: false,
          required: ["name", "target", "species", "forward", "reverse", "evidence", "source_url", "notes"],
          properties: {
            name: { type: "string", description: "Short pair name in the database's style" },
            target: { type: "string", description: "Gene/locus and region; coordinates and assembly when known" },
            species: { type: "string", description: "Binomial name, e.g. Mus musculus" },
            forward: { type: "string", description: "Forward primer 5'->3', uppercase DNA" },
            reverse: { type: "string", description: "Reverse primer 5'->3', uppercase DNA" },
            evidence: { type: "string", enum: ["published", "designed"] },
            source_url: { type: "string", description: "URL of the paper or database entry the sequences were read from; empty string for designed pairs" },
            notes: { type: "string", description: "Full notes field in the database's style" }
          }
        }
      }
    }
  }
};

const TOOLS = [
  { type: "web_search_20260209", name: "web_search", max_uses: 8 },
  { type: "web_fetch_20260209", name: "web_fetch", max_uses: 10 },
  PROPOSE_TOOL
];

function validPropose(input) {
  const keys = PROPOSE_TOOL.input_schema.properties.pairs.items.required;
  return input && Array.isArray(input.pairs) &&
    input.pairs.every(p => p && keys.every(k => typeof p[k] === "string") && ["published", "designed"].includes(p.evidence));
}

// One conversation. `ui` callbacks: onText(delta), onThinking(), onSources(list),
// onPairs(pairs) -> tool result object, onStatus(text), onDone(stopReason).
export class Assistant {
  constructor(apiKey, databaseTSV) {
    this.client = new Anthropic({ apiKey, dangerouslyAllowBrowser: true });
    // Frozen at conversation start so the cached prefix stays stable.
    this.system = [
      { type: "text", text: SYSTEM },
      { type: "text", text: "Current database (primers.tsv):\n" + databaseTSV, cache_control: { type: "ephemeral" } }
    ];
    this.messages = [];
    this.stream = null;
  }

  stop() { if (this.stream) this.stream.abort(); }

  // On an error or Stop, the unfinished exchange is dropped so the history stays valid.
  async send(text, ui) {
    const start = this.messages.length;
    try { await this.run(text, ui); }
    catch (e) { this.stream = null; this.messages.length = start; throw e; }
  }

  async run(text, ui) {
    this.messages.push({ role: "user", content: text });
    for (let turn = 0; turn < 20; turn++) {
      ui.onStatus(turn ? "Working…" : "Thinking…");
      this.stream = this.client.beta.messages.stream({
        model: MODEL,
        max_tokens: 64000,
        thinking: { type: "adaptive" },
        output_config: { effort: "high" },
        betas: ["server-side-fallback-2026-07-01"],
        fallbacks: "default",
        system: this.system,
        tools: TOOLS,
        messages: this.messages
      });
      this.stream.on("streamEvent", ev => {
        if (ev.type === "content_block_start") {
          const b = ev.content_block;
          if (b.type === "server_tool_use") ui.onStatus(b.name === "web_search" ? "Searching the web…" : "Reading a source…");
          else if (b.type === "tool_use") ui.onStatus("Preparing primer pairs…");
          else if (b.type === "thinking") ui.onStatus("Thinking…");
        } else if (ev.type === "content_block_delta" && ev.delta.type === "text_delta") {
          ui.onText(ev.delta.text);
        }
      });
      const msg = await this.stream.finalMessage();
      this.stream = null;
      // Append the whole content (thinking, server tool blocks, fallback blocks) unchanged.
      this.messages.push({ role: "assistant", content: msg.content });

      const sources = [];
      for (const b of msg.content) {
        if (b.type === "text" && b.citations) for (const c of b.citations) if (c.url) sources.push({ url: c.url, title: c.title || c.url });
      }
      if (sources.length) ui.onSources(sources);

      if (msg.stop_reason === "pause_turn") continue;
      if (msg.stop_reason === "tool_use") {
        const results = [];
        for (const b of msg.content) {
          if (b.type !== "tool_use") continue;
          if (b.name !== PROPOSE_TOOL.name || !validPropose(b.input)) {
            results.push({ type: "tool_result", tool_use_id: b.id, is_error: true,
              content: "INVALID_JSON: input did not match the propose_primer_pairs schema; call it again with complete pairs." });
            continue;
          }
          results.push({ type: "tool_result", tool_use_id: b.id, content: JSON.stringify(ui.onPairs(b.input.pairs)) });
        }
        this.messages.push({ role: "user", content: results });
        continue;
      }
      ui.onDone(msg.stop_reason);
      return;
    }
    ui.onDone("turn_limit");
  }
}
