# Aedes Local Adaptation Analysis

## ⚠️ WEB SEARCH: ONE TOOL PER QUERY — NO DUPLICATES ⚠️

**Pick ONE search tool per query. NEVER call both for the same question.**

| Need | Tool | Notes |
|------|------|-------|
| Quick answer | Perplexity `web_search` model=**sonar** | DEFAULT. $1/1K. **NEVER use sonar-pro unless user explicitly asks.** |
| Find paper URLs | Exa `research_papers` | Returns links + highlights |
| Find code examples | Exa `code_search` | GitHub/SO/docs |
| Fetch URL content | crawl4ai `fetch_page` | Free, saves to RAG. Exa `crawl` as fallback only. |
| Deep literature review | Perplexity `research_search` model=sonar | Upgrade to sonar-pro ONLY if sonar results are insufficient |

**RULES:**
- **sonar-pro costs 15x more than sonar.** Do NOT use it for general web searches.
- **Do NOT run Perplexity AND Exa for the same query.** That wastes money and tokens.
- **crawl4ai first** for fetching URLs (free). Exa crawl only if crawl4ai fails.

---

## Project Type: Research

Population genomics pipeline for studying local adaptation in *Aedes aegypti* populations.

---

## SESSION TRACKING (MANDATORY)

**Always use `/start` with explicit project argument. Never rely on auto-detect.**

```bash
/start aedes-local-adaptation session-name
```

**Workflow:**
1. Claude starts -> SessionStart hook captures session ID automatically
2. `/start aedes-local-adaptation <name>` -> locks project
3. Work -> `/wrapup` at end -> stores handoff
4. Resume: `claude --resume <session-id>` (shown by /pickup)

**Task ID prefix:** `ALA`

**Use MCP tools for session management:**
- `mcp__session-tracker__create_task` for task tracking
- `mcp__session-tracker__log_event` for documenting findings
- NEVER write YAML session files manually

---

## SHARED TOOLS (from tracking-hub)

```bash
# Search RAG for previous work
mcp__tracking-hub-rag__search_rag "local adaptation" --project aedes-local-adaptation

# Web search
mcp__perplexity__web_search "Aedes aegypti local adaptation genomics"

# Literature search
mcp__literature__search_pubmed "Aedes aegypti population genomics"
```

---

## CONTAINER

```bash
docker run -it ghcr.io/cosmelab/aedes-local-adaptation:latest
```

---

## RULES

1. **NEVER** add Claude, Anthropic, or AI attribution to commits/files
2. **NEVER** duplicate tools from tracking-hub
3. **ALWAYS** keep files under 100 lines for RAG
4. **ALWAYS** document progress via MCP session tracker
