# Use KilnBio from an assistant

This optional integration is in the source checkout. Use a clone of the revision that contains the `assistants` directory. Historical research release ZIPs contain the software used for their respective analyses. KilnBio runs its statistical code directly. This adapter uses the official MCP SDK and makes no model calls. The separate `evidence/pubmed_mcp.py` adapter queries PubMed.

## Codex on your computer

From the top of the checkout on Linux, macOS or WSL with Python 3.12:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r assistants/requirements.txt
codex mcp add kilnbio -- "$PWD/.venv/bin/python" "$PWD/assistants/kilnbio_mcp.py"
.venv/bin/python assistants/check_mcp.py
```

On Windows, from the top of the checkout in PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r assistants/requirements.txt
codex mcp add kilnbio -- (Resolve-Path .venv/Scripts/python.exe).Path (Resolve-Path assistants/kilnbio_mcp.py).Path
.\.venv\Scripts\python.exe assistants/check_mcp.py
```

Start a Codex session and ask: “Use KilnBio to show how I can compare RNA and protein from the same specimens.” For a calculation, try: “Two donors have RNA. One also has protein. Use `assay_reach` with synthetic identifiers to report how many are eligible in each analysis.”

| MCP tool | What it returns |
|---|---|
| `list_guides` | Checkout version and IDs for the shipped guides |
| `read_guide` | Installation, recipes, assay commands or study instructions for a listed ID |
| `assay_reach` | Counts of distinct donors with both, one or neither assay, plus base-only-policy and combined-policy reach |
| `compare_methods` | Per-feature positive, negative and missing calls, sign conflicts, and support among planned and tested methods for a common analysis context |

Enter small, authorized summary tables directly. The server accepts at most 500 pseudonymous donor rows for reach, or 24,000 characters per TSV for method comparison. Chat content and tool arguments can be visible to the assistant provider and server operator. Run larger tables, restricted datasets and sequencing jobs locally through the normal [KilnBio CLI](../CLI.txt). The MCP tools accept neither local file paths nor executable commands. The worked study preview writes a plan; completed output receipts document execution.

## ChatGPT

Run the same server with Streamable HTTP on the local loopback interface:

```bash
.venv/bin/python assistants/kilnbio_mcp.py --transport streamable-http --port 8000
```

From a second terminal, run `.venv/bin/python assistants/check_mcp.py --url http://127.0.0.1:8000/mcp`. In ChatGPT's Plugins settings, choose **Add custom MCP server**. Connect through a workspace's Secure MCP Tunnel to that local endpoint, or supply the HTTPS /mcp URL of an access-controlled deployment. When proxying to loopback, use the loopback upstream Host header and retain the SDK's transport checks. The server itself binds to loopback and implements no login. Restrict remote access to authorized users through your deployment. A GitHub clone supplies the code; a running service provides the connection. Client availability depends on the workspace. [OpenAI connection instructions](https://developers.openai.com/plugins/deploy/connect-chatgpt)

## Claude Code and Claude.ai

Claude Code discovers [the project skill](../.claude/skills/kilnbio/SKILL.md) when opened at the checkout. Invoke `/kilnbio` or describe a KilnBio analysis. The skill can use a local shell and the existing CLI. To connect the optional MCP tools from the same checkout on Linux or macOS, run `claude mcp add --transport stdio kilnbio -- "$PWD/.venv/bin/python" "$PWD/assistants/kilnbio_mcp.py"`.

For Claude.ai, download [kilnbio-skill.zip](../assistants/kilnbio-skill.zip) using GitHub's **Download raw file**, then upload that ZIP in Claude's Skills settings. The skill supplies workflow guidance by itself. Analysis execution requires a connected KilnBio server or an environment with KilnBio and the route-specific dependencies installed. [Anthropic skill instructions](https://code.claude.com/docs/en/skills)

The [license](../LICENSE) and [access terms](../ACCESS.md) apply to the integration and underlying code.
