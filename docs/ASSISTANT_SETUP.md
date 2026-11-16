# Use KilnBio from your assistant

Choose the route that matches what you want to do.

| I want to… | What I need |
|---|---|
| Ask Claude.ai which KilnBio workflow fits my data | The [Claude skill ZIP](https://raw.githubusercontent.com/varunorama/kilnbio-research/main/assistants/kilnbio-skill.zip) and a Claude account with Skills enabled. [Install the skill](#claudeai). |
| Have ChatGPT or Claude.ai look up the installed guides, count eligible donors, or compare methods | A running KilnBio MCP server that my account can reach. [Connect ChatGPT](#chatgpt) or [connect Claude.ai](#connect-claudeai-to-the-tools). |
| Use the skill or tools from a terminal | A local checkout. [Set up Claude Code or Codex](#claude-code-and-codex). |
| Use my existing editor or agent | Run the same server locally from [Google Antigravity](#google-antigravity), [Cursor or VS Code with GitHub Copilot](#cursor-vs-code-and-other-mcp-clients). |

The [MCP server](../assistants/kilnbio_mcp.py) and [skill](../.claude/skills/kilnbio/SKILL.md) ship in the public repository's default branch. The earlier research.3 manuscript archives predate the integration. Anyone can inspect the source and download it. Obtain written approval for your specific noncommercial research project before using or adapting KilnBio. Limit a deployment to approved users. [Access and permission requests](../ACCESS.md) · [License](../LICENSE)

## ChatGPT

1. Get the **HTTPS server address ending in `/mcp`** from the person running KilnBio for your project, or obtain the name or ID of your workspace's configured Secure MCP Tunnel. GitHub hosts the source code. Your project supplies the running service. If you are setting up the service yourself, start with [the operator instructions](#run-the-server-for-your-project).
2. In ChatGPT on the web, open [Plugins](https://chatgpt.com/plugins), choose **+ → Add custom MCP server**, and name it **KilnBio**. Under **Connection**, enter the HTTPS address, or choose **Tunnel** and select your project's tunnel.
3. Choose the authentication method supplied by the operator, review the permissions, and select **Create as a plugin**. Check that the four [tools listed below](#verify-a-live-connection) appear. Install the plugin.
4. Start a new chat, type `@`, select **KilnBio**, and try the [two-donor check](#verify-a-live-connection). If your workspace already provides KilnBio, start here. If the add control is missing, ask your workspace administrator whether custom MCP servers are enabled. [OpenAI connection guide](https://developers.openai.com/plugins/deploy/connect-chatgpt)

## Claude.ai

### Upload the skill for workflow guidance

1. Download the [KilnBio skill ZIP](https://raw.githubusercontent.com/varunorama/kilnbio-research/main/assistants/kilnbio-skill.zip). Keep it zipped. This small archive contains the KilnBio skill instructions and license.
2. Open **Customize → Skills → + → Create skill → Upload a skill**. Select the ZIP and turn **kilnbio** on. If Skills is unavailable, check **Settings → Capabilities → Create and edit files**. A managed account may need an owner to enable skills and code execution. [Anthropic skill setup](https://support.claude.com/en/articles/12512180-use-skills-in-claude)
3. In a new chat, try *“Use the KilnBio skill. I have RNA and protein measurements from the same donors. What identifiers do I need to check whether both came from the same specimen?”* The skill can walk you through a workflow and suggest local commands. Running a calculation in Claude.ai through KilnBio takes the separate connector below.

### Connect Claude.ai to the tools

1. Obtain the HTTPS `/mcp` address and sign-in instructions for your approved group's server. An individual account with custom connectors can open **Customize → Connectors → + Add → Add custom connector**. On a managed workspace, an authorized administrator adds the server under **Organization settings → Connectors → Add → Custom → Web**.
2. Name it **KilnBio**, enter the full endpoint, and follow the authentication settings supplied by the operator. After an administrator adds it, members open **Customize → Connectors → KilnBio → Connect**.
3. In a chat, use **+ → Connectors** to enable KilnBio. Try the [two-donor check](#verify-a-live-connection). The connector works independently of the optional skill. [Anthropic custom connector setup](https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp)

## Claude Code and Codex

Use Git to get the current public checkout and record the revision for your study. Run the following commands in a terminal.

```sh
git clone https://github.com/varunorama/kilnbio-research.git kilnbio-assistant
cd kilnbio-assistant
git rev-parse HEAD
```

For the **Claude Code skill**, start `claude` in this directory and type `/kilnbio` followed by your request. Claude Code loads the project skill from `.claude/skills/kilnbio/SKILL.md`. Commands that process your study use the local KilnBio installation and the dependencies for your assay. [Claude Code skills](https://code.claude.com/docs/en/skills) · [KilnBio installation](../README.md)

For **local MCP tools**, install the adapter with Python 3.12, run its check, then run the registration line for your client. The client starts the server when needed. Start a fresh client session in the checkout afterward.

**Linux, macOS or WSL**

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r assistants/requirements.txt
.venv/bin/python assistants/check_mcp.py

# Claude Code
claude mcp add --transport stdio kilnbio -- "$PWD/.venv/bin/python" "$PWD/assistants/kilnbio_mcp.py"
# Codex
codex mcp add kilnbio -- "$PWD/.venv/bin/python" "$PWD/assistants/kilnbio_mcp.py"
```

**Windows PowerShell**

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r assistants/requirements.txt
.\.venv\Scripts\python.exe assistants/check_mcp.py
$kilnbioPython = (Resolve-Path .venv/Scripts/python.exe).Path
$kilnbioServer = (Resolve-Path assistants/kilnbio_mcp.py).Path

# Claude Code
claude mcp add --transport stdio kilnbio -- $kilnbioPython $kilnbioServer
# Codex
codex mcp add kilnbio -- $kilnbioPython $kilnbioServer
```

In Claude Code, `/mcp` shows the connection. In a terminal with Codex installed, `codex mcp list` shows its registration. Use the [tool check](#verify-a-live-connection) in the client to verify that it can call the server. [Claude Code MCP reference](https://code.claude.com/docs/en/mcp)

## Google Antigravity

1. Follow the [terminal setup](#claude-code-and-codex) through the Python environment, dependency installation and self-check. Skip the Claude and Codex registration commands. For a local Antigravity session, the editor starts KilnBio when needed.
2. Open your MCP configuration. In the **Antigravity App**, go to **Customizations → Installed → MCP Servers** and use **Copy Config File Path** to locate the file. In the **Antigravity IDE**, open the agent side panel and select **… → MCP Servers → Manage MCP Servers → View raw config**. For the **Antigravity CLI**, use `.agents/mcp_config.json` in the workspace and open `/mcp` to reload or inspect the connection. [Google's MCP guide](https://antigravity.google/docs/mcp)
3. Add the `kilnbio` entry below under `mcpServers`, keeping any existing entries. Replace **both** example paths with absolute paths on the computer where the client runs. On Windows, forward slashes work in these JSON paths. On macOS or Linux, use the checkout's `.venv/bin/python` and `assistants/kilnbio_mcp.py`.

```json
{
  "mcpServers": {
    "kilnbio": {
      "command": "C:/research/kilnbio-assistant/.venv/Scripts/python.exe",
      "args": ["C:/research/kilnbio-assistant/assistants/kilnbio_mcp.py"]
    }
  }
}
```

4. Save, reload the server or restart the client, and try the [two-donor check](#verify-a-live-connection). If your team already hosts KilnBio, set `serverUrl` to its authenticated HTTPS `/mcp` address in place of `command` and `args`, then sign in as directed by the operator.

For workflow guidance alone, extract the [KilnBio skill ZIP](https://raw.githubusercontent.com/varunorama/kilnbio-research/main/assistants/kilnbio-skill.zip) into the approved project's `.agents/skills/` directory so the file is at `.agents/skills/kilnbio/SKILL.md`. The MCP setup above provides the live calculations. [Antigravity skills](https://antigravity.google/docs/skills)

## Cursor, VS Code and other MCP clients

| Client | How to connect |
|---|---|
| **Cursor desktop** | Use the [Antigravity example](#google-antigravity) in the project's `.cursor/mcp.json`. Add `"type": "stdio"` alongside `"command"` in the `kilnbio` entry. Open the client's MCP settings and check the tools. [Cursor MCP reference](https://cursor.com/docs/mcp) |
| **GitHub Copilot in VS Code** | In a local project, select Copilot and run **MCP: Add Server** from the Command Palette. Choose **Command (stdio)**. Enter the full path to the environment's Python followed by the full path to `assistants/kilnbio_mcp.py`, quoting each path if it contains spaces. Name the server `kilnbio` and choose the intended configuration scope. [VS Code connection guide](https://code.visualstudio.com/docs/agent-customization/mcp-servers) |
| **Another MCP client** | Check which transports it accepts. For local **stdio**, supply the same executable and single script argument. For **Streamable HTTP**, use the approved operator's reachable `/mcp` address and supported authentication method. A client with only a GitHub repository connector can read the source. Executing the tools requires MCP tool support. |

Confirm each connection by inspecting tool activity for the [two-donor check](#verify-a-live-connection). The adapter has been exercised with the MCP SDK using stdio and loopback HTTP. Client instructions follow each vendor's documentation.

## Picking up new pipelines

The server reads the command index, sequencing guide and recipe index from its installed checkout when you ask for them. When a maintainer documents a new workflow in one of those sources, `read_guide` can explain it across clients using that revision. Adding a standalone guide with its own ID requires a new entry in `GUIDES` in [the adapter](../assistants/kilnbio_mcp.py). A new calculation that users should execute **inside chat** needs an explicit MCP tool and a test. Full pipeline execution remains available in the KilnBio CLI with the workflow's own dependencies.

For updates within your project approval, run `git pull --ff-only` in a clean checkout tracking `main`, reinstall `assistants/requirements.txt` with the environment's Python, rerun `assistants/check_mcp.py`, and restart the server. Refresh the client's tools after a tool change. If you use a shared endpoint, its operator performs this update. Record the commit with `git rev-parse HEAD`. A Claude.ai skill upload or an extracted Antigravity skill is a snapshot; replace it when the skill instructions change. An analysis pinned to a commit continues to use that revision.

## Run the server for your project

This step is for the person providing a connection to a web client. Use the checkout, Python environment, dependency installation and local check from [the terminal instructions](#claude-code-and-codex). From the checkout, start the service and keep it running. In a second terminal at the same directory, run the HTTP check.

| Terminal | Linux, macOS or WSL | Windows PowerShell |
|---|---|---|
| Start the service | `.venv/bin/python assistants/kilnbio_mcp.py --transport streamable-http --port 8000` | `.\.venv\Scripts\python.exe assistants/kilnbio_mcp.py --transport streamable-http --port 8000` |
| Check the connection | `.venv/bin/python assistants/check_mcp.py --url http://127.0.0.1:8000/mcp` | `.\.venv\Scripts\python.exe assistants/check_mcp.py --url http://127.0.0.1:8000/mcp` |

The check prints `KilnBio MCP checks passed.` The address above is local to this machine. For ChatGPT, a workspace administrator can connect a [Secure MCP Tunnel](https://developers.openai.com/api/docs/guides/secure-mcp-tunnels) to the loopback service. For either web client, an operator can provide an authenticated HTTPS gateway with Streamable HTTP at `/mcp`. The adapter itself binds to loopback and implements no user authentication. Have the gateway enforce access for the approved project, forward the loopback upstream Host header and preserve the SDK's transport checks. Give users the resulting endpoint and sign-in instructions. Keep credentials out of chats.

## Verify a live connection

Ask *“Call KilnBio's `list_guides` and tell me which version it returns. Then call `read_guide` for `study` and show me the preview option.”* Check the tool activity and look for `--dry-run` in the returned guide.

Then try *“Call KilnBio's `assay_reach`. Donor D1 has RNA and protein. Donor D2 has RNA only. Use RNA as the base assay and report both counts.”* The tool result should include **base_reach = 2** and **combined_reach = 1**. The tool activity confirms the software actually ran. These counts describe donor availability. Same-specimen pairing needs specimen identifiers and the local pairing command.

| Tool | What you can do |
|---|---|
| `list_guides` and `read_guide` | Read versioned KilnBio instructions directly from the operator's checkout. |
| `assay_reach` | Count donors retained by one assay or by requiring both, with up to 500 distinct pseudonymous donor rows. |
| `compare_methods` | Compare directions, significant calls and missing results across methods for the same contrast using two small summary tables. |

Use synthetic examples or inputs approved for the chat provider and server operator. Keep restricted data and full sequencing runs in your approved compute environment. The standard [KilnBio command line](../CLI.txt) runs independently of these assistant integrations.
