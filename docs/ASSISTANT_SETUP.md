# Use KilnBio from your assistant

Choose the route that matches what you want to do.

| I want to… | What I need |
|---|---|
| Run a pipeline and get its output receipt | An approved project, a Linux/WSL or macOS execution host with Bash, the workflow dependencies, and a project-enabled connection. [Install and run the example](#install-and-run-a-pipeline). |
| Ask Claude.ai which KilnBio workflow fits my data | The [Claude skill ZIP](https://raw.githubusercontent.com/varunorama/kilnbio-research/main/assistants/kilnbio-skill.zip) and a Claude account with Skills enabled. [Install the skill](#claudeai). |
| Have ChatGPT or Claude.ai look up the installed guides, count eligible donors, or compare methods | A running KilnBio MCP server that my account can reach. [Connect ChatGPT](#chatgpt) or [connect Claude.ai](#connect-claudeai-to-the-tools). |
| Use the skill or tools from a terminal | A local checkout. [Set up Claude Code or Codex](#claude-code-and-codex). |
| Use my existing editor or agent | Run the same server locally from [Google Antigravity](#google-antigravity), [Cursor or VS Code with GitHub Copilot](#cursor-vs-code-and-other-mcp-clients). |
| Run on my lab's existing cloud VM | Connect a local MCP client over [SSH to AWS, Google Cloud or a lab Linux host](#use-your-labs-existing-compute). |
| Use KilnBio inside a scientific workbench | Configure a [custom connector in Claude Science](#claude-science). |

The [MCP server](../assistants/kilnbio_mcp.py) and [skill](../.claude/skills/kilnbio/SKILL.md) ship in the public repository's default branch. The earlier research.3 manuscript archives predate the integration. Anyone can inspect the source and download it. Obtain written approval for your specific noncommercial research project before using or adapting KilnBio. Limit a deployment to approved users. [Access and permission requests](../ACCESS.md) · [License](../LICENSE)

## Install and run a pipeline

**Local clients.** Install Git, Python 3.12 with venv/pip, and your chosen client. Codex and Claude Code require their CLI on PATH. From a terminal, clone this repository and enter the checkout:

```sh
git clone https://github.com/varunorama/kilnbio-research.git kilnbio-assistant
cd kilnbio-assistant
```

On Linux, macOS or inside WSL, run **one** command for the client installed in that environment:

| Client | Install the adapter and register the synthetic example |
|---|---|
| Codex | `python3.12 assistants/install_client.py codex --demo` |
| Claude Code | `python3.12 assistants/install_client.py claude --demo` |
| Google Antigravity | `python3.12 assistants/install_client.py antigravity --demo` |
| Other local MCP clients | `python3.12 assistants/install_client.py config --demo` prints the executable and arguments to register. |

For a native Windows client, run the installer from PowerShell with your WSL distro, for example `py -3.12 assistants/install_client.py codex --demo --wsl Ubuntu`. Substitute `claude` or `antigravity` for the other clients. Install Python 3.12 and Bash inside the distro first. Use `--linux-python /absolute/path/to/python3.12` if needed. The installer registers a Windows launch command that starts KilnBio inside WSL. Open the Claude Code project where you want its local registration to apply. Antigravity registration uses the user's MCP configuration.

The installer creates a Python environment in the checkout and copies a small fictional capture study to `~/.kilnbio/assistant-demo` on the execution host. Keep the checkout in place. For real analyses, substitute `--project /absolute/path/to/approved-project` for `--demo`. Install the selected workflow's dependencies on that host. A separate client registration can use `--name kilnbio-lab`.

For an existing scientific environment, configure the server process with the same `KILNBIO_PYTHON`, `KILNBIO_RSCRIPT`, R library and workflow runtime settings you use on the CLI. Without an explicit `KILNBIO_PYTHON`, Python stages use the adapter's interpreter. A FASTQ workflow also needs its reference inputs, Nextflow and container runtime.

**ChatGPT and Claude.ai.** A project operator sets up the [execution endpoint](#run-the-server-for-your-project) once. Approved users then connect it in [ChatGPT](#chatgpt) or [Claude.ai](#connect-claudeai-to-the-tools). The repository supplies the server software. The web clients need the project's endpoint or supported tunnel. The Claude skill ZIP alone supplies instructions.

**Already have the project's HTTPS endpoint?** A researcher can skip the clone and local runtime. In ChatGPT or Claude.ai, add the endpoint in the client as described below. For Codex, run `codex mcp add kilnbio --url https://YOUR-PROJECT-HOST/mcp`. For Claude Code, use `claude mcp add --transport http kilnbio https://YOUR-PROJECT-HOST/mcp`. In Antigravity, set the server's `serverUrl` to that same address. Replace the placeholder with the URL supplied by your operator and complete the project's sign-in. These are custom connections. Publishing to the default branch does not create a listing in a provider's marketplace.

**First run.** Reload the client's MCP connection or open a fresh session. Ask:

> List my KilnBio pipelines. Run `capture-demo` with run ID `first-panel`. Check its status until it finishes, then show the output directory and the files KilnBio verified.

Look for a `run_pipeline` call followed by `pipeline_status` returning `complete`. In the project, open `runs/first-panel/study/panel/summary.json`. The example's two overlapping target intervals cover **80 distinct bases**. Of **three** synthetic features, **one** lies fully inside. The receipt includes SHA-256 hashes for the summary, feature table and target BED. The assistant receives stage and artifact information; the measurements and full logs stay on the execution host.

To exercise the adapter independently of any assistant account, run `.venv/bin/python assistants/check_pipeline.py` on the execution host. The test runs the example, installs and launches a second two-stage manifest while connected, checks duplicate submissions and a missing input, and verifies the receipt over both stdio and loopback HTTP. It uses temporary synthetic projects.

## Use your lab's existing compute

An approved project can use its own Linux host and storage. The assistant subscription supplies the interface. The selected machine runs the KilnBio processes and stores their results. This setup uses an existing machine and leaves billing, approved inputs, storage permissions and compute limits with its owner.

1. **On the VM**, the operator installs the checkout, Python 3.12, Bash, workflow dependencies and project manifests as in the local setup. From the checkout, run `python3.12 assistants/install_client.py config --project /absolute/path/to/project` to prepare the adapter environment. Record the checkout and project paths. Input data must be accessible to the VM. For example, the operator can stage approved data from S3 or Cloud Storage using the VM's existing, scoped cloud identity.
2. **On the client**, configure an OpenSSH host alias such as `kilnbio-lab` and verify the host key through the lab's usual login procedure. Check that `ssh -T -o BatchMode=yes -o StrictHostKeyChecking=yes kilnbio-lab true` succeeds without prompting. EC2 researchers can use their existing Linux SSH access. Google Compute Engine users can use an OpenSSH configuration for an internal VM via Identity-Aware Proxy. [AWS EC2 SSH](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/connect-linux-inst-ssh.html) · [Google Compute Engine SSH with IAP](https://docs.cloud.google.com/compute/docs/connect/ssh-using-iap)
3. **Register the connection** using the checkout on a configuration machine. Substitute the VM's paths:

```sh
python assistants/install_client.py codex --name kilnbio-lab --ssh kilnbio-lab --remote-checkout /srv/kilnbio --project /srv/projects/approved-study
```

Choose `claude` for Claude Code, `antigravity` for Antigravity, or `config` to print the command and its argument list. A lab operator can give researchers those settings so each can enter them directly in a local-command client without a local KilnBio clone. Use the path to OpenSSH on the researcher's own machine. SSH performs authentication using that user's normal credentials. The client launches KilnBio on the VM and exchanges tool requests and receipts over the SSH connection. Files are read from and written to the VM. The lab can replace SSH access with a project-authenticated HTTPS gateway for [web-client access](#run-the-server-for-your-project).

This route uses an existing VM. A cloud account or storage bucket alone does not supply a running KilnBio worker. AWS Batch and Google Cloud Batch would require separate scheduler integration. The SSH transport was tested using a real local OpenSSH client and server with synthetic data. Live provider-specific IAM and network settings must be checked on the lab's host.

## Claude Science

Claude Science can provide the workbench for a KilnBio study. KilnBio supplies the versioned study manifest, assay and statistical code, and execution receipt. The same study can be run directly with the command line.

In Claude Science, open **Settings → Connectors → Add connector**. Choose **Local command** for a local or SSH-backed server and name it `kilnbio`. Copy the `command` and individual `args` produced by the installer's `config` mode into the corresponding fields. On Windows, use the full path to the executable, such as `C:/Windows/System32/wsl.exe`. Local commands run inside Claude Science's sandbox; the project must be accessible and permitted to write its run directory, and an SSH connector needs network access to the VM. Alternatively, choose **Remote** and enter your authorized project's HTTPS MCP URL. Managed organizations can control custom connector availability. Verify with the [synthetic capture study](#install-and-run-a-pipeline). These steps follow Anthropic's documentation; a signed-in Claude Science test is still pending. [Claude Science custom connectors](https://claude.com/docs/claude-science/custom-connectors)

## ChatGPT

1. Get the **HTTPS server address ending in `/mcp`** from the person running KilnBio for your project, or obtain the name or ID of your workspace's configured Secure MCP Tunnel. GitHub hosts the source code. Your project supplies the running service. If you are setting up the service yourself, start with [the operator instructions](#run-the-server-for-your-project).
2. In ChatGPT on the web, open [Plugins](https://chatgpt.com/plugins), choose **+ → Add custom MCP server**, and name it **KilnBio**. Under **Connection**, enter the HTTPS address, or choose **Tunnel** and select your project's tunnel.
3. Choose the authentication method supplied by the operator, review the permissions, and select **Create as a plugin**. Check the [tool list](#verify-a-live-connection). A project-enabled server also advertises `list_pipelines`, `run_pipeline`, and `pipeline_status`. Install the plugin.
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

On the execution host, place a reviewed [KilnBio study manifest](STUDY_RUNNER.md) at `PROJECT/pipelines/PIPELINE_ID/study.json` with its input files or configured paths. Use an ID of up to 64 lowercase letters, numbers, hyphens or underscores, starting with a letter or number. Test the recipe with the normal study runner. The next `list_pipelines` call discovers the manifest, and `run_pipeline` executes it through the same runner used by `kilnbio study`. No MCP code change, assistant reinstall or server restart is needed for a new manifest. A study that introduces a new underlying analysis command also needs a KilnBio version that supports that command and its runtime dependencies. [Example project](../examples/assistant-project/pipelines/capture-demo)

Keep manifests, input mappings and remote-service settings under project-operator control. All authorized users of the same endpoint share its catalog and run namespace. Use a separate endpoint, project directory and OS access boundary for each group with different permissions. Supply filenames and step IDs suitable for display in chat. Execution tools accept pipeline and run IDs. Study arguments and paths are taken from the operator's manifest.

Each distinct run ID creates a new directory. Repeating a request with the **same** pipeline and run ID returns the original run's receipt. Use the same ID after a dropped connection. A different pipeline needs a different ID. After a failure, inspect `runs/RUN_ID/launcher.log` and the stage logs on the host. Once the issue is fixed, use a new run ID or the documented [local study resume procedure](STUDY_RUNNER.md#preview-and-resume). Runs continue after a client disconnects while the host stays up.

For updates within your project approval, run `git pull --ff-only` in a clean checkout tracking `main`, reinstall `assistants/requirements.txt` with the environment's Python, rerun `assistants/check_mcp.py`, and restart the server. Refresh the client's tools after a tool change. If you use a shared endpoint, its operator performs this update. Record the commit with `git rev-parse HEAD`. A Claude.ai skill upload or an extracted Antigravity skill is a snapshot; replace it when the skill instructions change. An analysis pinned to a commit continues to use that revision.

## Run the server for your project

This step is for the person providing a connection to a web client. Use the checkout, Python environment, dependency installation and local check from [the terminal instructions](#claude-code-and-codex). From the checkout, start the service and keep it running. In a second terminal at the same directory, run the HTTP check.

| Terminal | Linux, macOS or WSL | Windows PowerShell |
|---|---|---|
| Start the service | `.venv/bin/python assistants/kilnbio_mcp.py --transport streamable-http --port 8000` | `.\.venv\Scripts\python.exe assistants/kilnbio_mcp.py --transport streamable-http --port 8000` |
| Check the connection | `.venv/bin/python assistants/check_mcp.py --url http://127.0.0.1:8000/mcp` | `.\.venv\Scripts\python.exe assistants/check_mcp.py --url http://127.0.0.1:8000/mcp` |

The check prints `KilnBio MCP checks passed.` To enable actual runs, start the server on the execution host with the project explicitly selected:

```bash
.venv/bin/python assistants/kilnbio_mcp.py --project /absolute/path/to/approved-project --transport streamable-http --port 8000
```

The service account needs read access to approved manifests and inputs, write access to `PROJECT/runs`, and the normal KilnBio workflow dependencies. The adapter binds to loopback and relies on its operator for authentication. On a shared host, isolate the loopback service from other users. For ChatGPT, an administrator can configure a [Secure MCP Tunnel](https://developers.openai.com/api/docs/guides/secure-mcp-tunnels) and limit the plugin to approved accounts. For either web client, place a project-authenticated HTTPS gateway in front of the private service and enforce membership for every request. Forward the loopback upstream Host header and preserve the SDK's transport checks. Give approved users the resulting endpoint and sign-in instructions. Keep credentials out of chats.

## Verify a live connection

Ask *“Call KilnBio's `list_guides` and tell me which version it returns. Then call `read_guide` for `study` and show me the preview option.”* Check the tool activity and look for `--dry-run` in the returned guide.

Then try *“Call KilnBio's `assay_reach`. Donor D1 has RNA and protein. Donor D2 has RNA only. Use RNA as the base assay and report both counts.”* The tool result should include **base_reach = 2** and **combined_reach = 1**. The tool activity confirms the software actually ran. These counts describe donor availability. Same-specimen pairing needs specimen identifiers and the local pairing command.

| Tool | What you can do |
|---|---|
| `list_guides` and `read_guide` | Read versioned KilnBio instructions directly from the operator's checkout. |
| `assay_reach` | Count donors retained by one assay or by requiring both, with up to 500 distinct pseudonymous donor rows. |
| `compare_methods` | Compare directions, significant calls and missing results across methods for the same contrast using two small summary tables. |
| `list_pipelines` | See the installed studies and their stage commands in the configured project. |
| `run_pipeline` | Launch a named project study with a stable run ID. |
| `pipeline_status` | Check completion, stage state, and relative locations and hashes of the declared outputs. |

Use synthetic examples or inputs approved for the chat provider and server operator. Keep restricted data and full sequencing runs in your approved compute environment. The standard [KilnBio command line](../CLI.txt) runs independently of these assistant integrations.
