"""Install the optional KilnBio MCP adapter and register an approved project."""
import argparse
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import venv


ROOT = Path(__file__).resolve().parents[1]


def prepare(args):
    if args.ssh:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.ssh):
            raise ValueError("--ssh expects a host alias from your OpenSSH configuration")
        if not args.remote_checkout or not args.project:
            raise ValueError("--ssh requires --remote-checkout and --project on the remote host")
        remote = PurePosixPath(args.remote_checkout)
        if not remote.is_absolute() or not PurePosixPath(args.project).is_absolute():
            raise ValueError("Use absolute Linux paths for the remote checkout and project")
        ssh = shutil.which("ssh")
        if not ssh:
            raise ValueError("Install OpenSSH on the client first")
        command = shlex.join([str(remote / ".venv/bin/python"), str(remote / "assistants/kilnbio_mcp.py"),
                              "--project", args.project])
        return {"command": ssh, "args": ["-T", "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=yes",
                                          args.ssh, command]}
    if args.remote_checkout:
        raise ValueError("--remote-checkout requires --ssh")
    if args.wsl:
        if os.name != "nt":
            raise ValueError("--wsl is for a Windows client connecting to a Linux execution host")
        prefix = ["wsl.exe", "-d", args.wsl, "--exec"]
        script = subprocess.check_output(prefix + ["wslpath", "-a", "-u", Path(__file__).resolve().as_posix()],
                                         encoding="utf-8").strip()
        project_args = ["--demo"] if args.demo else ["--project", args.project]
        config = json.loads(subprocess.check_output(
            prefix + [args.linux_python, script, "config", *project_args], text=True))
        launch = config["mcpServers"]["kilnbio"]
        return {"command": prefix[0], "args": prefix[1:] + [launch["command"], *launch["args"]]}
    if os.name != "posix" or sys.version_info[:2] < (3, 12) or not shutil.which("bash"):
        raise ValueError("Use Python 3.12+ and Bash on Linux/macOS, or select a Linux distro with --wsl.")
    project = Path(args.project).expanduser().resolve() if args.project else Path.home() / ".kilnbio/assistant-demo"
    if args.demo and not project.exists():
        shutil.copytree(ROOT / "examples/assistant-project/pipelines", project / "pipelines")
    if not (project / "pipelines").is_dir():
        raise ValueError("The project needs a pipelines directory. Use --demo for a synthetic example.")
    managed = os.environ.get("KILNBIO_MANAGED_LAUNCHER")
    if managed:
        if not Path(managed).is_file() or importlib.util.find_spec("mcp") is None:
            raise ValueError("Install the assistant profile along with your analysis profiles before registering a client.")
        # Keep R, Python and sequencing dependencies available when the client reconnects.
        return {"command": managed, "args": ["assistant-serve", "--project", str(project)]}
    python = ROOT / ".venv/bin/python"
    if not python.is_file():
        venv.EnvBuilder(with_pip=True).create(python.parents[1])
    subprocess.run([str(python), "-m", "pip", "install", "-q", "-r", str(ROOT / "assistants/requirements.txt")],
                   stdout=sys.stderr, check=True)
    return {"command": str(python), "args": [str(ROOT / "assistants/kilnbio_mcp.py"), "--project", str(project)]}


def install_antigravity(name, entry):
    path = Path.home() / ".gemini/config/mcp_config.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    # Concurrent installers must preserve one another's registrations.
    lock = path.with_suffix(".kilnbio-lock")
    lock.open("x").close()
    try:
        config = json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else {}
        if not isinstance(config, dict) or not isinstance(config.get("mcpServers", {}), dict):
            raise ValueError("Expected an MCP configuration object with an mcpServers map")
        servers = config.setdefault("mcpServers", {})
        if name in servers and servers[name] != entry:
            raise ValueError("A different server already uses this name. Select another --name.")
        servers[name] = entry
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False, encoding="utf-8") as f:
            json.dump(config, f, indent=2)
            f.write("\n")
        try:
            os.replace(f.name, path)
        finally:
            Path(f.name).unlink(missing_ok=True)
    finally:
        lock.unlink()
    print(f"Registered {name} in {path}. Reload the Antigravity MCP servers.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("client", choices=["codex", "claude", "antigravity", "config"])
    project = parser.add_mutually_exclusive_group(required=True)
    project.add_argument("--demo", action="store_true", help="Set up the supplied synthetic capture study")
    project.add_argument("--project", help="Approved project directory on the execution host")
    parser.add_argument("--name", default="kilnbio", help="Server registration name")
    host = parser.add_mutually_exclusive_group()
    host.add_argument("--wsl", metavar="DISTRO", help="Use this installed WSL distro from a Windows client")
    host.add_argument("--ssh", metavar="ALIAS", help="Use the lab's existing Linux VM via its OpenSSH host alias")
    parser.add_argument("--remote-checkout", help="Absolute path to the prepared KilnBio checkout on that VM")
    parser.add_argument("--linux-python", default="python3.12", help="Python executable inside the WSL distro")
    args = parser.parse_args()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", args.name):
        parser.error("Use a short lowercase server name")
    try:
        program = shutil.which(args.client) if args.client in {"codex", "claude"} else None
        if args.client in {"codex", "claude"} and not program:
            raise ValueError(f"Install the {args.client} CLI on this client first")
        entry = prepare(args)
        if args.client == "config":
            print(json.dumps({"mcpServers": {args.name: entry}}, indent=2))
        elif args.client == "antigravity":
            install_antigravity(args.name, entry)
        else:
            flags = ["--transport", "stdio"] if args.client == "claude" else []
            subprocess.run([program, "mcp", "add", *flags, args.name, "--", entry["command"], *entry["args"]],
                           check=True)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"KilnBio setup: {exc}\n")


if __name__ == "__main__":
    main()
