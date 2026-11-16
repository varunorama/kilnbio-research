"""Run operator-installed study manifests in one approved project (POSIX)."""
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading

from mcp import types
from mcp.server.mcpserver.exceptions import ToolError


STUDY = Path(__file__).resolve().parents[1] / "downstream" / "study.py"
ID = re.compile(r"[a-z0-9][a-z0-9_-]{0,63}")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def register(server, project):
    project = project.resolve(strict=True)
    catalog = project / "pipelines"
    runs = project / "runs"
    if not catalog.is_dir() or catalog.is_symlink() or runs.is_symlink():
        raise ValueError("Use a project with a local pipelines directory and a local runs directory")
    runs.mkdir(mode=0o700, exist_ok=True)

    def manifest(name):
        if not ID.fullmatch(name):
            raise ToolError("Invalid pipeline ID")
        folder = catalog / name
        path = folder / "study.json"
        if folder.is_symlink() or path.is_symlink() or not path.is_file():
            raise ToolError("Unknown pipeline. Choose an ID from list_pipelines.")
        return path

    def job_dir(run_id):
        if not ID.fullmatch(run_id):
            raise ToolError("Use a run ID of 1–64 lowercase letters, digits, hyphens or underscores, starting with a letter or digit.")
        job = runs / run_id
        if job.is_symlink():
            raise ToolError("Invalid run directory")
        return job

    read_only = types.ToolAnnotations(read_only_hint=True, destructive_hint=False,
                                      idempotent_hint=True, open_world_hint=False)

    @server.tool(annotations=read_only)
    def list_pipelines() -> dict[str, object]:
        """List this project's installed pipeline IDs, stages and remote-service setting."""
        entries = {}
        unavailable = []
        for folder in sorted(catalog.iterdir()):
            if not ID.fullmatch(folder.name) or folder.is_symlink() or not (folder / "study.json").is_file():
                continue
            try:
                cfg = read_json(manifest(folder.name))
                if cfg["schema_version"] != 1 or cfg["research_use"] is not True or not cfg["steps"]:
                    raise ValueError("Invalid study")
                entries[folder.name] = {
                    "stages": [{"id": step["id"], "command": step["command"]} for step in cfg["steps"]],
                    "allows_remote_scientific_services": cfg.get("allow_remote_scientific_services", False),
                }
            except (OSError, ValueError, TypeError, KeyError, ToolError):
                unavailable.append(folder.name)
        return {"pipelines": entries, "unavailable": unavailable}

    @server.tool(annotations=read_only)
    def pipeline_status(run_id: str) -> dict[str, object]:
        """Check execution and return stage status and completed artifact names and hashes."""
        job = job_dir(run_id)
        if not (job / "request.json").is_file():
            raise ToolError("Unknown run ID for this project")
        result = read_json(job / "request.json")
        # The runner inherits this flock. It survives an MCP restart and is released
        # by the OS on exit, even if the runner never wrote its final receipt.
        with (job / "active.lock").open("rb") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
                active = False
            except BlockingIOError:
                active = True
            statefile = job / "study" / "study-state.json"
            state = read_json(statefile) if statefile.is_file() else {}
        status = state.get("status", "starting")
        if not active and status not in {"complete", "failed"}:
            status = "interrupted" if state else "failed"
        result.update(status=status, steps={
            name: {key: rec[key] for key in ("status", "returncode", "outputs_sha256") if key in rec}
            for name, rec in state.get("steps", {}).items()
        }, output_directory=f"runs/{run_id}/study", launcher_log=f"runs/{run_id}/launcher.log")
        return result

    @server.tool(annotations=types.ToolAnnotations(read_only_hint=False, destructive_hint=True,
                                                   idempotent_hint=True, open_world_hint=True))
    def run_pipeline(pipeline_id: str, run_id: str) -> dict[str, object]:
        """Launch an installed study. Reuse the same project-unique run_id on retry.
        Poll pipeline_status for receipts. Inputs, outputs and commands are set on the host.
        """
        job = job_dir(run_id)
        with (runs / ".submit.lock").open("a+b") as submission:
            fcntl.flock(submission, fcntl.LOCK_EX)
            if job.exists():
                existing = pipeline_status(run_id)
                if existing["pipeline_id"] != pipeline_id:
                    raise ToolError("This run ID already belongs to a different pipeline")
                return existing
            config = manifest(pipeline_id)
            job.mkdir(mode=0o700)
            with (job / "active.lock").open("wb") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                request = {"pipeline_id": pipeline_id, "run_id": run_id}
                tmp = job / "request.tmp"
                tmp.write_text(json.dumps(request), encoding="utf-8")
                tmp.replace(job / "request.json")
                env = os.environ.copy()
                env.setdefault("KILNBIO_PYTHON", sys.executable)
                with (job / "launcher.log").open("wb") as log:
                    try:
                        process = subprocess.Popen(
                            [sys.executable, str(STUDY), str(config), str(job / "study")],
                            cwd=project, env=env, stdin=subprocess.DEVNULL, stdout=log,
                            stderr=subprocess.STDOUT, start_new_session=True, pass_fds=(lock.fileno(),),
                        )
                    except OSError as exc:
                        log.write(f"Could not start study runner: {exc}\n".encode())
                        raise ToolError(f"Launch failed. See runs/{run_id}/launcher.log on the host.") from exc
                threading.Thread(target=process.wait, daemon=True).start()
        return pipeline_status(run_id)
