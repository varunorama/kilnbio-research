"""Exercise real project runs through stdio MCP in a disposable synthetic project."""
import asyncio
import hashlib
import json
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile

from mcp.client import Client
from mcp.client.stdio import StdioServerParameters


ROOT = Path(__file__).resolve().parents[1]


async def result(client, name, **kwargs):
    reply = await client.call_tool(name, kwargs)
    assert not reply.is_error, reply.content
    return reply.structured_content


async def finish(client, run_id):
    for _ in range(300):
        receipt = await result(client, "pipeline_status", run_id=run_id)
        if receipt["status"] in {"complete", "failed", "interrupted"}:
            return receipt
        await asyncio.sleep(0.1)
    raise AssertionError("The synthetic study did not finish in time")


async def check(project, target):
    pipeline = project / "pipelines/capture-demo"
    async with Client(target) as client:
        assert set((await result(client, "list_pipelines"))["pipelines"]) == {"capture-demo"}
        await asyncio.gather(*(result(client, "run_pipeline", pipeline_id="capture-demo", run_id="first")
                               for _ in range(2)))
        assert (await client.call_tool("run_pipeline", {"pipeline_id": "../other", "run_id": "nope"})).is_error
        assert (await client.call_tool("pipeline_status", {"run_id": "../other"})).is_error
    # Simulate closing the assistant and reconnecting.
    async with Client(target) as client:
        receipt = await finish(client, "first")
        assert receipt["status"] == "complete", receipt
        output = project / receipt["output_directory"]
        summary = json.loads((output / "panel/summary.json").read_text())
        assert summary["target_union_bases"] == 80
        assert summary["features"]["tested"] == 3 and summary["features"]["fully_on_target"] == 1
        assert receipt["steps"]["panel"]["outputs_sha256"]["panel/features.tsv"] == hashlib.sha256(
            (output / "panel/features.tsv").read_bytes()).hexdigest()
        assert (await result(client, "run_pipeline", pipeline_id="capture-demo", run_id="first")) == receipt
        assert len(list((output / "logs").iterdir())) == 1

        # Add a two-stage study while the server is still running.
        second = project / "pipelines/two-panels"
        shutil.copytree(pipeline, second)
        manifest = second / "study.json"
        cfg = json.loads(manifest.read_text())
        step = dict(cfg["steps"][0], id="repeat", depends_on=["panel"])
        step["args"] = [s.replace("{{run}}/panel", "{{run}}/repeat") for s in step["args"]]
        step["outputs"] = [s.replace("panel/", "repeat/") for s in step["outputs"]]
        cfg["steps"].append(step)
        manifest.write_text(json.dumps(cfg))
        assert "two-panels" in (await result(client, "list_pipelines"))["pipelines"]
        assert (await client.call_tool("run_pipeline", {"pipeline_id": "two-panels", "run_id": "first"})).is_error
        await result(client, "run_pipeline", pipeline_id="two-panels", run_id="second")
        again = await finish(client, "second")
        assert again["status"] == "complete" and set(again["steps"]) == {"panel", "repeat"}
        (second / "sizes.tsv").unlink()
        await result(client, "run_pipeline", pipeline_id="two-panels", run_id="missing-input")
        assert (await finish(client, "missing-input"))["status"] == "failed"
        assert "Input missing" in (project / "runs/missing-input/launcher.log").read_text()
        print("MCP execution verified: 1 of 3 features fully on target, 80 target bases;"
              " two-stage discovery, duplicate submission, reconnection and failed input checks passed.")


async def exercise(project, transport):
    shutil.copytree(ROOT / "examples/assistant-project/pipelines", project / "pipelines")
    args = [str(ROOT / "assistants/kilnbio_mcp.py"), "--project", str(project)]
    if transport == "stdio":
        await check(project, StdioServerParameters(command=sys.executable, args=args))
        return
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    process = subprocess.Popen([sys.executable, *args, "--transport", "streamable-http", "--port", str(port)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(100):
            if process.poll() is not None:
                raise RuntimeError("Loopback server exited before initialization")
            try:
                with socket.create_connection(("127.0.0.1", port), timeout=0.1):
                    break
            except OSError:
                await asyncio.sleep(0.1)
        await check(project, f"http://127.0.0.1:{port}/mcp")
    finally:
        process.terminate()
        process.wait(timeout=10)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="kilnbio-mcp-") as temp:
        for transport in ("stdio", "http"):
            print(f"Testing {transport}", flush=True)
            asyncio.run(asyncio.wait_for(exercise(Path(temp) / f"{transport} example", transport), timeout=120))
