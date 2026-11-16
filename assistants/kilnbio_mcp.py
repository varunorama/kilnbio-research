"""MCP tools for the KilnBio research checkout. Run this file with the optional SDK."""
import argparse
import sys
from io import StringIO
from pathlib import Path
from typing import Annotated

from mcp import types
from mcp.server import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import BaseModel, ConfigDict, Field, StrictBool

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from downstream import assay_audit, method_agreement

VERSION = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
GUIDES = {
    "commands": "CLI.txt",
    "sequencing": "docs/USER_GUIDE.md",
    "integration": "docs/INTEGRATING_OMICS.md",
    "study": "docs/STUDY_RUNNER.md",
    "recipes": "recipes/README.md",
    "olink": "docs/MULTIOMIC_AND_OLINK.md",
    "capture": "docs/CAPTURE_PANELS.md",
    "exome": "docs/EXOME_AND_CH.md",
    "access": "ACCESS.md",
}
server = MCPServer("KilnBio", version=VERSION)
READ_ONLY = types.ToolAnnotations(read_only_hint=True, destructive_hint=False,
                                  idempotent_hint=True, open_world_hint=False)


@server.tool(annotations=READ_ONLY)
def list_guides() -> dict[str, object]:
    """List the guides and command index shipped in this KilnBio checkout."""
    return {"version": VERSION, "guides": GUIDES.copy()}


@server.tool(annotations=READ_ONLY)
def read_guide(guide: str) -> str:
    """Read a KilnBio guide by its ID from list_guides."""
    if guide not in GUIDES:
        raise ToolError("Unknown guide. Choose an ID returned by list_guides.")
    return (ROOT / GUIDES[guide]).read_text(encoding="utf-8")


class DonorAvailability(BaseModel):
    model_config = ConfigDict(extra="forbid")
    donor: Annotated[str, Field(min_length=1, max_length=128)]
    base_ok: StrictBool
    add_ok: StrictBool


@server.tool(annotations=READ_ONLY)
def assay_reach(roster: Annotated[list[DonorAvailability], Field(min_length=1, max_length=500)]) -> dict[str, object]:
    """Count donors eligible for the base assay and for both assays. Supply one row per
    distinct donor, using pseudonyms and explicit true/false availability. This counts
    availability; specimen pairing requires the separate local assay-audit pair command."""
    rows = [{"donor": r.donor, "base_ok": str(int(r.base_ok)), "add_ok": str(int(r.add_ok))}
            for r in roster]
    try:
        return assay_audit.reach(rows)
    except ValueError as error:
        raise ToolError(str(error)) from error


SmallTable = Annotated[str, Field(min_length=1, max_length=24000)]


@server.tool(annotations=READ_ONLY)
def compare_methods(matrix_tsv: SmallTable, manifest_tsv: SmallTable,
                    alpha: Annotated[float, Field(gt=0, lt=1)] = 0.05) -> dict[str, object]:
    """Compare signed effects and adjusted P values from methods for the same contrast.
    Pass tab-separated text with headers. Matrix: feature and METHOD_effect, METHOD_q
    for each method. Manifest: method, study, contrast, feature_family, feature_unit,
    assembly, channel. Blank or NA marks missing estimates. Support is descriptive."""
    def table(text):
        return method_agreement.read_table(StringIO(text.removeprefix("\ufeff"), newline=""))

    try:
        columns, rows = method_agreement.compare_tables(table(matrix_tsv), table(manifest_tsv), alpha)
    except ValueError as error:
        raise ToolError(str(error)) from error
    return {"columns": columns, "rows": rows, "version": VERSION}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport", choices=("stdio", "streamable-http"), default="stdio")
    parser.add_argument("--port", type=int, default=8000, help="Loopback HTTP port")
    args = parser.parse_args()
    if args.transport == "stdio":
        server.run(transport="stdio")
    else:
        server.run(transport="streamable-http", host="127.0.0.1", port=args.port,
                   stateless_http=True, json_response=True, max_request_body_size=262144)
