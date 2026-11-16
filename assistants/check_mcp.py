"""Exercise the optional KilnBio integration through the official MCP client."""
import argparse
import asyncio
from pathlib import Path
import sys

from mcp.client import Client
from mcp.client.stdio import StdioServerParameters


async def check(target):
    async with Client(target) as client:
        tools = (await client.list_tools()).tools
        assert {tool.name for tool in tools} == {
            "list_guides", "read_guide", "assay_reach", "compare_methods"
        }
        assert all(tool.annotations.read_only_hint for tool in tools)
        catalog = await client.call_tool("list_guides", {})
        assert "study" in catalog.structured_content["guides"]
        guide = await client.call_tool("read_guide", {"guide": "study"})
        assert "--dry-run" in guide.content[0].text
        assert (await client.call_tool("read_guide", {"guide": "../LICENSE"})).is_error

        one = {"donor": "example-1", "base_ok": True, "add_ok": True}
        two = {"donor": "example-2", "base_ok": True, "add_ok": False}
        reach = await client.call_tool("assay_reach", {"roster": [one, two]})
        assert reach.structured_content["base_reach"] == 2
        assert reach.structured_content["combined_reach"] == 1
        assert (await client.call_tool("assay_reach", {"roster": [one, one]})).is_error

        manifest = ("method\tstudy\tcontrast\tfeature_family\tfeature_unit\tassembly\tchannel\n"
                    "A\tstudy1\tcase-control\tgenes\tgene\tGRCh38\tRNA\n"
                    "B\tstudy1\tcase-control\tgenes\tgene\tGRCh38\tRNA\n")
        matrix = "feature\tA_effect\tA_q\tB_effect\tB_q\nX\t1\t0.01\t-1\t0.01\nY\t1\t0.01\tNA\tNA\n"
        result = await client.call_tool("compare_methods", {"matrix_tsv": matrix, "manifest_tsv": manifest})
        data = result.structured_content
        x, y = [dict(zip(data["columns"], row)) for row in data["rows"]]
        assert (x["status"], y["planned"], y["tested"]) == ("opposite_calls", 2, 1)
        mixed = manifest.replace("B\tstudy1", "B\tstudy2")
        assert (await client.call_tool("compare_methods", {"matrix_tsv": matrix, "manifest_tsv": mixed})).is_error


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", help="Test a running Streamable HTTP server; default is a new stdio process")
    args = parser.parse_args()
    target = args.url or StdioServerParameters(command=sys.executable, args=[str(Path(__file__).with_name("kilnbio_mcp.py").resolve())])
    asyncio.run(asyncio.wait_for(check(target), timeout=45))
    print("KilnBio MCP checks passed.")
