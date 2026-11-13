"""Read-only local PubMed MCP server using the official MCP Python SDK v2.

This is a KilnBio adapter to NCBI, not a server operated or endorsed by NCBI.
Never submit participant/client information as a search query.
"""
from mcp.server import MCPServer
from pubmed import search, fetch

server = MCPServer('KilnBio PubMed evidence')

@server.tool()
def pubmed_search(query: str, limit: int = 10, offset: int = 0) -> dict:
    """Search public PubMed records. Do not include private study or patient data."""
    return search(query, limit, offset)

@server.tool()
def pubmed_fetch(pmids: list[str], include_abstract: bool = False) -> dict:
    """Fetch authoritative citation metadata and indexed notice/retraction links.

    Indexing does not prove peer review, correctness, or the absence of future notices.
    """
    return {'records':fetch(pmids, include_abstract), 'scope':'PubMed metadata, not independent scientific approval'}

if __name__ == '__main__':
    server.run(transport='stdio')
