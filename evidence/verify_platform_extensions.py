"""Refresh public method metadata through the read-only PubMed MCP adapter."""
import asyncio
import datetime
import json
import sys
from pathlib import Path
from mcp import Client,StdioServerParameters

PAPERS={
    '38666213':('10.1093/nargab/lqae031','Sarek preprocessing, not independent exome accuracy'),
    '39251642':('10.1038/s41467-024-52302-9','Public CH screening methods; subset implemented, not all rescue rules'),
    '39631069':('10.1182/bloodadvances.2024014467','Candidate public CH/Olink/multi-omic demonstration, not a KilnBio result'),
    '18261164':('10.1111/j.1541-0420.2007.00984.x','Partial-conjunction hypothesis framework'),
    '33012899':('10.1080/01621459.2018.1554485','Exploratory ACAT global-null combination'),
    '25605792':('10.1093/nar/gkv007','limma differential modeling'),
    '24485249':('10.1186/gb-2014-15-2-r29','voom RNA count modeling'),
    '30858580':('10.1038/s41587-019-0054-x','GA4GH variant benchmarking design')}

async def main(out):
    async with Client(StdioServerParameters(command=sys.executable,args=[str(Path(__file__).with_name('pubmed_mcp.py'))])) as client:
        result=await client.call_tool('pubmed_fetch',{'pmids':list(PAPERS)})
        data=result.structured_content
        if data is None: data=json.loads(next(c.text for c in result.content if getattr(c,'type',None)=='text'))
    records=data['records']
    if {r['pmid'] for r in records}!=set(PAPERS): raise ValueError('Incomplete PubMed response')
    for row in records:
        doi,role=PAPERS[row['pmid']]
        if row['doi']!=doi: raise ValueError(f'Identifier mismatch for {row["pmid"]}: {row["doi"]}')
        row['role']=role
        print(row['pmid'],row['title'],row['doi'])
    receipt=dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        transport='Actual MCP stdio -> KilnBio PubMed adapter -> NCBI EFetch',
        scope='Bibliographic verification and indexed notices; not certification of methods or novelty',records=records)
    Path(out).write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    if len(sys.argv)!=2: raise SystemExit('Expected output JSON path')
    asyncio.run(main(sys.argv[1]))
