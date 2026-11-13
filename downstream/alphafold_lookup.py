"""Explicit AlphaFold Database lookup for research, not de novo structure inference."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys
import urllib.request
from datetime import datetime, timezone
from variant_annotations import digest, write_tsv

ACCESSION = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[1-9][0-9]*)?')


def map_protein(model, variant):
    if not variant:
        return 'structure_only', None
    m = re.fullmatch(r'([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])', variant)
    if not m or m[1] == m[3]:
        raise ValueError('protein_variant must be a one-letter missense change, e.g. R175H')
    position = int(m[2])
    sequence, start, end = model.get('sequence'), model.get('sequenceStart'), model.get('sequenceEnd')
    if not isinstance(sequence, str) or type(start) is not int or type(end) is not int or len(sequence) != end-start+1:
        return 'sequence_mapping_unavailable', None
    if not start <= position <= end:
        return 'outside_model_fragment', None
    if sequence[position-start] != m[1]:
        return 'reference_residue_mismatch', None
    return 'reference_residue_matches', position-start+1


def lookup(queries, output, execute=False, cache=None):
    if execute and cache:
        raise ValueError('Choose cached responses or explicit remote execution, not both')
    query = Path(queries)
    with query.open(encoding='utf-8-sig',newline='') as f:
        reader = csv.DictReader(f,delimiter='\t')
        if reader.fieldnames != ['uniprot_id','protein_variant']:
            raise ValueError('TSV needs uniprot_id and protein_variant; blank variant requests structure only')
        rows = list(reader)
    if not 1 <= len(rows) <= 1000:
        raise ValueError('Use 1..1000 protein rows')
    for row in rows:
        if not ACCESSION.fullmatch(row['uniprot_id'] or ''):
            raise ValueError('Invalid UniProt accession')
        map_protein({}, row['protein_variant'])
    if len({tuple(r.values()) for r in rows}) != len(rows):
        raise ValueError('Duplicate protein request')
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    receipt = {'status':'planned','research_use_only':True,'request_sha256':digest(query),
               'query_rows':len(rows),'recorded_utc':datetime.now(timezone.utc).isoformat(),
               'remote_execution':execute,'new_structure_inference':False,'sources':[]}
    def save():
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    requests = {r['uniprot_id']:'https://alphafold.ebi.ac.uk/api/prediction/'+r['uniprot_id'] for r in rows}
    (out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n',encoding='utf-8')
    save()
    if not execute and not cache:
        return receipt
    receipt['status'] = 'running'; save()
    models, results = {}, []
    try:
        for accession, url in requests.items():
            if cache:
                raw = (Path(cache)/(accession+'.json')).read_bytes()
            else:
                req = urllib.request.Request(url,headers={'User-Agent':'KilnBio research annotation'})
                with urllib.request.urlopen(req,timeout=60) as response:
                    raw = response.read(10*1024*1024+1)
                    if len(raw)>10*1024*1024:
                        raise ValueError('Unexpectedly large AlphaFold metadata response')
            records = json.loads(raw)
            if not isinstance(records,list) or any(not isinstance(r,dict) for r in records):
                raise ValueError('Unexpected AlphaFold response schema')
            local = out/(accession+'.json'); local.write_bytes(raw)
            models[accession] = [r for r in records if r.get('uniprotAccession') == accession]
            receipt['sources'].append({'uniprot_id':accession,'url':url,'sha256':digest(local),
                'cached_response':bool(cache),'matching_models':len(models[accession])})
            save()
        for row in rows:
            selected = models[row['uniprot_id']]
            if not selected:
                results.append({**row,'mapping_status':'no_matching_model'})
            for model in selected:
                status, residue = map_protein(model,row['protein_variant'])
                results.append({**row,'model_id':model.get('modelEntityId',''),
                    'model_version':model.get('latestVersion',''),'tool_used':model.get('toolUsed',''),
                    'mapping_status':status,'model_sequence_position':residue or '',
                    'sequence_start':model.get('sequenceStart',''),'sequence_end':model.get('sequenceEnd',''),
                    'mean_plddt':model.get('globalMetricValue',''), 'pdb_url':model.get('pdbUrl',''),
                    'cif_url':model.get('cifUrl','')})
        fields = ['uniprot_id','protein_variant','model_id','model_version','tool_used','mapping_status',
                  'model_sequence_position','sequence_start','sequence_end','mean_plddt','pdb_url','cif_url']
        write_tsv(out/'structures.tsv',results,fields)
        receipt.update(status='complete',output_rows=len(results))
        save()
    except Exception as exc:
        receipt.update(status='failed',error=str(exc)); save(); raise
    return receipt


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('proteins');p.add_argument('output');p.add_argument('--execute',action='store_true')
    p.add_argument('--cache',help='Directory with accession.json API responses for offline mapping')
    a=p.parse_args()
    try:
        print(json.dumps(lookup(a.proteins,a.output,a.execute,a.cache),indent=2))
    except (OSError,ValueError,KeyError) as exc:
        sys.exit(str(exc))
