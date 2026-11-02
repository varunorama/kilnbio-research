"""Resolve pinned file accessions to public source URLs and verify experimental roles."""
import argparse
import csv
import json
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def get_experiment(cache, accession):
    path=cache/(accession+'.json')
    if not path.exists():
        url=f'https://www.encodeproject.org/experiments/{accession}/?format=json'
        with urllib.request.urlopen(url,timeout=60) as source:
            data=json.load(source)
        path.write_text(json.dumps(data)+'\n')
    data=json.loads(path.read_text())
    assert data['accession']==accession and data['status']=='released'
    return data


def build(cache):
    cache.mkdir(parents=True,exist_ok=True)
    output=[]
    with (ROOT/'experiments.tsv').open() as f:
        experiments=list(csv.DictReader(f,delimiter='\t'))
    for row in experiments:
        treatment=get_experiment(cache,row['experiment'])
        assert treatment['biosample_ontology']['term_name']==row['cell']
        atac=row['assay']=='ATAC'
        if not atac:
            assert treatment['target']['label']==row['assay']
            assert row['control_experiment'] in [p['accession'] for p in treatment['possible_controls']]
        for key in ('rep1','rep2','controls','dnase_rep1','dnase_rep2'):
            if key.startswith('dnase') and not atac:
                continue
            for acc in filter(None,row[key].split(',')):
                control=key=='controls'
                dnase=key.startswith('dnase')
                exp=row['control_experiment'] if control else row['dnase_experiment'] if dnase else row['experiment']
                data=get_experiment(cache,exp)
                assert data['biosample_ontology']['term_name']==row['cell']
                file=next(f for f in data['files'] if f['accession']==acc)
                assert file['status']=='released' and file['assembly']=='GRCh38'
                if not control:
                    assert file['biological_replicates']==[int(key[-1])]
                if dnase:
                    assert data['assay_title']=='DNase-seq' and file['file_format']=='bed'
                    assert '.0.001.peaks.' in file['submitted_file_name']
                else:
                    assert file['output_type']=='alignments' and file['file_format']=='bam'
                    assert file['mapped_run_type']==('paired-ended' if atac else 'single-ended')
                if not atac and not control:
                    parents={f'/files/{i}/' for i in [acc,*row['controls'].split(',')]}
                    assert any(set(x['derived_from'])==parents for x in treatment['files']
                               if x.get('output_type')=='fold change over control' and x.get('assembly')=='GRCh38')
                depths = {q['usable_fragments'] for q in file.get('quality_metrics', [])
                          if isinstance(q, dict) and q.get('usable_fragments') is not None}
                warnings = {a['category'] for level in ('WARNING','NOT_COMPLIANT','ERROR')
                            for a in treatment.get('audit', {}).get(level, [])
                            if f'/files/{acc}/' in a.get('detail','')}
                output.append(dict(panel=row['experiment'],role=key,cell=row['cell'],assay=row['assay'],
                                   source_experiment=exp,file_accession=acc,assembly=file['assembly'],
                                   bytes=file['file_size'],md5=file['md5sum'],url=file['cloud_metadata']['url'],
                                   encode_usable_fragments=','.join(str(x) for x in sorted(depths)) or 'NA',
                                   treatment_audit_flags_naming_file='|'.join(sorted(warnings)) or 'none listed'))
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('metadata',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    records=build(args.metadata)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,records[0].keys(),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
