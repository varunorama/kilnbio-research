"""Cross-check the airway package, GEO sample records and ENA run records."""
import argparse,csv,re
from pathlib import Path

def build(soft,ena,metadata,out):
    geo={r.splitlines()[0]:r for r in re.split(r'(?m)^\^SAMPLE = ',soft.read_text(encoding='utf-8-sig'))[1:]}
    runs={r['run_accession']:r for r in csv.DictReader(ena.open(),delimiter='\t')}
    records=[]
    for m in csv.DictReader(metadata.open(encoding='utf-8-sig')):
        g=geo[m['SampleName']];e=runs[m['Run']]
        attr=dict(re.findall(r'^!Sample_characteristics_ch1 = ([^:\n]+): ([^\n\r]+)',g,re.M))
        condition={'Untreated':'untrt','Dexamethasone':'trt','Dex':'trt'}[attr['treatment']]
        if not (attr['cell line']==m['cell'] and condition==m['dex'] and m['albut']=='untrt' and
                m['Experiment']==e['experiment_accession'] and m['BioSample']==e['sample_accession'] and
                m['Experiment'] in g and m['BioSample'] in g and e['library_layout']=='PAIRED'):
            raise ValueError('GEO/ENA/airway metadata disagree: '+m['Run'])
        records.append(dict(sample=m['Run'],biological_id=m['cell'],condition=condition,geo_ercc_mix=attr['ercc_mix'],gsm=m['SampleName'],experiment=m['Experiment']))
    if len(records)!=8 or len({r['biological_id'] for r in records})!=4 or sum(r['geo_ercc_mix']=='1' for r in records)!=3:raise ValueError('Unexpected public design or spike assignment')
    with out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,records[0]);w.writeheader();w.writerows(records)
    print('Verified eight libraries, four donors, three documented ERCC mix-1 libraries')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['geo_soft','ena_tsv','airway_metadata','output']:p.add_argument(k,type=Path)
    a=p.parse_args();build(a.geo_soft,a.ena_tsv,a.airway_metadata,a.output)
