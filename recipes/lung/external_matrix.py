import csv,gzip,json
from pathlib import Path
import pandas as pd,numpy as np
b=Path(__file__).resolve().parent
o=b/'analysis';e=b/'external'
genes=set(pd.read_csv(o/'RNA.tsv',sep='\t',usecols=['gene']).gene)
with gzip.open(e/'GPL570.annot.gz','rt',encoding='utf-8',errors='replace') as f:
    for line in f:
        if line.startswith('!platform_table_begin'):break
    anno={}
    for row in csv.DictReader(f,delimiter='\t'):
        gene=(row['Gene symbol'] or '').strip()
        if gene in genes:anno[row['ID']]=gene
with gzip.open(e/'GSE31210_series_matrix.txt.gz','rt',encoding='utf-8') as f:
    data=[]; probes=[]; titles=accessions=tissue=status=None
    for row in csv.reader(f,delimiter='\t'):
        if not row:continue
        k,*v=row
        if k=='!series_matrix_table_begin':break
        if k=='!Sample_title':titles=v
        if k=='!Sample_geo_accession':accessions=v
        if k=='!Sample_characteristics_ch1' and any(p.startswith('tissue:') for p in v):tissue=v
        if k=='!Sample_characteristics_ch1' and any(p.startswith('gene alteration status:') for p in v):status=v
    assert all(isinstance(x,list) for x in (titles,accessions,tissue,status))
    header=next(csv.reader(f,delimiter='\t'))[1:]
    assert header==accessions
    for row in csv.reader(f,delimiter='\t'):
        if row and row[0] in anno:
            probes.append(row[0]);data.append([float(x) if x not in ('null','NA','') else np.nan for x in row[1:]])
meta=pd.DataFrame({'sample':accessions,'patient_label':titles,'tissue':tissue,'status':status}).set_index('sample')
ext=meta[meta.tissue.eq('tissue: primary lung tumor')].copy()
v={'gene alteration status: EGFR mutation +':1,'gene alteration status: EGFR/KRAS/ALK -':0,'gene alteration status: KRAS mutation +':0,'gene alteration status: ALK-fusion +':0}
assert set(ext.status)==set(v) and ext.patient_label.is_unique
ext['egfr']=ext.status.map(v)
x=pd.DataFrame(data,index=probes,columns=accessions)[ext.index]
assert x.min().min()>0 and x.stack().quantile(.99)>100
expr=np.log2(x).groupby([anno[i] for i in x.index]).median()
eligible=expr.std(axis=1)>0
for g in (0,1):eligible &= expr[ext[ext.egfr==g].index].notna().mean(axis=1)>=0.8
expr[eligible].rename_axis('gene').to_csv(o/'external_expression.tsv',sep='\t',na_rep='NA')
ext.to_csv(o/'external_patients.tsv',sep='\t')
pd.DataFrame(sorted(anno.items()),columns=['probe','symbol']).to_csv(o/'external_probe_map.tsv',sep='\t',index=False)
print('tumors',len(ext),'groups',ext.egfr.value_counts().to_dict(),'eligible genes',int(eligible.sum()),'probes',len(probes))


