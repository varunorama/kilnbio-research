"""Check within-library ERCC Mix 1 dose ordering; never scale endogenous genes here."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SOURCE='https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
EXPECTED='5151b394332477301851404d37b907cb6f27c9aaa2569e288531a839d6ca2b3b'

def evaluate(table,genes,analysis,out):
    if hashlib.sha256(table.read_bytes()).hexdigest()!=EXPECTED:raise ValueError('Concentration source differs from frozen manufacturer file; review first')
    mix=pd.read_csv(table,sep='\t').set_index('ERCC ID')['concentration in Mix 1 (attomoles/ul)']
    if len(mix)!=92 or not mix.index.is_unique or not (mix>0).all():raise ValueError('Incomplete or invalid concentration table')
    qc=pd.read_csv(analysis/'samples.csv',dtype=str)
    if len(qc)!=8 or not qc['sample'].is_unique or (qc.geo_ercc_mix=='1').sum()!=3:raise ValueError('Unexpected airway design')
    counts=pd.read_csv(analysis/'featurecounts/counts.tsv',sep='\t').set_index('feature_id')
    lengths=pd.read_csv(analysis/'featurecounts/featureCounts.tsv',sep='\t',comment='#',usecols=['Geneid','Length']).set_index('Geneid').Length
    idx=mix.index
    if not (lengths.loc[idx]>0).all():raise ValueError('Missing annotated ERCC lengths')
    per=[];summary=[]
    for s,flag in zip(qc['sample'],qc.geo_ercc_mix):
        if flag not in ['1','-']:raise ValueError('Unexpected ERCC assignment')
        rsem=pd.read_csv(genes/(s+'.genes.results'),sep='\t').set_index('gene_id')
        for index in [counts.index,rsem.index,lengths.index]:
            if set(index[index.str.startswith('ERCC-')])!=set(idx) or not index.is_unique:raise ValueError('Measured ERCC IDs do not match the reference table')
        tpm=rsem.loc[idx,'TPM'];fc=counts.loc[idx,s];rpk=1000*fc/lengths.loc[idx]
        m={'sample':s,'declared_mix':flag,'rsem_detected':int((tpm>0).sum()),'fc_detected':int((fc>0).sum())}
        for k,v in [('rsem_TPM',tpm),('featureCounts_per_kb',rpk)]:
            if not np.isfinite(v).all() or (v<0).any():raise ValueError('Invalid control abundance')
            m[k+'_spearman_all92']=float(spearmanr(mix,v).statistic) if flag=='1' and v.nunique()>1 else None
        summary.append(m)
        for gene in idx:
            per.append(dict(sample=s,declared_mix=flag,ERCC_ID=gene,stock_Mix1_attomoles_per_ul=mix[gene],
                RSEM_TPM=tpm[gene],featureCounts_fragments=fc[gene],featureCounts_per_kb=rpk[gene]))
    out.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(per).to_csv(out/'controls.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'by_sample.tsv',sep='\t',index=False,na_rep='NA')
    (out/'provenance.json').write_text(json.dumps({'source':SOURCE,'sha256':EXPECTED,
        'scope':'Secondary control: within-sample Spearman correlation of all 92 measured controls including zeros against stock molarity in declared mix-1-positive libraries only. FeatureCounts is length-scaled; RSEM uses its reported TPM. Stock molarity is not final per-cell abundance. No post hoc pass cutoff or endogenous normalization.'},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['concentrations','full_rsem_folder','analysis','out']:p.add_argument(k,type=Path)
    a=p.parse_args();evaluate(a.concentrations,a.full_rsem_folder,a.analysis,a.out)
