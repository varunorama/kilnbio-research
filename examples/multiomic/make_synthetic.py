"""Make a small fictional genotype/RNA/Olink/methylation study for learning the CLI.

These are simulated measurements, not an exome accuracy or biological benchmark.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

def make(out):
    out=Path(out)
    if out.exists(): raise ValueError('Choose a new directory')
    out.mkdir(parents=True)
    rng=np.random.default_rng(391); n=60; f=200
    ids=[f'D{i:02}' for i in range(n)]; genes=[f'SIM_G{i:03}' for i in range(f)]
    carrier=rng.integers(0,2,n); age=rng.normal(50,10,n)
    donors=pd.DataFrame(dict(donor=ids,carrier=carrier,age=age))
    donors.to_csv(out/'donors.tsv',sep='\t',index=False)
    pd.DataFrame(dict(specimen=ids,donor=ids)).to_csv(out/'specimens.tsv',sep='\t',index=False)
    pd.DataFrame(dict(sample=ids,carrier=carrier)).to_csv(out/'olink_samples.tsv',sep='\t',index=False)
    effect=np.zeros(f); effect[:8]=1.5
    # Protein, RNA and methylation share a contrast but retain their own scales.
    protein=7+rng.normal(size=(f,n)) + effect[:,None]*carrier
    mean=60*np.exp(effect[:,None]*carrier)
    counts=rng.negative_binomial(10,10/(10+mean))
    mvalue=rng.normal(size=(f,n))-effect[:,None]*carrier
    beta=1/(1+2**(-mvalue))
    for name,x,features in [('rna',counts,genes),('methylation',beta,[f'cg_sim_{i:03}' for i in range(f)])]:
        pd.DataFrame(x,index=features,columns=ids).rename_axis('feature').to_csv(out/(name+'.tsv'),sep='\t')
    pd.DataFrame(dict(feature=[f'cg_sim_{i:03}' for i in range(f)],gene=genes)).to_csv(out/'methylation_map.tsv',sep='\t',index=False)
    oid=[f'OID_SIM_{i:03}' for i in range(f)]
    pd.DataFrame(dict(feature=oid,gene=genes)).to_csv(out/'olink_map.tsv',sep='\t',index=False)
    rows=[dict(SampleID=s,OlinkID=oid[i],NPX=protein[i,j],SampleType='SAMPLE',AssayType='assay',SampleQC='PASS',AssayQC='PASS',PlateID='SIM_P1',WellID=f'A{j+1}') for i in range(f) for j,s in enumerate(ids)]
    pd.DataFrame(rows).to_csv(out/'olink_export.tsv',sep='\t',index=False)
    layers=[dict(name='rna',matrix='rna.tsv',scale='counts',feature_map='identity',orientation=1),
            dict(name='olink',matrix='olink/eligible_npx.tsv',scale='log2',feature_map='olink_map.tsv',orientation=1),
            dict(name='methylation',matrix='methylation.tsv',scale='beta',feature_map='methylation_map.tsv',orientation=-1)]
    for layer in layers: layer['specimens']='specimens.tsv'
    cfg=dict(donors='donors.tsv',estimand='association',trait='carrier',covariates=['age'],cohort_policy='matched',layers=layers)
    (out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (out/'README.txt').write_text('Synthetic teaching data only. First 8 genes have simulated associations.\nOrientation -1 was chosen before simulation for this illustrative methylation layer; it is not a general biological rule.\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out');make(p.parse_args().out)
