import json,os,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
def run(*a):
    print('RUN',a,flush=True);subprocess.run([str(x) for x in a],check=True)
run(sys.executable,root/'examples/multiomic/make_study_demo.py',out/'matrix_inputs')
run('bash',root/'run.sh','study',out/'matrix_inputs/study.json',out/'matrix_run')
inp=out/'inputs';inp.mkdir();rng=np.random.default_rng(3381);samples=[f'S{i}' for i in range(8)]
meta=pd.DataFrame({'sample':samples,'biological_id':samples,'condition':['A']*4+['B']*4});meta.to_csv(inp/'samples.csv',index=False)
counts=rng.negative_binomial(8,.12,(1200,8));counts[:60,4:]*=3
x=pd.DataFrame(counts,columns=samples);x.insert(0,'feature_id',[f'SYN{i}' for i in range(1200)]);x.to_csv(inp/'counts.tsv',sep='\t',index=False)
for engine in ['edger','deseq2']:
    cfg=dict(samples='samples.csv',numerator='B',denominator='A',assay='rnaseq',engine=engine,count_source='raw_counts',counts='counts.tsv');cp=inp/(engine+'.json');cp.write_text(json.dumps(cfg));run('bash',root/'run.sh','de',cp,out/engine)
manifest=[]
for i,s in enumerate(samples):
    n=rng.integers(20,50,400);fraction=np.full(400,.45);fraction[:30]=.25 if i<4 else .8
    pd.DataFrame({'chr':['chr1']*400,'pos':np.arange(400)*20+1000,'N':n,'X':rng.binomial(n,fraction)}).to_csv(inp/(s+'.tsv'),sep='\t',index=False);manifest.append({'sample':s,'counts':s+'.tsv'})
pd.DataFrame(manifest).to_csv(inp/'methylation.csv',index=False)
cfg=dict(samples='samples.csv',numerator='B',denominator='A',engine='dss',modification='modC',assembly='GRCh38',coordinates='1-based-cytosine',measurement_reviewed=True,methylation_manifest='methylation.csv',smoothing=False)
cp=inp/'dm.json';cp.write_text(json.dumps(cfg));run('bash',root/'run.sh','dm',cp,out/'dss')
print('ALL REAL FIT REPORTS COMPLETE',flush=True)
