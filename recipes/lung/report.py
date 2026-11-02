import json,subprocess,sys,hashlib
from pathlib import Path
import numpy as np,pandas as pd
root=Path(__file__).resolve().parent;o=root/'analysis';repo=Path(sys.argv[1]).resolve()
read=lambda n:pd.read_csv(o/n,sep='\t')
f=read('frozen.tsv');test=read('external_tests.tsv');boot=read('bootstrap.tsv')
lock=(o/'frozen.sha256').read_text(encoding='utf-8-sig').strip().lower()
assert hashlib.sha256((o/'frozen.tsv').read_bytes()).hexdigest()==lock
f=f.merge(test[['gene','effect','q']],how='left',on='gene',validate='many_to_one')
f['testable']=f.q.notna(); f['confirmed']=f.testable & (f.q<=.05) & (f.direction*f.effect>0)
f.to_csv(o/'target_decisions.tsv',sep='\t',index=False)
audit=o/'evidence';audit.mkdir(exist_ok=True)
pd.DataFrame([['RNA','GSE31210']],columns=['assay','cohort']).to_csv(audit/'endpoints.tsv',sep='\t',index=False)
for arm,frame in f.groupby('arm'):
    a=audit/arm;a.mkdir(exist_ok=True)
    nom=frame[['gene','direction']].copy();nom['assay']='RNA';nom.to_csv(a/'nominations.tsv',sep='\t',index=False)
    obs=frame[['gene','effect','q','testable']].copy();obs['testable']=obs.testable.astype(int);obs['assay']='RNA';obs['cohort']='GSE31210';obs.to_csv(a/'observations.tsv',sep='\t',index=False,na_rep='')
    subprocess.run([sys.executable,str(repo/'downstream/candidate_evidence.py'),str(a/'nominations.tsv'),str(audit/'endpoints.tsv'),str(a/'observations.tsv'),'--out',str(a/'counts.json')],check=True)
    actual={c['gene'] for c in json.loads((a/'counts.json').read_text())['candidates'] if c['endpoints'][0]['status']=='confirmed'}
    assert actual==set(frame.loc[frame.confirmed,'gene'])
count=f.groupby('arm').agg(nominated=('gene','size'),testable=('testable','sum'),confirmed=('confirmed','sum')).to_dict('index')
comparisons={}
for label,arm1,arm2 in [('restriction','A_RNA_all','B_RNA_paired'),('added_protein','B_RNA_paired','C_joint'),('both','A_RNA_all','C_joint'),('source_id_vs_audited','standard_exact','C_joint'),('planted_donor_join_vs_clean','C_joint','planted_donor_join'),('planted_withhold_vs_clean','C_joint','planted_withhold')]:
    a=f[f.arm==arm1].set_index('gene');b=f[f.arm==arm2].set_index('gene')
    common=sorted(set(a.index)&set(b.index))
    reversals=[g for g in common if a.at[g,'direction']!=b.at[g,'direction']]
    changed=lambda df,other:[{'gene':g,'direction':int(df.at[g,'direction']),'confirmed':bool(df.at[g,'confirmed'])} for g in sorted(set(df.index)-set(other.index))]
    comparisons[label]=dict(from_arm=arm1,to_arm=arm2,shared=len(common),direction_reversed=reversals,leaving=changed(a,b),entering=changed(b,a),yield_difference=int(b.confirmed.sum()-a.confirmed.sum()),bootstrap_central95=np.quantile(boot[arm2]-boot[arm1],[.025,.975],method='inverted_cdf').tolist())
result=dict(frozen_sha256=lock,external_family=len(test),arms=count,comparisons=comparisons)
(o/'comparison.json').write_text(json.dumps(result,indent=2))
for k,d in comparisons.items():
    print(k,d['shared'],[v['gene'] for v in d['leaving']],[v['gene'] for v in d['entering']],d['direction_reversed'],d['yield_difference'],d['bootstrap_central95'])
print(count)

