"""Prepare real genomic BAM/quantifications and evaluate paired airway controls.

Use `prepare` on Linux with Docker after successful KilnBio STAR/RSEM execution;
`fit` may run on another OS with the prepared outputs and KilnBio's R environment.
"""
import argparse,csv,json,math,re,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def resolve_strand(text):
    def fraction(label):
        match=re.search(re.escape(label)+r':\s*([0-9.]+)',text)
        if not match:raise ValueError('Unexpected paired RSeQC report: '+label)
        return float(match[1])
    if 'PairEnd Data' not in text:raise ValueError('Expected paired RNA library')
    forward=fraction('"1++,1--,2+-,2-+"');reverse=fraction('"1+-,1-+,2++,2--"')
    if forward+reverse<.2:raise ValueError('Insufficient strand-informative fraction')
    share=forward/(forward+reverse)
    if share>=.8:strand='forward'
    elif share<=.2:strand='reverse'
    elif abs(2*share-1)<=.1:strand='unstranded'
    else:raise ValueError(f'Ambiguous strand, forward share {share:.3f}')
    return strand,forward,reverse

def prepare(run,gtf,samples,out):
    out.mkdir(parents=True,exist_ok=False);shutil.copy2(samples,out/'samples.csv')
    metadata=list(csv.DictReader(samples.open()));ids=[r['sample'] for r in metadata];mix={r['sample']:r['geo_ercc_mix'] for r in metadata}
    if len(ids)!=8 or len(set(ids))!=8:raise ValueError('Expected eight airway libraries')
    quant=run/'results/star_rsem';endogenous=out/'rsem_endogenous';endogenous.mkdir()
    rsem_cmd={}
    for cmd in (run/'work').glob('*/*/.command.sh'):
        t=cmd.read_text(errors='replace')
        if not re.search(r'(?m)^\s*rsem-calculate-expression\s',t):continue
        for sample in ids:
            if re.search(r'\b'+sample+r'\b',t):
                m=re.search(r'--strandedness\s+(forward|reverse|none)',t)
                rsem_cmd[sample]=m[1] if m and m[1]!='none' else 'unstranded'
    audit=[]
    with (out/'libraries.csv').open('w',newline='') as lf,(out/'rsem.csv').open('w',newline='') as rf:
        lw=csv.writer(lf);rw=csv.writer(rf);lw.writerow(['sample','bam','strandedness']);rw.writerow(['sample','genes_results'])
        for s in ids:
            bam=quant/(s+'.markdup.sorted.bam');gene=quant/(s+'.genes.results')
            strand,fw,rv=resolve_strand((quant/'rseqc/infer_experiment'/(s+'.infer_experiment.txt')).read_text())
            if not bam.is_file() or strand!=rsem_cmd.get(s):raise ValueError(f'BAM absent or RSEM/RSeQC strand disagreement: {s}, {strand}, {rsem_cmd.get(s)}')
            lw.writerow([s,bam.resolve(),strand]);rw.writerow([s,'rsem_endogenous/'+gene.name])
            ercc=[]
            with gene.open() as src,(endogenous/gene.name).open('w',newline='') as dst:
                rows=csv.DictReader(src,delimiter='\t');writer=csv.DictWriter(dst,rows.fieldnames,delimiter='\t',lineterminator='\n');writer.writeheader()
                for row in rows:
                    if row['gene_id'].startswith('ERCC-'):ercc.append(float(row['expected_count']))
                    else:writer.writerow(row)
            if len(ercc)!=92 or (mix[s]=='1' and sum(ercc)<=0):raise ValueError('ERCC annotation/documented-positive quantification check failed')
            log=(quant/'log'/(s+'.Log.final.out')).read_text()
            q={a.strip():b.strip() for line in log.splitlines() if '|' in line for a,b in [line.split('|',1)]}
            audit.append(dict(sample=s,geo_ercc_mix=mix[s],strand=strand,forward=fw,reverse=rv,rsem_strand=rsem_cmd[s],star_input_pairs=int(q['Number of input reads']),
              star_unique_percent=float(q['Uniquely mapped reads %'].strip('%')),ercc_features=len(ercc),ercc_detected=sum(v>0 for v in ercc),ercc_expected_count=sum(ercc)))
    cfg=out/'featurecounts.json';cfg.write_text(json.dumps(dict(gtf=str(gtf.resolve()),libraries='libraries.csv',layout='paired',threads=4,runtime='docker'),indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'downstream/featurecounts.py'),str(cfg),str(out/'featurecounts')],check=True)
    sums={s:0 for s in ids};ec={s:0 for s in ids}
    with (out/'featurecounts/counts.tsv').open() as src,(out/'counts_endogenous.tsv').open('w',newline='') as dst:
        rows=csv.DictReader(src,delimiter='\t');w=csv.DictWriter(dst,rows.fieldnames,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:
            if row['feature_id'].startswith('ERCC-'):
                for s in ids:ec[s]+=int(row[s])
            else:
                w.writerow(row)
                for s in ids:sums[s]+=int(row[s])
    for r in audit:r.update(featurecounts_endogenous=sums[r['sample']],featurecounts_ercc=ec[r['sample']])
    with (out/'qc.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,audit[0].keys(),delimiter='\t');w.writeheader();w.writerows(audit)
    (out/'prepare.complete').write_text('Raw quantification/strand/count checks complete\n')

def fit(out,baseline):
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr
    if not (out/'prepare.complete').is_file():raise ValueError('Input preparation incomplete')
    samples=pd.read_csv(out/'samples.csv',dtype=str)
    cfg=dict(samples='samples.csv',assay='rnaseq',denominator='untrt',numerator='trt',group='condition',covariates=[dict(name='biological_id',type='factor')],min_count=10,min_samples=4,fdr=.05)
    checks=[];effects={};per_donor=[]
    for route,engine,source in [('RSEM_DESeq2','deseq2','rsem_tximport'),('featureCounts_DESeq2','deseq2','raw_counts'),('featureCounts_edgeR','edger','raw_counts')]:
        c=dict(cfg,engine=engine,count_source=source);c.update({'quant_manifest':'rsem.csv'} if source=='rsem_tximport' else {'counts':'counts_endogenous.tsv'})
        config=out/(route+'.json');config.write_text(json.dumps(c,indent=2)+'\n',encoding='utf-8')
        folder=out/route
        if not (folder/'provenance/status.json').exists():subprocess.run([sys.executable,str(ROOT/'downstream/analysis.py'),'de',str(config),str(folder)],check=True)
        if not json.loads((folder/'provenance/status.json').read_text())['complete']:raise ValueError('Incomplete fit: '+route)
        d=pd.read_csv(folder/'results/all_features.tsv',sep='\t').set_index('feature_id')
        lfc='log2FoldChange' if engine=='deseq2' else 'logFC';q='padj' if engine=='deseq2' else 'FDR';effects[route]=d[lfc]
        n=pd.read_csv(folder/'results'/('normalized_counts.tsv' if engine=='deseq2' else 'normalized_cpm.tsv'),sep='\t').set_index('feature_id');gene='ENSG00000103196'
        signs=[]
        for donor,group in samples.groupby('biological_id'):
            values={r.condition:float(n.loc[gene,r['sample']]) for _,r in group.iterrows()}
            signs.append(values['trt']>values['untrt']);per_donor.append(dict(route=route,donor=donor,**values))
        checks.append(dict(route=route,log2FC=float(d.loc[gene,lfc]),FDR=float(d.loc[gene,q]),donors_up=sum(signs),donors=len(signs),tested=len(d)))
    for engine,col in [('deseq2','log2FoldChange'),('edger','logFC')]:
        effects['released_counts_'+engine]=pd.read_csv(baseline/engine/'results/all_features.tsv',sep='\t').set_index('feature_id')[col]
    correlations=[]
    for i,(a,left) in enumerate(effects.items()):
        for b,right in list(effects.items())[i+1:]:
            pair=pd.concat([left,right],axis=1,join='inner').replace([np.inf,-np.inf],np.nan).dropna()
            correlations.append(dict(left=a,right=b,matched_features=len(pair),spearman=float(spearmanr(pair.iloc[:,0],pair.iloc[:,1]).statistic)))
    pd.DataFrame(per_donor).to_csv(out/'CRISPLD2_donors.tsv',sep='\t',index=False)
    record={'scope':'Five million complete pairs per donor-condition; eight public libraries, four donors. Orphan FASTQs excluded from matched branch comparison.','primary':'CRISPLD2 induction under dexamethasone','results':checks,'effect_concordance':correlations,'biological_control_pass':all(r['log2FC']>0 and r['FDR']<.05 and r['donors_up']==4 for r in checks),'empirical_FDR_calibrated':False}
    (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(json.dumps(checks))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='action',required=True)
    a=sp.add_parser('prepare');b=sp.add_parser('fit')
    for k in ['run','gtf','samples','output']:a.add_argument(k,type=Path)
    b.add_argument('output',type=Path);b.add_argument('released_count_benchmark',type=Path)
    v=p.parse_args()
    if v.action=='prepare':prepare(v.run,v.gtf,v.samples,v.output)
    else:fit(v.output,v.released_count_benchmark)
