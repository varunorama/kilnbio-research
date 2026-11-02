"""Audit completed public biological examples. Writes only to the study directory.

python summarize.py <biological-validation directory>
"""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def de_record(folder, ids):
    records=[]
    for engine in ('deseq2','edger'):
        run=folder/engine
        assert json.loads((run/'provenance/status.json').read_text())['complete']
        d=pd.read_csv(run/'results/all_features.tsv',sep='\t').set_index('feature_id')
        normalized=pd.read_csv(run/'results'/('normalized_counts.tsv' if engine=='deseq2' else 'normalized_cpm.tsv'),sep='\t').set_index('feature_id')
        samples=pd.read_csv(folder/'inputs/samples.csv' if (folder/'inputs/samples.csv').exists() else folder/'samples.csv',dtype=str)
        cfg=json.loads((run/'provenance/config.json').read_text())
        for gene,symbol in ids.items():
            row=d.loc[gene];lfc=row['log2FoldChange' if engine=='deseq2' else 'logFC'];q=row['padj' if engine=='deseq2' else 'FDR']
            signs=[]
            for donor,sub in samples.groupby('biological_id'):
                lo=sub[sub.condition==cfg['denominator']].iloc[0]['sample'];hi=sub[sub.condition==cfg['numerator']].iloc[0]['sample']
                signs.append(bool(normalized.loc[gene,hi]>normalized.loc[gene,lo]))
            records.append(dict(engine=engine,gene=gene,symbol=symbol,log2_fold_change=float(lfc),FDR=float(q),
                                donors_up=sum(signs),donors=len(signs),tested_features=len(d)))
    return records

def main(base):
    result={'scope':'Public positive-control recovery and real-input execution; not independent FDR or clinical calibration'}
    result['airway']=de_record(base/'airway',{'ENSG00000103196':'CRISPLD2'})
    mono=base/'kang/pseudobulk/celltype_002'
    result['kang']=de_record(mono,{'ENSG00000185745':'IFIT1','ENSG00000187608':'ISG15','ENSG00000157601':'MX1'})
    qc=json.loads((base/'kang/qc/provenance/run.json').read_text())
    result['kang_qc']={k:qc[k] for k in ('cells_input','cells_retained','genes_retained')}
    chrom=base/'chromatin'
    atac=pd.read_csv(chrom/'ATAC/de/results/all_features.tsv',sep='\t').set_index('feature_id')
    dhs=pd.read_csv(chrom/'dnase_support.tsv',sep='\t').set_index('feature_id')
    joined=atac.join(dhs);records=[]
    for cell,other,sign in [('K562','GM12878',1),('GM12878','K562',-1)]:
        all_ref=(dhs[cell+'_r1']==1)&(dhs[cell+'_r2']==1)&(dhs[other+'_r1']==0)&(dhs[other+'_r2']==0)
        selected=(joined[cell+'_r1']==1)&(joined[cell+'_r2']==1)&(joined[other+'_r1']==0)&(joined[other+'_r2']==0)
        values=joined.loc[selected]
        records.append(dict(cell=cell,reference_tiles=int(all_ref.sum()),tested_tiles=len(values),
                            direction_agreement=float((sign*values.logFC>0).mean()),
                            direction_and_FDR05=int(((sign*values.logFC>0)&(values.FDR<.05)).sum()),
                            median_logFC=float(values.logFC.median())))
    result['atac_dnase']=records
    result['chromatin_tests']={};result['cross_assay_exploratory']=[]
    for assay in ('ATAC','H3K27me3','H3K36me3'):
        r=chrom/assay/'de';assert json.loads((r/'provenance/status.json').read_text())['complete']
        d=pd.read_csv(r/'results/all_features.tsv',sep='\t').set_index('feature_id')
        result['chromatin_tests'][assay]={'tested':len(d),'FDR05':int((d.FDR<.05).sum())}
        if assay!='ATAC':
            z=atac[['logFC']].rename(columns={'logFC':'ATAC'}).join(d[['logFC']].rename(columns={'logFC':assay}),how='inner')
            result['cross_assay_exploratory'].append({'assay':assay,'tiles':len(z),'spearman_with_ATAC':float(spearmanr(z.ATAC,z[assay]).statistic)})
    control=pd.read_csv(chrom/'input_controls/counts.tsv',sep='\t').set_index('feature_id')
    result['chip_input_adjusted_exploratory']=[]
    for row in json.loads((chrom/'input_controls/lineage.json').read_text()):
        counts=pd.read_csv(chrom/row['assay']/'counts.tsv',sep='\t').set_index('feature_id')
        ip=counts[[row['rep1'],row['rep2']]].sum(axis=1);bg=control[row['controls'].split(',')].sum(axis=1)
        enrichment=np.log2((ip+.5)/(ip.sum()+.5*len(ip)))-np.log2((bg+.5)/(bg.sum()+.5*len(bg)))
        active=(dhs[row['cell']+'_r1']==1)&(dhs[row['cell']+'_r2']==1)
        absent=(dhs[row['cell']+'_r1']==0)&(dhs[row['cell']+'_r2']==0)
        result['chip_input_adjusted_exploratory'].append({'assay':row['assay'],'cell':row['cell'],
            'DHS_minus_absent_median_log2_enrichment':float(enrichment[active].median()-enrichment[absent].median()),
            'scope':'Pooled replicate descriptive 10kb bins; input-adjusted enrichment is not spike-in normalization or a differential p-value'})
    colon=base/'colon';d=pd.read_csv(colon/'dss/results/all_loci.tsv',sep='\t')
    assert json.loads((colon/'dss/provenance/status.json').read_text())['complete']
    inp=[pd.read_csv(colon/f'inputs/{s}.tsv',sep='\t').set_index(['chr','pos']) for s in ['C1','C2','C3','N1','N2','N3']]
    cov=pd.concat([t.N.rename(i) for i,t in enumerate(inp)],axis=1);keep=(cov>=5).all(axis=1)
    differences=[]
    for i in range(3):
        delta=(inp[i].X/inp[i].N-inp[i+3].X/inp[i+3].N)[keep]
        differences.append({'pair':i+1,'median_cancer_minus_normal':float(delta.median()),'mean_cancer_minus_normal':float(delta.mean()),'CpGs':len(delta)})
    result['colon']={'tested_loci':len(d),'BH05_loci':int((d.FDR_BH<.05).sum()),
        'median_unadjusted_delta':float(d.unadjusted_delta_fraction.median()),'paired_common_coverage':differences,
        'candidate_regions':len(pd.read_csv(colon/'dss/results/DMR_candidates.tsv',sep='\t')),
        'region_FDR_validated':False,'coordinate_audit':json.loads((colon/'reference_audit.json').read_text())}
    result['ont']={channel:json.loads((base/'ont'/channel/'status.json').read_text())['sites_written'] for channel in ('5mC','5hmC')}
    result['ont_scope']='Ten real ONT control reads; not human tissue truth or a current Dorado benchmark. PacBio HG002 fixture excluded.'
    (base/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
