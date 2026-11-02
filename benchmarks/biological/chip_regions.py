"""Freeze hg19 regions; combine two matched PC9 spike-in campaigns after QC.

The `freeze` step is independent of treatment signal. `assemble` requires
completed per-campaign alignments. DE tests use biological cultures, campaign
blocking and raw integer counts. Promoter summaries are exploratory, not the
original qPCR amplicons. No empirical FDR truth set is available.
"""
import argparse,collections,csv,gzip,hashlib,json,re,subprocess,sys,urllib.request
from array import array
from pathlib import Path
from chip_spikein import IMAGES,R1,SAMPLES
ROOT=Path(__file__).resolve().parents[2]
URLS={'hg19.chrom.sizes':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
      'gap.txt.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz',
      'blacklist.bed.gz':'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg19-blacklist.v2.bed.gz'}

def write_table(path,rows,delimiter='\t'):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,rows[0].keys(),delimiter=delimiter,lineterminator='\n');w.writeheader();w.writerows(rows)

def freeze(gtf,out):
    out.mkdir(parents=True,exist_ok=False);hashes={}
    for name,url in URLS.items():
        data=urllib.request.urlopen(url,timeout=90).read();(out/name).write_bytes(data);hashes[name]=hashlib.sha256(data).hexdigest()
    chroms=['chr'+str(i) for i in range(1,23)]
    sizes={c:int(n) for c,n in (l.split() for l in (out/'hg19.chrom.sizes').read_text().splitlines()) if c in chroms}
    if len(sizes)!=22:raise ValueError('Incomplete hg19 autosome list')
    masked=set()
    for name,col in [('gap.txt.gz',1),('blacklist.bed.gz',0)]:
        with gzip.open(out/name,'rt') as f:
            for line in f:
                c,a,b=line.split()[col:col+3]
                if c in sizes:
                    a,b=int(a),int(b)
                    if not 0<=a<b<=sizes[c]:raise ValueError('Invalid mask interval')
                    masked.update((c,i) for i in range(a//10000,(b-1)//10000+1))
    n=0
    with (out/'regions.tsv').open('w') as f:
        f.write('feature_id\tchrom\tstart\tend\n')
        for c in chroms:
            for a in range(0,sizes[c],10000):
                if (c,a//10000) in masked:continue
                b=min(a+10000,sizes[c]);f.write(f'{c}:{a}-{b}\t{c}\t{a}\t{b}\n');n+=1
    genes={'MYT1':{},'CCND2':{}}
    with gtf.open() as f:
        for line in f:
            if 'gene_name "MYT1"' not in line and 'gene_name "CCND2"' not in line:continue
            v=line.rstrip().split('\t')
            if len(v)!=9 or v[2]!='exon' or 'chr'+v[0] not in sizes:continue
            sym=re.search(r'gene_name "([^"]+)"',v[8])[1];tx=re.search(r'transcript_id "([^"]+)"',v[8])[1]
            p=int(v[3])-1 if v[6]=='+' else int(v[4])-1
            old=genes[sym].get(tx,('chr'+v[0],p))[1]
            genes[sym][tx]=('chr'+v[0],min(old,p) if v[6]=='+' else max(old,p))
    for gene,tx in genes.items():
        if not tx:raise ValueError('Promoter annotation absent: '+gene)
        regions=[]
        for c,a,b in sorted((c,max(0,p-2000),min(sizes[c],p+2001)) for c,p in tx.values()):
            if regions and regions[-1][0]==c and a<=regions[-1][2]:regions[-1][2]=max(b,regions[-1][2])
            else:regions.append([c,a,b])
        (out/(gene+'.bed')).write_text(''.join(f'h_{c}\t{a}\t{b}\n' for c,a,b in regions))
    rec={'urls':URLS,'source_sha256':hashes,'regions_sha256':hashlib.sha256((out/'regions.tsv').read_bytes()).hexdigest(),'regions':n,
         'coordinate_system':'hg19 0-based half-open','excluded_tiles':len(masked),'promoters':'Union of +/- 2 kb around all annotated Ensembl 75 transcript TSS; distinct from exact qPCR amplicons',
         'blacklist_citation':'https://doi.org/10.1038/s41598-019-45839-z'}
    (out/'freeze.json').write_text(json.dumps(rec,indent=2)+'\n');print(n,'frozen regions')

def assemble(ref,r1,r2,out):
    for folder in [r1,r2]:
        if not (folder/'benchmark.complete').is_file():raise ValueError('Incomplete campaign: '+str(folder))
    out.mkdir(parents=True,exist_ok=False)
    regions=list(csv.DictReader((ref/'regions.tsv').open(),delimiter='\t'));index={(r['chrom'],int(r['start'])):i for i,r in enumerate(regions)}
    meta=[(sample,mark,cond,camp,folder) for camp,folder,ss in [('R1',r1,R1),('R2',r2,SAMPLES)] for sample,(mark,cond) in ss.items()]
    samples={s:json.loads((folder/s/'counts.json').read_text()) for s,mark,cond,camp,folder in meta}
    ratios=[];promoters=[]
    for mark in ['H3K27me3','H3K4me3']:
        group=[r for r in meta if r[1]==mark];dest=out/mark;dest.mkdir()
        design=[dict(sample=s,biological_id=s,condition=cond,campaign=camp) for s,m,cond,camp,folder in group]
        write_table(dest/'samples.csv',design,',')
        rows=[dict(sample=s,normalization_group=m,target_reads=samples[s]['primary']['human'],reference_reads=samples[s]['primary']['fly'],target_amount=30000,reference_amount=750,amount_unit='ng_chromatin',reference_material='OSS_H2Av') for s,m,c,camp,f in group]
        write_table(dest/'spike_counts.csv',rows,',')
        subprocess.run([sys.executable,str(ROOT/'downstream/spikein.py'),str(dest/'spike_counts.csv'),str(dest/'factors'),'--group',mark],check=True)
        columns=[]
        for s,m,c,camp,folder in group:
            col=array('I',[0])*len(regions)
            with (folder/s/'human_10kb_leftmost_starts.tsv').open() as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    i=index.get((r['chrom'][2:],int(r['start'])))
                    if i is not None:col[i]=int(r['reads'])
            columns.append(col)
            for gene in ['MYT1','CCND2']:
                cmd=['docker','run','--rm','--network','none','-v',f'{folder}:/bam:ro','-v',f'{ref}:/regions:ro',IMAGES['samtools'],'samtools','view','-c','-M','-F','3844','-q','30','-e','[NM]<=2','-L',f'/regions/{gene}.bed',f'/bam/{s}/marked.bam']
                reads=int(subprocess.check_output(cmd,text=True).strip())
                promoters.append(dict(sample=s,mark=m,condition=c,campaign=camp,gene=gene,reads=reads,per_million_fly=1e6*reads/samples[s]['primary']['fly'],per_million_human=1e6*reads/samples[s]['primary']['human']))
        with (dest/'counts.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['feature_id',*[s[0] for s in group]])
            for i,r in enumerate(regions):w.writerow([r['feature_id'],*[col[i] for col in columns]])
        for method in ['TMM','spikein']:
            cfg=dict(samples='samples.csv',counts='counts.tsv',assay='chipseq',count_source='raw_counts',engine='edger',group='condition',denominator='DMSO',numerator='EZH2i',covariates=[dict(name='campaign',type='factor')],feature_universe='external_fixed',normalization_reviewed=True,min_count=10,min_samples=2,fdr=.05)
            if method=='spikein':cfg['external_size_factors']='factors/size_factors.csv'
            (dest/(method+'.json')).write_text(json.dumps(cfg,indent=2)+'\n')
        for campaign in ['R1','R2']:
            g={c:s for s,m,c,camp,folder in group if camp==campaign}
            for method in ['primary','duplicates_retained','mapq26']:
                a=samples[g['DMSO']][method];b=samples[g['EZH2i']][method]
                ratio=(b['human']/b['fly'])/(a['human']/a['fly'])
                ratios.append(dict(campaign=campaign,mark=mark,method=method,dmso_human=a['human'],dmso_fly=a['fly'],inhibitor_human=b['human'],inhibitor_fly=b['fly'],ratio=ratio))
    write_table(out/'global_ratios.tsv',ratios);write_table(out/'promoter_counts.tsv',promoters)
    (out/'assembly.complete').write_text('Two campaigns assembled; reference-normalized and ordinary edgeR fits can now run.\n')

def fit(out):
    """Compare the same regions, samples and design with two denominators."""
    import pandas as pd
    if not (out/'assembly.complete').is_file(): raise ValueError('Campaign assembly is incomplete')
    details=[]
    for mark in ['H3K27me3','H3K4me3']:
        folder=out/mark; tested=None
        for method in ['TMM','spikein']:
            dest=folder/('de_'+method); status=dest/'provenance/status.json'
            if not status.is_file():
                subprocess.run([sys.executable,str(ROOT/'downstream/analysis.py'),'de',
                                str(folder/(method+'.json')),str(dest)],check=True)
            if not json.loads(status.read_text())['complete']: raise ValueError('Incomplete fit: '+str(dest))
            tab=pd.read_csv(dest/'results/all_features.tsv',sep='\t')
            ids=frozenset(tab.feature_id)
            if len(ids)!=len(tab) or (tested is not None and tested!=ids):
                raise ValueError('Features differ across fits or are duplicated')
            tested=ids
            details.append(dict(mark=mark,method=method,tested=len(tab),nominal_BH05=int((tab.FDR<.05).sum()),
                                median_logFC=float(tab.logFC.median()),q25=float(tab.logFC.quantile(.25)),q75=float(tab.logFC.quantile(.75))))
    ratios=pd.read_csv(out/'global_ratios.tsv',sep='\t')
    primary=ratios[ratios.method=='primary'].pivot(index='campaign',columns='mark',values='ratio')
    if set(primary.index)!={'R1','R2'} or set(primary.columns)!={'H3K27me3','H3K4me3'} or primary.isna().any().any():
        raise ValueError('Both verified campaigns and both marks are required')
    passed=(primary.H3K27me3<.5)&primary.H3K4me3.between(.5,2)&((primary.H3K27me3/primary.H3K4me3)<.5)
    result={'scope':'Raw-read reanalysis of the two published PC9 campaigns, five million randomly sampled reads/library; not an independent cohort',
            'global':ratios.to_dict('records'),'regional':details,'biological_control_pass':bool(passed.all()),
            'residual_df_per_mark':1,'empirical_FDR_calibrated':False}
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='action',required=True)
    a=sp.add_parser('freeze');a.add_argument('gtf',type=Path);a.add_argument('out',type=Path)
    b=sp.add_parser('assemble')
    for k in ['regions','campaign1','campaign2','out']:b.add_argument(k,type=Path)
    c=sp.add_parser('fit');c.add_argument('out',type=Path)
    v=p.parse_args()
    if v.action=='freeze':freeze(v.gtf,v.out)
    elif v.action=='assemble':assemble(v.regions.resolve(),v.campaign1.resolve(),v.campaign2.resolve(),v.out.resolve())
    else:fit(v.out.resolve())
