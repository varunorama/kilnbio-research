from pathlib import Path
import csv, hashlib,json,re,sys,subprocess
import numpy as np,pandas as pd
root=Path(__file__).resolve().parent
data=root/'discovery'
out=root/'analysis'
out.mkdir(exist_ok=True)
repo=Path(sys.argv[1]).resolve()
def save(df,name,index=False): df.to_csv(out/name,sep='\t',index=index,na_rep='NA')
def audit(mode,inputs,options,name,expect_failure=False):
    path=out/(name+'.json')
    if path.exists():raise FileExistsError(path)
    cmd=[sys.executable,str(repo/'downstream/assay_audit.py'),mode,*[str(out/i) for i in inputs],*options,'--out',str(path)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    if expect_failure:
        assert p.returncode and 'conflicting donors or tissues' in p.stderr
        (out/(name+'.log')).write_text(p.stderr);return
    p.check_returncode()
    return json.loads(path.read_text())
clinical=pd.read_csv(data/'data_clinical_sample.txt',sep='\t',comment='#',dtype=str)
clinical=clinical[clinical.SAMPLE_ID.str.fullmatch(r'TCGA-[^-]+-[^-]+-01')].copy()
assert clinical.SAMPLE_ID.is_unique and clinical.PATIENT_ID.is_unique
donors=clinical.set_index('SAMPLE_ID').PATIENT_ID.to_dict()
profile=set(next(l.split(':',1)[1] for l in (data/'case_lists/cases_sequenced.txt').read_text().splitlines() if l.startswith('case_list_ids:')).split())
classes={'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Translation_Start_Site','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins'}
mutated=set()
egfr_classes={}
for chunk in pd.read_csv(data/'data_mutations.txt',sep='\t',usecols=['Hugo_Symbol','Variant_Classification','Tumor_Sample_Barcode'],dtype=str,chunksize=50000):
    x=chunk[chunk.Hugo_Symbol.eq('EGFR')]
    for v,n in x.Variant_Classification.value_counts().items():egfr_classes[v]=egfr_classes.get(v,0)+int(n)
    mutated.update(x.loc[x.Variant_Classification.isin(classes),'Tumor_Sample_Barcode'].str.slice(0,15))
rna=pd.read_csv(data/'data_mrna_seq_v2_rsem.txt',sep='\t',index_col='Hugo_Symbol').drop(columns='Entrez_Gene_Id')
rna=rna.loc[rna.index.notna() & ~rna.index.duplicated(keep=False)]
prot=pd.read_csv(data/'data_rppa.txt',sep='\t',index_col=0)
records=[]
for k in prot.index:
    gene,probe=k.split('|',1);gene=gene.strip()
    valid=(re.fullmatch(r'[A-Za-z0-9-]+',gene) and gene!='EGFR' and not re.search(r'(?:_|-)p[STY]\d|cleav|acetyl|phosph|methyl|ubiquit|mutant',probe,re.I))
    records.append((k,gene,probe,bool(valid)))
mapping=pd.DataFrame(records,columns=['source','gene','probe','total_single_candidate'])
keep=mapping[mapping.total_single_candidate]
keep=keep[~keep.gene.duplicated(keep=False)]
prot=prot.loc[keep.source];prot.index=keep.gene
primary=set(donors)
with_truth=primary & profile
all_rna=sorted(with_truth & set(rna.columns))
rna_only=sorted(set(all_rna)-set(prot.columns))
paired=sorted(set(all_rna)&set(prot.columns))
labels=pd.DataFrame({'specimen':all_rna,'donor':[donors[s] for s in all_rna],'egfr':[int(s in mutated) for s in all_rna],'base_ok':1,'add_ok':[int(s in paired) for s in all_rna]})
save(labels,'roster.tsv')
meta=labels.set_index('specimen')
obs=[]
for s in all_rna:
    for assay in (['RNA','protein'] if s in paired else ['RNA']):obs.append((donors[s],s,'primary',assay,assay+'_'+s))
obs=pd.DataFrame(obs,columns=['donor','specimen','tissue','assay','measurement'])
save(obs,'source_manifest.tsv')
opts=['--base','RNA','--add','protein','--tissue','primary']
r=audit('pair',['roster.tsv','source_manifest.tsv'],opts,'clean_pairing')
assert set(x['specimen'] for x in r['unique_pairs'])==set(paired) and not r['review_required']
audit('reach',['roster.tsv'],[],'reach')
mixed=obs.copy()
swaps=[min(s for s in paired if int(s in mutated)==i) for i in (0,1)]
for s,t in zip(swaps,reversed(swaps)):
    mixed.loc[(mixed.specimen==s)&(mixed.assay=='protein'),'donor']=donors[t]
save(mixed,'planted_manifest.tsv')
audit('pair',['roster.tsv','planted_manifest.tsv'],opts,'planted_pairing',True)
conflicts=mixed.groupby('specimen').donor.nunique()
flagged=sorted(conflicts[conflicts>1].index)
assert set(flagged)==set(swaps)
save(mixed[mixed.specimen.isin(flagged)],'conflicting_records.tsv')
save(mixed[~mixed.specimen.isin(flagged)],'withheld_manifest.tsv')
withheld=audit('pair',['roster.tsv','withheld_manifest.tsv'],opts,'withheld_pairing')
assert {x['specimen'] for x in withheld['unique_pairs']}==set(paired)-set(flagged)
# Counts on each native scale; comparisons use the same set of genes.
rall=np.log2(1+rna[all_rna]); pr=prot[paired]
genes=sorted(set(rall.index)&set(pr.index))
eligible=pd.Series(True,index=genes)
for mat in (rall,rall[paired],pr):
    x=mat.loc[genes]
    eligible &= x.std(axis=1)>0
    for g in [0,1]:
        cols=[s for s in x.columns if int(s in mutated)==g]
        eligible &= x[cols].notna().mean(axis=1)>=0.8
genes=eligible[eligible].index
rall=rall.loc[genes].rename_axis('gene')
pr=pr.loc[genes].rename_axis('gene')
save(rall,'RNA.tsv',True);save(pr,'protein.tsv',True)
save(mapping,'antibodies.tsv')
summary=dict(public_primary=len(primary),mutation_profiled_primary=len(with_truth),missing_RNA_among_profiled=len(with_truth-set(all_rna)),all_rna=len(all_rna),matched=len(paired),RNA_only=len(rna_only),mutation_class_counts=egfr_classes,groups=labels.groupby(['egfr','add_ok']).size().unstack(fill_value=0).to_dict(),common_genes=len(genes),conflicts_clean=0,planted=swaps,identified_from_conflicting_keys=flagged,software_failure_message=(out/'planted_pairing.log').read_text().strip(),same_sample_ID_proves_same_vial=False)
(out/'design_counts.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary))

