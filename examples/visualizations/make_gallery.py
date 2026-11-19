"""Render the public gallery from fictional, deterministic output fixtures.

These are renderer examples, not fitted biological findings or accuracy benchmarks.
Run: python examples/visualizations/make_gallery.py NEW_OUTPUT [--publish docs/gallery]
"""
import argparse
import math
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'downstream'))
from auto_reports import render
from plot_design import panel_guide


def fixtures(out):
    rng=np.random.default_rng(18026);samples=[f'S{i+1}' for i in range(8)];genes=[f'SYN{i+1:03d}' for i in range(80)]
    counts=rng.negative_binomial(15,.2,(80,8));counts[:10,4:]*=3
    mat=pd.DataFrame(counts,columns=samples);mat.insert(0,'feature',genes)
    meta=pd.DataFrame({'sample':samples,'group':['control']*4+['treated']*4,'batch':['A','B']*4})
    effects=rng.normal(0,.35,80);effects[:12]+=np.tile([-1.8,1.8],6)
    se=rng.uniform(.2,.6,80);pvals=np.array([math.erfc(abs(x)/s/math.sqrt(2)) for x,s in zip(effects,se)])
    order=np.argsort(pvals);qvals=np.empty(80);qvals[order]=np.minimum(1,np.minimum.accumulate((pvals[order]*80/np.arange(1,81))[::-1])[::-1])
    stats=pd.DataFrame({'feature':genes,'effect':effects,'p':pvals,'q':qvals,'baseMean':counts.mean(axis=1)})
    stats.loc[::13,['effect','p','q']]=np.nan
    def folder(n):p=out/n;p.mkdir();return p
    def tsv(t,p):t.to_csv(p,sep='\t',index=False)
    def js(o,p):p.write_text(json.dumps(o,indent=2)+'\n')
    inputs=folder('inputs');tsv(mat,inputs/'counts.tsv');tsv(meta,inputs/'samples.tsv')
    jobs=[]
    def job(c,p,a=None):jobs.append((c,[str(x) for x in (a or [])],p))
    p=folder('differential');r=p/'results';r.mkdir();tsv(stats,r/'all_features.tsv');x=mat.astype({c:float for c in samples});x.iloc[:,1:]=np.log2(counts+1);tsv(x,r/'display_transformed.tsv');tsv(meta,r/'display_samples.tsv');tsv(pd.DataFrame({'sample':samples,'library_size':counts.sum(axis=0),'normalization_factor':rng.uniform(.9,1.1,8)}),r/'sample_diagnostics.tsv');job('de',p)
    p=folder('methylation');r=p/'results';r.mkdir();s=stats.drop(columns='baseMean').copy();s['chr']=['chr1']*40+['chr2']*40;s['pos']=np.arange(80)*10+1000;s=s.rename(columns={'effect':'unadjusted_delta_fraction'});s['unadjusted_delta_fraction']/=10;tsv(s,r/'all_loci.tsv');x=mat.astype({c:float for c in samples});x['feature']=['chr1:'+str(i*10+1000) for i in range(80)];x.iloc[:,1:]=rng.uniform(0,1,(80,8));x.iloc[3,2]=np.nan;tsv(x,r/'display_fractions.tsv');depth=mat.copy();depth['feature']=x.feature;tsv(depth,r/'display_depth.tsv');tsv(meta,r/'display_samples.tsv');job('dm',p)
    p=folder('method_comparison');t=pd.DataFrame({'feature':genes})
    for n in ['DESeq2','edgeR','limma_voom']:
        t[n+'_effect']=stats.effect+rng.normal(0,.15,80);t[n+'_q']=np.clip(stats.q+rng.uniform(-.06,.06,80),0,1)
    tsv(t,p/'features.tsv');tsv(pd.DataFrame({'feature':genes,'eligibility':['included']*70+['filtered']*10}),p/'eligibility.tsv');job('differential-compare',p,['rna-counts',inputs/'counts.tsv',inputs/'samples.tsv'])
    p=folder('olink');x=mat.astype({c:float for c in samples});x.iloc[:,1:]=rng.normal(7,1.4,(80,8));x.iloc[0:8,3]=np.nan;tsv(x,p/'npx.tsv');tsv(meta.rename(columns={'sample':'SampleID'}),p/'samples.tsv');tsv(pd.DataFrame({'SampleQC':['PASS']*7+['WARN'],'AssayQC':['PASS','PASS','WARN','PASS']*2,'above_LOD':[True]*6+[False]*2}),p/'measurement_qc.tsv');tsv(pd.DataFrame({'cramers_v':[0,.3,.8]}),p/'aliasing.tsv');js({'samples':8,'assays':80,'scale':'NPX log2','data':'synthetic renderer fixture'},p/'qc.json');job('olink-import',p)
    p=folder('modifications');tsv(pd.DataFrame({'chr':['chr1']*80,'pos':np.arange(80)*10+1000,'N':[20]*80,'X':rng.binomial(20,.6,80)}),p/'counts.tsv');js({'complete':True,'chemistry':'synthetic modC counts','coordinates':'1-based-cytosine','rows_written':80},p/'status.json');job('import-mods',p)
    p=folder('capture');(p/'targets.bed').write_text('chr1\t100\t400\nchr1\t1000\t1400\nchr2\t500\t1200\n');tsv(pd.DataFrame({'feature':genes,'fully_on_target':[True]*60+[False]*20}),p/'features.tsv');js({'assembly':'synthetic','target_rows':3,'target_union_bases':1400,'alignment_depth':{'mean_depth':31.2,'fraction_at_min_depth':.91,'min_depth':20}},p/'summary.json');job('capture-panel',p)
    p=folder('genomic');tsv(pd.DataFrame([{'event':g,'source_sample':s,'state':['altered','reference','unknown'][(i+j)%3]} for i,g in enumerate(['GENE_A','GENE_B','GENE_C']) for j,s in enumerate(samples)]),p/'calls.tsv');job('genomic-panel',p/'calls.tsv')
    p=folder('CH');tsv(pd.DataFrame({'sample':samples,'ch_status':['candidate_only']*3+['unknown']*5}),p/'samples.tsv');tsv(pd.DataFrame({'gene':['GENE_A','GENE_B']*4,'vaf':rng.uniform(.02,.25,8),'dp':rng.integers(50,400,8),'alt_count':rng.integers(10,30,8),'candidate':[True]*3+[False]*5,'reasons':['review']*8}),p/'variants.tsv');job('ch-screen',p)
    p=folder('annotations');tsv(pd.DataFrame({n+'_status':['matched']*5+['not_found']*2+['multiple'] for n in ['clinvar','alphamissense','vep']}),p/'variants.tsv');job('variant-annotations',p)
    p=folder('structures');tsv(pd.DataFrame({'uniprot_id':['SYN_P1','SYN_P2','SYN_P3'],'protein_variant':['p.A12V','p.G51S','p.R86H'],'mapping_status':['mapped','mapped','mapped'],'tool_used':['database_lookup']*3,'mean_plddt':[90,72,81],'sequence_start':[1]*3,'sequence_end':[100,150,120],'model_sequence_position':[12,51,86]}),p/'structures.tsv');job('alphafold',p)
    p=folder('alphagenome');tsv(pd.DataFrame({'variant_id':['SYN1','SYN2']*4,'output_type':['RNA_SEQ']*4+['ATAC']*4,'raw_score':rng.normal(0,1,8),'quantile_score':rng.uniform(0,1,8)}),p/'scores.tsv');js({'status':'complete','source':'synthetic predictions for renderer testing'},p/'receipt.json');job('alphagenome',p)
    p=folder('clocks');tsv(pd.DataFrame({'sample':samples,'synthetic_age_years':rng.normal(50,10,8),'synthetic_pace':rng.normal(1,.1,8)}),p/'scores.tsv');tsv(pd.DataFrame({'model':['synthetic_age','synthetic_pace'],'required':[80,60],'present':[80,55],'complete_profiles':[8,5],'scored_profiles':[8,5]}),p/'coverage.tsv');job('clock-panel',p)
    p=folder('prediction');pred=pd.DataFrame({'donor':samples,'observed':np.arange(8)+20,'RNA':np.arange(8)+20+rng.normal(0,2,8),'joint':np.arange(8)+20+rng.normal(0,1,8)});tsv(pred,p/'predictions.tsv');tsv(pd.DataFrame({'model':['joint'],'improvement_MAE':[.5],'lower':[-.1],'upper':[1.1]}),p/'paired_differences.tsv');job('compare-assays',p)
    p=folder('spikein');tsv(pd.DataFrame({'sample':samples,'target_reads':counts.sum(axis=0)*100,'reference_reads':rng.integers(1000,3000,8),'size_factor':rng.uniform(.5,2,8)}),p/'diagnostics.tsv');job('spikein-factors',p)
    p=folder('audit');js({'donors':8,'patterns':{'both':5,'base_only':2,'add_only':1,'neither':0},'base_reach':7,'combined_reach':5,'both_scored':{'n':5,'base_mae':2,'policy_mae':1.5,'gain':.5,'interval':[-.2,1.2]}},p/'reach.json');job('assay-audit',p/'reach.json')
    js({'nominated_genes':3,'alpha':.05,'candidates':[{'gene':g,'endpoints':[{'assay':a,'cohort':'synthetic validation','status':['confirmed','not_confirmed','unmeasured'][(i+j)%3]} for j,a in enumerate(['RNA','protein','ATAC'])]} for i,g in enumerate(['SYN001','SYN002','SYN003'])]},p/'candidates.json');job('candidate-evidence',p/'candidates.json')
    p=folder('graft');tsv(pd.DataFrame({'lot':['A','B','C'],'predictor':[.2,.5,.8],'median_outcome':[20,30,35],'enrolled_recipients':[4,5,4],'measured_recipients':[3,5,4]}),p/'lots.tsv');job('graft-readouts',p)
    p=folder('vep');(p/'annotated.vcf').write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|IMPACT|BIOTYPE">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\tPASS\tCSQ=G|missense_variant|MODERATE|protein_coding\nchr2\t200\t.\tC\tT\t.\tPASS\tCSQ=T|intron_variant|MODIFIER|protein_coding\n');job('vep-annotate',p)
    p=folder('decisions');js({'method':'rank_assays','result':{'features':['SYN1','SYN2','SYN3'],'arms':{'RNA':{'ranking':['SYN1','SYN2','SYN3']},'protein':{'ranking':['SYN3','SYN1','SYN2']},'combined':{'ranking':['SYN1','SYN3','SYN2']}},'transitions':[{'from':'RNA','to':'combined','rank_change':{'SYN1':0,'SYN2':1,'SYN3':-1}}]}},p/'ranks.json');job('decision-methods',p/'ranks.json')
    p=folder('bladder_style');tsv(pd.DataFrame([{'arm':arm,'gene':g,'rank':i+1} for arm,gs in [('A',genes[:20]),('B',genes[2:20]+genes[:2]),('C',genes[1:20]+genes[:1])] for i,g in enumerate(gs)]),p/'rankings.tsv');tsv(pd.DataFrame({'gene':genes[:3],'status':['retained','entered','exited']}),p/'changes.tsv');js({'budget':15,'common_genes':20,'data':'fictional rankings, not the bladder cohort'},p/'report.json');job('bladder-example',p)
    # An existing embedding is provided only to exercise the renderer, not to claim a fitted cluster result.
    import anndata
    p=folder('singlecell');n=160;obs=pd.DataFrame({'sample':np.repeat(samples,20),'leiden':pd.Categorical(np.tile(['0','1','2','3'],40)),'total_counts':rng.integers(500,5000,n),'n_genes_by_counts':rng.integers(100,1000,n),'pct_counts_mt':rng.uniform(0,20,n),'qc_pass':[True]*n},index=[f'C{i}' for i in range(n)])
    ad=anndata.AnnData(X=np.log1p(rng.poisson(2,(n,30))).astype('float32'),obs=obs,var=pd.DataFrame({'highly_variable':[True]*30},index=genes[:30]));ad.obsm['X_umap']=rng.normal(size=(n,2))+np.asarray(obs.leiden.astype(int))[:,None]*2;ad.write_h5ad(p/'filtered.h5ad');tsv(obs.reset_index(names='cell_id'),p/'cell_qc.tsv');job('scrna',p)
    return jobs


def publish(reports,dest):
    from PIL import Image,ImageOps,ImageDraw,ImageFont
    from matplotlib import font_manager
    dest.mkdir(parents=True,exist_ok=False);items=[]
    for cmd,report in reports:
        receipt=json.loads((report/'visualization.json').read_text())
        for panel in receipt['panels']:
            name=cmd+'__'+panel['name'];src=report/'figures'/(panel['name']+'.png');dst=dest/(name+'.png');shutil.copyfile(src,dst)
            data=report/'plotted_data'/(panel['name']+'.tsv');shutil.copyfile(data,dest/(name+'.tsv'))
            items.append(dict(command=cmd,name=panel['name'],image=dst.name,data=name+'.tsv',caption=panel['caption'],**panel_guide(panel['name'])))
    (dest/'manifest.json').write_text(json.dumps({'data':'synthetic renderer fixtures; not biological validation','generator':'examples/visualizations/make_gallery.py','panels':items},indent=2)+'\n')
    text=['# KilnBio visualization gallery','','Every example below uses **fictional data**. These are the same renderers called automatically by the listed commands. They demonstrate output design, not biological validation. Each image links to its full resolution; each TSV contains the plotted values.','','Run `python examples/visualizations/make_gallery.py NEW_OUTPUT --with-enrichment` in the singlecell environment to regenerate these examples and their offline HTML/PDF reports. For real data, run the analysis command normally. See [automatic reports](../VISUALIZATIONS.md) for required inputs, command aliases, output locations, display limits, and native upstream reports.','','The matrix-QC designs shown for `de` and Olink are also used by `multiomic-layers`, `featurecounts` and donor-level `pseudobulk`. Method-comparison panels also serve `method-agreement`; modification panels also serve completed `ont` outputs; score panels serve compatible methylation and other assay clocks. A new fitted model is never inferred from a similar-looking chart.','']
    commands=list(dict.fromkeys(x['command'] for x in items))
    text+=['## Browse by analysis','',' · '.join(f'[`{cmd}`](#{cmd})' for cmd in commands),'']
    for cmd in commands:
        text+=['## '+cmd,'']
        for x in [v for v in items if v['command']==cmd]:
            text+=['### '+x['title'],'','**Why it matters.** '+x['why'],'',f"[![{x['title']}]({x['image']})]({x['image']})",'','**How to read it.** '+x['read'],'','<details><summary>Method and display details</summary>','',x['caption'],'','</details>','',f"[Full-size figure]({x['image']}) · [Plotted values (TSV)]({x['data']})",'']
    (dest/'README.md').write_text('\n'.join(text)+'\n')
    selections=[('assays',['de','dm','olink-import','scrna']),('integration',['differential-compare','compare-assays','assay-audit','candidate-evidence']),('platforms',['capture-panel','genomic-panel','alphafold','clock-panel'])]
    for title,commands in selections:
        chosen=[]
        for cmd in commands:
            choices=[v for v in items if v['command']==cmd and v['name'] not in ['audit_metrics']]
            preferred={'de':'all_features_volcano','dm':'all_loci_manhattan','olink-import':'NPX_feature_heatmap','scrna':'cell_expression_dotplot','differential-compare':'method_calls','capture-panel':'panel_intervals'}
            chosen.append(next((x for x in choices if x['name']==preferred.get(cmd)),choices[-1]))
        sheet=Image.new('RGB',(2000,1550),'#EEF2F6');draw=ImageDraw.Draw(sheet)
        regular=font_manager.findfont('DejaVu Sans');bold=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
        fonts={s:ImageFont.truetype(regular,s) for s in [21,24]};heading=ImageFont.truetype(bold,44);subhead=ImageFont.truetype(bold,26)
        titles={'assays':'Understand each assay','integration':'Connect evidence. Inspect decisions.','platforms':'Add molecular context'}
        draw.text((52,30),'KILNBIO  /  VISUAL ANALYSIS',font=fonts[21],fill='#286FA3')
        draw.text((52,69),titles[title],font=heading,fill='#203447')
        draw.text((52,130),'Automatic reports  ·  Inspectable values  ·  Fictional demonstration data',font=fonts[24],fill='#637587')
        for i,item in enumerate(chosen):
            x=42+(i%2)*980;y=195+(i//2)*655
            draw.rounded_rectangle((x,y,x+936,y+625),radius=20,fill='white',outline='#D9E2EB',width=2)
            draw.text((x+24,y+20),f'{i+1:02d}  '+item['title'],font=subhead,fill='#203447')
            draw.text((x+24,y+62),item['command']+'  /  SYNTHETIC EXAMPLE',font=fonts[21],fill='#637587')
            im=Image.open(dest/item['image']).convert('RGB');im=ImageOps.contain(im,(900,500));sheet.paste(im,(x+(936-im.width)//2,y+112+(500-im.height)//2))
        sheet.save(dest/(title+'-overview.png'))


def integration_examples(out):
    import importlib.util
    import multiomic_meta,pathway_compare
    p=out/'joint_inputs';p.mkdir();layers=[]
    for i,n in enumerate(['RNA','protein','ATAC','methylation']):
        genes=[f'SYN{k+1:02d}' for k in range(24)]
        effect=np.array([1 if k%3 else -1 for k in range(24)],float);pv=np.array([10**(-5.5+.35*k+.18*i) if k<12 else .03+.03*(k-12)+.02*i for k in range(24)])
        effect[8:12]*=(-1 if i in [1,3] else 1)
        for k in [6+i,15+i,19+i]:effect[k]=np.nan;pv[k]=np.nan
        effect[22]=np.nan;pv[22]=np.nan;effect[23]=0;pv[23]=1
        pd.DataFrame({'feature':genes,'effect':effect,'p':pv}).to_csv(p/(n+'.tsv'),sep='\t',index=False)
        layers.append(dict(layer=n,results=n+'.tsv',feature_map='identity',orientation=1))
    pd.DataFrame(layers).to_csv(p/'manifest.tsv',sep='\t',index=False)
    multiomic_meta.run(SimpleNamespace(manifest=str(p/'manifest.tsv'),out=str(out/'joint'),min_layers=2,alpha=.05))
    spec=importlib.util.spec_from_file_location('demo',ROOT/'examples/pathways/make_demo.py');demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo);demo.create(out/'pathway_inputs');pathway_compare.run(out/'pathway_inputs/compare.json',out/'pathway_report')
    return [('multiomic-meta',out/'joint'),('pathway-compare',out/'pathway_report')]


def enrichment_example(out):
    p=out/'enrichment_inputs';p.mkdir();rng=np.random.default_rng(105);ranks=rng.normal(0,1,200);ranks[:20]+=2;ranks[20:40]-=2
    pd.DataFrame({'gene_id':[f'G{i}' for i in range(1,201)],'statistic':ranks}).to_csv(p/'ranks.tsv',sep='\t',index=False)
    sets={'SYNTHETIC_POSITIVE':[f'G{i}' for i in range(1,21)],'SYNTHETIC_NEGATIVE':[f'G{i}' for i in range(21,41)],'SYNTHETIC_MIXED':[f'G{i}' for i in range(5,150,5)]}
    (p/'sets.gmt').write_text(''.join(n+'\tsynthetic\t'+'\t'.join(gs)+'\n' for n,gs in sets.items()))
    cfg=dict(ranks='ranks.tsv',gmt='sets.gmt',species='human',gene_id_type='synthetic',rank_definition='other_reviewed_signed_statistic',annotation_reviewed=True)
    (p/'config.json').write_text(json.dumps(cfg));subprocess.run(['bash',str(ROOT/'run.sh'),'enrich',str(p/'config.json'),str(out/'enrichment')],check=True)
    return [('enrich',out/'enrichment')]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('--publish',type=Path);p.add_argument('--with-enrichment',action='store_true',help='Also run installed R/fgsea on fictional rankings');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    reports=[]
    for cmd,args,output in fixtures(a.output):
        print('Rendering '+cmd,flush=True);reports.append((cmd,render(cmd,args,output)))
    reports.extend(integration_examples(a.output))
    if a.with_enrichment:reports.extend(enrichment_example(a.output))
    if a.publish:publish(reports,a.publish)
    (a.output/'gallery_runs.json').write_text(json.dumps([{'command':c,'report':str(r)} for c,r in reports],indent=2)+'\n')


if __name__=='__main__':main()
