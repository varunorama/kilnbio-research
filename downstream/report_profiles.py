"""Reports for the typed outputs of KilnBio commands. No analysis refitting."""
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from analysis_visualizations import Builder, methods, loci, predictions, record_table, intervals, numbers
from visual_reports import canvas, matrix, empty, BLUE, ORANGE, MISSING
from plot_design import PALETTE


def option(args,flag,default=None):
    for i,a in enumerate(args):
        if a==flag and i+1<len(args):return args[i+1]
        if a.startswith(flag+'='):return a.split('=',1)[1]
    return default


def categorical_columns(b,t,cols,prefix):
    if t is not None:
        for c in cols:
            if c in t:b.categories(prefix+'_'+c,t[c],c.replace('_',' '))


def binary_colors(fig,ax):
    cmap=ListedColormap(['#fffaf2',ORANGE]);cmap.set_bad(MISSING)
    if ax.images:ax.images[-1].set_cmap(cmap);ax.images[-1].set_clim(0,1);fig.axes[-1].set_yticks([0,1])


def interval_estimates(b,t,prefix='comparison'):
    if t is None or t.empty:return
    lower=next((c for c in ['paired_lower','lower','ci_lower','lower_95'] if c in t),None)
    upper=next((c for c in ['paired_upper','upper','ci_upper','upper_95'] if c in t),None)
    value=next((c for c in ['reduction_vs_baseline','improvement_MAE','gain','estimate'] if c in t),None)
    if not all([lower,upper,value]):return
    d=t.copy();x=numbers(d[value]);lo=numbers(d[lower]);hi=numbers(d[upper]);ok=x.notna()&lo.notna()&hi.notna()
    d=d[ok];x=x[ok];lo=lo[ok];hi=hi[ok]
    label=next((c for c in ['model','comparison','policy'] if c in d),d.columns[0])
    fig,ax=canvas(len(d));yy=np.arange(len(d));ax.hlines(yy,lo,hi,color=BLUE);ax.scatter(x,yy,color=ORANGE);ax.axvline(0,color='#555',ls='--');ax.set_yticks(yy,d[label].astype(str));ax.set_xlabel(value);ax.set_title('Recorded comparison and uncertainty')
    b.add(prefix+'_intervals',fig,'Intervals are read from the completed analysis, not recomputed by the visualization. Their inferential scope is the scope of the source analysis.',d)


def audit(b,path,kind):
    obj=json.loads(Path(path).read_text());record_table(b,path,'audit_metrics')
    if kind=='candidate-evidence':
        rows=[]
        for entry in obj.get('candidates',[]):
            for e in entry.get('endpoints',[]):
                if isinstance(e,dict):rows.append(dict(gene=entry.get('gene',''),endpoint=str(e.get('assay',''))+' / '+str(e.get('cohort','')),status=e.get('status','unmeasured')))
        if rows:
            t=pd.DataFrame(rows);codes={'confirmed':1,'not_confirmed':0,'unmeasured':np.nan}
            t['display']=t.status.map(codes);p=t.pivot_table(index='gene',columns='endpoint',values='display',aggfunc='first',dropna=False)
            fig,ax=canvas(min(len(p),30),12);matrix(fig,ax,p.head(30),p.index[:30],p.columns,'Confirmation status (1 confirmed; 0 not confirmed)',1)
            binary_colors(fig,ax)
            b.add('candidate_confirmation',fig,'Measured endpoint confirmations remain separate from perturbation evidence. Missing evidence is not a negative experiment.',t)
    patterns=obj.get('patterns')
    if isinstance(patterns,dict):
        t=pd.DataFrame([{'pattern':k,'donors':v} for k,v in patterns.items() if isinstance(v,(int,float))])
        if len(t):
            fig,ax=canvas(len(t));ax.barh(t.pattern,t.donors,color=BLUE);ax.set_xlabel('Donors');ax.set_title('Assay availability patterns')
            b.add('assay_reach',fig,'Donor counts for the recorded assay-availability patterns; denominators come from this roster.',t)
    groups=[]
    for key,val in obj.items():
        if isinstance(val,dict) and 'gain' in val:
            rec=dict(policy=key,**val)
            if isinstance(val.get('interval'),list) and len(val['interval'])==2:rec.update(lower=val['interval'][0],upper=val['interval'][1])
            groups.append(rec)
    if groups:interval_estimates(b,pd.DataFrame(groups),'policy')


def singlecell(b,root):
    qc=b.load(root/'cell_qc.tsv')
    if qc is not None:
        b.numeric('cell_qc',qc,['total_counts','n_genes_by_counts','pct_counts_mt'],'Single-cell QC distributions')
        categorical_columns(b,qc,['qc_pass','sample'],'cells')
        if {'total_counts','n_genes_by_counts'}<=set(qc):
            fig,ax=canvas(10);ax.scatter(qc.total_counts,qc.n_genes_by_counts,s=3,alpha=.4,c=np.where(qc.get('qc_pass',pd.Series(True,index=qc.index)).astype(bool),BLUE,ORANGE));ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=1);ax.set_xlabel('Total counts');ax.set_ylabel('Detected genes');ax.set_title('Cell complexity and QC')
            b.add('cell_complexity',fig,'Each point is a cell read from cell_qc.tsv; orange indicates a cell excluded by the recorded QC flags.',qc)
    path=root/'filtered.h5ad'
    if not path.exists():b.skip('cell_embedding','No fitted single-cell AnnData output');return
    import anndata
    ad=anndata.read_h5ad(path,backed='r')
    try:
        n=min(ad.n_obs,10000);obs=ad.obs.iloc[:n].copy();obs.index=obs.index.astype(str)
        if 'X_umap' in ad.obsm:
            xy=np.asarray(ad.obsm['X_umap'][:n]);fields=[c for c in ['leiden','sample','total_counts','pct_counts_mt'] if c in obs]
            for c in fields:
                fig,ax=canvas(10);s=obs[c]
                if pd.api.types.is_numeric_dtype(s):im=ax.scatter(xy[:,0],xy[:,1],c=numbers(s),s=3,cmap='viridis');fig.colorbar(im,ax=ax,label=c)
                else:
                    cats=s.astype(str)
                    for i,name in enumerate(sorted(cats.unique())):
                        k=(cats==name).to_numpy();ax.scatter(xy[k,0],xy[k,1],s=8,label=f'{name} (n={int(k.sum()):,})',color=PALETTE[i%len(PALETTE)],alpha=.7,edgecolors='none')
                    if cats.nunique()<=20:ax.legend(fontsize=7,markerscale=3)
                ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_title('Recorded UMAP / '+c)
                b.add('umap_'+c,fig,'Uses the stored embedding; no re-clustering. First 10,000 retained cells displayed. Clusters are exploratory groups, not automatically validated cell types.',pd.DataFrame({'cell':obs.index,'UMAP1':xy[:,0],'UMAP2':xy[:,1],c:s.to_numpy()}))
        if {'sample','leiden'}<=set(obs):
            ct=pd.crosstab(obs['sample'],obs['leiden'],normalize='index');fig,ax=canvas(len(ct),10);ct.plot.bar(stacked=True,ax=ax,color=[PALETTE[i%len(PALETTE)] for i in range(len(ct.columns))]);ax.set_ylabel('Cell fraction');ax.set_title('Cluster composition by sample');ax.legend(title='Recorded cluster',loc='upper left',bbox_to_anchor=(1.02,1))
            b.add('cell_composition',fig,'Descriptive fractions among the first 10,000 retained cells. Donors, not individual cells, are the replication unit for comparative inference.',ct.reset_index())
        if 'leiden' in obs and ad.n_vars:
            # Use declared HVGs if available; do not label these as tested marker genes.
            ix=np.flatnonzero(ad.var.get('highly_variable',pd.Series(True,index=ad.var_names)).to_numpy())[:20]
            if not len(ix):ix=np.arange(min(20,ad.n_vars))
            xx=ad[:n,ix].X;xx=xx.toarray() if hasattr(xx,'toarray') else np.asarray(xx)
            rows=[]
            for group in sorted(obs.leiden.astype(str).unique()):
                take=(obs.leiden.astype(str)==group).to_numpy()
                for j,gene in enumerate(ad.var_names[ix]):rows.append(dict(cluster=group,gene=str(gene),mean_log_expression=float(np.mean(xx[take,j])),fraction_nonzero=float(np.mean(xx[take,j]>0))))
            t=pd.DataFrame(rows);gs=list(t.gene.drop_duplicates());cs=list(t.cluster.drop_duplicates());fig,ax=canvas(len(gs),max(7.5,min(14,4+len(cs)*.55)))
            for j in range(len(gs)):
                if j%2==0:ax.axhspan(j-.5,j+.5,color='#F3F6FA',zorder=0)
            im=ax.scatter(t.cluster.map({x:i for i,x in enumerate(cs)}),t.gene.map({x:i for i,x in enumerate(gs)}),s=180*t.fraction_nonzero,c=t.mean_log_expression,cmap='viridis',edgecolors='white',linewidths=.6,zorder=3)
            fig.colorbar(im,ax=ax,label='Mean log-normalized expression',shrink=.65,pad=.04)
            cn=obs.leiden.astype(str).value_counts();ax.set_xticks(range(len(cs)),[f'{c}\nn={cn[c]:,}' for c in cs]);ax.set_yticks(range(len(gs)),gs)
            ax.set_xlim(-.6,len(cs)-.4);ax.set_ylim(len(gs)-.5,-.5);ax.set_xlabel('Recorded Leiden cluster');ax.set_title('Expression and detection by cluster')
            handles=[ax.scatter([],[],s=180*k,color='#637587',edgecolors='white',label=f'{k:.0%}') for k in [.25,.5,.75,1]]
            ax.legend(handles=handles,title='Cells with a nonzero measurement',ncol=4,loc='upper center',bbox_to_anchor=(.5,-.13),columnspacing=1.4)
            b.add('cell_expression_dotplot',fig,'Dot area is the nonzero fraction; color is mean stored log-normalized expression. First 20 declared variable genes (or genes if flags absent). This is descriptive, not a differential marker test.',t)
    finally:ad.file.close()


def state_matrix(b,t,prefix):
    if t is None:return
    state=next((c for c in ['state','ch_status','genomic_state'] if c in t),None)
    if state:
        b.categories(prefix+'_states',t[state],'Recorded genomic states')
        sample=next((c for c in ['sample','source_sample','specimen','donor'] if c in t),None)
        feature=next((c for c in ['event','gene','variant'] if c in t),None)
        if sample and feature:
            d=t.copy();d['display']=d[state].map({'altered':1,'reference':0,'unknown':np.nan})
            if not d.duplicated([sample,feature]).any():
                h=d.pivot(index=feature,columns=sample,values='display').iloc[:30,:60];fig,ax=canvas(len(h),12);matrix(fig,ax,h,h.index,h.columns,'Altered (1), reference (0), unknown (gray)',1)
                binary_colors(fig,ax)
                b.add(prefix+'_oncoprint',fig,'States come from reviewed calls. Missing or unknown entries are never converted to reference; this view does not infer clinical CH status.',d)


def proteins(b,t):
    if t is None:return
    categorical_columns(b,t,['mapping_status','tool_used'],'structure')
    b.numeric('structure_confidence',t,['mean_plddt'],'Recorded model confidence')
    cols={'model_sequence_position','sequence_start','sequence_end'}
    if cols<=set(t):
        d=t.head(30);fig,ax=canvas(len(d),12)
        for i,row in enumerate(d.to_dict('records')):
            lo=numbers(pd.Series([row['sequence_start']])).iloc[0];hi=numbers(pd.Series([row['sequence_end']])).iloc[0];pos=numbers(pd.Series([row['model_sequence_position']])).iloc[0]
            if np.isfinite(lo) and np.isfinite(hi):ax.hlines(i,lo,hi,color=MISSING,lw=5)
            if np.isfinite(pos):ax.scatter(pos,i,color=ORANGE)
        labels=d.get('protein_variant',d.iloc[:,0]).astype(str);ax.set_yticks(range(len(d)),labels);ax.set_xlabel('Sequence coordinate supplied by the lookup');ax.set_title('Variant positions in available model spans')
        b.add('protein_positions',fig,'Model span and mapped position from the lookup receipt, not a rendered protein structure or a pathogenicity prediction.',d)


def upstream(b,root):
    # Link existing FastQC/MultiQC/Nextflow reports without recomputing diagnostics.
    found=[];beds=[];seen=0
    for directory,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=[d for d in dirs if d not in {'work','.git','figures','plotted_data'} and 'visualization' not in d]
        for name in sorted(files):
            seen+=1;p=Path(directory)/name
            if p.suffix.lower() in ['.html','.pdf'] and any(s in p.name.lower() for s in ['multiqc','fastqc','report','qc','enrichment','diagnostic','voom']):found.append(p.relative_to(root).as_posix())
            if p.suffix.lower() in ['.bed','.narrowpeak','.broadpeak'] and len(beds)<3:beds.append(p)
        if len(found)>=200 or seen>=20000:break
    b.report.table('upstream_reports',pd.DataFrame({'relative_path':found}))
    if found:
        b.report.notes.append('Existing upstream reports (paths relative to the analysis output): '+', '.join(found))
        b.report.links.extend((name,'../'+name) for name in found)
        fig,ax=canvas(min(len(found),20),13);ax.axis('off');ax.text(0,1,'\n'.join(found[:20]),va='top',fontsize=9,transform=ax.transAxes);ax.set_title('Upstream QC reports')
        b.add('upstream_reports',fig,'Native FastQC/MultiQC/model reports retain their original algorithms and detailed plots. Paths are indexed here; no new QC score is invented.',pd.DataFrame({'relative_path':found}))
    else:b.skip('upstream_reports','No native QC HTML/PDF found')
    for i,p in enumerate(beds):intervals(b,p,'intervals_'+str(i+1)+'_'+p.stem)


def vep(b,path):
    if not path.exists():b.skip('consequences','No completed annotated VCF');return
    rows=[];fields=[]
    with path.open() as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ') and 'Format: ' in line:fields=line.split('Format: ',1)[1].split('"',1)[0].split('|')
            if line.startswith('#'):continue
            parts=line.rstrip().split('\t')
            if len(parts)<8:continue
            info=dict(v.split('=',1) for v in parts[7].split(';') if '=' in v)
            for value in info.get('CSQ','').split(','):
                if fields and value:rows.append(dict(zip(fields,value.split('|'))))
            if len(rows)>=100000:break
    t=pd.DataFrame(rows[:100000]);categorical_columns(b,t,['Consequence','IMPACT','BIOTYPE'],'vep')
    b.report.notes.append('VEP distributions count transcript-consequence records, not unique people or pathogenic variants; first 100,000 records maximum.')


def nanopore(b,root):
    path=root/'cpg.bed';t=b.load(path,header=None)
    if t is not None and t.shape[1]>=13:
        for channel,sub in t.groupby(3):
            # modkit bedMethyl: start=0-based; valid coverage col10; modified count col12.
            c=pd.DataFrame({'chr':sub[0],'pos':numbers(sub[1])+1,'N':numbers(sub[9]),'X':numbers(sub[11])})
            loci(b,c,'modification_'+str(channel))
    record_table(b,root/'status.json','nanopore_status')


def coverage_bar(b,t):
    if t is None or not {'required','present'}<=set(t):return
    d=t.copy();d['fraction_present']=numbers(d.present)/numbers(d.required).where(numbers(d.required)>0)
    labels=d['model'].astype(str) if 'model' in d else d.index.astype(str)
    fig,ax=canvas(len(d));ax.barh(labels,d.fraction_present,color=BLUE);ax.set_xlim(0,1);ax.set_xlabel('Required features present / required features');ax.set_title('Biomarker feature coverage')
    b.add('clock_coverage',fig,'Feature presence is distinct from complete profiles and successful scoring. Coverage does not establish transferability to a new tissue or assay.',d)


def methylation_locus(b,root):
    fractions=b.load(root/'display_fractions.tsv');depth=b.load(root/'display_depth.tsv')
    if fractions is None or depth is None:return
    ids=fractions.iloc[:,0].astype(str).str.rsplit(':',n=1,expand=True)
    if ids.shape[1]!=2:return
    pos=numbers(ids[1]);chrom=ids[0].iloc[0];start=pos.iloc[0];sel=(ids[0]==chrom)&pos.between(start,start+2000)
    fig,(ax,dd)=plt.subplots(2,1,figsize=(12,7),layout='constrained',sharex=True);rows=[]
    for i,s in enumerate(fractions.columns[1:13]):
        x=pos[sel];y=numbers(fractions.loc[sel,s]);n=numbers(depth.loc[sel,s]);color=plt.get_cmap('tab20')(i)
        ax.plot(x,y,'.-',lw=.5,ms=3,color=color,label=s);dd.plot(x,n,'.-',lw=.5,ms=3,color=color)
        rows.extend(dict(sample=s,chrom=chrom,pos=p,fraction=f,depth=c) for p,f,c in zip(x,y,n))
    ax.set_ylabel('Modified / total');ax.set_ylim(0,1);dd.set_ylabel('Depth');dd.set_xlabel(f'{chrom}: one-based cytosine position');ax.legend(fontsize=7,ncol=4);ax.set_title('Per-sample modification and coverage')
    b.add('methylation_sample_locus',fig,'First 2 kb of retained covered loci and first 12 samples. Lines connect measured positions for display; they are not interpolation, smoothing or a DMR call. Missing values remain gaps.',pd.DataFrame(rows))


def clock_age(b,root,args):
    path=option(args,'--metadata')
    if not path:return
    meta=b.load(path);scores=b.load(root/'scores.tsv');coverage=b.load(root/'coverage.tsv')
    if meta is None or scores is None or coverage is None or not {'sample','age'}<=set(meta):return
    if not {'model','endpoint'}<=set(coverage):return
    names=coverage.loc[coverage.endpoint.astype(str).str.lower().isin(['age','dnamage','age (years)','years']),'model']
    merged=scores.merge(meta[['sample','age']],on='sample',validate='one_to_one',how='left')
    for name in names:
        if name not in merged:continue
        x=numbers(merged.age);y=numbers(merged[name]);ok=x.notna()&y.notna();fig,ax=canvas(10)
        ax.scatter(x[ok],y[ok],color=BLUE)
        if ok.any():lo=min(x[ok].min(),y[ok].min());hi=max(x[ok].max(),y[ok].max());ax.plot([lo,hi],[lo,hi],ls='--',color='#555')
        ax.set_xlabel('Chronological age (years)');ax.set_ylabel(name+' (years)');ax.set_title('Age-score comparison')
        b.add('age_'+name,fig,'Only models explicitly declaring an age endpoint are compared with chronological age. This is a descriptive score check, not held-out validation unless the study design establishes independence.',merged[['sample','age',name]])


def decision(b,path):
    obj=json.loads(Path(path).read_text());r=obj.get('result',{});record_table(b,path,'method_metrics')
    if not isinstance(r,dict):return
    if {'features','effect','p'}<=set(r):b.stats(pd.DataFrame({k:r[k] for k in ['features','effect','p']}),'fitted_contrast',feature='features')
    if isinstance(r.get('status'),dict):b.categories('validation_status',pd.Series(r['status']),'Independent endpoint status')
    for col in ['phase','replicate_fractions']:
        if isinstance(r.get(col),list):
            if col=='phase':b.categories('donor_partitions',pd.Series(r[col]),'Donor partition assignments')
            else:b.numeric('element_fractions',pd.DataFrame({col:r[col]}),[col],'Recorded element methylation by replicate')
    if {'features','arms'}<=set(r):
        d=pd.DataFrame({k:pd.Series(np.arange(1,len(v['ranking'])+1),index=v['ranking']) for k,v in r['arms'].items()});show=d.loc[d.min(axis=1).sort_values().head(30).index]
        fig,ax=canvas(len(show),12);matrix(fig,ax,show,show.index,show.columns,'Rank (smaller is higher)');ax.set_title('Ranks under each declared assay choice')
        b.add('assay_rank_matrix',fig,'Same declared feature family across the arms. Rank is descriptive; it is not a p-value or a measure of independent validation.',d.reset_index(names='feature'))
    for i,transition in enumerate(r.get('transitions',[])):
        d=pd.DataFrame(list(transition.get('rank_change',{}).items()),columns=['feature','rank_change']);show=d.reindex(d.rank_change.abs().sort_values(ascending=False).head(30).index)
        fig,ax=canvas(len(show));ax.barh(show.feature,show.rank_change,color=np.where(show.rank_change>0,BLUE,ORANGE));ax.set_xlabel('After rank - before rank');ax.set_title(str(transition.get('from'))+' to '+str(transition.get('to')))
        b.add('rank_transition_'+str(i+1),fig,'Negative change means a better numerical rank after the specified transition; same feature family, without a causal attribution.',d)


def build(command,args,output,dest):
    root=Path(output);b=Builder(dest,command);isdir=root.is_dir();t=None
    if command in ['de','dm']:
        r=root/'results';meta=b.load(r/'display_samples.tsv')
        for file in ['all_features.tsv','all_loci.tsv','all_regions.tsv']:
            t=b.load(r/file)
            if t is not None:b.stats(t,file[:-4])
        for file,scale in [('display_transformed.tsv','model transformation'),('display_fractions.tsv','modified / total fraction')]:
            if (r/file).exists():b.matrix_qc(b.load_matrix(r/file),file[:-4],meta,scale)
        t=b.load(r/'sample_diagnostics.tsv')
        if t is not None:b.numeric('library_diagnostics',t,['library_size','normalization_factor','effective_library_size','size_factor'],'Library normalization diagnostics')
        t=b.load(r/'display_depth.tsv')
        if t is not None:b.numeric('methylation_depth',pd.DataFrame({'depth':t.iloc[:,1:].to_numpy().ravel()}),['depth'],'Depth at displayed eligible loci')
        scope=b.load(r/'display_scope.tsv')
        if scope is not None:b.report.notes.append('Methylation display subset: '+scope.to_json(orient='records'))
        if command=='dm':methylation_locus(b,r)
        upstream(b,root)
    elif command=='multiomic-layers':
        manifest=b.load(root/'manifest.tsv')
        if manifest is not None:
            for row in manifest.itertuples():
                b.stats(b.load(root/row.results),str(row.layer));b.matrix_qc(b.load_matrix(root/(row.layer+'_display.tsv')),str(row.layer),b.load(root/(row.layer+'_samples.tsv')))
        upstream(b,root)
    elif command=='differential-compare':
        methods(b,b.load(root/'features.tsv'));b.matrix_qc(b.load_matrix(args[1]),'input',b.load(args[2]),'counts' if args[0]=='rna-counts' else 'supplied log2 scale')
        categorical_columns(b,b.load(root/'eligibility.tsv'),['eligibility'],'feature')
    elif command=='method-agreement':
        methods(b,b.load(args[0]),float(option(args,'--alpha','.05')));categorical_columns(b,b.load(root),['status'],'feature')
    elif command=='olink-import':
        b.matrix_qc(b.load_matrix(root/'npx.tsv'),'NPX',b.load(root/'samples.tsv'),'NPX (log2)')
        categorical_columns(b,b.load(root/'measurement_qc.tsv'),['SampleQC','AssayQC','above_LOD'],'measurement')
        a=b.load(root/'aliasing.tsv')
        if a is not None:b.numeric('metadata_association',a,['cramers_v'],'Recorded metadata association strength')
        record_table(b,root/'qc.json','olink_receipt')
    elif command=='scrna':singlecell(b,root)
    elif command=='pseudobulk':
        for p in sorted(root.glob('celltype_*/counts.tsv')):
            b.matrix_qc(b.load_matrix(p),p.parent.name,b.load(p.parent/'samples.csv'),'counts')
        t=b.load(root/'excluded_groups.tsv')
        if t is not None:b.numeric('excluded_groups',t,['cells'],'Cells in excluded donor/cell-type groups')
    elif command=='featurecounts':b.matrix_qc(b.load_matrix(root/'counts.tsv'),'counts',scale='counts')
    elif command=='import-mods':loci(b,b.load(root/'counts.tsv'));record_table(b,root/'status.json','chemistry_and_scope')
    elif command in ['compare-assays','compare-predictions']:
        predictions(b,b.load(root/'predictions.tsv') if command=='compare-assays' else b.load(args[0]))
        for f in ['paired_differences.tsv','paired_mae.tsv']:interval_estimates(b,b.load(root/f))
    elif command=='decision-methods':decision(b,root)
    elif command in ['assay-audit','candidate-evidence','epidiff-coverage']:
        audit(b,root,command)
    elif command=='capture-panel':
        intervals(b,root/'targets.bed','panel');categorical_columns(b,b.load(root/'features.tsv'),['fully_on_target'],'features');record_table(b,root/'summary.json','panel_coverage')
    elif command in ['genomic-panel','ch-screen','ch-import']:
        if not isdir:state_matrix(b,b.load(root),'genomic')
        else:
            state_matrix(b,b.load(root/'samples.tsv'),'CH')
            for f in ['variants.tsv','screen.tsv']:
                t=b.load(root/f)
                if t is not None:
                    b.numeric(f[:-4]+'_alleles',t,['vaf','dp','alt_count'],'Variant allele support');categorical_columns(b,t,['gene','candidate','reasons','evidence','review'],f[:-4])
    elif command=='genomic-associate':b.stats(b.load(root/'association.tsv'),'genomic_association');state_matrix(b,b.load(root/'linkage.tsv'),'linked_assay')
    elif command=='variant-annotations':
        t=b.load(root/'variants.tsv');categorical_columns(b,t,['clinvar_status','alphamissense_status','vep_status'],'annotation')
    elif command=='alphafold':proteins(b,b.load(root/'structures.tsv'))
    elif command=='vep-annotate':vep(b,root/'annotated.vcf')
    elif command=='ont':nanopore(b,root)
    elif command=='alphagenome':
        t=b.load(root/'scores.tsv')
        if t is not None:
            categorical_columns(b,t,['output_type','biosample_name','variant_id'],'prediction')
            keys=[c for c in ['output_type','variant_scorer'] if c in t]
            groups=t.groupby(keys,dropna=False,sort=False) if keys else [('declared scoring output',t)]
            for i,(label,group) in enumerate(groups):
                if i>=12:b.skip('additional_regulatory_scores','More than 12 scoring strata; inspect the source table');break
                b.numeric('regulatory_predictions_'+str(i+1),group,['raw_score','quantile_score','score'],'Provider scores: '+str(label))
        record_table(b,root/'receipt.json','prediction_receipt')
    elif command in ['clock-panel','clocks','assay-clocks','erythroid-stage']:
        t=b.load(root/'scores.tsv') if isdir else b.load(root)
        if t is not None:
            cols=[c for c in t.columns[1:] if numbers(t[c]).notna().any()]
            for c in cols:b.numeric('score_'+c,t,[c],'Recorded biomarker score: '+c)
            if len(cols)>1:
                corr=t[cols].apply(numbers).corr(min_periods=3);fig,ax=canvas(len(cols));matrix(fig,ax,corr,cols,cols,'Pairwise Pearson correlation',1)
                b.add('clock_correlation',fig,'Pairwise complete profiles (minimum three). Scores retain different endpoints; correlation does not make them interchangeable.',corr.reset_index(names='model'))
        for file in ['coverage.tsv','sample_qc.tsv']:
            tt=b.load(root/file) if isdir else None
            if tt is not None:b.numeric(file[:-4],tt,[c for c in ['required','present','complete_profiles','scored_profiles','fraction_missing','missing_fraction'] if c in tt],file[:-4].replace('_',' '))
            if file=='coverage.tsv':coverage_bar(b,tt)
        if command=='erythroid-stage' and isdir:
            pp=b.load(root/'predictions.tsv')
            if pp is not None:predictions(b,pp.rename(columns={'stage':'observed','estimate':'predicted','assay':'model'}))
        b.report.notes.append('Each biomarker retains its own endpoint and units. Score correlation is descriptive; pace, risk and age estimates are not interchangeable.')
        if command=='clock-panel':clock_age(b,root,args)
    elif command=='spikein-factors':
        t=b.load(root/'diagnostics.tsv')
        if t is not None:
            b.numeric('spikein_counts',t,['target_reads','reference_reads','size_factor'],'Spike-in reference diagnostics')
            if {'target_reads','reference_reads'}<=set(t):
                fig,ax=canvas(10);ax.scatter(t.reference_reads,t.target_reads,color=BLUE);ax.set_xlabel('Reference reads');ax.set_ylabel('Target reads');ax.set_title('Target versus reference exposure');b.add('spikein_exposure',fig,'Recorded target and reference counts. Reference exposure assumptions are specified by the analysis, not inferred from this correlation.',t)
    elif command=='graft-readouts':
        t=b.load(root/'lots.tsv')
        if t is not None:
            b.numeric('graft_readouts',t,['predictor','median_outcome','enrolled_recipients','measured_recipients'],'Lot-level readouts');b.report.notes.append('Recipient measurements are nested within lots/donors; these distributions do not treat recipients as independent lots.')
    elif command=='bladder-example':
        from bladder_visualizations import bladder_report
        bladder_report(b,root)
        record_table(b,root/'report.json','recipe_denominators')
    else:upstream(b,root)
    b.finish()
    return b
