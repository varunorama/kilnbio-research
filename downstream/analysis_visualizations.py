"""Typed, bounded visualizations of completed analyses; no inferential refitting."""
import hashlib
import html
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from visual_reports import Report, canvas, matrix, empty, strength, BLUE, ORANGE, MISSING
from plot_design import MUTED, PALETTE

MAX_ROWS=100000


def read(path, limit=MAX_ROWS, **kwargs):
    path=Path(path)
    if not path.is_file():return None
    try:
        table=pd.read_csv(path,sep=',' if path.suffix=='.csv' else '\t',nrows=limit+1,**kwargs)
    except pd.errors.EmptyDataError:return pd.DataFrame()
    table.attrs['truncated']=len(table)>limit
    return table.iloc[:limit].copy()


def safe(value):return re.sub(r'[^A-Za-z0-9_-]','_',str(value))[:100]


def numbers(series):return pd.to_numeric(series,errors='coerce').replace([np.inf,-np.inf],np.nan)


class Builder:
    def __init__(self,out,command):
        self.report=Report(out,'KilnBio '+command+' report',[
            'Automatically generated from this completed analysis. Figures are descriptive unless they display statistics already computed by the analysis.',
            'Missing measurements remain missing. Display limits and unavailable panels are recorded in panel_status.tsv. No LLM calls or remote services are used by reporting.'])
        self.status=[];self.sources=[];self.command=command

    def load(self,path,**kwargs):
        t=read(path,**kwargs)
        if t is not None:
            p=Path(path);stat=p.stat()
            # A bounded digest avoids scanning multi-GB outputs a second time.
            with p.open('rb') as f:sha=hashlib.sha256(f.read(8*1024*1024)).hexdigest()
            self.sources.append(dict(file=p.name,size_bytes=stat.st_size,sha256_first_8MiB=sha,full_hash=stat.st_size<=8*1024*1024,rows_read=len(t),truncated=t.attrs.get('truncated',False)))
        return t

    def skip(self,name,reason):self.status.append(dict(panel=name,status='unavailable',reason=reason))

    def load_matrix(self,path):
        path=Path(path)
        if not path.is_file():return None
        columns=pd.read_csv(path,sep=',' if path.suffix=='.csv' else '\t',nrows=0).columns.tolist()
        t=self.load(path,usecols=columns[:301])
        self.sources[-1].update(columns_available=len(columns),columns_read=min(301,len(columns)))
        return t

    def add(self,name,fig,caption,data):
        name=safe(name)
        self.report.add(name,fig,caption,data)
        self.status.append(dict(panel=name,status='generated',reason=caption))

    def finish(self):
        if not self.report.panels:
            fig,ax=canvas(4);empty(ax,'No measured data available for the requested panels.\nSee panel_status.tsv and the upstream receipts.')
            self.report.add('availability',fig,'This is an availability notice, not a biological result.',pd.DataFrame(self.status))
        self.report.table('panel_status',pd.DataFrame(self.status,columns=['panel','status','reason']))
        self.report.finish(dict(command=self.command,max_table_rows=MAX_ROWS,input_files=self.sources,panel_status=self.status))

    def categories(self,name,series,title,caption='Counts of records; not a significance test.'):
        if series is None or not len(series):self.skip(name,'No categorical records');return
        counts=series.fillna('Missing').astype(str).value_counts(sort=True).rename_axis('category').reset_index(name='records')
        fig,ax=canvas(min(len(counts),25));shown=counts.head(25)
        ax.barh(shown.category,shown.records,color=BLUE);ax.invert_yaxis();ax.set_xlabel('Records');ax.set_title(title)
        self.add(name,fig,caption+' Largest 25 categories shown; TSV includes all categories read.',counts)

    def numeric(self,name,frame,columns,title):
        cols=[c for c in columns if c in frame and numbers(frame[c]).notna().any()]
        if not cols:self.skip(name,'No finite values in expected columns: '+', '.join(columns));return
        fig,axes=plt.subplots(1,len(cols),figsize=(min(18,4*len(cols)),4),layout='constrained',squeeze=False)
        data=frame[cols].apply(numbers)
        for ax,col in zip(axes[0],cols):
            ax.hist(data[col].dropna(),bins=30,color=BLUE,alpha=.85);ax.set_xlabel(col);ax.set_ylabel('Records')
        fig.suptitle(title)
        self.add(name,fig,'Separate scales preserve each quantity. Histograms use the available finite values; no missing values are imputed. At most 100,000 input rows are read.',data)

    def stats(self,table,name='differential',effect=None,p=None,q=None,feature=None,abundance=None):
        if table is None or table.empty:self.skip(name,'No differential-result rows');return
        def choose(given,candidates):return given or next((c for c in candidates if c in table),None)
        effect=choose(effect,['log2FoldChange','logFC','effect','median_difference','unadjusted_delta_fraction','stat'])
        p=choose(p,['pvalue','PValue','p','pvals','pval']);q=choose(q,['padj','FDR','q','FDR_BH','qval','adj.P.Val'])
        feature=choose(feature,['feature_id','feature','gene','probe']);abundance=choose(abundance,['baseMean','logCPM','AveExpr','mean'])
        if effect is None or p is None:self.skip(name,'Requires an explicitly recognized effect and raw p-value column');return
        d=table.copy();d['display_effect']=numbers(d[effect]);d['display_p']=numbers(d[p]);d['display_q']=numbers(d[q]) if q else np.nan
        d['label']=d[feature].astype(str) if feature else d.index.astype(str)
        if ((d.display_p.dropna()<0)|(d.display_p.dropna()>1)).any():raise ValueError('Invalid p-value in report input')
        valid=d.display_effect.notna()&d.display_p.notna();sig=d.display_q.le(.05)
        fig=plt.figure(figsize=(11,6.4),layout='constrained');gs=fig.add_gridspec(2,1,height_ratios=[1,4],hspace=.04)
        marginal=fig.add_subplot(gs[0]);ax=fig.add_subplot(gs[1],sharex=marginal)
        groups=[(valid&~sig,'Other / q unavailable','#B7C3CE'),(valid&sig&(d.display_effect<0),'Negative · q ≤ 0.05',BLUE),(valid&sig&(d.display_effect>=0),'Nonnegative · q ≤ 0.05',ORANGE)]
        for mask,label,color in groups:
            ax.scatter(d.loc[mask,'display_effect'],strength(d.loc[mask,'display_p']),s=30 if color!='#B7C3CE' else 16,c=color,alpha=.85 if color!='#B7C3CE' else .55,edgecolors='white',linewidths=.35,label=f'{label} (n={int(mask.sum()):,})',zorder=3)
        if valid.any():marginal.hist(d.loc[valid,'display_effect'],bins=35,color='#C6D7E4',edgecolor='white',lw=.4)
        marginal.set_ylabel('Count',fontsize=9);marginal.tick_params(labelbottom=False);marginal.spines['bottom'].set_visible(False)
        marginal.set_title(f'Effect and evidence  |  {int(valid.sum()):,} plotted / {len(d):,} rows',fontsize=14)
        ax.axvline(0,color='#A1B0BE',lw=.8,ls='--');ax.set_ylim(bottom=0);ax.margins(x=.12,y=.15)
        for positive in [False,True]:
            label_rows=d[valid&sig&(d.display_effect.ge(0) if positive else d.display_effect.lt(0))].sort_values(['display_p','label'],kind='stable').head(3)
            for i,(_,row) in enumerate(label_rows.iterrows()):
                ax.annotate(str(row.label),xy=(row.display_effect,float(strength([row.display_p])[0])),xytext=(.97 if positive else .03,.93-.13*i),textcoords='axes fraction',ha='right' if positive else 'left',fontsize=9,color=MUTED,arrowprops=dict(arrowstyle='-',color='#9BABBA',lw=.6),bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2))
        ax.set_xlabel(effect);ax.set_ylabel('−log₁₀(raw p), display cap 16');ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=9)
        self.add(name+'_volcano',fig,'Blue/orange mark negative/nonnegative effects with supplied adjusted p <= 0.05; gray includes other points and unavailable q. The marginal histogram shows the plotted effect distribution. At most three qualifying features per direction are labeled by raw p then identifier. No horizontal FDR line is inferred on the raw-p axis. The analysis discovery threshold and result table remain authoritative.',d)
        if abundance:
            a=numbers(d[abundance]);display=np.log2(a+1) if abundance=='baseMean' else a
            fig,ax=canvas(10);ax.scatter(display,d.display_effect,s=8,c=np.where(sig,ORANGE,BLUE),alpha=.5);ax.axhline(0,color='#666',lw=.7);ax.set_xlabel('log2(baseMean+1)' if abundance=='baseMean' else abundance);ax.set_ylabel(effect);ax.set_title('Abundance versus effect')
            dd=d.copy();dd['display_abundance']=display;self.add(name+'_ma',fig,'MA/MD-style descriptive diagnostic; original model values are not altered.',dd)
        pp=np.sort(d.display_p.dropna().to_numpy());pp=pp[(pp>=0)&(pp<=1)]
        fig,ax=canvas(10)
        if len(pp):
            exp=-np.log10((np.arange(1,len(pp)+1)-.5)/len(pp));obs=strength(pp);ax.scatter(exp,obs,s=6,color=BLUE);ax.plot([0,max(exp)],[0,max(exp)],color='#555',ls='--');ax.set_xlabel('Expected -log10(uniform p)');ax.set_ylabel('Observed -log10(p), capped at 16')
        else:empty(ax,'No tested p-values')
        ax.set_title('P-value QQ diagnostic')
        self.add(name+'_qq',fig,'A diagnostic comparison with uniform p-values. Alternatives, dependence, discretization and confounding can all affect this plot; it is not a calibration proof.',pd.DataFrame({'p_sorted':pp,'expected_uniform_quantile':(np.arange(len(pp))+.5)/len(pp) if len(pp) else []}))
        chrom=next((c for c in ['chr','chrom','seqnames'] if c in d),None);pos=next((c for c in ['pos','start'] if c in d),None)
        if chrom and pos:
            fig,ax=canvas(10,12);offset=0;ticks=[];labels=[];coords=[]
            spans=[float(numbers(t[pos]).max()-numbers(t[pos]).min()) for _,t in d.groupby(chrom) if numbers(t[pos]).notna().any()]
            gap=max(1,sum(spans)*.035)
            for i,(c,sub) in enumerate(d.groupby(chrom,sort=True)):
                positions=numbers(sub[pos]);ok=positions.notna()&sub.display_p.notna();xx=positions[ok]-positions[ok].min()+offset
                ax.scatter(xx,strength(sub.loc[ok,'display_p']),s=16,color=BLUE if i%2==0 else ORANGE,alpha=.8)
                if ok.any():ticks.append((xx.min()+xx.max())/2);labels.append(str(c));offset=float(xx.max())+gap
                coords.extend(dict(feature=k,chrom=c,position=v,display_x=x,raw_p=pv) for k,v,x,pv in zip(sub.loc[ok,'label'],positions[ok],xx,sub.loc[ok,'display_p']))
            ax.set_xticks(ticks,labels,rotation=0);ax.set_ylabel('−log₁₀(raw p), display cap 16');ax.set_xlabel('Observed spans, concatenated by contig');ax.set_title('Where tested evidence occurs')
            self.add(name+'_manhattan',fig,'Observed contig spans are concatenated lexically with gaps equal to 3.5% of their summed spans (minimum one coordinate unit); unobserved contig prefixes are trimmed for display. Coordinates in the TSV retain the source convention. This is not a full-genome coverage map. No genome-wide threshold is inferred.',pd.DataFrame(coords))

    def matrix_qc(self,frame,prefix='matrix',metadata=None,scale='already transformed'):
        if frame is None or frame.empty or len(frame.columns)<2:self.skip(prefix+'_matrix','No feature-by-sample matrix');return
        ids=frame.iloc[:,0].astype(str);x=frame.iloc[:,1:301].apply(numbers);x.index=ids
        if ids.duplicated().any():raise ValueError('Duplicate feature identifiers in report matrix')
        if scale=='counts':
            if (x<0).any().any():raise ValueError('Negative counts in display input')
            lib=x.sum();x=np.log2(x.div(lib.where(lib>0),axis=1)*1e6+1)
            scale='log2(CPM+1), visualization only'
        if scale=='nonnegative':x=np.log2(x+1);scale='log2(value+1), visualization only'
        missing=x.isna().mean().rename_axis('sample').reset_index(name='fraction_missing')
        fig,ax=canvas(8,12);ax.bar(np.arange(len(missing)),missing.fraction_missing,color=BLUE);ax.set_xticks(range(len(missing)),missing['sample'],rotation=90,fontsize=6);ax.set_ylabel('Fraction missing');ax.set_title(prefix+' missing measurements')
        self.add(prefix+'_missingness',fig,'Fraction of missing values among the feature rows read. Zero values remain observed measurements.',missing)
        selected=x.columns[:60];fig,ax=canvas(9,12)
        vals=[x[c].dropna().to_numpy() for c in selected];boxes=ax.boxplot(vals,showfliers=False)
        for box in boxes['boxes']:box.set_facecolor('#D8E6F0');box.set_linewidth(.9)
        ax.set_xticks(range(1,len(selected)+1),selected);ax.tick_params(axis='x',labelrotation=90,labelsize=6);ax.set_ylabel(scale);ax.set_title(prefix+' distributions')
        self.add(prefix+'_distributions',fig,'First 60 sample columns displayed, up to 300 analyzed. Whiskers use the standard 1.5-IQR rule; outliers are hidden in the figure but retained in plotted data.',x[selected].reset_index(names='feature'))
        complete=x.dropna();complete=complete.loc[complete.std(axis=1)>0]
        if len(complete)<2 or x.shape[1]<2:self.skip(prefix+'_pca','Fewer than two complete variable features or samples');return
        # PCA/heatmaps are display calculations, not transformations supplied to inference.
        chosen=complete.var(axis=1).sort_values(ascending=False,kind='stable').head(2000).index
        z=complete.loc[chosen].to_numpy(copy=True).T;z-=z.mean(axis=0)
        u,s,_=np.linalg.svd(z,full_matrices=False);pcs=u[:,:2]*s[:2];var=s*s/(s*s).sum()
        pdata=pd.DataFrame({'sample':x.columns,'PC1':pcs[:,0],'PC2':pcs[:,1]})
        if metadata is not None and len(metadata):
            key=next((c for c in ['sample','specimen','SampleID'] if c in metadata),metadata.columns[0]);meta=metadata.rename(columns={key:'sample'}).copy();meta['sample']=meta['sample'].astype(str)
            if meta['sample'].duplicated().any():raise ValueError('Ambiguous sample metadata in report')
            pdata=pdata.merge(meta,on='sample',how='left',validate='one_to_one')
        fields=[c for c in pdata if c not in ['sample','PC1','PC2'] and pdata[c].nunique(dropna=True)>1][:4] or [None]
        fig,axes=plt.subplots(1,len(fields),figsize=(6*len(fields),5),layout='constrained',squeeze=False)
        for ax,field in zip(axes[0],fields):
            if field:
                cats=pdata[field].fillna('Missing').astype(str);unique=sorted(cats.unique())
                for j,c in enumerate(unique):
                    ok=cats==c;ax.scatter(pdata.loc[ok,'PC1'],pdata.loc[ok,'PC2'],s=65,label=c,color=PALETTE[j%len(PALETTE)],edgecolors='white',linewidths=.8)
                if len(unique)<=12:ax.legend(fontsize=7)
            else:ax.scatter(pdata.PC1,pdata.PC2,color=BLUE)
            ax.set_title(str(field or 'samples'));ax.set_xlabel(f'PC1 ({var[0]:.1%})');ax.set_ylabel(f'PC2 ({var[1]:.1%})')
        self.add(prefix+'_pca',fig,'Unscaled PCA on at most 2,000 variable, complete features among the rows read. Missing features are omitted from PCA; samples are not imputed. Separation is descriptive, not proof of a biological or batch effect.',pdata)
        corr=complete[selected].corr();fig,ax=canvas(min(len(corr),35),10);matrix(fig,ax,corr,corr.index,corr.columns,'Pearson correlation',1);ax.set_title(prefix+' sample correlation')
        self.add(prefix+'_correlation',fig,'Pearson correlation on complete variable features in the displayed matrix. Constant columns have undefined correlation.',corr.reset_index(names='sample'))
        top=complete.var(axis=1).sort_values(ascending=False,kind='stable').head(30).index;heat=complete.loc[top,selected];heat=heat.sub(heat.mean(axis=1),axis=0).div(heat.std(axis=1),axis=0)
        fig,ax=canvas(len(top),12);matrix(fig,ax,heat,heat.index,heat.columns,'Within-feature z score');ax.set_title(prefix+' variable-feature heatmap')
        self.add(prefix+'_feature_heatmap',fig,'Thirty highest-variance complete features among rows read; within-feature z scores are only for display and are never used as new assay measurements.',heat.reset_index(names='feature'))
        top=top[:4];fig,axes=plt.subplots(1,len(top),figsize=(5*len(top),4),layout='constrained',squeeze=False);points=[]
        group=next((c for c in ['condition','group','trait','genomic_state'] if c in pdata),None)
        for ax,gene in zip(axes[0],top):
            values=complete.loc[gene];labels=pdata[group].fillna('Missing').astype(str) if group else pd.Series('All samples',index=pdata.index)
            groups=sorted(labels.unique());loc={g:i for i,g in enumerate(groups)}
            for g in groups:
                ix=np.flatnonzero(labels==g);ax.scatter(np.full(len(ix),loc[g])+np.linspace(-.12,.12,len(ix)),values.iloc[ix],s=18,color=BLUE)
                points.extend(dict(feature=gene,sample=values.index[i],group=g,value=values.iloc[i]) for i in ix)
            ax.set_xticks(range(len(groups)),groups,rotation=30,ha='right');ax.set_title(gene);ax.set_ylabel(scale)
        self.add(prefix+'_feature_values',fig,'Individual measurements for four variable features, selected for visualization rather than significance. Points represent samples; repeated donor measurements require a paired statistical model.',pd.DataFrame(points))


def methods(b,frame,alpha=.05):
    if frame is None or frame.empty:b.skip('method_agreement','No method matrix');return
    names=[c[:-2] for c in frame if c.endswith('_q') and c[:-2]+'_effect' in frame]
    if len(names)<2:b.skip('method_agreement','Need at least two named effect/q pairs');return
    calls=pd.DataFrame({n:np.where(numbers(frame[n+'_q']).isna(),np.nan,np.where(numbers(frame[n+'_q'])<=alpha,np.sign(numbers(frame[n+'_effect'])),0)) for n in names})
    calls.index=frame.iloc[:,0].astype(str);order=calls.abs().sum(axis=1).sort_values(ascending=False,kind='stable').head(30).index
    fig,ax=canvas(len(order));matrix(fig,ax,calls.loc[order],order,names,f'Call direction at q <= {alpha:g}',1);ax.set_title('Agreement across methods')
    b.add('method_calls',fig,'Signed calls from supplied adjusted values. Gray denotes untestable; zero denotes no call. Agreement is descriptive, not an additional FDR procedure.',calls.reset_index(names='feature'))
    patterns=calls.value_counts(dropna=False).rename('features').reset_index();shown=patterns.head(25)
    fig,(ax,count)=plt.subplots(1,2,figsize=(11,max(4.8,2+.34*len(shown))),width_ratios=[max(3,len(names)),2],layout='constrained')
    matrix(fig,ax,shown[names],range(1,len(shown)+1),names,f'Call at q ≤ {alpha:g}: −1 / 0 / +1',1)
    ax.set_title('Method-call combinations');ax.set_ylabel('Pattern')
    count.barh(range(len(shown)),shown.features,color=BLUE,height=.55);count.set_ylim(len(shown)-.5,-.5);count.set_yticks([]);count.set_xlabel('Features');count.set_title('Exclusive counts');count.margins(x=.2)
    for i,k in enumerate(shown.features):count.text(k+.03*max(shown.features.max(),1),i,str(k),va='center',fontsize=9)
    b.add('method_intersections',fig,'Each row is an exclusive combination of signed calls at the declared adjusted-p threshold. Zero is no call, hatched gray is untestable; up to 25 most frequent patterns are shown. The TSV contains all combinations. This is descriptive robustness, not a new FDR procedure.',patterns)
    pairs=[]
    for i,a in enumerate(names):
        for bb in names[i+1:]:
            x=numbers(frame[a+'_effect']);y=numbers(frame[bb+'_effect']);ok=x.notna()&y.notna()
            pairs.append((a,bb,x,y,ok))
    fig,axes=plt.subplots(1,min(6,len(pairs)),figsize=(6*min(6,len(pairs)),5),layout='constrained',squeeze=False);data=[]
    for ax,(a,bb,x,y,ok) in zip(axes[0],pairs):
        ax.scatter(x[ok],y[ok],s=8,alpha=.5,color=BLUE);ax.axhline(0,color='#777',lw=.5);ax.axvline(0,color='#777',lw=.5);ax.set_xlabel(a+' effect');ax.set_ylabel(bb+' effect');ax.set_title(f'{int(ok.sum())} comparable features')
        data.extend(dict(feature=f,method_A=a,method_B=bb,effect_A=u,effect_B=v) for f,u,v in zip(frame.iloc[:,0][ok],x[ok],y[ok]))
    b.add('method_effects',fig,'First six method pairs; values retain each method\'s effect definition. Shared direction or ranking is not proof that estimands or numerical scales are identical.',pd.DataFrame(data))


def loci(b,counts,name='methylation',regions=None):
    if counts is None or counts.empty:b.skip(name,'No covered loci');return
    if not {'N','X'}<=set(counts):b.skip(name,'Requires actual modified X and total N counts');return
    t=counts.copy();t['N']=numbers(t.N);t['X']=numbers(t.X);t['fraction']=t.X/t.N.where(t.N>0)
    b.numeric(name+'_counts',t,['N','fraction'],'Covered-site depth and modification fraction')
    if 'chr' in t:b.categories(name+'_contigs',t.chr,'Covered sites by contig','Counts apply to rows read, not whole-genome coverage when the input exceeds the report limit.')
    if not {'chr','pos'}<=set(t):return
    c=str(t.chr.iloc[0]);start=int(t.pos.iloc[0]);end=start+2000
    if regions is not None and len(regions):
        cc=next((k for k in ['chr','chrom','seqnames'] if k in regions),None)
        if cc and 'start' in regions and 'end' in regions:c=str(regions[cc].iloc[0]);start=int(regions.start.iloc[0]);end=int(regions.end.iloc[0])
    sub=t[(t.chr.astype(str)==c)&numbers(t.pos).between(start,end)].copy();fig,(ax,depth)=plt.subplots(2,1,figsize=(11,6),sharex=True,layout='constrained')
    if len(sub):ax.scatter(sub.pos,sub.fraction,c=BLUE,s=20);depth.vlines(sub.pos,0,sub.N,color=BLUE)
    else:empty(ax,'No loaded covered sites overlap the displayed interval')
    ax.set_ylim(0,1);ax.set_ylabel('Modified / total');depth.set_ylabel('Depth');depth.set_xlabel(f'{c}: one-based position');ax.set_title('Modification locus evidence')
    b.add(name+'_locus',fig,f'Display interval {c}:{start}-{end}; absent sites remain gaps. The fraction uses actual X/N and does not resolve modification channels beyond the declared chemistry.',sub)


def predictions(b,frame):
    if frame is None or frame.empty or 'observed' not in frame:b.skip('predictions','No held-out predictions with observed outcomes');return
    if {'model','predicted'}<=set(frame):long=frame.copy()
    else:
        labels=[c for c in frame if c not in ['donor','sample','fold','observed','partition'] and pd.api.types.is_numeric_dtype(frame[c])]
        long=frame.melt(id_vars=[c for c in ['donor','observed'] if c in frame],value_vars=labels,var_name='model',value_name='predicted')
    models=list(long.model.drop_duplicates())[:6];fig,axes=plt.subplots(1,max(1,len(models)),figsize=(5*max(1,len(models)),5),layout='constrained',squeeze=False)
    for ax,n in zip(axes[0],models):
        s=long[long.model==n];x=numbers(s.observed);y=numbers(s.predicted);ok=x.notna()&y.notna();ax.scatter(x[ok],y[ok],s=20,color=BLUE,alpha=.6)
        if ok.any():low=min(x[ok].min(),y[ok].min());hi=max(x[ok].max(),y[ok].max());ax.plot([low,hi],[low,hi],ls='--',color='#555')
        ax.set_title(str(n));ax.set_xlabel('Observed');ax.set_ylabel('Predicted')
    b.add('heldout_predictions',fig,'Supplied held-out or out-of-fold predictions; the report does not train or select a model. The dashed line is identity, not a fitted calibration line.',long)
    long['residual']=numbers(long.predicted)-numbers(long.observed);fig,ax=canvas(10)
    for n in models:
        s=long[long.model==n];ax.scatter(s.predicted,s.residual,label=n,s=12,alpha=.5)
    ax.axhline(0,color='#555',ls='--');ax.set_xlabel('Predicted');ax.set_ylabel('Predicted - observed');ax.legend(fontsize=8);ax.set_title('Prediction residuals')
    b.add('prediction_residuals',fig,'Residual pattern by supplied model, without refitting. Units are those of the outcome.',long)


def record_table(b,path,prefix='summary'):
    """Show mixed-unit audit metrics as a table, never as an invalid common-axis bar chart."""
    if not Path(path).is_file():b.skip(prefix,'No summary receipt');return
    obj=json.loads(Path(path).read_text(encoding='utf-8'));rows=[]
    def walk(v,key=''):
        if isinstance(v,dict):
            for k,x in v.items():
                if k not in ['input_sha256','inputs_sha256','source_sha256','inputs','argv','config','resources','variants','path']:walk(x,key+'.'+k if key else k)
        elif isinstance(v,(str,int,float,bool)) or v is None:rows.append(dict(metric=key,value=str(v)))
    walk(obj);table=pd.DataFrame(rows)
    if table.empty:b.skip(prefix,'No scalar summary fields');return
    shown=table.head(24).copy();shown['value']=shown.value.str.slice(0,100)
    fig,ax=canvas(len(shown),13);ax.axis('off');grid=ax.table(cellText=shown.values,colLabels=['Recorded metric','Value'],loc='center',colWidths=[.45,.55]);grid.auto_set_font_size(False);grid.set_fontsize(8);grid.scale(1,1.45);ax.set_title('Recorded analysis metrics (original units)')
    b.add(prefix,fig,'Receipt values retain their individual units; they are not compared on a shared numeric axis. Up to 24 rows displayed; TSV includes all scalar entries.',table)


def intervals(b,path,prefix='regions'):
    t=b.load(path,header=None,comment='#')
    if t is None or t.empty or t.shape[1]<3:b.skip(prefix,'No interval rows');return
    t=t.iloc[:,:min(10,t.shape[1])];t=t.rename(columns={0:'chrom',1:'start',2:'end'});t['width']=numbers(t.end)-numbers(t.start)
    if (t.width<=0).any():raise ValueError('Invalid BED interval in report')
    b.numeric(prefix+'_width',t,['width'],'Genomic interval width (bp)');b.categories(prefix+'_contigs',t.chrom,'Intervals by contig')
    sub=t.head(30);fig,ax=canvas(len(sub),12)
    for i,row in enumerate(sub.itertuples()):ax.plot([row.start,row.end],[i,i],lw=3,color=BLUE)
    ax.set_yticks(range(len(sub)),[str(c) for c in sub.chrom]);ax.set_xlabel('Source contig coordinate (0-based, half-open)');ax.set_title('First 30 genomic intervals (separate rows)')
    b.add(prefix+'_intervals',fig,'Each row is a separate interval; coordinates on different contigs are not a shared genomic location. This is interval geometry, not evidence of enhancer–gene linkage.',sub)
