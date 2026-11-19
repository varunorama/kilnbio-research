"""Offline scientific reports; figures never change tests or selection rules."""
import base64
import hashlib
import html
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from plot_design import BLUE, ORANGE, MISSING, CMAP, INK, MUTED, TEAL, finish_figure, panel_guide

def strength(p):
    """Cap display only; original p/q values remain in tables."""
    return -np.log10(np.clip(np.asarray(p, float), 1e-16, 1))


class Report:
    def __init__(self, out, title, notes):
        self.out = Path(out)
        (self.out/'figures').mkdir()
        (self.out/'plotted_data').mkdir()
        self.title, self.notes, self.panels = title, notes, []
        self.links = []
        self.pdf = PdfPages(self.out/'visualizations.pdf', metadata={
            'Title': title, 'CreationDate': None, 'ModDate': None})

    def table(self, name, frame):
        frame.to_csv(self.out/'plotted_data'/f'{name}.tsv', sep='\t', index=False, na_rep='NA')

    def add(self, name, fig, caption, frame):
        self.table(name, frame)
        try:
            finish_figure(fig)
            fig.savefig(self.out/'figures'/f'{name}.png', dpi=200, bbox_inches='tight')
            fig.savefig(self.out/'figures'/f'{name}.svg', bbox_inches='tight', metadata={'Date': None})
            self.pdf.savefig(fig, bbox_inches='tight')
        finally:
            plt.close(fig)
        self.panels.append((name, caption))

    def finish(self, details):
        self.pdf.close()
        body = ['<!doctype html><html lang="en"><meta charset="utf-8">',
                '<meta name="viewport" content="width=device-width,initial-scale=1">',
                '<style>*{box-sizing:border-box}html{scroll-behavior:smooth}body,main,header,section{min-width:0}p,h1,h2,h3,a{overflow-wrap:anywhere}body{font:16px/1.65 system-ui,sans-serif;margin:0;color:#203447;background:#eef2f6}main{max-width:1220px;margin:auto;padding:48px 32px}header,section{background:white;border:1px solid #dfe6ed;border-radius:16px;padding:32px;margin-bottom:28px;box-shadow:0 4px 18px #20344705}header{border-top:5px solid #ce713b}.eyebrow{font-size:12px;font-weight:750;letter-spacing:.16em;color:#286fa3;text-transform:uppercase}h1{font-size:32px;line-height:1.2;max-width:900px;margin:12px 0 24px}h2{font-size:24px;line-height:1.3;margin:10px 0 18px}.guide{display:grid;grid-template-columns:1fr 1fr;gap:28px;background:#f4f7fa;border-radius:10px;padding:24px;margin:20px 0}.guide h3{font-size:12px;letter-spacing:.08em;text-transform:uppercase;margin:0 0 8px;color:#637587}.guide p{margin:0}img{display:block;max-width:100%;height:auto;margin:24px auto}a{color:#21638d;text-underline-offset:3px}.downloads{display:flex;flex-wrap:wrap;gap:12px;margin-top:24px}.downloads a{border:1px solid #cbd6e1;border-radius:7px;padding:7px 14px;font-size:14px;text-decoration:none}nav{columns:2;column-gap:36px;margin:24px 0}nav a{display:block;padding:6px 0;break-inside:avoid}details{font-size:14px;color:#53687a}summary{cursor:pointer;font-weight:600}footer{color:#637587;font-size:13px;padding:0 8px 24px}@media(max-width:720px){main{padding:20px 12px}header,section{padding:20px}.guide{grid-template-columns:1fr}nav{columns:1}h1{font-size:26px}}@media print{body{background:white}main{padding:0}section{box-shadow:none;break-inside:avoid}nav{display:none}.downloads{display:none}}</style>',
                '<title>'+html.escape(self.title)+'</title><body><main><header><div class="eyebrow">KilnBio / Reproducible analysis</div><h1>'+html.escape(self.title)+'</h1>']
        body.append('<p>Inspect the measurements, compare the evidence, and trace what supports each follow-up decision. Every figure explains why it is shown and provides the exact plotted values.</p><details><summary>Analysis scope, thresholds and display rules</summary>')
        body.extend('<p>'+html.escape(n)+'</p>' for n in self.notes)
        body.append('</details>')
        for label,target in self.links:
            body.append('<p><a href="'+html.escape(target,quote=True)+'">'+html.escape(label)+'</a></p>')
        body.append('<div class="downloads"><a href="visualizations.pdf">Download all figures · PDF</a><a href="visualization.json">Report provenance</a></div><nav aria-label="Report panels">')
        for i,(name, _) in enumerate(self.panels,1):
            body.append('<a href="#'+html.escape(name,quote=True)+'">'+f'{i:02d} · '+html.escape(panel_guide(name)['title'])+'</a>')
        body.append('</nav></header>')
        for i,(name, caption) in enumerate(self.panels,1):
            guide=panel_guide(name)
            encoded = base64.b64encode((self.out/'figures'/f'{name}.png').read_bytes()).decode('ascii')
            body.append('<section id="'+html.escape(name,quote=True)+'"><div class="eyebrow">Figure '+f'{i:02d}'+'</div><h2>'+html.escape(guide['title'])+'</h2><div class="guide"><div><h3>Why this matters</h3><p>'+html.escape(guide['why'])+'</p></div><div><h3>How to read it</h3><p>'+html.escape(guide['read'])+'</p></div></div><img loading="lazy" alt="'+html.escape(caption,quote=True)+'" src="data:image/png;base64,'+encoded+'"><details><summary>Method, display rules and scope</summary><p>'+html.escape(caption)+'</p></details><div class="downloads"><a href="figures/'+name+'.svg">Vector SVG</a><a href="figures/'+name+'.png">PNG image</a><a href="plotted_data/'+name+'.tsv">Plotted values · TSV</a></div></section>')
        body.append('<footer>Generated locally by KilnBio. Figures and explanatory text use deterministic code; no remote service or language model is called.</footer></main></body></html>')
        (self.out/'report.html').write_text('\n'.join(body), encoding='utf-8')
        record = dict(schema_version=1, title=self.title, rendering='offline; no LLM or external requests',
                      matplotlib=matplotlib.__version__, numpy=np.__version__, pandas=pd.__version__,
                      display_p_floor=1e-16, details=details,
                      design_version=2, panels=[dict(name=n,caption=c,**panel_guide(n)) for n,c in self.panels],
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        record['design_source_sha256']=hashlib.sha256(Path(__file__).with_name('plot_design.py').read_bytes()).hexdigest()
        record['artifact_sha256'] = {str(p.relative_to(self.out)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(self.out.rglob('*')) if p.is_file() and (p.parent.name in {'figures','plotted_data'} or p.name in {'report.html','visualizations.pdf'})}
        (self.out/'visualization.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n', encoding='utf-8')


def canvas(rows=12, width=10):
    return plt.subplots(figsize=(width, max(4.6, 2.1 + rows*.31)), layout='constrained')


def empty(ax, text):
    ax.axis('off'); ax.text(.5,.5,text,ha='center',va='center',wrap=True,transform=ax.transAxes)


def matrix(fig, ax, values, rows, columns, label, limit=None, colorbar=True):
    values = np.asarray(values, float)
    if not values.size:
        empty(ax,'No eligible rows for this view.'); return
    finite = np.abs(values[np.isfinite(values)])
    limit = limit or max(1, float(finite.max()) if len(finite) else 1)
    im = ax.imshow(np.ma.masked_invalid(values), aspect='auto', cmap=CMAP, vmin=-limit, vmax=limit, interpolation='nearest')
    ax.set_yticks(range(len(rows)), rows, fontsize=10 if len(rows)<=30 else 8)
    ax.set_xticks(range(len(columns)), columns, rotation=0 if len(columns)<=5 else 45, ha='center' if len(columns)<=5 else 'right', fontsize=10 if len(columns)<=15 else 7)
    ax.set_xticks(np.arange(-.5,len(columns),1),minor=True);ax.set_yticks(np.arange(-.5,len(rows),1),minor=True)
    ax.grid(which='minor',color='white',lw=1.4);ax.tick_params(which='minor',bottom=False,left=False)
    for spine in ax.spines.values():spine.set_visible(False)
    if values.size<=600:
        for y,x in zip(*np.where(~np.isfinite(values))):
            ax.add_patch(Rectangle((x-.5,y-.5),1,1,facecolor='none',edgecolor='#A4B2BF',hatch='///',linewidth=0))
    if values.size<=80:
        discrete=np.isin(values[np.isfinite(values)],[-1,0,1]).all()
        for y,x in zip(*np.where(np.isfinite(values))):
            v=values[y,x];cell_label=f'{v:.0f}' if discrete else f'{v:.2f}' if limit==1 else f'{v:.1f}'
            ax.text(x,y,cell_label,ha='center',va='center',fontsize=10,color='white' if abs(v)>.65*limit else INK)
    if colorbar:fig.colorbar(im, ax=ax, label=label, shrink=.75, pad=.03)
    return im


def evidence_board(values, rows, columns, q, measured, alpha, title, subtitle, scale):
    """Aligned evidence, coverage and joint-q views of the same displayed rows."""
    n=len(rows);fig=plt.figure(figsize=(13,max(5.6,2.8+.32*n)),layout='constrained')
    gs=fig.add_gridspec(1,3,width_ratios=[max(3,len(columns)),1.15,2.2],wspace=.14)
    ax=fig.add_subplot(gs[0]);cov=fig.add_subplot(gs[1]);rank=fig.add_subplot(gs[2])
    im=matrix(fig,ax,values,rows,columns,scale,colorbar=False)
    ax.set_title('A  /  Assay evidence',fontsize=12)
    if im is not None:fig.colorbar(im,ax=ax,orientation='horizontal',label=scale,shrink=.85,pad=.035,fraction=.035,aspect=32)
    yy=np.arange(n);cov.barh(yy,measured,height=.56,color=TEAL)
    cov.set_ylim(n-.5,-.5);cov.set_yticks([]);cov.set_xlim(0,len(columns)+.6)
    cov.set_xticks(range(0,len(columns)+1,max(1,len(columns)//4)));cov.set_xlabel('Measured layers')
    for y,k in enumerate(measured):cov.text(k+.12,y,str(int(k)),va='center',fontsize=9)
    cov.set_title('B  /  Coverage',fontsize=12)
    cov.text(0,-.13,'Hatched cell:\nno testable evidence',transform=cov.transAxes,fontsize=9,color=MUTED,va='top')
    q=np.asarray(q,float);valid=np.isfinite(q);xx=strength(q[valid]);colors=np.where(q[valid]<=alpha,ORANGE,BLUE)
    rank.hlines(yy[valid],0,xx,color='#CDDAE5',lw=2.5);rank.scatter(xx,yy[valid],c=colors,s=45,zorder=3)
    rank.axvline(-np.log10(alpha),ls='--',color=ORANGE,lw=1)
    rank.set_ylim(n-.5,-.5);rank.set_yticks([]);rank.set_xlim(-.12,max(2,-np.log10(alpha)+.7,float(max(xx,default=0))+.5))
    rank.set_xlabel('−log₁₀(joint BY q)');rank.set_title('C  /  Joint evidence',fontsize=12)
    rank.text(.98,.02,f'Dashed line: q = {alpha:g}',ha='right',va='bottom',transform=rank.transAxes,fontsize=9,color=MUTED)
    for y in yy[~valid]:rank.text(.05,y,'Untestable',va='center',color=MUTED,fontsize=9)
    fig.suptitle(title+'\n'+subtitle,fontsize=16,fontweight='bold',ha='left',x=.02)
    return fig


def multiomic_report(out, result, manifest, summary, circos_config=None):
    # Import lazily to avoid changing the statistical module's public interface.
    from multiomic_meta import directional_pc
    layers = list(manifest.layer)
    result = result.sort_values('gene').reset_index(drop=True)
    p = result[['p_'+n for n in layers]].to_numpy(float)
    up = result[['p_up_'+n for n in layers]].to_numpy(float)
    down = result[['p_down_'+n for n in layers]].to_numpy(float)
    signs = np.sign(down-up)
    top = result.sort_values(['p_pc','gene'], na_position='last').head(30).index
    notes = [f"{len(result)} declared genes; {summary['testable_genes']} testable; {summary['convergent_genes']} selected at BY q <= {summary['alpha']}. At least {summary['min_layers']} layers required.",
             'Plots describe the supplied tests. Color represents prespecified orientation-adjusted evidence direction, not a common effect size or pathway activation. Gray means untestable/missing; pale means weak evidence.',
             'Top panels display at most 30 genes, ordered by joint raw p then identifier. All genes remain in the statistical family. Display scores are capped at 16; TSVs retain original values.',
             'Layer orientations: '+', '.join(f'{r.layer}={r.orientation}' for r in manifest.itertuples())]
    report = Report(out, 'KilnBio multi-omic evidence report', notes)
    long = pd.DataFrame([dict(gene=result.gene[i],layer=n,p=p[i,j],direction=signs[i,j],display_score=signs[i,j]*strength([p[i,j]])[0])
                         for i in top for j,n in enumerate(layers)])
    long=long.merge(result[['gene','layers_measured','q_pc']],on='gene',how='left',validate='many_to_one')
    fig=evidence_board(signs[top]*strength(p[top]),result.gene.iloc[top],layers,result.q_pc.iloc[top],result.layers_measured.iloc[top],summary['alpha'],
                       'Connect the evidence before choosing a target',
                       f"{len(result):,} declared genes  ·  {summary['testable_genes']:,} testable  ·  {summary['convergent_genes']:,} selected  ·  {len(top)} displayed",
                       'Orientation-adjusted −log₁₀(layer p)')
    report.add('evidence_heatmap',fig,'Rows are aligned across evidence, measured-layer counts and joint BY q. Direction is defined by the declared orientation. Per-layer p values include feature-to-gene aggregation. Hatched gray cells are untestable; pale cells have weak evidence. Colors retain the same direction scale across all displayed assays.',long)

    ranked=result.loc[top].copy();fig,ax=canvas(len(top))
    valid=ranked.q_pc.notna()
    ax.hlines(np.flatnonzero(valid),0,strength(ranked.loc[valid,'q_pc']),color='#D7E1EA',lw=2)
    ax.scatter(strength(ranked.loc[valid,'q_pc']),np.flatnonzero(valid),s=25+20*ranked.loc[valid,'layers_measured'],c=np.where(ranked.loc[valid,'class']=='convergent',ORANGE,BLUE),edgecolors='white',linewidths=.8,zorder=3)
    ax.set_yticks(range(len(ranked)),ranked.gene,fontsize=9);ax.invert_yaxis()
    ax.axvline(-np.log10(summary['alpha']),color='#666',linestyle='--',label='BY threshold')
    ax.set_xlabel('−log₁₀(joint BY q), display cap 16');ax.set_title('Joint evidence ranking')
    for k in sorted(ranked.layers_measured.unique()):ax.scatter([],[],s=25+20*k,color=MUTED,label=f'{k} measured layers')
    ax.legend(fontsize=9,loc='upper left',bbox_to_anchor=(1.02,1))
    report.add('combined_ranking',fig,'Orange marks selected genes; point area increases with measured layer count. Untestable rows have no point. Ranking uses raw joint p; the horizontal scale uses BY q.',ranked)

    concordance=[]
    for a in range(len(layers)):
        for b in range(a+1,len(layers)):
            both=np.isfinite(p[:,a])&np.isfinite(p[:,b]);nonzero=both&(signs[:,a]!=0)&(signs[:,b]!=0)
            same=int((nonzero&(signs[:,a]==signs[:,b])).sum());opp=int(nonzero.sum())-same
            concordance.append(dict(pair=layers[a]+' / '+layers[b],both_measured=int(both.sum()),same=same,opposite=opp,neutral=int(both.sum()-nonzero.sum()),not_both=int(len(result)-both.sum())))
    conc=pd.DataFrame(concordance);fig,ax=canvas(len(conc));left=np.zeros(len(conc))
    for key,color in [('same',BLUE),('opposite',ORANGE),('neutral','#e5d4a2'),('not_both',MISSING)]:
        ax.barh(conc.pair,conc[key],left=left,label=key.replace('_',' '),color=color);left+=conc[key]
    ax.set_xlabel('Declared genes (all evidence strengths)');ax.set_title('Pairwise directional concordance');ax.legend(fontsize=9,loc='upper left',bbox_to_anchor=(1.02,1))
    report.add('concordance',fig,'Counts include all declared genes, without a significance filter. Same direction is not independent replication. Neutral denotes tied directional evidence; not-both includes unmeasured/untestable genes.',conc)

    missing=pd.DataFrame(np.isfinite(p).astype(int),columns=layers);patterns=missing.value_counts().rename('genes').reset_index().sort_values(['genes',*layers],ascending=[False,*([True]*len(layers))])
    shown=patterns.head(25);fig,(ax,counts)=plt.subplots(1,2,figsize=(11,max(4.8,2+.34*len(shown))),width_ratios=[len(layers),2],layout='constrained')
    for i,(_,row) in enumerate(shown.iterrows()):
        present=np.flatnonzero([row[n] for n in layers]);ax.scatter(range(len(layers)),[i]*len(layers),s=55,color=MISSING)
        if len(present):ax.plot([present.min(),present.max()],[i,i],color=BLUE,lw=2);ax.scatter(present,[i]*len(present),s=65,color=BLUE,zorder=3)
    ax.set_xticks(range(len(layers)),layers);ax.set_yticks([]);ax.set_ylim(len(shown)-.5,-.5);ax.set_title('Measured assay combinations')
    counts.barh(range(len(shown)),shown.genes,height=.55,color=TEAL);counts.set_ylim(len(shown)-.5,-.5);counts.set_yticks([]);counts.set_xlabel('Genes');counts.set_title('Exclusive counts');counts.margins(x=.2)
    for i,k in enumerate(shown.genes):counts.text(k+.03*max(shown.genes.max(),1),i,str(k),va='center',fontsize=10)
    report.add('measurement_patterns',fig,'Exclusive measurement patterns preserve the distinction between absence of evidence and evidence of no association. The TSV contains every observed pattern.',patterns)

    rows=[];eligibility=[];r=summary['min_layers'];full=result.p_pc.to_numpy(float)
    for j,name in enumerate(layers):
        if len(layers)-1<r:
            eligibility.append(dict(removed=name,status='not_estimable_fewer_layers_than_r',common_genes=0,lost_testability=int(np.isfinite(full).sum())))
            continue
        reduced,_=directional_pc(np.delete(up,j,axis=1),np.delete(down,j,axis=1),r)
        common=np.isfinite(full)&np.isfinite(reduced)
        fullrank=pd.Series(full[common]).rank(method='average').to_numpy();reducedrank=pd.Series(reduced[common]).rank(method='average').to_numpy()
        eligibility.append(dict(removed=name,status='computed',common_genes=int(common.sum()),lost_testability=int((np.isfinite(full)&~np.isfinite(reduced)).sum())))
        for i,rf,rr in zip(np.flatnonzero(common),fullrank,reducedrank):
            rows.append(dict(gene=result.gene[i],removed=name,p_full=full[i],p_without=reduced[i],rank_full_common=rf,rank_without_common=rr,rank_gain_with_assay=rr-rf))
    changes=pd.DataFrame(rows,columns=['gene','removed','p_full','p_without','rank_full_common','rank_without_common','rank_gain_with_assay'])
    report.table('rank_eligibility',pd.DataFrame(eligibility));fig,ax=canvas(18)
    if changes.empty:empty(ax,'Leave-one-assay-out comparison is not estimable\nwith the required layer count held fixed.\nSee rank_eligibility.tsv.')
    else:
        scores=changes.groupby('gene').rank_gain_with_assay.apply(lambda x:x.abs().max()).sort_values(ascending=False,kind='stable').head(20)
        mat=changes.pivot(index='gene',columns='removed',values='rank_gain_with_assay').reindex(index=scores.index,columns=layers)
        matrix(fig,ax,mat,mat.index,mat.columns,'Rank without assay - rank with assay')
        ax.set_title('Rank sensitivity on common testable genes')
    report.add('rank_sensitivity',fig,'Each column removes one assay and holds r fixed. Average ranks handle ties within the same common testable gene set. Positive values indicate a better rank with that assay. This is sensitivity of supplied evidence, not measured predictive gain; eligibility losses are recorded separately.',changes)
    details=dict(top_n=30,min_layers=r,rank_ties='average',rank_eligibility=eligibility,alpha=summary['alpha'])
    if circos_config is not None:
        from circos_report import add_circos
        details['circos']=add_circos(report,circos_config,allowed_layers=layers)
    report.finish(details)


def enrichment_report(out, cfg):
    out=Path(out);table=pd.read_csv(out/'all_pathways.tsv',sep='\t');top=table.sort_values(['padj','pathway']).head(25).copy()
    report=Report(out,'KilnBio gene-set enrichment report',[
        'Rank-based fgsea enrichment from the complete supplied signed gene ranking and the declared local GMT.',
        f"{len(table)} pathway rows; {int((table.padj<cfg['fdr']).sum())} below FDR {cfg['fdr']}. Labels are gene-set names, not proven activation states.",
        'The enrichment-curves.pdf file contains running enrichment scores computed by fgsea. Leading-edge membership and set overlap are descriptive.'])
    fig,ax=canvas(len(top));valid=top.NES.notna()&top.padj.notna()
    nes_limit=max(1,float(top.loc[valid,'NES'].abs().max())) if valid.any() else 1
    scatter=ax.scatter(top.loc[valid,'NES'],np.flatnonzero(valid),s=25+strength(top.loc[valid,'padj'])*20,c=top.loc[valid,'NES'],cmap=CMAP,vmin=-nes_limit,vmax=nes_limit,edgecolors='white',linewidths=.6,zorder=3)
    ax.set_yticks(range(len(top)),top.pathway,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#666',lw=.7);ax.set_xlabel('Normalized enrichment score (NES)');ax.set_title('Top pathways by adjusted p value')
    for exponent in [1,2,4]:ax.scatter([],[],s=25+20*exponent,c=MUTED,label=f'q = {10.**(-exponent):g}')
    ax.legend(title='Adjusted evidence · point area',loc='upper left',bbox_to_anchor=(1.02,1))
    report.add('pathway_ranking',fig,'Point area encodes -log10(FDR q), capped at 16; position and color encode NES. Negative NES means enrichment toward the negative end of the supplied ranking.',top)
    edge={str(row.pathway):set(str(row.leadingEdge).split(';'))-{'','nan'} for row in top.itertuples()}
    freq={gene:sum(gene in s for s in edge.values()) for gene in set().union(*edge.values())} if edge else {}
    genes=sorted(freq,key=lambda g:(-freq[g],g))[:35]
    mat=np.array([[int(g in edge[n]) for g in genes] for n in edge])
    fig,ax=canvas(len(edge),12)
    if len(genes):
        ax.imshow(mat,aspect='auto',cmap='Blues',vmin=0,vmax=1,interpolation='nearest');ax.set_yticks(range(len(edge)),list(edge),fontsize=8);ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=8)
    else:empty(ax,'No leading-edge genes returned.')
    ax.set_title('Leading-edge membership: top 35 shared genes')
    long=pd.DataFrame([dict(pathway=n,gene=g,leading_edge=1) for n,s in edge.items() for g in sorted(s)],columns=['pathway','gene','leading_edge'])
    report.add('leading_edge',fig,'Blue indicates leading-edge membership; white is absence from that leading edge, not a negative gene effect. The TSV includes all leading-edge genes in the displayed pathways.',long)
    report.links.append(('Native fgsea curves (vector PDF)','enrichment-curves.pdf'))
    for path in sorted(out.glob('enrichment_curve_*.png')):
        fig,ax=plt.subplots(figsize=(9,6),layout='constrained');ax.imshow(plt.imread(path));ax.axis('off')
        report.add(path.stem,fig,'Native fgsea plotEnrichment output. The TSV records the full ranked universe and set membership; the separate enrichment-curves.pdf retains native vector curves.',pd.read_csv(path.with_suffix('.tsv'),sep='\t'))
        shutil.copyfile(path,out/'figures'/path.name)
        shutil.copyfile(path.with_suffix('.svg'),out/'figures'/path.with_suffix('.svg').name)
    report.finish(dict(top_n=25,leading_edge_max_genes=35,alpha=cfg['fdr'],gmt=Path(cfg['gmt']).name))
