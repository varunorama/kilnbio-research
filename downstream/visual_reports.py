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
import numpy as np
import pandas as pd

BLUE = '#2166ac'
ORANGE = '#b35806'
MISSING = '#bdbdbd'
CMAP = LinearSegmentedColormap.from_list('evidence', [BLUE, '#fffaf2', ORANGE])
CMAP.set_bad(MISSING)


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
            fig.savefig(self.out/'figures'/f'{name}.png', dpi=160, bbox_inches='tight')
            fig.savefig(self.out/'figures'/f'{name}.svg', bbox_inches='tight', metadata={'Date': None})
            self.pdf.savefig(fig, bbox_inches='tight')
        finally:
            plt.close(fig)
        self.panels.append((name, caption))

    def finish(self, details):
        self.pdf.close()
        body = ['<!doctype html><html lang="en"><meta charset="utf-8">',
                '<meta name="viewport" content="width=device-width,initial-scale=1">',
                '<style>body{font:16px system-ui,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#182b3a;background:#fafbfc}section{background:white;padding:1.3rem;margin:1rem 0;border:1px solid #dce2e7;border-radius:8px}img{max-width:100%;height:auto}p{line-height:1.55}a{color:#2166ac}h1,h2{line-height:1.2}</style>',
                '<title>'+html.escape(self.title)+'</title><h1>'+html.escape(self.title)+'</h1>']
        body.extend('<p>'+html.escape(n)+'</p>' for n in self.notes)
        for label,target in self.links:
            body.append('<p><a href="'+html.escape(target,quote=True)+'">'+html.escape(label)+'</a></p>')
        body.append('<p><a href="visualizations.pdf">All figures (PDF)</a> | <a href="visualization.json">Report provenance</a></p>')
        for name, caption in self.panels:
            encoded = base64.b64encode((self.out/'figures'/f'{name}.png').read_bytes()).decode('ascii')
            body.append('<section><h2>'+html.escape(name.replace('_',' ').capitalize())+'</h2><img alt="'+html.escape(caption,quote=True)+'" src="data:image/png;base64,'+encoded+'"><p>'+html.escape(caption)+'</p><a href="figures/'+name+'.svg">Vector SVG</a> | <a href="plotted_data/'+name+'.tsv">Plotted data (TSV)</a></section>')
        body.append('</html>')
        (self.out/'report.html').write_text('\n'.join(body), encoding='utf-8')
        record = dict(schema_version=1, title=self.title, rendering='offline; no LLM or external requests',
                      matplotlib=matplotlib.__version__, numpy=np.__version__, pandas=pd.__version__,
                      display_p_floor=1e-16, details=details,
                      panels=[dict(name=n,caption=c) for n,c in self.panels],
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        record['artifact_sha256'] = {str(p.relative_to(self.out)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(self.out.rglob('*')) if p.is_file() and (p.parent.name in {'figures','plotted_data'} or p.name in {'report.html','visualizations.pdf'})}
        (self.out/'visualization.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n', encoding='utf-8')


def canvas(rows=12, width=10):
    return plt.subplots(figsize=(width, max(4, 1.7 + rows*.29)), layout='constrained')


def empty(ax, text):
    ax.axis('off'); ax.text(.5,.5,text,ha='center',va='center',wrap=True,transform=ax.transAxes)


def matrix(fig, ax, values, rows, columns, label, limit=None):
    values = np.asarray(values, float)
    if not values.size:
        empty(ax,'No eligible rows for this view.'); return
    finite = np.abs(values[np.isfinite(values)])
    limit = limit or max(1, float(finite.max()) if len(finite) else 1)
    im = ax.imshow(np.ma.masked_invalid(values), aspect='auto', cmap=CMAP, vmin=-limit, vmax=limit, interpolation='nearest')
    ax.set_yticks(range(len(rows)), rows, fontsize=9)
    ax.set_xticks(range(len(columns)), columns, rotation=35, ha='right', fontsize=9)
    fig.colorbar(im, ax=ax, label=label, shrink=.8)


def multiomic_report(out, result, manifest, summary):
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
    fig,ax=canvas(len(top)); matrix(fig,ax,signs[top]*strength(p[top]),result.gene.iloc[top],layers,'Signed -log10(layer p)')
    ax.set_title('Gene-by-assay evidence; gray = missing')
    report.add('evidence_heatmap',fig,'Direction is defined by the declared orientation. Per-layer p values include feature-to-gene aggregation; color is descriptive and does not identify individually significant layers.',long)

    ranked=result.loc[top].copy();fig,ax=canvas(len(top))
    valid=ranked.q_pc.notna()
    ax.scatter(strength(ranked.loc[valid,'q_pc']),np.flatnonzero(valid),s=25+20*ranked.loc[valid,'layers_measured'],c=np.where(ranked.loc[valid,'class']=='convergent',ORANGE,BLUE))
    ax.set_yticks(range(len(ranked)),ranked.gene,fontsize=9);ax.invert_yaxis()
    ax.axvline(-np.log10(summary['alpha']),color='#666',linestyle='--',label='BY threshold')
    ax.set_xlabel('-log10(joint BY q), display cap 16');ax.set_title('Joint evidence ranking');ax.legend(fontsize=8)
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
    ax.set_xlabel('Declared genes (all evidence strengths)');ax.set_title('Pairwise directional concordance');ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.13),ncol=4)
    report.add('concordance',fig,'Counts include all declared genes, without a significance filter. Same direction is not independent replication. Neutral denotes tied directional evidence; not-both includes unmeasured/untestable genes.',conc)

    missing=pd.DataFrame(np.isfinite(p).astype(int),columns=layers);patterns=missing.value_counts().rename('genes').reset_index().sort_values(['genes',*layers],ascending=[False,*([True]*len(layers))])
    shown=patterns.head(25);fig,ax=canvas(len(shown));labels=[' '.join(f'{n}:{"yes" if row[n] else "NA"}' for n in layers) for _,row in shown.iterrows()]
    ax.barh(labels,shown.genes,color=BLUE);ax.invert_yaxis();ax.set_xlabel('Genes');ax.set_title('Assay coverage patterns (largest 25)')
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
    report.finish(dict(top_n=30,min_layers=r,rank_ties='average',rank_eligibility=eligibility,alpha=summary['alpha']))


def enrichment_report(out, cfg):
    out=Path(out);table=pd.read_csv(out/'all_pathways.tsv',sep='\t');top=table.sort_values(['padj','pathway']).head(25).copy()
    report=Report(out,'KilnBio gene-set enrichment report',[
        'Rank-based fgsea enrichment from the complete supplied signed gene ranking and the declared local GMT.',
        f"{len(table)} pathway rows; {int((table.padj<cfg['fdr']).sum())} below FDR {cfg['fdr']}. Labels are gene-set names, not proven activation states.",
        'The enrichment-curves.pdf file contains running enrichment scores computed by fgsea. Leading-edge membership and set overlap are descriptive.'])
    fig,ax=canvas(len(top));valid=top.NES.notna()&top.padj.notna()
    scatter=ax.scatter(top.loc[valid,'NES'],np.flatnonzero(valid),s=25+strength(top.loc[valid,'padj'])*20,c=top.loc[valid,'NES'],cmap=CMAP)
    ax.set_yticks(range(len(top)),top.pathway,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#666',lw=.7);ax.set_xlabel('Normalized enrichment score (NES)');ax.set_title('Top pathways by adjusted p value')
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
