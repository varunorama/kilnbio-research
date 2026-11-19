"""Compare full pathway test families across assays; generate an offline report.

The optional joint test is nondirectional Bonferroni partial conjunction followed
by BY. Missing tests are assigned p=1 within the fixed declared assay family.
It combines valid pathway p-values, not gene ranks, q-values or effect sizes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile

import numpy as np
import pandas as pd
from multiomic_meta import adjust
from visual_reports import Report, canvas, matrix, empty, strength, BLUE, ORANGE


def read_gmt(path):
    sets={}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        fields=line.split('\t')
        if len(fields)<3 or not fields[0] or fields[0] in sets:
            raise ValueError('GMT requires unique nonempty pathway names, description and gene IDs')
        genes=set(fields[2:])
        if '' in genes:raise ValueError('Empty GMT gene identifier')
        sets[fields[0]]=genes
    if not sets:raise ValueError('Empty GMT')
    return sets


def partial_conjunction(p, r):
    """Fixed m assays, untestable p=1: dependence-safe for valid raw p-values."""
    p=np.asarray(p,float)
    if p.ndim!=2 or not 2<=r<=p.shape[1]:raise ValueError('Invalid partial-conjunction r')
    valid=np.isfinite(p)
    if ((p[valid]<=0)|(p[valid]>1)).any():raise ValueError('Raw p must be in (0,1]')
    combined=np.minimum(1,(p.shape[1]-r+1)*np.sort(np.where(valid,p,1),axis=1)[:,r-1])
    combined[valid.sum(axis=1)<r]=np.nan
    q=adjust(np.where(np.isfinite(combined),combined,1))
    q[~np.isfinite(combined)]=np.nan
    return combined,q


def from_enrichment(directory, gmt, cfg):
    """Read a completed local enrich run without requiring hand-edited tables."""
    directory=Path(directory)
    source=directory/'config.resolved.json';status=directory/'status.json'
    if json.loads(status.read_text()).get('complete') is not True:
        raise ValueError('Enrichment run is not complete')
    meta=json.loads(source.read_text());ranks=Path(meta['ranks']);original_gmt=Path(meta['gmt'])
    if hashlib.sha256(original_gmt.read_bytes()).digest()!=hashlib.sha256(gmt.read_bytes()).digest():
        raise ValueError('All enrichment runs must use the identical frozen GMT')
    if meta['species']!=cfg['species'] or meta['gene_id_type']!=cfg['gene_id_type']:
        raise ValueError('Species/identifier mismatch across enrichment runs')
    table_path=directory/'all_pathways.tsv'
    table=pd.read_csv(table_path,sep='\t',dtype={'pathway':str,'leadingEdge':str})
    sets=read_gmt(gmt)
    if table.pathway.duplicated().any() or not set(table.pathway)<=set(sets):raise ValueError('Unknown or duplicate enrichment pathways')
    table=table.set_index('pathway').reindex(sorted(sets)).reset_index()
    table=table[['pathway','pval','NES','leadingEdge']].rename(columns={'pval':'p','NES':'score','leadingEdge':'contributing_genes'})
    background=pd.read_csv(ranks,sep='\t',dtype={'gene_id':str})[['gene_id']]
    return table,background,[source,status,table_path,ranks,original_gmt]


def run(config, output):
    config=Path(config).resolve();out=Path(output).resolve()
    if out.exists():raise ValueError('Choose a new output directory')
    cfg=json.loads(config.read_text(encoding='utf-8-sig'))
    allowed={'gmt','collection','species','gene_id_type','annotation_reviewed','contrast','layers','min_assays','alpha','joint_test'}
    if set(cfg)-allowed:raise ValueError('Unknown pathway comparison key')
    if cfg.get('annotation_reviewed') is not True or not all(isinstance(cfg.get(k),str) and cfg[k].strip() for k in ['collection','species','gene_id_type','contrast']):
        raise ValueError('Declare collection/version, species, gene ID type and contrast; review annotations')
    layers=cfg.get('layers',[]);names=[x.get('layer','') for x in layers]
    if len(layers)<2 or len(set(names))!=len(names) or not all(re.fullmatch('[A-Za-z][A-Za-z0-9_]*',n) for n in names):
        raise ValueError('At least two unique simple layer names required')
    r=cfg.get('min_assays',2);alpha=cfg.get('alpha',.05);joint=cfg.get('joint_test',False)
    if type(r) is not int or not 2<=r<=len(layers) or type(alpha) not in (int,float) or not 0<alpha<1 or type(joint) is not bool:
        raise ValueError('Invalid min_assays, alpha or joint_test')
    gmt=(config.parent/cfg['gmt']).resolve();sets=read_gmt(gmt);pathways=sorted(sets)
    paths=[config,gmt];evidence=[];reviews=[]
    for layer in layers:
        if set(layer)=={'layer','enrichment_dir'}:
            tab,background,extra=from_enrichment((config.parent/layer['enrichment_dir']).resolve(),gmt,cfg)
            paths.extend(extra)
            layer=dict(layer,method='fgseaMultilevel from completed KilnBio enrich run',score_kind='signed')
        else:
            if set(layer)!={'layer','results','background','method','score_kind','complete_family_reviewed'} or layer['complete_family_reviewed'] is not True:
                raise ValueError('Each layer needs enrichment_dir OR results/background/method/score_kind/complete_family_reviewed=true')
            result_path=(config.parent/layer['results']).resolve();background_path=(config.parent/layer['background']).resolve();paths.extend([result_path,background_path])
            background=pd.read_csv(background_path,sep='\t',dtype=str,keep_default_na=False)
            tab=pd.read_csv(result_path,sep='\t',dtype={'pathway':str,'contributing_genes':str})
        if layer['score_kind'] not in ('signed','unsigned') or not isinstance(layer['method'],str) or not layer['method'].strip():
            raise ValueError('Declare upstream method and signed/unsigned score_kind')
        if list(background)!=['gene_id'] or background.gene_id.isna().any() or background.gene_id.duplicated().any() or not background.gene_id.str.strip().astype(bool).all():
            raise ValueError('Background needs unique nonempty gene_id rows')
        universe=set(background.gene_id)
        if not universe:raise ValueError('Empty assay background')
        if set(tab)!={'pathway','p','score','contributing_genes'} or tab.pathway.isna().any() or tab.pathway.duplicated().any() or set(tab.pathway)!=set(pathways):
            raise ValueError('Results must contain every GMT pathway once: pathway,p,score,contributing_genes. Use NA for untestable rows; significant-only tables are rejected.')
        tab=tab.set_index('pathway').reindex(pathways)
        for col in ['p','score']:tab[col]=pd.to_numeric(tab[col],errors='raise')
        if np.isinf(tab[['p','score']]).any().any() or ((tab.p.dropna()<=0)|(tab.p.dropna()>1)).any():raise ValueError('Finite scores and raw p in (0,1], or NA, required')
        if tab.loc[tab.p.notna(),'score'].isna().any():raise ValueError('Tested rows need a descriptive score')
        if layer['score_kind']=='signed' and ((tab.score==0)&(tab.p<1)).any():raise ValueError('Signed p<1 rows require a nonzero descriptive score')
        if layer['score_kind']=='unsigned' and (tab.score.dropna()<0).any():raise ValueError('Unsigned scores cannot be negative')
        q=adjust(tab.p.fillna(1));q[tab.p.isna()]=np.nan
        for i,(pathway,row) in enumerate(tab.iterrows()):
            overlap=sets[pathway]&universe
            genes=set(str(row.contributing_genes).split(';')) if pd.notna(row.contributing_genes) else set()
            genes.discard('')
            if not genes<=overlap:raise ValueError('Contributing genes must belong to the pathway and measured background')
            if pd.notna(row.p) and not overlap:raise ValueError('A tested pathway must overlap the measured background')
            if pd.isna(row.p) and genes:raise ValueError('Untestable rows cannot have contributing genes')
            evidence.append(dict(pathway=pathway,layer=layer['layer'],p=row.p,q_within=q[i],score=row.score,score_kind=layer['score_kind'],genes_in_background=len(overlap),contributing_genes=';'.join(sorted(genes))))
        reviews.append(dict(layer=layer['layer'],method=layer['method'],score_kind=layer['score_kind'],background_genes=len(universe),tested_pathways=int(tab.p.notna().sum())))
    evidence=pd.DataFrame(evidence)
    pmat=evidence.pivot(index='pathway',columns='layer',values='p').reindex(index=pathways,columns=names)
    qmat=evidence.pivot(index='pathway',columns='layer',values='q_within').reindex(index=pathways,columns=names)
    result=pd.DataFrame({'pathway':pathways,'assays_tested':pmat.notna().sum(axis=1).to_numpy(),'assays_passing_within_BY':(qmat<=alpha).sum(axis=1).to_numpy()})
    if joint:result['p_pc'],result['q_pc']=partial_conjunction(pmat.to_numpy(),r)
    result['min_layer_p']=pmat.min(axis=1).to_numpy() # Display sorting only, not an adjusted joint test.
    order=['p_pc','pathway'] if joint else ['min_layer_p','pathway']
    result=result.sort_values(order,na_position='last')
    summary=dict(collection=cfg['collection'],contrast=cfg['contrast'],pathway_family=len(pathways),assays=reviews,min_assays=r,alpha=alpha,joint_test=joint,
                 joint_method='nondirectional Bonferroni partial conjunction, fixed declared m with missing p=1; BY across full GMT' if joint else 'descriptive comparison; no joint selection',
                 selected_pathways=int((result.q_pc<=alpha).sum()) if joint else None,
                 interpretation='Evidence in at least r assays under their respective upstream pathway null hypotheses; different genes may contribute in each assay.')
    out.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.pathway-',dir=out.parent))
    try:
        result.to_csv(stage/'pathways.tsv',sep='\t',index=False,na_rep='NA');evidence.to_csv(stage/'assay_evidence.tsv',sep='\t',index=False,na_rep='NA')
        (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
        (stage/'provenance.json').write_text(json.dumps({'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
        report=Report(stage,'KilnBio cross-assay pathway report',[
            f"Collection: {cfg['collection']}. Contrast: {cfg['contrast']}. {len(pathways)} pathways in the fixed family.",
            summary['joint_method'],summary['interpretation'],
            'Gray means untestable. Signed test scores describe rank/effect direction under the upstream method; unsigned enrichment has no activation direction. Pathway support can arise from different genes.'])
        selected=result.head(25).pathway.tolist()
        view=evidence.set_index(['pathway','layer']).reindex(pd.MultiIndex.from_product([selected,names],names=['pathway','layer'])).reset_index()
        view['display_score']=strength(view.p)*np.where(view.score_kind=='signed',np.sign(view.score),1)
        fig,ax=canvas(len(selected));matrix(fig,ax,view.display_score.to_numpy().reshape(len(selected),len(names)),selected,names,'Signed -log10(raw p); unsigned = positive')
        ax.set_title('Pathway evidence across assays')
        report.add('pathway_evidence',fig,'Top 25 by joint p when enabled, otherwise by minimum layer p for display only. Unsigned enrichment is plotted positive regardless of mechanism. All raw p and descriptive scores remain in the TSV.',view)
        fig,ax=canvas(len(selected));coverage=view.genes_in_background.to_numpy().reshape(len(selected),len(names))
        im=ax.imshow(coverage,aspect='auto',cmap='Blues');ax.set_yticks(range(len(selected)),selected,fontsize=8);ax.set_xticks(range(len(names)),names,rotation=30,ha='right');fig.colorbar(im,ax=ax,label='Measured genes in pathway');ax.set_title('Pathway coverage by assay')
        report.add('pathway_coverage',fig,'Gene coverage uses the assay-specific measured background. A large overlap does not itself establish enrichment.',view[['pathway','layer','genes_in_background']])
        support=[]
        for pathway in selected:
            subset=view[view.pathway==pathway]
            gene_layers={}
            for row in subset.itertuples():
                for gene in str(row.contributing_genes).split(';'):
                    if gene:gene_layers.setdefault(gene,[]).append(row.layer)
            for gene,supporters in sorted(gene_layers.items()):support.append(dict(pathway=pathway,gene=gene,assays=len(supporters),layers=';'.join(supporters)))
        support=pd.DataFrame(support,columns=['pathway','gene','assays','layers']);report.table('contributing_genes',support)
        counts=pd.DataFrame([dict(pathway=n,one_assay=int(((support.pathway==n)&(support.assays==1)).sum()),multiple_assays=int(((support.pathway==n)&(support.assays>1)).sum())) for n in selected])
        fig,ax=canvas(len(selected));ax.barh(counts.pathway,counts.one_assay,label='One assay',color=BLUE);ax.barh(counts.pathway,counts.multiple_assays,left=counts.one_assay,label='Multiple assays',color=ORANGE);ax.invert_yaxis();ax.set_xlabel('Distinct contributing genes');ax.set_title('Which genes support pathway agreement?');ax.legend(fontsize=8)
        report.add('shared_gene_support',fig,'Genes supplied by upstream analyses, such as leading-edge genes or ORA hits, are counted once per pathway. Empty contributions mean the upstream test did not supply gene-level support.',counts)
        # Overlap network concerns the annotation collection, not inferred regulation.
        edges=[];fig,ax=canvas(18,12);theta=np.linspace(0,2*np.pi,len(selected),endpoint=False);positions=np.column_stack([np.cos(theta),np.sin(theta)])
        for i,a in enumerate(selected):
            for j in range(i+1,len(selected)):
                b=selected[j];intersection=len(sets[a]&sets[b]);union=len(sets[a]|sets[b]);jac=intersection/union if union else 0
                if jac>=.25:
                    edges.append(dict(pathway_a=a,pathway_b=b,jaccard=jac,shared_genes=intersection));ax.plot(positions[[i,j],0],positions[[i,j],1],color='#bdbdbd',lw=1+3*jac,zorder=1)
        ax.scatter(positions[:,0],positions[:,1],color=BLUE,s=80,zorder=2)
        for i,n in enumerate(selected):ax.annotate(str(i+1),positions[i],ha='center',va='center',color='white',fontsize=7)
        ax.text(1.25,1,'\n'.join(f'{i+1}. {n}' for i,n in enumerate(selected)),fontsize=8,va='top')
        ax.set_xlim(-1.15,3.2);ax.set_ylim(-1.2,1.2)
        ax.set_aspect('equal');ax.axis('off');ax.set_title('Pathway overlap: numbered key in plotted-data TSV')
        edgeframe=pd.DataFrame(edges,columns=['pathway_a','pathway_b','jaccard','shared_genes'])
        report.table('pathway_network_key',pd.DataFrame({'node':range(1,len(selected)+1),'pathway':selected}))
        report.add('pathway_overlap',fig,'Edges join pathways with full-GMT gene-set Jaccard overlap >= 0.25. This shows annotation redundancy, not molecular interactions or independent biological confirmations.',edgeframe)
        report.finish(dict(top_n=25,network_jaccard=.25,joint_test=joint,min_assays=r,alpha=alpha))
        if out.exists():raise ValueError('Output appeared during execution')
        stage.rename(out)
    finally:
        if stage.exists():shutil.rmtree(stage)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('config');parser.add_argument('output');args=parser.parse_args()
    try:print(json.dumps(run(args.config,args.output)))
    except (ValueError,KeyError,OSError) as error:parser.exit(2,str(error)+'\n')


if __name__=='__main__':main()
