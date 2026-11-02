"""10x gene-expression MEX -> preserved raw counts, QC, exploratory Scanpy, pseudobulk.

Cell clusters are exploratory labels, never automatic cell-type annotations.
"""
import argparse, csv, gzip, hashlib, json, platform, sys, os
from importlib.metadata import version
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from analysis import need, rows, file_at, digest

def load_mex(folder):
    import anndata as ad
    import numpy as np
    from scipy import io, sparse
    def locate(names):
        for name in names:
            p=folder/name
            if p.is_file():return p
        raise ValueError(f'Missing MEX file in {folder}: {names}')
    matrix=locate(['matrix.mtx.gz','matrix.mtx'])
    features=locate(['features.tsv.gz','features.tsv','genes.tsv.gz','genes.tsv'])
    barcodes=locate(['barcodes.tsv.gz','barcodes.tsv'])
    def read(path):
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rt') as f:return list(csv.reader(f,delimiter='\t'))
    feat=read(features);bc=read(barcodes)
    opener=gzip.open if matrix.suffix=='.gz' else open
    with opener(matrix,'rb') as f:x=sparse.csr_matrix(io.mmread(f).T)
    need(x.shape==(len(bc),len(feat)),'MEX dimensions do not match barcodes/features')
    need(np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all(),'10x matrix must contain raw integer counts')
    need(len({f[0] for f in feat})==len(feat),'Duplicate feature IDs')
    need(len({b[0] for b in bc})==len(bc),'Duplicate cell barcodes within a library')
    a=ad.AnnData(x.astype('int64'))
    a.var_names=[f[0] for f in feat];a.var['gene_symbol']=[f[1] if len(f)>1 else f[0] for f in feat]
    a.var['feature_type']=[f[2] if len(f)>2 else 'Gene Expression' for f in feat]
    a.obs_names=[b[0] for b in bc]
    a=a[:,a.var['feature_type']=='Gene Expression'].copy()
    return a,[matrix,features,barcodes]

def process(config,output):
    os.environ.setdefault('MPLBACKEND','Agg')
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    cfg=json.loads(config.read_text(encoding='utf-8-sig'));out=output.resolve()
    need(isinstance(cfg,dict) and not set(cfg)-{'libraries','cell_metadata','assembly','matrix_stage','min_genes','max_genes','max_percent_mito','mitochondrial_prefix','doublet_detection','leiden_resolution'},'Unknown single-cell configuration keys')
    need(not out.exists(),'Use a new output directory')
    need(cfg.get('assembly'),'Assembly/reference annotation must be explicit')
    need(cfg.get('matrix_stage') in ('cell_called','raw_droplets'),'Declare raw_droplets or cell_called')
    need(cfg['matrix_stage']=='cell_called','Initial downstream module requires a cell-called matrix; run upstream empty-droplet/cell calling first')
    for k in ('min_genes','max_genes','max_percent_mito','mitochondrial_prefix','doublet_detection'):
        need(k in cfg,f'Explicit QC setting required: {k}')
    need(type(cfg['min_genes']) is int and type(cfg['max_genes']) is int and 0<=cfg['min_genes']<cfg['max_genes'],'Invalid gene thresholds')
    need(0<=cfg['max_percent_mito']<=100 and isinstance(cfg['mitochondrial_prefix'],str) and cfg['mitochondrial_prefix'],'Invalid mitochondrial QC settings')
    need(cfg['doublet_detection'] in ('scrublet','none'),'Doublet detection must be scrublet or explicit none')
    manifest=file_at(cfg['libraries'],config.parent);header,data=rows(manifest)
    need({'sample','matrix_dir'}<=set(header),'Libraries CSV needs sample,matrix_dir')
    need(len({s['sample'] for s in data})==len(data),'One row per 10x library; merge lanes upstream')
    matrices=[];hashes={str(manifest):digest(manifest),str(config):digest(config)}
    feature_ids=None
    for row in data:
        need(row['sample'] and ':' not in row['sample'],'Library ID must be nonempty without colon')
        folder=(manifest.parent/row['matrix_dir']).resolve()
        a,files=load_mex(folder)
        if feature_ids is None:feature_ids=list(a.var_names)
        need(set(a.var_names)==set(feature_ids),'Libraries must share the same gene annotation; no silent outer join/zero filling')
        a=a[:,feature_ids].copy()
        if matrices:
            need(a.var['gene_symbol'].equals(matrices[0].var['gene_symbol']),'Gene IDs have inconsistent symbols/annotation between libraries')
        a.obs['barcode']=a.obs_names
        a.obs_names=[row['sample']+':'+str(b) for b in a.obs_names]
        for key,value in row.items():
            if key!='matrix_dir':a.obs[key]=value
        for key in ('biological_id','condition'):
            if key not in a.obs:a.obs[key]='unknown'
        a.var['mt']=a.var['gene_symbol'].str.startswith(cfg['mitochondrial_prefix']).to_numpy()
        sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
        if cfg['doublet_detection']=='scrublet':
            need(a.n_obs>=100,'Scrublet needs an adequately sized library; use a reviewed alternative for tiny inputs')
            sc.pp.scrublet(a,random_state=20261031,n_prin_comps=min(30,a.n_vars-2,a.n_obs-2))
        else:
            a.obs['doublet_score']=np.nan;a.obs['predicted_doublet']=False
        a.obs['qc_pass']=(a.obs.n_genes_by_counts>=cfg['min_genes']) & (a.obs.n_genes_by_counts<=cfg['max_genes']) & (a.obs.pct_counts_mt<=cfg['max_percent_mito']) & ~a.obs.predicted_doublet
        matrices.append(a)
        hashes.update({str(p):digest(p) for p in files})
    combined=ad.concat(matrices,join='inner',merge='same',index_unique=None)
    if 'cell_metadata' in cfg:
        metadata_path=file_at(cfg['cell_metadata'],config.parent)
        metadata=pd.read_csv(metadata_path,sep='\t',dtype=str,keep_default_na=False)
        need('cell_id' in metadata and 'biological_id' in metadata,'Cell metadata requires cell_id and biological_id')
        need(not metadata.cell_id.duplicated().any() and set(metadata.cell_id)==set(combined.obs_names),
             'Cell metadata must match every input cell exactly, using sample:barcode IDs')
        metadata=metadata.set_index('cell_id').loc[combined.obs_names]
        need(not metadata.biological_id.str.lower().isin(['','unknown','na','nan']).any(),'Cell donors cannot be missing or unknown')
        for field in metadata:
            if field in combined.obs:
                current=combined.obs[field].astype(str)
                compatible=(current==metadata[field])
                if field=='biological_id':compatible=compatible | (current=='unknown')
                need(compatible.all(),f'Cell metadata conflicts with library/QC metadata: {field}')
            combined.obs[field]=metadata[field].to_numpy()
        hashes[str(metadata_path)]=digest(metadata_path)
    combined.layers['counts']=combined.X.copy()
    combined.uns['kilnbio']={'assembly':cfg['assembly'],'matrix_stage':cfg['matrix_stage'],
                             'doublet_detection':cfg['doublet_detection'],'seed':20261031}
    # Explicit float workspace prevents integer truncation during normalization.
    filtered=combined[combined.obs.qc_pass].copy();filtered.X=filtered.X.astype('float64')
    need(filtered.n_obs>=4 and filtered.n_vars>=4,'Too few retained cells/genes for exploration')
    out.mkdir(parents=True);(out/'provenance').mkdir()
    combined.obs.to_csv(out/'cell_qc.tsv',sep='\t',index_label='cell_id')
    combined.write_h5ad(out/'all_cells_raw.h5ad')
    # Raw counts retain ALL genes; variable-gene selection affects exploration only.
    sc.pp.normalize_total(filtered,target_sum=1e4);sc.pp.log1p(filtered)
    detected=np.asarray((filtered.layers['counts']>0).sum(axis=0)).ravel()>0
    explore=filtered[:,detected].copy()
    need(explore.n_vars>=4,'Too few expressed genes for exploration')
    sc.pp.highly_variable_genes(explore,n_top_genes=min(2000,explore.n_vars),flavor='seurat')
    if int(explore.var.highly_variable.sum())<3:explore.var['highly_variable']=True
    sc.pp.pca(explore,n_comps=min(30,explore.n_obs-2,int(explore.var.highly_variable.sum())-1),random_state=20261031)
    sc.pp.neighbors(explore,n_neighbors=min(15,explore.n_obs-1),random_state=20261031)
    sc.tl.umap(explore,random_state=20261031)
    sc.tl.leiden(explore,resolution=float(cfg.get('leiden_resolution',.5)),random_state=20261031,flavor='igraph',n_iterations=2,directed=False)
    filtered.obsm['X_pca']=explore.obsm['X_pca'];filtered.obsm['X_umap']=explore.obsm['X_umap']
    filtered.obs['leiden']=explore.obs.leiden.astype('str')
    filtered.uns['kilnbio']['cell_types']='not automatically annotated; Leiden labels are exploratory'
    filtered.write_h5ad(out/'filtered.h5ad')
    filtered.obs.to_csv(out/'cells.tsv',sep='\t',index_label='cell_id')
    sc.settings.figdir=out
    figure=sc.pl.umap(explore,color=['leiden','sample'],show=False,return_fig=True)
    figure.savefig(out/'umap_exploration.png',bbox_inches='tight',dpi=150)
    record={'input_sha256':hashes,'config':cfg,'python':platform.python_version(),'scanpy':version('scanpy'),
            'cells_input':combined.n_obs,'cells_retained':filtered.n_obs,'genes_retained':filtered.n_vars,
            'source_sha256':digest(Path(__file__)),
            'limitations':'Exploratory embeddings; no automatic cell typing, integration or biological DE inference'}
    (out/'provenance/run.json').write_text(json.dumps(record,indent=2)+'\n')
    return record

def pseudobulk(config,output):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    cfg=json.loads(config.read_text(encoding='utf-8-sig'));out=output.resolve()
    need(isinstance(cfg,dict) and not set(cfg)-{'h5ad','annotations','annotations_reviewed','covariates','min_cells'},'Unknown pseudobulk configuration keys')
    need(not out.exists(),'Use a new output directory')
    h5=file_at(cfg['h5ad'],config.parent);annotation=file_at(cfg['annotations'],config.parent)
    need(cfg.get('annotations_reviewed') is True,'Cell-type assignments require explicit review')
    a=ad.read_h5ad(h5);need('counts' in a.layers,'Raw counts layer is required; normalized X is not used')
    meta=pd.read_csv(annotation,sep='\t',dtype=str,keep_default_na=False)
    need(list(meta.columns)==['cell_id','cell_type'] and not meta.cell_id.duplicated().any(),'Annotations need unique cell_id,cell_type columns')
    need(set(meta.cell_id)==set(a.obs_names),'Annotations must match all retained cells exactly')
    a.obs['cell_type']=meta.set_index('cell_id').loc[a.obs_names,'cell_type']
    fields=['biological_id','condition','cell_type']+cfg.get('covariates',[])
    need(len(fields)==len(set(fields)) and all(f in a.obs for f in fields),'Missing/repeated pseudobulk metadata')
    need(not a.obs[fields].isna().any().any(),'Missing pseudobulk metadata')
    need(not a.obs[fields].astype(str).apply(lambda x:x.str.lower().isin(['','unknown','na','nan'])).any().any(),'Unknown donor/group/cell-type labels cannot become independent replicates')
    raw=sparse.csr_matrix(a.layers['counts'])
    need(np.isfinite(raw.data).all() and (raw.data>=0).all() and (raw.data==np.floor(raw.data)).all(),'Pseudobulk requires raw integer counts')
    minimum=cfg.get('min_cells',20);need(type(minimum) is int and minimum>=1,'min_cells must be positive integer')
    groups=a.obs.groupby(['cell_type','biological_id','condition'],observed=True,sort=True).indices
    # Check everything before creating output; one count matrix per reviewed cell type.
    bytype={};exclusions=[]
    for (ctype,donor,condition),indices in groups.items():
        sub=a.obs.iloc[indices]
        for field in cfg.get('covariates',[]):need(sub[field].nunique()==1,'Inconsistent covariates within donor/condition; resolve technical libraries explicitly')
        if len(indices)<minimum:
            exclusions.append(dict(cell_type=ctype,biological_id=donor,condition=condition,cells=len(indices)));continue
        bytype.setdefault(str(ctype),[]).append((indices,sub.iloc[0]))
    need(bytype,'No pseudobulk groups pass the cell threshold')
    out.mkdir(parents=True);mapping=[]
    for i,(ctype,items) in enumerate(sorted(bytype.items()),1):
        folder=out/f'celltype_{i:03d}';folder.mkdir()
        counts={};records=[]
        for j,(indices,row) in enumerate(items,1):
            sid=f'PB{j:04d}';counts[sid]=np.asarray(raw[indices,:].sum(axis=0)).ravel().astype('int64')
            rec={'sample':sid,**{f:str(row[f]) for f in fields},'cells':len(indices)};records.append(rec)
        pd.DataFrame(counts,index=a.var_names).to_csv(folder/'counts.tsv',sep='\t',index_label='feature_id')
        pd.DataFrame(records).to_csv(folder/'samples.csv',index=False)
        mapping.append({'directory':folder.name,'cell_type':ctype,'pseudobulk_samples':len(items)})
    (out/'cell_types.json').write_text(json.dumps(mapping,indent=2)+'\n')
    pd.DataFrame(exclusions,columns=['cell_type','biological_id','condition','cells']).to_csv(out/'excluded_groups.tsv',sep='\t',index=False)
    (out/'provenance.json').write_text(json.dumps({'config':cfg,'input_sha256':{str(h5):digest(h5),str(annotation):digest(annotation)},'source_sha256':digest(Path(__file__))},indent=2)+'\n')
    return mapping

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['scrna','pseudobulk']);p.add_argument('config');p.add_argument('output')
    args=p.parse_args()
    try:
        result=(process if args.mode=='scrna' else pseudobulk)(Path(args.config).resolve(),Path(args.output))
        print(json.dumps(result,indent=2))
    except (ValueError,KeyError,OSError,ImportError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()
