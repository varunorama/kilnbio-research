"""Complete the Kang MEX example after QC; keeps author cell types explicit.

python kang_pseudobulk.py <study directory containing qc/filtered.h5ad>
KILNBIO_RSCRIPT and KILNBIO_R_LIBRARY select the reviewed R installation.
"""
import json, subprocess, sys
from pathlib import Path
import anndata as ad
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]

def main(base):
    a=ad.read_h5ad(base/'qc/filtered.h5ad')
    author=pd.read_csv(base/'inputs/cell_metadata.tsv',sep='\t',dtype=str).set_index('cell_id')
    assert a.obs.biological_id.astype(str).equals(author.loc[a.obs_names,'biological_id'].rename('biological_id'))
    author.loc[a.obs_names,['cell_type']].to_csv(base/'annotations.tsv',sep='\t',index_label='cell_id')
    cfg={'h5ad':'qc/filtered.h5ad','annotations':'annotations.tsv','annotations_reviewed':True,'min_cells':20}
    (base/'pseudobulk.json').write_text(json.dumps(cfg,indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'singlecell/analyze.py'),'pseudobulk',str(base/'pseudobulk.json'),str(base/'pseudobulk')],check=True)
    mapping=json.loads((base/'pseudobulk/cell_types.json').read_text())
    for entry in mapping:
        if entry['cell_type']!='CD14+ Monocytes':continue
        folder=base/'pseudobulk'/entry['directory']
        samples=pd.read_csv(folder/'samples.csv',dtype=str)
        assert samples.biological_id.nunique()==8 and len(samples)==16
        for engine in ('edger','deseq2'):
            cfg={'samples':'samples.csv','counts':'counts.tsv','assay':'rnaseq','engine':engine,
                 'count_source':'raw_counts','denominator':'ctrl','numerator':'stim','group':'condition',
                 'covariates':[{'name':'biological_id','type':'factor'}],'min_count':10,'min_samples':8,'fdr':.05}
            (folder/f'{engine}.json').write_text(json.dumps(cfg,indent=2)+'\n')
            subprocess.run([sys.executable,str(ROOT/'downstream/analysis.py'),'de',str(folder/f'{engine}.json'),str(folder/engine)],check=True)
        (base/'primary_cell_type.json').write_text(json.dumps(entry,indent=2)+'\n')

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
