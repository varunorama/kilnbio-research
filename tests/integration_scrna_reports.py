import json,os,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import io,sparse
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
rng=np.random.default_rng(291);manifest=[]
for i in range(6):
    p=out/f'library{i}';p.mkdir();x=rng.poisson(2,(120,30));x[:10]+=i%2
    io.mmwrite(p/'matrix.mtx',sparse.coo_matrix(x));(p/'features.tsv').write_text(''.join(f'G{j}\t'+('MT-' if j<3 else '')+f'G{j}\tGene Expression\n' for j in range(120)));(p/'barcodes.tsv').write_text(''.join(f'CELL{j}\n' for j in range(30)))
    manifest.append(dict(sample=f'S{i}',matrix_dir=p.name,biological_id=f'D{i}',condition='A' if i<3 else 'B'))
pd.DataFrame(manifest).to_csv(out/'libraries.csv',index=False)
cfg=dict(libraries='libraries.csv',assembly='synthetic-human-compatible-IDs',matrix_stage='cell_called',min_genes=10,max_genes=200,max_percent_mito=30,mitochondrial_prefix='MT-',doublet_detection='none')
(out/'scrna.json').write_text(json.dumps(cfg))
def run(*a):subprocess.run(['bash',str(root/'run.sh'),*[str(x) for x in a]],check=True)
run('scrna',out/'scrna.json',out/'scrna')
cells=pd.read_csv(out/'scrna/cells.tsv',sep='\t');pd.DataFrame({'cell_id':cells.cell_id,'cell_type':['SYNTHETIC_TYPE']*len(cells)}).to_csv(out/'annotations.tsv',sep='\t',index=False)
(out/'pb.json').write_text(json.dumps(dict(h5ad='scrna/filtered.h5ad',annotations='annotations.tsv',annotations_reviewed=True,min_cells=10)))
run('pseudobulk',out/'pb.json',out/'pseudobulk')
print('SINGLE CELL AND PSEUDOBULK REPORTS COMPLETE',flush=True)
