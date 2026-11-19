"""Create fictional mapped omics for the offline circular report."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def create(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(117)
    sizes=pd.DataFrame({'chrom':['chr1','chr2','chr3','chr4'], 'length':[1000000,850000,700000,600000]})
    sizes.to_csv(out/'chrom_sizes.tsv',sep='\t',index=False)
    specs=[]
    for i,(name,unit,scale) in enumerate([('RNA','log2FC',3),('protein','log2FC',2),('ATAC','log2FC',3),('methylation','delta fraction',.5),('ChIP','log2FC',2)]):
        rows=[]
        for chromosome,length in sizes.itertuples(index=False):
            for j,start in enumerate(np.linspace(0,length-1000,120,dtype=int)):
                value=scale*np.clip(.6*np.sin(j/13+i*.7)+rng.normal(0,.15),-1,1)
                rows.append(dict(chrom=chromosome,start=int(start),end=int(start+1000),value=float(value) if j%31 else np.nan,feature=f'SYN_{name}_{chromosome}_{j}'))
        pd.DataFrame(rows).to_csv(out/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
        specs.append(dict(name=name,file=name+'.tsv',genome_build='synthetic-demo-v1',value_label=unit,value_min=-scale,value_max=scale))
    links=[]
    for i in range(12):
        c1=sizes.chrom.iloc[i%4];c2=sizes.chrom.iloc[(i+1)%4]
        s1=int(70000+i*25000);s2=int(120000+i*23000)
        links.append(dict(chrom1=c1,start1=s1,end1=s1+1000,chrom2=c2,start2=s2,end2=s2+1000,relation='fictional feature link'))
    pd.DataFrame(links).to_csv(out/'links.tsv',sep='\t',index=False)
    cfg=dict(schema_version=1,genome_build='synthetic-demo-v1',coordinate_system='0-based-half-open',chrom_sizes='chrom_sizes.tsv',tracks=specs,links=dict(file='links.tsv',genome_build='synthetic-demo-v1'))
    (out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (out/'README.txt').write_text('Fictional coordinates, values and links. This tests rendering, not biological validation.\n')
    return out/'config.json'


if __name__=='__main__':print(create(sys.argv[1]))
