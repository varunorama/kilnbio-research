"""Count the exact ENCODE input-control lineage in the same fixed tile universe."""
import csv,json,subprocess,sys
from pathlib import Path
from chromatin_counts import ROOT,IMAGE,publish_matrix

def main(inputs,base):
    rows=list(csv.DictReader((ROOT/'benchmarks/chromatin/experiments.tsv').open(),delimiter='\t'))
    selected=[r for r in rows if r['assay']!='ATAC']
    ids=sorted({i for r in selected for i in r['controls'].split(',')})
    out=base/'input_controls';out.mkdir(exist_ok=False)
    bams=[str(inputs/(i+'.bam')) for i in ids]
    cmd=['docker','run','--rm','--network','none','-v',f'{inputs}:{inputs}:ro','-v',f'{base}:{base}',
         IMAGE,'featureCounts','-F','SAF','-s','0','-T','4','-Q','30','--primary','--ignoreDup',
         '-a',str(base/'reference/regions.saf'),'-o',str(out/'featureCounts.tsv'),*bams]
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (out/'featureCounts.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
    publish_matrix(out,ids,bams)
    (out/'lineage.json').write_text(json.dumps(selected,indent=2)+'\n')

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
