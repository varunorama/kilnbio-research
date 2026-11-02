"""Audit legacy chr21 coordinate identity by exact CpG agreement to hg18/hg19.

This is an explicit sequence-based inference, not missing author metadata.
python colon_reference.py <colon study folder> <reference cache>
"""
import csv,gzip,hashlib,json,sys,urllib.request
from pathlib import Path

def main(study,cache):
    cache.mkdir(parents=True,exist_ok=True)
    positions=[int(r['pos']) for r in csv.DictReader((study/'inputs/positions.tsv').open(),delimiter='\t')]
    record={}
    for build in ('hg18','hg19'):
        url=f'https://hgdownload.soe.ucsc.edu/goldenPath/{build}/chromosomes/chr21.fa.gz'
        path=cache/(build+'.chr21.fa.gz')
        if not path.exists():
            b=urllib.request.urlopen(url,timeout=120).read()
            if b[:2]!=b'\x1f\x8b':raise ValueError('Not gzip sequence')
            path.write_bytes(b)
        with gzip.open(path,'rt') as handle:sequence=''.join(l.strip() for l in handle if not l.startswith('>')).upper()
        hits=sum(sequence[p-1:p+1]=='CG' for p in positions)
        record[build]={'url':url,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                       'matching_CpGs':hits,'total_sites':len(positions),'fraction':hits/len(positions)}
    selected=max(record,key=lambda x:record[x]['fraction'])
    if record[selected]['fraction']<.999:raise ValueError('Reference unresolved; do not assign assembly')
    record['inferred_assembly']=selected
    record['limitation']='Native data object lacks assembly metadata; inference from exact one-based forward-C CpG positions, not an author declaration. No feature annotation transfer.'
    (study/'reference_audit.json').write_text(json.dumps(record,indent=2)+'\n')
    cfg={'samples':'inputs/samples.csv','methylation_manifest':'inputs/methylation.csv',
         'engine':'dss','denominator':'normal','numerator':'cancer','group':'condition',
         'covariates':[{'name':'biological_id','type':'factor'}],'fdr':.05,'min_coverage':5,
         'modification':'modC','assembly':selected+' (sequence-audited inference)',
         'coordinates':'1-based-cytosine','measurement_reviewed':True,'smoothing':False}
    (study/'dss.json').write_text(json.dumps(cfg,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
