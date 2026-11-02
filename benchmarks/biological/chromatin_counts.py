"""Full-autosome ENCODE BAM-start control experiment; run on Linux with Docker.

Usage: python chromatin_counts.py <cached ENCODE input folder> <new output>
Uses nonoverlapping 10 kb tiles fixed independently of ChIP/ATAC signal, and
excludes every tile intersecting assembly gaps or the ENCODE blacklist.
"""
import csv, gzip, hashlib, json, subprocess, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'downstream'))
from featurecounts import IMAGE, publish_matrix

URLS = {
    'chrom.sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'gap.txt.gz': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz',
    'blacklist.bed.gz': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz',
}
WIDTH = 10000


def main(inputs, out):
    out.mkdir(parents=True, exist_ok=False)
    ref = out / 'reference'; ref.mkdir()
    provenance = {'source': URLS, 'sha256': {}, 'tile_bp': WIDTH}
    for name, url in URLS.items():
        data = urllib.request.urlopen(url, timeout=90).read()
        if not data: raise ValueError('Empty reference download')
        (ref/name).write_bytes(data)
        provenance['sha256'][name] = hashlib.sha256(data).hexdigest()
    sizes = dict(line.split() for line in (ref/'chrom.sizes').read_text().splitlines())
    sizes = {f'chr{i}': int(sizes[f'chr{i}']) for i in range(1,23)}
    excluded = set()
    for name, offset in [('gap.txt.gz',1),('blacklist.bed.gz',0)]:
        with gzip.open(ref/name,'rt') as handle:
            for line in handle:
                chrom,start,end = line.split()[offset:offset+3]
                if chrom in sizes:
                    excluded.update((chrom,i) for i in range(int(start)//WIDTH,(int(end)-1)//WIDTH+1))
    regions = [(c,i*WIDTH,min((i+1)*WIDTH,n)) for c,n in sizes.items()
               for i in range((n+WIDTH-1)//WIDTH) if (c,i) not in excluded]
    saf = ref/'regions.saf'
    with saf.open('w') as handle:
        handle.write('GeneID\tChr\tStart\tEnd\tStrand\n')
        for c,a,b in regions: handle.write(f'{c}:{a}-{b}\t{c}\t{a+1}\t{b}\t.\n')
    provenance['eligible_tiles'] = len(regions)
    source = list(csv.DictReader((ROOT/'benchmarks/chromatin/files.tsv').open(),delimiter='\t'))
    experiments = list(csv.DictReader((ROOT/'benchmarks/chromatin/experiments.tsv').open(),delimiter='\t'))
    checks = {}
    for row in source:
        accession=row['file_accession']; path=inputs/(accession+('.bam' if row['role'].startswith(('rep','control')) else '.bed.gz'))
        if accession in checks: continue
        if path.stat().st_size != int(row['bytes']): raise ValueError(f'Size mismatch: {path}')
        h=hashlib.md5()
        with path.open('rb') as handle:
            for block in iter(lambda:handle.read(8*1024*1024),b''):h.update(block)
        if h.hexdigest()!=row['md5']:raise ValueError(f'MD5 mismatch: {path}')
        checks[accession]={'md5':h.hexdigest(),'bytes':path.stat().st_size,'url':row['url']}
    provenance['inputs']=checks
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # Independent DNase peak support is evaluated at the same fixed tile scale.
    support={}
    for row in experiments:
        if row['assay']!='ATAC':continue
        for rep in (1,2):
            found=set()
            with gzip.open(inputs/(row[f'dnase_rep{rep}']+'.bed.gz'),'rt') as handle:
                for line in handle:
                    c,a,b=line.split()[:3]
                    if c in sizes:found.update((c,i) for i in range(int(a)//WIDTH,(int(b)-1)//WIDTH+1))
            support[f'{row["cell"]}_r{rep}']=found
    with (out/'dnase_support.tsv').open('w') as handle:
        handle.write('feature_id\t'+'\t'.join(support)+'\n')
        for c,a,b in regions:
            handle.write(f'{c}:{a}-{b}\t'+'\t'.join(str(int((c,a//WIDTH) in v)) for v in support.values())+'\n')
    for assay in ('H3K27me3','H3K36me3','ATAC'):
        folder=out/assay;folder.mkdir()
        selected=[r for r in experiments if r['assay']==assay]
        ids=[r[f'rep{i}'] for r in selected for i in (1,2)]
        bams=[str(inputs/(i+'.bam')) for i in ids]
        cmd=['docker','run','--rm','--network','none','-v',f'{inputs}:{inputs}:ro',
             '-v',f'{out}:{out}',IMAGE,'featureCounts','-F','SAF','-s','0','-T','4',
             '-Q','30','--primary','--ignoreDup','-a',str(saf),'-o',str(folder/'featureCounts.tsv')]
        if assay=='ATAC':cmd += ['-p','--countReadPairs','-B','-C']
        cmd += bams
        (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        with (folder/'featureCounts.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
        publish_matrix(folder,ids,bams)
        with (folder/'samples.csv').open('w',newline='') as handle:
            writer=csv.writer(handle);writer.writerow(['sample','biological_id','condition'])
            for row in selected:
                for rep in (1,2):writer.writerow([row[f'rep{rep}'],row['cell']+f'_culture{rep}',row['cell']])
        cfg={'assay':'atacseq' if assay=='ATAC' else 'chipseq','engine':'edger','count_source':'raw_counts',
             'samples':'samples.csv','counts':'counts.tsv','denominator':'GM12878','numerator':'K562',
             'group':'condition','covariates':[],'fdr':.05,'min_count':20,'min_samples':2,
             'feature_universe':'external_fixed','normalization_reviewed':True}
        (folder/'de.json').write_text(json.dumps(cfg,indent=2)+'\n')
        print('Counted',assay,flush=True)

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
