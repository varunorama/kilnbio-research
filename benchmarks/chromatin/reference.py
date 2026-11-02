"""Prepare the GRCh38 chromosome panel, blacklist and assembly-gap mask."""
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

CHROMS = ('chr1','chr19','chr21')
URLS = {
    'chrom.sizes':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'gap.txt.gz':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz',
    'blacklist.bed.gz':'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz',
}


def prepare(folder):
    folder.mkdir(parents=True,exist_ok=True)
    hashes = {}
    for name, url in URLS.items():
        data = urllib.request.urlopen(url,timeout=60).read()
        (folder/name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    sizes = {c:int(n) for c,n in (l.split() for l in (folder/'chrom.sizes').read_text().splitlines()) if c in CHROMS}
    assert len(sizes)==len(CHROMS)
    (folder/'panel.sizes').write_text(''.join(f'{c}\t{sizes[c]}\n' for c in CHROMS))
    (folder/'panel.bed').write_text(''.join(f'{c}\t0\t{sizes[c]}\n' for c in CHROMS))
    masks = {c:[] for c in CHROMS}
    for name, offset in [('blacklist.bed.gz',0),('gap.txt.gz',1)]:
        with gzip.open(folder/name,'rt') as stream:
            for line in stream:
                chrom,start,end = line.split()[offset:offset+3]
                if chrom in masks:
                    masks[chrom].append((int(start),int(end)))
    eligible = sum(sizes.values())
    with (folder/'mask.bed').open('w') as out:
        for chrom in CHROMS:
            merged=[]
            for start,end in sorted(masks[chrom]):
                if merged and start <= merged[-1][1]:
                    merged[-1][1]=max(merged[-1][1],end)
                else:
                    merged.append([start,end])
            for start,end in merged:
                assert 0 <= start < end <= sizes[chrom]
                out.write(f'{chrom}\t{start}\t{end}\n')
                eligible -= end-start
    (folder/'eligible.txt').write_text(str(eligible)+'\n')
    (folder/'references.json').write_text(json.dumps({'source':URLS,'sha256':hashes,'eligible_bp':eligible},indent=2)+'\n')


if __name__=='__main__':
    prepare(Path(sys.argv[1]))
