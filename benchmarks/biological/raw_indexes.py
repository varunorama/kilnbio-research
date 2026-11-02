"""Build measured, full references for the public raw-read examples (Docker)."""
import argparse
import gzip
import hashlib
import json
import subprocess
import tarfile
from pathlib import Path
from airway_reference import sha256
from chip_spikein import IMAGES

STAR='community.wave.seqera.io/library/htslib_samtools_star_gawk@sha256:4a468118dbd7491a69bf9813c68233afa8558d1f3380fd8cab03e0e3d3135190'
RSEM='community.wave.seqera.io/library/rsem_star@sha256:5f807a3a7e9850209445f02dbcbe64778d516ca68cf5a778f49a5b5bc61fea0a'

def run(ref,label,memory,image,argv):
    command=['docker','run','--rm','--network','none','--memory',memory,'-v',f'{ref}:/ref',image,*argv]
    (ref/(label+'.command.json')).write_text(json.dumps(command,indent=2)+'\n')
    with (ref/(label+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)

def audit_rna(ref):
    with (ref/'star/transcriptInfo.tab').open() as f:
        declared=int(f.readline());star_names=[line.split('\t')[0] for line in f];star=set(star_names)
    with (ref/'rsem/genome.transcripts.fa').open() as f:
        names=[line[1:].split()[0] for line in f if line.startswith('>')]
    rsem=set(names)
    audit={'star_declared':declared,'star_unique':len(star),'rsem_unique':len(rsem),
           'star_only':sorted(star-rsem),'rsem_only':sorted(rsem-star)}
    if not declared or declared!=len(star_names) or declared!=len(star) or len(names)!=len(rsem) or star!=rsem:
        raise ValueError('Nonunique or discordant transcript identifiers; do not quantify')
    return audit

def airway(ref):
    audit=json.loads((ref/'reference-audit.json').read_text())
    for name in ['genome.fa','genome.gtf']:
        if sha256(ref/name)!=audit['output_sha256'][name]:raise ValueError('Reference changed: '+name)
    if any((ref/k).exists() for k in ['star','rsem','indexes.complete']):raise ValueError('Use a freshly prepared reference directory')
    (ref/'star').mkdir();(ref/'rsem').mkdir()
    run(ref,'rsem','10g',RSEM,['rsem-prepare-reference','--gtf','/ref/genome.gtf','--num-threads','4','/ref/genome.fa','/ref/rsem/genome'])
    run(ref,'star','10g',STAR,['STAR','--runMode','genomeGenerate','--genomeDir','/ref/star','--genomeFastaFiles','/ref/genome.fa',
        '--sjdbGTFfile','/ref/genome.gtf','--runThreadN','4','--limitGenomeGenerateRAM','9500000000','--genomeSAsparseD','8',
        '--genomeSAindexNbases','12','--sjdbOverhang','62','--outFileNamePrefix','/ref/star/'])
    (ref/'transcript-set-audit.json').write_text(json.dumps(audit_rna(ref),indent=2)+'\n')
    (ref/'indexes.complete').write_text('Both indexes built and transcript sets matched.\n')

def pc9(human,fly,ref):
    expected=[(human,'806c02398f5ac5da8ffd6da2d1d5d1a9'),(fly,'8699f76606a48f79fc4c0215f731c328')]
    for path,want in expected:
        h=hashlib.md5()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
        if h.hexdigest()!=want:raise ValueError('File disagrees with UCSC checksum: '+str(path))
    ref.mkdir(parents=True,exist_ok=False);names=set()
    with (ref/'genome.fa').open('wb') as dst:
        def append(src,prefix):
            for line in src:
                if line.startswith(b'>'):
                    name=prefix+line[1:].split()[0]
                    if name in names:raise ValueError('Duplicate contig: '+repr(name))
                    names.add(name);line=b'>'+name+b'\n'
                dst.write(line)
        with gzip.open(human,'rb') as src:append(src,b'h_')
        with tarfile.open(fly,'r:gz') as tar:
            for member in tar:
                if member.isfile() and member.name.endswith('.fa'):
                    with tar.extractfile(member) as src:append(src,b'd_')
    (ref/'reference-audit.json').write_text(json.dumps({'source_md5':{p.name:m for p,m in expected},
        'combined_sha256':sha256(ref/'genome.fa'),'contigs':sorted(n.decode() for n in names)},indent=2)+'\n')
    run(ref,'index','6g',IMAGES['bwa'],['bwa','index','-a','bwtsw','-b','100000000','/ref/genome.fa'])
    (ref/'index.complete').write_text('BWA composite index built.\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='mode',required=True)
    a=sp.add_parser('airway');a.add_argument('reference',type=Path)
    c=sp.add_parser('pc9')
    for k in ['hg19_fa_gz','dm3_chromFa_tar_gz','reference']:c.add_argument(k,type=Path)
    a=p.parse_args()
    if a.mode=='airway':airway(a.reference.resolve())
    else:pc9(a.hg19_fa_gz,a.dm3_chromFa_tar_gz,a.reference.resolve())
