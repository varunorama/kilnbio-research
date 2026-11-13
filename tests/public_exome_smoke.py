"""Stage tiny public nf-core fixtures for a real exome-adapter smoke run.

Not an exome accuracy or low-VAF CH benchmark. One read pair file is used once.
Run with Python 3; next run the printed launcher command on Linux/WSL2.
"""
import argparse
import hashlib
import json
import shlex
import urllib.request
from pathlib import Path

COMMIT='3e548a3b34994380db1602506e0721fadbad6de1'
BASE=f'https://raw.githubusercontent.com/nf-core/test-datasets/{COMMIT}/data/genomics/homo_sapiens/'
FILES={
    'fasta':'genome/genome.fasta','fasta_fai':'genome/genome.fasta.fai',
    'dbsnp':'genome/vcf/dbsnp_146.hg38.vcf.gz','dbsnp_tbi':'genome/vcf/dbsnp_146.hg38.vcf.gz.tbi',
    'known_indels':'genome/vcf/mills_and_1000G.indels.vcf.gz','known_indels_tbi':'genome/vcf/mills_and_1000G.indels.vcf.gz.tbi',
    'fastq_1':'illumina/fastq/test_1.fastq.gz','fastq_2':'illumina/fastq/test_2.fastq.gz'}

def stage(out):
    out=Path(out).resolve()
    if out.exists(): raise ValueError('Choose a new directory')
    if any(c.isspace() for c in str(out)): raise ValueError('Use a path without spaces')
    out.mkdir(parents=True)
    receipt=[]; paths={}
    for key,rel in FILES.items():
        path=out/Path(rel).name
        with urllib.request.urlopen(BASE+rel,timeout=90) as src, path.open('wb') as dst:
            size=0
            while True:
                block=src.read(1024**2)
                if not block: break
                size+=len(block)
                if size>100*1024**2: raise ValueError('Fixture unexpectedly exceeds 100 MiB')
                dst.write(block)
        receipt.append(dict(url=BASE+rel,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        paths[key]=str(path)
    bed=out/'targets.bed'
    contigs=[line.split('\t')[:2] for line in Path(paths['fasta_fai']).read_text().splitlines()]
    bed.write_text(''.join(f'{c}\t0\t{n}\n' for c,n in contigs if c.removeprefix('chr').isdigit() and 1<=int(c.removeprefix('chr'))<=22))
    params={k:v for k,v in paths.items() if not k.startswith('fastq')}
    # This fixture produces an empty germline callset. GATK tranche fitting rejects
    # empty VCFs; explicitly omit that stage in this smoke, never claim filtered accuracy.
    params.update(wes=True,aligner='bwa-mem2',tools='haplotypecaller',igenomes_ignore=True,intervals=str(bed),split_fastq=0,skip_tools='haplotypecaller_filter')
    (out/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    (out/'samples.csv').write_text('patient,sex,status,sample,lane,fastq_1,fastq_2\n'+f'test,XX,0,test,L1,{paths["fastq_1"]},{paths["fastq_2"]}\n')
    (out/'source_receipt.json').write_text(json.dumps(dict(scope='Tiny unfiltered caller/QC smoke only; not biological, post-filter or variant-accuracy validation',commit=COMMIT,files=receipt),indent=2)+'\n')
    cmd=['bash',str(Path(__file__).resolve().parents[1]/'run.sh'),'exome',str(out/'samples.csv'),str(out/'params.json'),str(out/'run'),'--cpus','2','--memory-gb','10']
    print(shlex.join(cmd))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out');stage(p.parse_args().out)
