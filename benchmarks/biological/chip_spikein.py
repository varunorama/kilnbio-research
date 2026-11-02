"""GSE64243 fixed-depth raw ChIP spike-in benchmark; no unreplicated DE tests.

Inputs: sample_fastq.py outputs for the four prespecified runs; full composite
BWA index named genome.fa with h_ (hg19) and d_ (dm3) contig prefixes. See
RAW_PROTOCOL.md. This dedicated benchmark does not invent input controls or
bypass the ordinary ChIP adapter's matched-input requirement.
"""
import argparse, collections, csv, hashlib, json, subprocess, sys, time
from pathlib import Path

SAMPLES={'SRR2532668':('H3K27me3','DMSO'),'SRR2532669':('H3K27me3','EZH2i'),
         'SRR2532670':('H3K4me3','DMSO'),'SRR2532667':('H3K4me3','EZH2i')}
R1={'SRR1721586':('H3K27me3','DMSO'),'SRR1721587':('H3K27me3','EZH2i'),
    'SRR1721588':('H3K4me3','DMSO'),'SRR1721589':('H3K4me3','EZH2i')}
IMAGES={
'trim':'community.wave.seqera.io/library/trim-galore@sha256:51c7ccb48048e4b08bcbc60796d9be8e41cda648dba02e8710f71c1b74f008d8',
'fastqc':'quay.io/biocontainers/fastqc@sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900',
'bwa':'quay.io/biocontainers/bwa@sha256:a8ea7d74457b55395e35f1e4200e3ac4a44e3f7c4b361d7628da86cec2133e03',
'samtools':'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'}

def docker(image,inputs,reference,out,*args,stdin=False):
    return ['docker','run','--rm',*(['-i'] if stdin else []),'--network','none','--memory','5g',
            '-v',f'{inputs}:/reads:ro','-v',f'{reference}:/ref:ro','-v',f'{out}:/out',IMAGES[image],*map(str,args)]

def run(cmd,folder,name,stdout_file=None):
    marker=folder/(name+'.complete')
    if marker.exists():return
    (folder/(name+'.command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
    with (folder/(name+'.log')).open('wb') as log:
        if stdout_file:
            with stdout_file.open('wb') as target:subprocess.run(cmd,check=True,stdout=target,stderr=log)
        else:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
    marker.write_text('Completed\n')

def count_bam(cmd,folder):
    output=folder/'counts.json'
    if output.exists():return json.loads(output.read_text())
    stats={k:collections.Counter() for k in ['primary','duplicates_retained','mapq26']}
    bins=collections.Counter();lengths={}
    (folder/'count.command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    with (folder/'count.log').open('wb') as log:
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,text=True)
        for line in proc.stdout:
            if line.startswith('@'):
                if line.startswith('@SQ\t'):
                    fields=dict(v.rstrip().split(':',1) for v in line.split('\t')[1:])
                    lengths[fields['SN']]=int(fields['LN'])
                continue
            row=line.split('\t');flag=int(row[1]);chrom=row[2];mq=int(row[4]);pos=int(row[3])-1
            if chrom not in lengths or not 0<=pos<lengths[chrom]:raise ValueError('Read outside the declared reference')
            species='human' if chrom.startswith('h_') else 'fly' if chrom.startswith('d_') else 'unexpected'
            if species=='unexpected':raise ValueError('Composite genome species prefix missing')
            if not flag&1024:stats['mapq26'][species]+=1
            if mq<30:continue
            stats['duplicates_retained'][species]+=1
            if flag&1024:continue
            stats['primary'][species]+=1
            if species=='human':bins[(chrom,pos//10000)]+=1
        if proc.wait():raise RuntimeError('samtools count stream failed')
    with (folder/'human_10kb_leftmost_starts.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['chrom','start','end','reads'])
        for (chrom,b),count in sorted(bins.items()):w.writerow([chrom,b*10000,min((b+1)*10000,lengths[chrom]),count])
    result={k:dict(v) for k,v in stats.items()}
    temp=output.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(output)
    return result

def main(reads,ref,out,campaign='R2',after=None):
    samples=R1 if campaign=='R1' else SAMPLES
    reads,ref,out=reads.resolve(),ref.resolve(),out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    plan={'samples':samples,'inputs':str(reads),'reference':str(ref),'images':IMAGES,'reads_per_sample':5000000,
          'seed':74129,'mapq':30,'max_mismatches':2,'layout':'single','biological_replicates_per_condition_mark':1}
    planfile=out/'plan.json'
    if planfile.exists() and json.loads(planfile.read_text())!=json.loads(json.dumps(plan)):raise ValueError('Existing plan differs')
    planfile.write_text(json.dumps(plan,indent=2)+'\n')
    for sample in samples:
        folder=out/sample;folder.mkdir(exist_ok=True)
        audit=reads/sample/'audit.json'
        while not audit.exists():time.sleep(30)
        a=json.loads(audit.read_text())
        if a['outputs'][0]['reads']!=5000000 or a['plan']['seed']!=74129:raise ValueError('Sampling protocol mismatch')
        raw=f'/reads/{sample}/{sample}.fastq.gz'
        run(docker('trim',reads,ref,out,'trim_galore','--cores','2','--illumina','--quality','20','--length','20','--gzip','--output_dir',f'/out/{sample}',raw),folder,'trim')
        trimmed=f'/out/{sample}/{sample}_trimmed.fq.gz'
        run(docker('fastqc',reads,ref,out,'fastqc','--threads','2','--outdir',f'/out/{sample}',raw,trimmed),folder,'fastqc')
    while not (ref/'index.complete').exists():time.sleep(30)
    if after:
        while not after.is_file():time.sleep(30)
    reports={}
    for sample in samples:
        folder=out/sample;trimmed=f'/out/{sample}/{sample}_trimmed.fq.gz'
        print('Aligning',sample,flush=True)
        run(docker('bwa',reads,ref,out,'bwa','aln','-t','4','-n','2','/ref/genome.fa',trimmed),folder,'aln',folder/'reads.sai')
        if not (folder/'sort.complete').exists():
            align=docker('bwa',reads,ref,out,'bwa','samse','/ref/genome.fa',f'/out/{sample}/reads.sai',trimmed)
            sort=docker('samtools',reads,ref,out,'samtools','sort','-@','3','-m','512M','-o',f'/out/{sample}/sorted.bam','-',stdin=True)
            (folder/'sort.commands.json').write_text(json.dumps([align,sort],indent=2)+'\n')
            with (folder/'samse.log').open('wb') as alog,(folder/'sort.log').open('wb') as slog:
                p1=subprocess.Popen(align,stdout=subprocess.PIPE,stderr=alog)
                p2=subprocess.Popen(sort,stdin=p1.stdout,stderr=slog);p1.stdout.close()
                r2=p2.wait();r1=p1.wait()
            if r1 or r2:raise RuntimeError(f'samse/sort failure {r1}/{r2}')
            (folder/'sort.complete').write_text('Completed\n')
        run(docker('samtools',reads,ref,out,'samtools','markdup','-@','3','-s','-f',f'/out/{sample}/duplicates.txt',f'/out/{sample}/sorted.bam',f'/out/{sample}/marked.bam'),folder,'markdup')
        run(docker('samtools',reads,ref,out,'samtools','index','-@','3',f'/out/{sample}/marked.bam'),folder,'index')
        run(docker('samtools',reads,ref,out,'samtools','quickcheck','-v',f'/out/{sample}/marked.bam'),folder,'quickcheck')
        run(docker('samtools',reads,ref,out,'samtools','flagstat','-@','3',f'/out/{sample}/marked.bam'),folder,'flagstat')
        cmd=docker('samtools',reads,ref,out,'samtools','view','-h','-F','2820','-q','26','-e','[NM]<=2',f'/out/{sample}/marked.bam')
        reports[sample]=count_bam(cmd,folder)
    # Each mark uses its own reference normalization; no cross-antibody pooling.
    for sensitivity in ['primary','duplicates_retained','mapq26']:
        table=out/(sensitivity+'.csv')
        with table.open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['sample','normalization_group','target_reads','reference_reads','target_amount','reference_amount','amount_unit','reference_material'])
            for sample,(mark,condition) in samples.items():
                c=reports[sample][sensitivity]
                w.writerow([sample,mark,c.get('human',0),c.get('fly',0),30000,750,'ng_chromatin','OSS_H2Av'])
        for mark in ['H3K27me3','H3K4me3']:
            dest=out/(sensitivity+'_'+mark)
            if not dest.exists():subprocess.run([sys.executable,str(Path(__file__).resolve().parents[2]/'downstream/spikein.py'),str(table),str(dest),'--group',mark],check=True)
    (out/'benchmark.complete').write_text('Raw-read processing and descriptive reference normalization completed; biological control evaluation is separate.\n')
    print(out,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('sampled_reads',type=Path);p.add_argument('reference',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--campaign',choices=('R1','R2'),default='R2');p.add_argument('--after',type=Path,help='Optional completion marker for scheduling on a memory-limited host')
    a=p.parse_args();main(a.sampled_reads,a.reference,a.output,a.campaign,a.after)
