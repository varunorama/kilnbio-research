"""Stage fixed-depth public benchmark inputs from verified ENA downloads."""
import argparse
import csv
import json
from pathlib import Path
from sample_fastq import sample

HERE=Path(__file__).resolve().parent

def stage(fastq,study,seqtk):
    with (HERE/'public_fastq_manifest.tsv').open() as f:files=list(csv.DictReader(f,delimiter='\t'))
    groups={}
    for f in files:
        raw=fastq/f['file'];receipt=json.loads(raw.with_name(raw.name+'.verified.json').read_text())
        if not (receipt['verified'] and raw.stat().st_size==int(f['bytes'])==receipt['bytes']
                and receipt['md5']==f['md5'] and receipt['url']==f['url']):
            raise ValueError('Input does not match ENA receipt: '+str(raw))
        groups.setdefault(f['run'],[]).append(raw)
    orphans=[]
    for run,paths in sorted(groups.items()):
        if run.startswith('SRR1039'):
            mates=[fastq/(run+s) for s in ('_1.fastq.gz','_2.fastq.gz')]
            if not all(p in paths for p in mates):raise ValueError('Incomplete paired input: '+run)
            orphans.extend(str(p) for p in paths if p not in mates);paths=mates
        elif len(paths)!=1:raise ValueError('Expected single-end PC9: '+run)
        sample(paths,study/'sampled-5m'/run,seqtk,5000000,74129)
    (study/'excluded_orphans.json').write_text(json.dumps(orphans,indent=2)+'\n')

def rna_plan(study,ref):
    if not (ref/'indexes.complete').is_file():raise ValueError('Human reference/index audit incomplete')
    with (HERE/'airway_design.csv').open() as f:libs=list(csv.DictReader(f))
    out=study/'rna-inputs';out.mkdir(exist_ok=False)
    with (out/'samples.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['sample','fastq_1','fastq_2','strandedness'])
        for r in libs:
            s=r['sample'];folder=study/'sampled-5m'/s;a=json.loads((folder/'audit.json').read_text())
            reads=a['outputs']
            if len(reads)!=2 or any(x['reads']!=5000000 for x in reads) or reads[0]['read_id_sha256']!=reads[1]['read_id_sha256'] or a['plan']['seed']!=74129:raise ValueError('Wrong pair selection: '+s)
            w.writerow([s,folder/(s+'_1.fastq.gz'),folder/(s+'_2.fastq.gz'),'auto'])
    params=dict(fasta=str(ref/'genome.fa'),gtf=str(ref/'genome.gtf'),star_index=str(ref/'star'),rsem_index=str(ref/'rsem'),transcript_fasta=str(ref/'rsem/genome.transcripts.fa'),
        igenomes_ignore=True,aligner='star_rsem',star_ignore_sjdbgtf=True,skip_gtf_filter=True,
        save_align_intermeds=True,save_reference=False,skip_deseq2_qc=True,skip_stringtie=True,skip_bigwig=True,skip_qualimap=True,
        skip_dupradar=True,skip_biotype_qc=True,rseqc_modules='bam_stat,infer_experiment,read_distribution',publish_dir_mode='copy')
    (out/'params.json').write_text(json.dumps(params,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='mode',required=True)
    a=sp.add_parser('sample');a.add_argument('fastq',type=Path);a.add_argument('study',type=Path);a.add_argument('seqtk',type=Path)
    b=sp.add_parser('rna-plan');b.add_argument('study',type=Path);b.add_argument('reference',type=Path)
    a=p.parse_args()
    if a.mode=='sample':stage(a.fastq.resolve(),a.study.resolve(),a.seqtk.resolve())
    else:rna_plan(a.study.resolve(),a.reference.resolve())
