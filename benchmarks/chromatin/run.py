"""Run the specified caller benchmark in one experiment directory per accession (Linux)."""
import argparse
import concurrent.futures
import csv
import json
import random
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MACS2 = 'quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3'
MACS3 = 'quay.io/biocontainers/macs3:3.1.0--py312h9639c01_0'


def count_lines(path):
    with path.open('rb') as file:
        return sum(1 for _ in file)


def sample(files, output, total, target, seed, excluded=None):
    rng = random.Random(seed)
    remaining = total
    with output.open('wb') as chosen:
        for f in files:
            with f.open('rb') as src:
                for line in src:
                    if rng.randrange(remaining) < target:
                        chosen.write(line)
                        target -= 1
                    elif excluded:
                        excluded.write(line)
                    remaining -= 1
    assert remaining == target == 0


def balanced(base, row):
    out = base/row['experiment']
    out.mkdir(exist_ok=True)
    if (out/'depths.json').is_file():
        return
    controls = list(filter(None,row['controls'].split(',')))
    ids = [row['rep1'],row['rep2'],*controls]
    paths = {i:base/'tags'/(i+'.bed') for i in ids}
    for acc in ids:
        subprocess.run(['bash',str(ROOT/'extract.sh'),str(base),acc,'ATAC' if row['assay']=='ATAC' else 'ChIP'],check=True)
    counts = {i:count_lines(p) for i,p in paths.items()}
    target = min(counts[row['rep1']],counts[row['rep2']])
    if target <= 0:
        raise ValueError('Empty replicate')
    for rep in (1,2):
        acc=row[f'rep{rep}']
        sample([paths[acc]],out/f'rep{rep}.bed', counts[acc],target,20261101+rep)
    depths = {'eligible':counts,'balanced_treatment':target}
    if controls:
        total = sum(counts[i] for i in controls)
        with (out/'control.test.bed').open('wb') as test:
            sample([paths[i] for i in controls],out/'control.train.bed',total,total//2,20261101,test)
        depths.update(control_training=total//2,control_evaluation=total-total//2)
    (out/'depths.json').write_text(json.dumps(depths,indent=2)+'\n')


def execute(base, row, rep, caller, jar):
    exp = row['experiment']
    out = base/exp/f'r{rep}-{caller}'
    out.mkdir(exist_ok=True)
    status = out/'status.json'
    treatment=f'/b/{exp}/rep{rep}.bed'
    dest=f'/b/{exp}/{out.name}'
    docker=['docker','run','--rm','--network','none','--cpus','3','-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1','-v',f'{base}:/b','-w',dest]
    if caller=='Omnipeak':
        (out/'cache').mkdir(exist_ok=True)
        cmd=['java','--add-modules=jdk.incubator.vector','-Xmx4g','-jar',str(jar),'analyze',
             '-t',str(base/exp/f'rep{rep}.bed'),'-c',str(base/exp/'control.train.bed'),
             '--cs',str(base/'ref/panel.sizes'),'--bl',str(base/'ref/mask.bed'),'--fr','200','--kd','--threads','2',
             '-w',str(out/'cache'),'-p',str(out/'calls.broadPeak')]
    elif caller=='HMMRATAC':
        cmd=docker+[MACS3,'macs3','hmmratac','-i',treatment,'-f','BEDPE','--randomSeed','20261101',
                    '--threads','2','-n','calls','--outdir',dest]
    else:
        binary=caller.lower()
        cmd=docker+[MACS2 if caller=='MACS2' else MACS3,binary,'callpeak','-t',treatment,'-n','calls',
                    '--outdir',dest,'-g',(base/'ref/eligible.txt').read_text().strip(),'--keep-dup','all','-q','0.05']
        if row['assay']=='ATAC':
            cmd+=['-f','BEDPE']
        else:
            cmd+=['-f','BED','-c',f'/b/{exp}/control.train.bed','--nomodel','--extsize','200','--broad','--broad-cutoff','0.1']
    if status.exists():
        previous=json.loads(status.read_text())
        if previous['exit']==0:
            if previous['command']!=cmd:
                raise ValueError(f'Command changed for {out}. Use a fresh benchmark directory.')
            return
    start=time.monotonic()
    with (out/'tool.log').open('w') as log:
        proc=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT)
    status.write_text(json.dumps({'command':cmd,'exit':proc.returncode,'seconds':round(time.monotonic()-start,2)},indent=2)+'\n')
    if proc.returncode:
        raise RuntimeError(f'See {out}/tool.log')
    print(exp,rep,caller,'complete',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('experiment')
    p.add_argument('omnipeak_jar',type=Path)
    p.add_argument('--rep-jobs',type=int,choices=(1,2),default=1,
                   help='Number of independent replicate jobs to run at once (default 1)')
    args=p.parse_args()
    with (ROOT/'experiments.tsv').open() as f:
        row=next(r for r in csv.DictReader(f,delimiter='\t') if r['experiment']==args.experiment)
    base=args.directory.resolve()
    balanced(base,row)
    def run_replicate(rep):
        for caller in (['MACS2','MACS3','HMMRATAC'] if row['assay']=='ATAC' else ['MACS3','Omnipeak']):
            execute(base,row,rep,caller,args.omnipeak_jar.resolve())
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.rep_jobs) as pool:
        list(pool.map(run_replicate, (1,2)))
