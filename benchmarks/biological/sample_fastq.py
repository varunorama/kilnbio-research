"""Reproducible fixed-depth FASTQ sampling with pair and complete-file audits.

Requires seqtk 1.5 and pigz. Uses seqtk's low-memory two-pass sampling, identical
seeds for paired inputs, and verifies every output record and mate identifier.
"""
import argparse, gzip, hashlib, json, subprocess
from collections import Counter
from pathlib import Path

def audit(path):
    ids=hashlib.sha256();lengths=Counter();count=0
    with gzip.open(path,'rb') as f:
        while True:
            header=f.readline()
            if not header:break
            seq=f.readline().rstrip(b'\r\n');plus=f.readline();qual=f.readline().rstrip(b'\r\n')
            if not header.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq)!=len(qual):
                raise ValueError(f'Malformed FASTQ record {count+1}: {path}')
            name=header.split()[0][1:]
            if name.endswith((b'/1',b'/2')):name=name[:-2]
            ids.update(name+b'\n');lengths[len(seq)]+=1;count+=1
    return {'reads':count,'read_id_sha256':ids.hexdigest(),'length_counts':dict(lengths)}

def sample(files,out,seqtk,n,seed):
    out.mkdir(parents=True,exist_ok=True)
    plan={'inputs':[str(f.resolve()) for f in files],'reads_per_file':n,'seed':seed,'seqtk':str(seqtk.resolve()),'method':'seqtk sample -2, paired identical seed; complete output ID audit'}
    frozen=out/'plan.json'
    if frozen.exists() and json.loads(frozen.read_text())!=plan:raise ValueError('Existing sampling plan differs')
    frozen.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    records=[]
    for source in files:
        final=out/source.name
        if not final.exists():
            print('Sampling',source.name,flush=True)
            with (out/(source.name+'.partial')).open('wb') as dst, (out/(source.name+'.log')).open('wb') as log:
                first=subprocess.Popen([str(seqtk),'sample','-2',f'-s{seed}',str(source),str(n)],stdout=subprocess.PIPE,stderr=log)
                second=subprocess.Popen(['pigz','-p','2','-c','-n'],stdin=first.stdout,stdout=dst,stderr=log)
                first.stdout.close();ret2=second.wait();ret1=first.wait()
            if ret1 or ret2:raise RuntimeError(f'Sampling failed: seqtk={ret1}, pigz={ret2}')
            (out/(source.name+'.partial')).replace(final)
        result=audit(final)
        if result['reads']!=n:raise ValueError(f'Unexpected sampled count: {result}')
        result['file']=str(final);records.append(result)
    if len(records)==2 and records[0]['read_id_sha256']!=records[1]['read_id_sha256']:
        raise ValueError('Sampled mates differ; reject before alignment')
    (out/'audit.json').write_text(json.dumps({'plan':plan,'outputs':records},indent=2)+'\n',encoding='utf-8')
    print('Verified sample',out,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path);p.add_argument('fastq',type=Path,nargs='+')
    p.add_argument('--seqtk',type=Path,required=True);p.add_argument('--reads',type=int,default=5000000);p.add_argument('--seed',type=int,default=74129)
    a=p.parse_args()
    if len(a.fastq) not in (1,2) or a.reads<1:p.error('Supply one single-end or two paired files and a positive depth')
    sample(a.fastq,a.output,a.seqtk,a.reads,a.seed)
