"""Real upstream ONT control reads -> KilnBio runner -> mC and hmC count imports.

Requires Modkit 0.6.4 and Samtools 1.24 on PATH; no Dorado command is executed.
python ont_controls.py <input cache> <new study directory>
"""
import hashlib,json,subprocess,sys,urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]

def main(cache,out):
    if out.exists():raise ValueError('Use a new study directory')
    cache.mkdir(parents=True,exist_ok=True)
    source=json.loads(Path(__file__).with_name('ont_sources.json').read_text())
    for name,record in source['files'].items():
        path=cache/name
        if not path.exists():
            payload=urllib.request.urlopen(record['url'],timeout=90).read()
            if hashlib.sha256(payload).hexdigest()!=record['sha256']:raise ValueError('Download checksum mismatch')
            path.write_bytes(payload)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('Cached checksum mismatch')
    out.mkdir(parents=True)
    (out/'sources.json').write_text(json.dumps(source,indent=2)+'\n')
    cfg={'input_stage':'aligned_modbam','input':str(cache/'bc_anchored_10_reads.sorted.bam'),
         'reference':str(cache/'CGI_ladder_3.6kb_ref.fa'),'assembly':'ONT synthetic oligonucleotide control reference',
         'channels':['5mC','5hmC'],'models_reviewed':True,'duplicate_policy':'retain_pcr_free',
         'threads':2,'filter_threshold':.8,'versions':{'dorado':'2.1.2','modkit':'0.6.4','samtools':'1.24'}}
    (out/'ont.json').write_text(json.dumps(cfg,indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'longread/ont.py'),str(out/'ont.json'),str(out/'runner')],check=True)
    for channel in ('5mC','5hmC'):
        cfg={'format':'modkit_bedmethyl','input':'runner/cpg.bed','reference':str(cache/'CGI_ladder_3.6kb_ref.fa'),
             'assembly':'ONT synthetic oligonucleotide control reference','modification':channel,'chemistry':'nanopore',
             'pipeline_version':'0.6.4','measurement_reviewed':True}
        p=out/(channel+'.json');p.write_text(json.dumps(cfg,indent=2)+'\n')
        subprocess.run([sys.executable,str(ROOT/'modifications/import_counts.py'),str(p),str(out/channel)],check=True)

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
