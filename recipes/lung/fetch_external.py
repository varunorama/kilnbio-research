from pathlib import Path
from hashlib import sha256
import urllib.request, gzip, json,os,time
root=Path(__file__).resolve().parent
pred=root/'analysis/frozen.tsv'
lock=(root/'analysis/frozen.sha256').read_text(encoding='utf-8-sig').strip().lower()
assert lock==sha256(pred.read_bytes()).hexdigest()
out=root/'external'
out.mkdir(exist_ok=True)
base='https://ftp.ncbi.nlm.nih.gov/geo/'
files={'GSE31210_series_matrix.txt.gz':base+'series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz','GPL570.annot.gz':base+'platforms/GPLnnn/GPL570/annot/GPL570.annot.gz'}
info={}
for name,url in files.items():
    dest=out/name
    if not dest.exists():
        for i in range(3):
            try:
                with urllib.request.urlopen(url,timeout=180) as r,dest.with_suffix('.part').open('wb') as w:
                    assert int(r.headers.get('Content-Length',0)) < 150_000_000, 'Source exceeds study file cap'
                    while buf:=r.read(1024*1024):w.write(buf)
                with gzip.open(dest.with_suffix('.part'),'rb') as g:
                    while g.read(1024*1024):pass
                os.replace(dest.with_suffix('.part'),dest);break
            except OSError:
                if i==2:raise
                time.sleep(2**i)
    info[name]={'url':url,'sha256':sha256(dest.read_bytes()).hexdigest(),'bytes':dest.stat().st_size}
    print(name,info[name]['bytes'],flush=True)
(out/'sources.json').write_text(json.dumps({'frozen_sha256_before_access':lock,'files':info},indent=2))


