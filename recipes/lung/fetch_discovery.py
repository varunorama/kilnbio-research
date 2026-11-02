from pathlib import Path
import urllib.request as http, hashlib,json,time, configparser
root=Path(__file__).resolve().parent
data=root/'discovery'
sha='860b63a7f2d23540c94e25dd35302e133c29f48b'
data.mkdir(exist_ok=True)
(data/'revision.json').write_text(json.dumps({'sha':sha}))
base=f'https://raw.githubusercontent.com/cBioPortal/datahub/{sha}/'
study='public/luad_tcga_pan_can_atlas_2018/'
def request(url,post=None):
    headers={'User-Agent':'KilnBio-public-validation','Accept':'application/vnd.git-lfs+json' if post else '*/*'}
    if post:headers['Content-Type']='application/vnd.git-lfs+json'
    for i in range(4):
        try:
            with http.urlopen(http.Request(url,data=json.dumps(post).encode() if post else None,headers=headers),timeout=180) as response:return response.read()
        except OSError:
            if i==3:raise
            time.sleep(2**i)
config=configparser.ConfigParser()
config.read_string(request(base+'.lfsconfig').decode())
service=config['lfs']['url'].strip('/').strip()+'/objects/batch'
files=['data_clinical_sample.txt','case_lists/cases_sequenced.txt','data_mutations.txt','meta_rppa.txt','data_rppa.txt','meta_mrna_seq_v2_rsem.txt','data_mrna_seq_v2_rsem.txt']
manifest={}
for name in files:
    url=base+study+name
    payload=request(url)
    if payload.startswith(b'version https://git-lfs.github.com'):
        lines=payload.decode().splitlines()
        oid=lines[1].split(':')[1]; size=int(lines[2].split()[1])
        target=data/name
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==oid:
            payload=target.read_bytes()
        else:
            item=json.loads(request(service,{'operation':'download','transfers':['basic'],'objects':[{'oid':oid,'size':size}]}))['objects'][0]
            payload=request(item['actions']['download']['href'])
        assert len(payload)==size and hashlib.sha256(payload).hexdigest()==oid,name
    target=data/name
    target.parent.mkdir(exist_ok=True)
    target.write_bytes(payload)
    manifest[name]={'url':url,'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}
    print(name,len(payload),flush=True)
(data/'sources.json').write_text(json.dumps(manifest,indent=2))

