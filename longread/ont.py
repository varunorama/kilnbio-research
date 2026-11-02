"""Experimental native ONT DNA CpG path. Preserve MM/ML tags; pin local models."""
import argparse, json, os, re, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from analysis import need,file_at,digest

def check_platform(header):
    platforms={field[3:].upper() for line in header.splitlines() if line.startswith('@RG\t')
               for field in line.split('\t') if field.startswith('PL:')}
    need(not platforms-{'ONT'},'ONT route received a BAM declaring a different sequencing platform: '+','.join(sorted(platforms)))
    return {'declared_platforms':sorted(platforms),'ONT_declared':platforms=={'ONT'},
            'note':'No platform declaration is not evidence of ONT origin; review source provenance' if not platforms else 'Header platform checked'}

def plan(config,output):
    cfg=json.loads(config.read_text(encoding='utf-8-sig'));out=output.resolve()
    allowed={'input_stage','input','reference','assembly','channels','models_reviewed','basecall_model','modification_model','device','threads','filter_threshold','duplicate_policy','versions'}
    need(isinstance(cfg,dict) and not set(cfg)-allowed,'Unknown ONT keys')
    need(not out.exists(),'Use a new output directory')
    need(cfg.get('input_stage') in ('pod5','unaligned_modbam','aligned_modbam'),'Declare ONT input stage')
    need(cfg.get('assembly') and cfg.get('models_reviewed') is True,'Declare assembly and review pore/basecall/modification model compatibility')
    need(cfg.get('channels') in (['5mC'],['5mC','5hmC']),'Declare measured channels in order; 5hmC cannot be inferred from a 5mC-only model')
    need(cfg.get('duplicate_policy')=='retain_pcr_free','Initial ONT path retains PCR-free molecules; PCR-amplified inputs require another reviewed policy')
    threads=cfg.get('threads',4);need(type(threads) is int and 1<=threads<=256,'Invalid threads')
    threshold=cfg.get('filter_threshold');need(type(threshold) in (int,float) and 0<threshold<1,'Explicit modkit probability filter threshold required')
    versions=cfg.get('versions');need(isinstance(versions,dict) and set(versions)=={'dorado','modkit','samtools'},'Pin dorado/modkit/samtools versions')
    need(all(isinstance(v,str) and re.fullmatch(r'\d+\.\d+(?:\.\d+)?',v) for v in versions.values()),'Use exact numeric tool versions')
    fasta=file_at(cfg['reference'],config.parent);fai=file_at(str(fasta)+'.fai',config.parent)
    source=(config.parent/cfg['input']).resolve();inputs=[config,fasta,fai]
    need(source.exists(),'Missing ONT input')
    commands=[];calls=source
    if cfg['input_stage']=='pod5':
        files=sorted(source.rglob('*.pod5')) if source.is_dir() else [source]
        need(files and all(f.is_file() and f.suffix=='.pod5' for f in files),'POD5 input required')
        inputs+=files
        models=[]
        for key in ('basecall_model','modification_model'):
            folder=(config.parent/cfg[key]).resolve();need(folder.is_dir(),f'Local model directory required: {key}')
            fs=sorted(f for f in folder.rglob('*') if f.is_file());need(fs,'Empty model directory');inputs+=fs;models.append(folder)
        calls=out/'calls.mod.bam'
        commands.append({'argv':['dorado','basecaller',str(models[0]),str(source),'--modified-bases-models',str(models[1]),'--modified-bases-threshold','0','--device',cfg.get('device','cpu'),'--recursive'],'stdout':str(calls)})
    else:need(source.is_file() and source.suffix=='.bam','A modification-tagged BAM is required');inputs.append(source)
    if cfg['input_stage']!='aligned_modbam':
        aligned=out/'aligned.unsorted.bam'
        commands.append({'argv':['dorado','aligner',str(fasta),str(calls)],'stdout':str(aligned)})
    else:aligned=source
    sorted_bam=out/'aligned.sorted.mod.bam'
    commands += [
        {'argv':['samtools','sort','-@',str(threads),'-o',str(sorted_bam),str(aligned)]},
        {'argv':['samtools','index','-@',str(threads),str(sorted_bam)]},
        {'argv':['samtools','quickcheck','-v',str(sorted_bam)]},
        {'argv':['samtools','flagstat',str(sorted_bam)],'stdout':str(out/'alignment.flagstat.txt')},
        {'argv':['modkit','modbam','check-tags',str(sorted_bam),'--head','100'],'stdout':str(out/'modification_tags.txt')},
        {'argv':['modkit','pileup',str(sorted_bam),str(out/'cpg.bed'),'--reference',str(fasta),'--modified-bases',*cfg['channels'],'--cpg','--combine-strands','--filter-threshold',str(threshold),'--threads',str(threads)]}]
    return cfg,out,{'config':cfg,'commands':commands,'input_sha256':{str(f):digest(f) for f in inputs},'source_sha256':digest(Path(__file__)),
                    'status':'Experimental recipe; chemistry-matched biological benchmark and tag/QC review required'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('config');p.add_argument('output');p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    try:
        cfg,out,record=plan(Path(a.config).resolve(),Path(a.output))
        if a.dry_run:print(json.dumps(record,indent=2));return
        # Refuse version drift instead of silently taking whatever is on PATH.
        actual={}
        needed={c['argv'][0] for c in record['commands']}
        for name in sorted(needed):
            need(shutil.which(name),f'Tool not installed: {name}')
            flag='--version' if name!='samtools' else '--version-only'
            r=subprocess.run([name,flag],capture_output=True,text=True,check=True)
            text=r.stdout+r.stderr;need(re.search(r'(?<![0-9.])'+re.escape(cfg['versions'][name])+r'(?![0-9.])',text),f'{name} version differs from pinned configuration');actual[name]=text.strip()
        if cfg['input_stage']!='pod5':
            source=(Path(a.config).resolve().parent/cfg['input']).resolve()
            header=subprocess.run(['samtools','view','-H',str(source)],capture_output=True,text=True,check=True).stdout
            record['platform_check']=check_platform(header)
        out.mkdir(parents=True);record['actual_versions']=actual
        (out/'plan.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        try:
            for i,step in enumerate(record['commands'],1):
                with (out/f'step_{i:02d}.log').open('wb') as log:
                    if 'stdout' in step:
                        with Path(step['stdout']).open('wb') as dst:subprocess.run(step['argv'],stdout=dst,stderr=log,check=True)
                    else:subprocess.run(step['argv'],stdout=log,stderr=subprocess.STDOUT,check=True)
            (out/'status.json').write_text(json.dumps({'complete':True,'QC_review_required':True})+'\n')
        except Exception:
            (out/'status.json').write_text('{"complete":false}\n');raise
        print(f'ONT commands completed: {out}; review tags, coverage and controls before interpretation')
    except (ValueError,KeyError,OSError,subprocess.CalledProcessError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()
