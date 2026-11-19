"""Run a declared KilnBio research study in dependency order, without shell code."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from variant_annotations import digest

ROOT = Path(__file__).resolve().parents[1]
COMMANDS = {'rnaseq','scrnaseq','atacseq','chipseq','wgbs','rrbs','exome','ch-somatic','pathway-compare',
    'olink-import','multiomic-layers','multiomic-meta','variant-annotations','vep-annotate','alphafold','alphagenome',
    'ch-screen','ch-import','capture-panel','genomic-panel','genomic-associate','method-agreement',
    'differential-compare','assay-audit','candidate-evidence','compare-assays','compare-predictions',
    'clock-panel','featurecounts','spikein-factors','de','dm','enrich','scrna','pseudobulk',
    'import-mods','clocks','assay-clocks','epidiff-coverage','erythroid-stage','graft-readouts','ont','decision-methods'}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def save(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def source_digest():
    paths = [ROOT/'run.sh',ROOT/'prepare.py',ROOT/'versions.json',
             ROOT/'benchmarks/public_modeling/paired_predictions.py']
    for folder in ['downstream','singlecell','modifications','aging','longread','schemas']:
        paths.extend(p for p in (ROOT/folder).rglob('*') if p.suffix in {'.py','.R','.json'} and '__pycache__' not in p.parts)
    return hashlib.sha256(json.dumps({str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths))},sort_keys=True).encode()).hexdigest()


def plan(config, run_dir):
    cp = Path(config).resolve(); run_dir = Path(run_dir).resolve()
    cfg = json.loads(cp.read_text(encoding='utf-8'))
    if set(cfg) - {'schema_version','research_use','inputs','configs','steps','allow_remote_scientific_services'}:
        raise ValueError('Unknown study configuration field')
    if cfg.get('schema_version') != 1 or cfg.get('research_use') is not True:
        raise ValueError('Require schema_version=1 and research_use=true')
    if type(cfg.get('allow_remote_scientific_services',False)) is not bool:
        raise ValueError('allow_remote_scientific_services must be boolean')
    def expand(value):
        if isinstance(value,str):
            value=value.replace('{{study}}',str(cp.parent)).replace('{{run}}',str(run_dir))
            if '{{' in value or '}}' in value:
                raise ValueError('Unknown placeholder; use {{study}} and {{run}}')
            return value
        if isinstance(value,list):return [expand(v) for v in value]
        if isinstance(value,dict):return {k:expand(v) for k,v in value.items()}
        return value
    if not isinstance(cfg.get('inputs'),list) or not cfg['inputs']:
        raise ValueError('List every external study input in inputs for change detection')
    inputs={}
    for value in cfg['inputs']:
        p=(cp.parent/Path(expand(value))).resolve()
        if not p.is_file():raise ValueError('Input missing: '+str(p))
        inputs[str(p)]=digest(p)
    templates={}
    for name,path in cfg.get('configs',{}).items():
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*',name):raise ValueError('Unsafe config name')
        p=(cp.parent/Path(path)).resolve()
        templates[name]=expand(json.loads(p.read_text(encoding='utf-8')))
        inputs[str(p)]=digest(p)
    steps=cfg.get('steps')
    if not isinstance(steps,list) or not steps:raise ValueError('Provide at least one step')
    by_id={}
    for original in steps:
        if set(original) != {'id','command','args','depends_on','outputs'}:
            raise ValueError('Each step needs id, command, args, depends_on, outputs')
        step=expand(original); ident=step['id']
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*',ident) or ident in by_id:
            raise ValueError('Unique safe step IDs required')
        if step['command'] not in COMMANDS or not isinstance(step['args'],list) or not all(isinstance(a,str) for a in step['args']):
            raise ValueError('Use an allowed KilnBio command and string argument list')
        if step['command'] in {'alphafold','alphagenome'} and '--execute' in step['args'] and not cfg.get('allow_remote_scientific_services',False):
            raise ValueError('Remote scientific service requires explicit allow_remote_scientific_services=true')
        if '--dry-run' in step['args']:
            raise ValueError('Study steps must execute; use study --dry-run for a whole-plan preview')
        if not isinstance(step['depends_on'],list) or len(set(step['depends_on']))!=len(step['depends_on']):
            raise ValueError('depends_on must list unique step IDs')
        if not isinstance(step['outputs'],list) or not step['outputs']:
            raise ValueError('Declare expected output files')
        for output in step['outputs']:
            p=Path(output)
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!=ident or len(p.parts)<2:
                raise ValueError('Outputs must lie under their step ID: step_id/file')
        by_id[ident]=step
    ordered=[];remaining=dict(by_id)
    while remaining:
        ready=[s for s in remaining.values() if set(s['depends_on']).issubset({v['id'] for v in ordered})]
        if not ready:raise ValueError('Dependency cycle or unknown dependency')
        for s in ready:ordered.append(s);del remaining[s['id']]
    result={'schema_version':1,'research_use_only':True,'config_path':str(cp),'config_sha256':digest(cp),
        'study_dir':str(cp.parent),'run_dir':str(run_dir),'inputs_sha256':inputs,
        'source_sha256':source_digest(),'configs':templates,'steps':ordered,
        'runtime_selection':{k:os.environ.get(k,'') for k in ['KILNBIO_PYTHON','KILNBIO_RSCRIPT','KILNBIO_R_LIBRARY']},
        'allow_remote_scientific_services':cfg.get('allow_remote_scientific_services',False)}
    result['plan_sha256']=hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()
    return result


def execute(config, output, dry_run=False, resume=False):
    p=plan(config,output);out=Path(output).resolve()
    if out.exists() and not resume:raise ValueError('Output exists; use --resume or a new directory')
    out.mkdir(parents=True,exist_ok=True)
    oldfile=out/'study-state.json';planfile=out/'study-plan.json'
    if resume:
        if not planfile.is_file() or not oldfile.is_file():raise ValueError('No study state to resume')
        oldplan=json.loads(planfile.read_text(encoding='utf-8'))
        if oldplan['plan_sha256']!=p['plan_sha256']:
            raise ValueError('Study inputs, configuration, runtime selection or source changed; start a new run')
        state=json.loads(oldfile.read_text(encoding='utf-8'))
    else:
        save(planfile,p)
        state={'status':'planned','created_utc':stamp(),'research_use_only':True,'steps':{}}
        save(oldfile,state)
    if dry_run:return state
    lock=out/'.study.lock'
    try:
        handle=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError:
        raise ValueError('Study lock exists. Confirm no runner is active before removing the stale .study.lock')
    with os.fdopen(handle,'w') as f:f.write(str(os.getpid()))
    try:
        (out/'configs').mkdir(exist_ok=True);(out/'logs').mkdir(exist_ok=True)
        for name,value in p['configs'].items():save(out/'configs'/(name+'.json'),value)
        for step in p['steps']:
            ident=step['id'];record=state['steps'].get(ident,{})
            if record.get('status')=='complete':
                for name,sha in record['outputs_sha256'].items():
                    path=out/name
                    if not path.is_file() or digest(path)!=sha:
                        raise ValueError('Completed output changed: '+name+'; start a new run')
                continue
            # Interrupted or failed stage outputs are preserved, never overwritten/deleted.
            stage=out/ident
            if stage.exists():
                resolved=stage.resolve()
                if resolved.parent!=out or stage.is_symlink():raise ValueError('Unsafe stage path')
                archive=out/'attempts'/ident
                archive.mkdir(parents=True,exist_ok=True)
                target=archive/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
                stage.rename(target)
            attempt=int(record.get('attempt',0))+1
            log=out/'logs'/(ident+f'.{attempt}.log')
            argv=['bash',str(ROOT/'run.sh'),step['command'],*step['args']]
            from auto_reports import output_for, SUPPORTED
            if step['command'] in SUPPORTED:
                target=output_for(step['command'],step['args'])
                if target is not None and any(target==(out/n).resolve() for n in step['outputs']):
                    # File-producing commands need their step directory; directory-producing
                    # analyses deliberately require a new destination and create it themselves.
                    target.parent.mkdir(parents=True,exist_ok=True)
            print('Running '+ident+' ('+step['command']+'); log: '+str(log),flush=True)
            record={'status':'running','attempt':attempt,'started_utc':stamp(),'argv':argv,'log':str(log.relative_to(out))}
            state['steps'][ident]=record;state['status']='running';save(oldfile,state)
            with log.open('w',encoding='utf-8') as f:
                result=subprocess.run(argv,cwd=p['study_dir'],stdout=f,stderr=subprocess.STDOUT)
            record['returncode']=result.returncode
            if result.returncode:raise ValueError('Step failed: '+ident+'; see '+str(log))
            outputs={}
            for name in step['outputs']:
                path=out/name
                if not path.is_file() or out not in path.resolve().parents:
                    raise ValueError('Missing or external output: '+name)
                outputs[name]=digest(path)
            # Include automatic report assets even when the user's step lists only tables.
            reports=set((out/ident).glob('visualization*/visualization.json'))
            for name in step['outputs']:
                path=out/name
                reports.update(path.parent.glob('visualization*/visualization.json'))
                reports.update(path.parent.glob(path.name+'.visualization*/visualization.json'))
            for receipt in sorted(reports):
                if out not in receipt.resolve().parents:raise ValueError('External report receipt')
                for asset in receipt.parent.rglob('*'):
                    if asset.is_file():
                        if out not in asset.resolve().parents:raise ValueError('External report asset')
                        outputs[str(asset.relative_to(out))]=digest(asset)
            record.update(status='complete',finished_utc=stamp(),outputs_sha256=outputs)
            save(oldfile,state)
            print('Completed '+ident,flush=True)
        state.update(status='complete',finished_utc=stamp());save(oldfile,state)
    except BaseException as exc:
        state.update(status='failed',error=str(exc),updated_utc=stamp())
        for r in state['steps'].values():
            if r['status']=='running':r['status']='failed'
        save(oldfile,state);raise
    finally:
        lock.unlink()
    return state


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config');parser.add_argument('output')
    parser.add_argument('--dry-run',action='store_true');parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    try:
        state=execute(args.config,args.output,args.dry_run,args.resume)
        print('Study '+state['status']+'. See '+str(Path(args.output)/'study-state.json'))
    except (OSError,ValueError,KeyError) as exc:parser.exit(2,str(exc)+'\n')
