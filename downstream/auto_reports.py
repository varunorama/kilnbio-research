"""Run an analysis and transactionally publish its offline visual report."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
POSITION={'variant-annotations':1,'alphafold':1,'vep-annotate':1,'alphagenome':1,
 'multiomic-layers':1,'genomic-associate':5,'differential-compare':5,'compare-assays':3,
 'featurecounts':1,'spikein-factors':1,'de':1,'dm':1,'scrna':1,'pseudobulk':1,
 'import-mods':1,'ont':1,**{c:2 for c in ['wgbs','rrbs','atacseq','chipseq','rnaseq','scrnaseq','exome','ch-somatic']}}
FLAG={'olink-import','ch-screen','ch-import','capture-panel','genomic-panel','method-agreement',
 'assay-audit','candidate-evidence','compare-predictions','clock-panel','clocks','assay-clocks',
 'epidiff-coverage','erythroid-stage','graft-readouts','bladder-example','decision-methods'}
SUPPORTED=set(POSITION)|FLAG


def option(args,flag):
    for i,a in enumerate(args):
        if a==flag and i+1<len(args):return args[i+1]
        if a.startswith(flag+'='):return a.split('=',1)[1]


def output_for(command,args):
    if any(a in args for a in ['-h','--help','--dry-run','--list-features']):return None
    if command in ['clocks','assay-clocks'] and (not args or args[0]!='score'):return None
    if command=='epidiff-coverage' and args and args[0]=='make-map':return None
    if command in ['alphagenome','vep-annotate'] and '--execute' not in args:return None
    if command=='alphafold' and '--execute' not in args and option(args,'--cache') is None:return None
    value=option(args,'--out') if command in FLAG else args[POSITION[command]] if command in POSITION and len(args)>POSITION[command] else None
    return Path(value).resolve() if value else None


def report_location(output):return output/'visualization' if output.is_dir() else Path(str(output)+'.visualization')


def render(command,args,output):
    from report_profiles import build
    dest=report_location(output)
    # Preserve previous reports on resumed upstream runs. Never silently serve stale plots.
    if dest.exists():
        i=2
        while dest.with_name(dest.name+f'.{i}').exists():i+=1
        dest=dest.with_name(dest.name+f'.{i}')
    dest.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.'+dest.name+'.pending-',dir=dest.parent))
    try:
        build(command,args,output,stage)
        receipt=stage/'visualization.json';obj=json.loads(receipt.read_text())
        source_files=[Path(__file__),ROOT/'downstream/report_profiles.py',ROOT/'downstream/analysis_visualizations.py']
        if command=='bladder-example':source_files.append(ROOT/'downstream/bladder_visualizations.py')
        obj['report_source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
        receipt.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
        stage.rename(dest)
        print('Automatic report: '+str(dest/'report.html'),flush=True)
        return dest
    except BaseException as exc:
        (stage/'report_failure.json').write_text(json.dumps({'status':'failed','error':str(exc)},indent=2)+'\n')
        raise RuntimeError(f'Analysis outputs are preserved; report failed in {stage}: {exc}') from exc


def main():
    argv=sys.argv[1:];command=argv[0];args=argv[1:];output=output_for(command,args)
    if output is not None:
        try:import numpy,pandas,matplotlib
        except ImportError as exc:raise SystemExit('Automatic reports need NumPy, pandas and Matplotlib. Run bash setup.sh from this checkout, or install environments/report-requirements.txt. '+str(exc))
    env=dict(os.environ,KILNBIO_REPORT_WRAPPED=command)
    code=subprocess.call(['bash',str(ROOT/'run.sh'),*argv],env=env)
    if code:return code
    if output is not None:
        if not output.exists():raise RuntimeError('Completed command did not create declared output: '+str(output))
        render(command,args,output)
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except Exception as exc:print('ERROR: '+str(exc),file=sys.stderr);sys.exit(2)
