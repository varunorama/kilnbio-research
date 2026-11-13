"""Plan or explicitly run local Ensembl VEP 116.3 for research annotation."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from variant_annotations import digest

VERSION = '116.3'


def prepare(config, output):
    cp = Path(config).resolve()
    cfg = json.loads(cp.read_text(encoding='utf-8'))
    if set(cfg) != {'input_vcf','assembly','fasta','cache_dir','forks'}:
        raise ValueError('Provide input_vcf, assembly, fasta, cache_dir, forks')
    if cfg['assembly'] not in {'GRCh37','GRCh38'} or type(cfg['forks']) is not int or not 1 <= cfg['forks'] <= 32:
        raise ValueError('Use GRCh37/GRCh38 and 1..32 forks')
    def path(value):
        return (cp.parent/Path(value).expanduser()).resolve()
    vcf, fasta, cache = [path(cfg[k]) for k in ('input_vcf','fasta','cache_dir')]
    info = cache/'homo_sapiens'/('116_'+cfg['assembly'])/'info.txt'
    for p in (vcf,fasta,info):
        if not p.is_file():
            raise ValueError('Required input/cache metadata missing: '+str(p))
    out = Path(output).resolve()
    argv = ['vep','--offline','--cache','--cache_version','116','--dir_cache',str(cache),
            '--species','homo_sapiens','--assembly',cfg['assembly'],'--fasta',str(fasta),
            '--check_ref','--format','vcf','--vcf','--input_file',str(vcf),
            '--output_file',str(out/'annotated.vcf'),'--warning_file',str(out/'warnings.txt'),
            '--symbol','--protein','--uniprot','--canonical','--mane','--hgvs','--allele_number',
            '--no_stats','--fork',str(cfg['forks'])]
    return {'status':'planned','research_use_only':True,'vep_required_version':VERSION,
            'assembly':cfg['assembly'],'argv':argv,
            'input_sha256':digest(vcf),'fasta_sha256':digest(fasta),
            'cache_info_sha256':digest(info),'cache_dir':str(cache),
            'cache_hash_scope':'info.txt only; archive the exact downloaded cache and checksum separately',
            'all_transcripts_retained':True,'network_annotation':False}


def run(config, output, execute=False):
    plan = prepare(config, output)
    out = Path(output);out.mkdir(parents=True,exist_ok=False)
    def save():
        (out/'receipt.json').write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    save()
    if not execute:
        return plan
    try:
        help_result = subprocess.run(['vep','--help'],capture_output=True,text=True,check=True)
        version_text = help_result.stdout+'\n'+help_result.stderr
        (out/'vep-version.txt').write_text(version_text,encoding='utf-8')
        if not re.search(r'ensembl-vep\s*[:=]?\s*116\.3\b', version_text, re.I):
            raise ValueError('Installed VEP must report ensembl-vep 116.3 in --help')
        plan['status'] = 'running';save()
        with (out/'execution.log').open('w',encoding='utf-8') as log:
            subprocess.run(plan['argv'],stdout=log,stderr=subprocess.STDOUT,check=True)
        if not (out/'annotated.vcf').is_file():
            raise ValueError('VEP did not create annotated.vcf')
        plan.update(status='complete',output_sha256=digest(out/'annotated.vcf'))
        save()
    except Exception as exc:
        plan.update(status='failed',error=str(exc));save();raise
    return plan


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('config');p.add_argument('output');p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    try:
        print(json.dumps(run(a.config,a.output,a.execute),indent=2))
    except (OSError,ValueError,subprocess.CalledProcessError) as exc:
        sys.exit(str(exc))
