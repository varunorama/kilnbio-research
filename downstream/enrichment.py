"""Rank-based gene-set enrichment with explicit identifiers and a local GMT."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from analysis import need,rows,file_at,digest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('config');p.add_argument('output');a=p.parse_args()
    try:
        config=Path(a.config).resolve();cfg=json.loads(config.read_text(encoding='utf-8-sig'));out=Path(a.output).resolve()
        keys={'ranks','gmt','species','gene_id_type','rank_definition','annotation_reviewed','min_size','max_size','fdr'}
        need(isinstance(cfg,dict) and not set(cfg)-keys,'Unknown enrichment keys')
        need(not out.exists(),'Use a new output directory')
        need(cfg.get('species')=='human' and cfg.get('gene_id_type') and cfg.get('annotation_reviewed') is True,'Review species and exact matching gene identifier types')
        need(cfg.get('rank_definition') in ('wald_statistic','signed_sqrt_QLF','other_reviewed_signed_statistic'),'Declare a signed ranking statistic; include the full tested gene universe')
        ranks=file_at(cfg['ranks'],config.parent);gmt=file_at(cfg['gmt'],config.parent)
        h,data=rows(ranks,'\t');need(h==['gene_id','statistic'],'Ranks need gene_id,statistic columns')
        need(len({r['gene_id'] for r in data})==len(data) and all(r['gene_id'] for r in data),'One rank per nonempty gene ID; resolve transcript/gene mappings first')
        import math
        values=[float(r['statistic']) for r in data]
        need(all(math.isfinite(x) for x in values) and min(values)<0<max(values),'Finite signed ranks with both directions required')
        cfg.setdefault('min_size',15);cfg.setdefault('max_size',500);cfg.setdefault('fdr',.05)
        need(type(cfg['min_size']) is int and type(cfg['max_size']) is int and 1<=cfg['min_size']<=cfg['max_size'],'Invalid gene-set sizes')
        need(type(cfg['fdr']) in (int,float) and 0<cfg['fdr']<1,'Invalid FDR')
        cfg.update(ranks=str(ranks),gmt=str(gmt))
        rscript=os.environ.get('KILNBIO_RSCRIPT') or shutil.which('Rscript');need(rscript,'Rscript not installed')
        out.mkdir(parents=True);resolved=out/'config.resolved.json';resolved.write_text(json.dumps(cfg,indent=2)+'\n')
        record={'input_sha256':{str(f):digest(f) for f in (config,ranks,gmt)},'source_sha256':{f.name:digest(f) for f in (Path(__file__),Path(__file__).with_suffix('.R'),Path(__file__).with_name('visual_reports.py'))},'ties':len(values)-len(set(values)),'seed':20261031}
        (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
        with (out/'R.log').open('w') as log:r=subprocess.run([rscript,'--vanilla',str(Path(__file__).with_suffix('.R')),str(resolved),str(out)],stdout=log,stderr=subprocess.STDOUT)
        (out/'status.json').write_text(json.dumps({'complete':False,'exit_code':r.returncode,'stage':'statistics' if r.returncode else 'visualization'})+'\n')
        need(r.returncode==0,f'Enrichment failed; inspect {out / "R.log"}')
        from visual_reports import enrichment_report
        enrichment_report(out,cfg)
        (out/'status.json').write_text(json.dumps({'complete':True,'exit_code':0,'stage':'complete'})+'\n')
        print(f'Enrichment completed: {out}')
    except (ValueError,KeyError,OSError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()
