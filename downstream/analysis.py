"""Design-aware differential analysis. Validates before calling documented R methods."""
import argparse, csv, gzip, hashlib, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
SAFE = re.compile(r'[A-Za-z][A-Za-z0-9_]*\Z')

def need(ok, message):
    if not ok: raise ValueError(message)

def rows(path, delimiter=','):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        r = csv.DictReader(handle, delimiter=delimiter)
        need(r.fieldnames and len(r.fieldnames)==len(set(r.fieldnames)), 'Missing or duplicate headers')
        out = list(r)
        need(out and all(None not in x and all(v is not None for v in x.values()) for x in out), 'Empty/malformed table')
        return r.fieldnames, out

def file_at(value, base):
    need(isinstance(value,str) and value and '://' not in value, 'Use staged local input files')
    p = Path(value).expanduser()
    p = (base/p).resolve() if not p.is_absolute() else p.resolve()
    need(p.is_file() and p.stat().st_size>0, f'Missing or empty file: {p}')
    return p

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def rank(matrix):
    a = [list(map(float,r)) for r in matrix]; pivot = 0
    for col in range(len(a[0])):
        i = max(range(pivot,len(a)), key=lambda i:abs(a[i][col]), default=pivot)
        if pivot >= len(a) or abs(a[i][col]) < 1e-10: continue
        a[pivot],a[i] = a[i],a[pivot]
        scale = a[pivot][col]; a[pivot] = [v/scale for v in a[pivot]]
        for j in range(pivot+1,len(a)):
            s = a[j][col]; a[j] = [u-s*v for u,v in zip(a[j],a[pivot])]
        pivot += 1
        if pivot == len(a): break
    return pivot

def validate_design(cfg, base):
    table = file_at(cfg['samples'], base)
    header, samples = rows(table)
    group = cfg.get('group','condition')
    need(SAFE.fullmatch(group) and group not in ('sample','biological_id'), 'Invalid group column')
    need({'sample','biological_id',group} <= set(header), 'Samples need sample, biological_id and group columns')
    ids = [s['sample'] for s in samples]
    need(len(ids)==len(set(ids)) and all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s) for s in ids), 'Unique safe sample IDs required')
    need(all(s['biological_id'].strip() and s[group].strip() for s in samples), 'Biological IDs and groups cannot be missing')
    den, num = cfg['denominator'], cfg['numerator']
    need(isinstance(den,str) and isinstance(num,str) and den!=num, 'Different numerator and denominator required')
    groups = sorted({s[group] for s in samples})
    need(den in groups and num in groups, 'Contrast levels absent from samples')
    for level in groups:
        biological = [s['biological_id'] for s in samples if s[group]==level]
        need(len(biological)==len(set(biological)), 'Merge technical replicates; repeated donor/group rows are unsupported')
        need(len(set(biological))>=2, 'At least two independent biological replicates per group required; three or more recommended')
    cov = cfg.get('covariates',[])
    need(isinstance(cov,list) and all(set(c)=={'name','type'} and c['type'] in ('factor','numeric') for c in cov), 'Covariates need name/type declarations')
    names = [c['name'] for c in cov]
    need(len(names)==len(set(names)) and group not in names and all(SAFE.fullmatch(n) and n in header for n in names), 'Invalid/repeated covariate names')
    repeated = len({s['biological_id'] for s in samples}) < len(samples)
    need(not repeated or {'name':'biological_id','type':'factor'} in cov, 'Repeated donors require biological_id as a factor covariate')
    matrix = [[1.] for _ in samples]
    for term in cov + [{'name':group,'type':'factor'}]:
        values = [s[term['name']] for s in samples]
        need(all(v.strip() and v not in ('NA','NaN','nan') for v in values), 'Missing design values are not silently dropped')
        if term['type']=='numeric':
            values = [float(v) for v in values]
            need(all(math.isfinite(v) for v in values), 'Nonfinite numeric covariate')
            for row,v in zip(matrix,values): row.append(v)
        else:
            levels = sorted(set(values))
            if term['name']==group: levels = [den] + [x for x in levels if x!=den]
            need(len(levels)>=2, 'Remove constant covariates')
            for row,v in zip(matrix,values): row.extend(float(v==level) for level in levels[1:])
    ncol = len(matrix[0])
    need(rank(matrix)==ncol, 'Design is rank deficient: group and batch/donor/covariates may be confounded')
    need(len(samples)>ncol, 'No residual degrees of freedom')
    return table, samples, group, matrix

def check_rsem(files):
    # tximport's gene-level RSEM reader combines rows by position.
    reference_ids = None
    zero_length_genes = set()
    for path in files:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
            table = csv.DictReader(handle, delimiter='\t')
            required = {'gene_id', 'transcript_id(s)', 'expected_count', 'effective_length', 'TPM'}
            need(required <= set(table.fieldnames or []) and len(set(table.fieldnames)) == len(table.fieldnames), 'Expected per-sample RSEM genes.results')
            ids = []
            for row in table:
                need(None not in row and all(row[k] for k in required), 'Incomplete RSEM row')
                count, length, tpm = (float(row[k]) for k in ('expected_count', 'effective_length', 'TPM'))
                need(all(math.isfinite(v) and v >= 0 for v in (count, length, tpm)), 'Invalid RSEM estimate')
                # Capture unusable offsets before tximport floors RSEM lengths to one.
                if length == 0:
                    zero_length_genes.add(row['gene_id'])
                ids.append(row['gene_id'])
            need(ids and all(ids) and len(set(ids)) == len(ids), 'Empty or duplicate RSEM gene IDs')
            need(reference_ids is None or ids == reference_ids, 'RSEM gene IDs/order differ across samples; check the reference')
            reference_ids = ids
    return sorted(zero_length_genes)

def prepare(mode, config, output):
    path = Path(config).resolve(); cfg = json.loads(path.read_text(encoding='utf-8-sig'))
    need(isinstance(cfg,dict), 'Configuration must be an object')
    common={'samples','group','denominator','numerator','covariates','fdr','min_count','min_samples'}
    specific=({'assay','engine','count_source','counts','quant_manifest','tx2gene','feature_universe','normalization_reviewed','external_size_factors'} if mode=='de' else
              {'engine','modification','assembly','coordinates','measurement_reviewed','methylation_manifest','smoothing','min_coverage','dmr_p_threshold','dmr_min_length','dmr_min_cpg','dmr_merge_distance','dmr_fraction_significant','permutation_reviewed','max_permutations','region_cutoff'})
    need(not (set(cfg)-common-specific),f'Unknown configuration keys: {sorted(set(cfg)-common-specific)}')
    table, samples, group, matrix = validate_design(cfg,path.parent)
    cfg.update(samples=str(table),group=group,mode=mode)
    for key, default, low, high in [('fdr',.05,0,1),('min_count',10,0,float('inf')),('min_samples',2,1,len(samples))]:
        value = cfg.setdefault(key,default)
        need(type(value) in (int,float) and math.isfinite(value) and low<=value<=high, f'Invalid {key}')
    need(0<cfg['fdr']<1, 'fdr must be between 0 and 1')
    need(type(cfg['min_samples']) is int, 'min_samples must be an integer')
    inputs = {str(path),str(table)}
    if mode=='de':
        need(cfg.get('assay') in ('rnaseq','atacseq','chipseq','pseudobulk'), 'Declare assay')
        need(cfg.get('engine') in ('deseq2','edger'), 'Choose deseq2 or edger')
        source = cfg.get('count_source')
        need(source in ('raw_counts','salmon_tximport','rsem_tximport'), 'Use raw integer counts or a Salmon/RSEM tximport manifest')
        if 'external_size_factors' in cfg:
            need(source=='raw_counts','External size factors currently require raw counts; tximport length offsets cannot be replaced')
            need(cfg.get('normalization_reviewed') is True,'Review experimental reference assumptions before supplying external size factors')
            sf=file_at(cfg['external_size_factors'],path.parent)
            h,dat=rows(sf)
            need(h==['sample','size_factor'] and len(dat)==len(samples) and {r['sample'] for r in dat}=={s['sample'] for s in samples},'Size factors must match each sample exactly once')
            values=[float(r['size_factor']) for r in dat]
            need(all(math.isfinite(v) and v>0 for v in values),'Size factors must be finite and positive')
            cfg['external_size_factors']=str(sf);inputs.add(str(sf))
        if cfg['assay'] in ('atacseq','chipseq'):
            need(cfg.get('feature_universe') in ('external_fixed','label_blind_consensus'), 'Declare a fixed/label-blind chromatin feature universe')
            need(cfg.get('normalization_reviewed') is True, 'Review normalization assumptions and spike-ins/global shifts before chromatin testing')
            need(source=='raw_counts','Chromatin requires raw region counts')
        if source=='raw_counts':
            p = file_at(cfg['counts'],path.parent); cfg['counts']=str(p); inputs.add(str(p))
            header, data = rows(p,'\t')
            ids = [s['sample'] for s in samples]
            need(header[0]=='feature_id' and set(header[1:])==set(ids), 'Counts must have feature_id then exactly the sample IDs')
            seen=set()
            for r in data:
                need(r['feature_id'] and r['feature_id'] not in seen, 'Feature IDs must be unique and nonempty'); seen.add(r['feature_id'])
                values=[float(r[s]) for s in ids]
                need(all(math.isfinite(v) and 0<=v<=2147483647 and v.is_integer() for v in values), 'Counts must be finite nonnegative integers, not normalized or silently rounded')
        else:
            need(cfg['assay']=='rnaseq' and cfg['engine']=='deseq2', 'Estimated-count imports require RNA/DESeq2')
            column = 'genes_results' if source == 'rsem_tximport' else 'quant_sf'
            p = file_at(cfg['quant_manifest'],path.parent); header,data=rows(p)
            need(header==['sample',column] and len(data)==len(samples) and {r['sample'] for r in data}=={s['sample'] for s in samples}, f'Quant manifest must match samples exactly and have sample,{column}')
            cfg['quant_files']={r['sample']:str(file_at(r[column],p.parent)) for r in data}
            files = list(cfg['quant_files'].values())
            need(len(set(files))==len(files), 'Each sample requires a separate quantification file')
            inputs.update([str(p), *files])
            if source == 'rsem_tximport':
                need('tx2gene' not in cfg, 'RSEM genes.results already contains gene-level estimates')
                cfg['rsem_zero_length_genes'] = check_rsem(files)
            else:
                tx = file_at(cfg['tx2gene'],path.parent); h,t=rows(tx,'\t')
                need(h==['transcript_id','gene_id'] and len({r['transcript_id'] for r in t})==len(t), 'tx2gene needs unique transcript_id and gene_id columns')
                cfg['tx2gene']=str(tx); inputs.add(str(tx))
    else:
        need(cfg.get('engine') in ('dss','dmrseq'), 'Choose dss (loci/candidate regions) or dmrseq (region inference)')
        if cfg['engine']=='dmrseq':
            need(len({s[group] for s in samples})==2,'Initial dmrseq preset requires two groups')
            need(len({s['biological_id'] for s in samples})==len(samples),'Initial dmrseq preset requires independent subjects, not repeated measures')
            need(cfg.get('permutation_reviewed') is True,'Review exchangeability and permutation validity; covariate adjustment alone does not guarantee it')
            cfg.setdefault('max_permutations',100);cfg.setdefault('region_cutoff',.1)
            need(type(cfg['max_permutations']) is int and cfg['max_permutations']>=10,'At least ten requested permutations; assess attainable resolution')
            need(type(cfg['region_cutoff']) in (int,float) and 0<cfg['region_cutoff']<1,'Invalid region discovery cutoff')
        need(cfg.get('modification') in ('modC','5mC','5hmC'), 'Modification channel must be explicit')
        need(isinstance(cfg.get('assembly'),str) and cfg['assembly'].strip(), 'Explicit reference assembly required')
        need(cfg.get('coordinates')=='1-based-cytosine', 'Canonical methylation coordinates are 1-based cytosine positions, strand-collapsed CpGs')
        need(cfg.get('measurement_reviewed') is True, 'Review chemistry, conversion/calling errors and channel comparability before DM')
        p=file_at(cfg['methylation_manifest'],path.parent); h,data=rows(p)
        need(h==['sample','counts'] and len(data)==len(samples) and {r['sample'] for r in data}=={s['sample'] for s in samples}, 'Methylation manifest must match samples exactly')
        cfg['methylation_files']={r['sample']:str(file_at(r['counts'],p.parent)) for r in data}
        inputs.update([str(p),*cfg['methylation_files'].values()])
        for f in cfg['methylation_files'].values():
            sidecar=Path(f).parent/'status.json'
            if sidecar.exists():
                record=json.loads(sidecar.read_text(encoding='utf-8'))
                need(record.get('complete') is True,'Methylation importer did not complete')
                if 'config' in record:
                    need(record['config'].get('assembly')==cfg['assembly'] and record['config'].get('modification')==cfg['modification'],'Imported methylation assembly/channel disagrees with study config')
                inputs.add(str(sidecar))
        need(type(cfg.get('smoothing',False)) is bool, 'smoothing must be boolean')
        cfg.setdefault('smoothing',False)
        cfg.setdefault('min_coverage',5);cfg.setdefault('dmr_p_threshold',1e-4)
        cfg.setdefault('dmr_min_length',50);cfg.setdefault('dmr_min_cpg',3);cfg.setdefault('dmr_merge_distance',50)
        cfg.setdefault('dmr_fraction_significant',.5)
        need(type(cfg['min_coverage']) is int and cfg['min_coverage']>=1, 'min_coverage must be positive integer')
        need(0<cfg['dmr_p_threshold']<1 and 0<cfg['dmr_fraction_significant']<=1, 'Invalid DMR thresholds')
        for k in ('dmr_min_length','dmr_min_cpg','dmr_merge_distance'):
            need(type(cfg[k]) is int and cfg[k]>=1, 'Invalid DMR region parameter')
    output = Path(output).resolve()
    need(not output.exists(), 'Use a new output directory; statistical runs are immutable')
    hashes={p:digest(p) for p in sorted(inputs)}
    return cfg, {'mode':mode,'config':cfg,'input_sha256':hashes,'design_matrix':matrix,
                 'source_sha256':{p.name:digest(p) for p in (ROOT/'analysis.py',ROOT/'differential.R')},
                 'validation_scope':'Input/design validation, not biological study approval'}, output

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['de','dm'])
    p.add_argument('config');p.add_argument('output');p.add_argument('--dry-run',action='store_true')
    args=p.parse_args()
    try:
        cfg,plan,out=prepare(args.mode,args.config,args.output)
        if args.dry_run:
            print(json.dumps(plan,indent=2));return
        rscript=os.environ.get('KILNBIO_RSCRIPT') or shutil.which('Rscript')
        need(rscript,'Rscript not found; see environments/setup_r.R')
        (out/'provenance').mkdir(parents=True);(out/'results').mkdir()
        (out/'provenance/plan.json').write_text(json.dumps(plan,indent=2)+'\n')
        resolved=out/'provenance/config.json';resolved.write_text(json.dumps(cfg,indent=2)+'\n')
        with (out/'provenance/R.log').open('w') as log:
            result=subprocess.run([rscript,'--vanilla',str(ROOT/'differential.R'),str(resolved),str(out)],stdout=log,stderr=subprocess.STDOUT)
        (out/'provenance/status.json').write_text(json.dumps({'exit_code':result.returncode,'complete':result.returncode==0})+'\n')
        if result.returncode: raise RuntimeError(f'R analysis failed; inspect {out / "provenance/R.log"}')
        print(f'Completed {args.mode}: {out}; assess QC and model diagnostics before interpretation')
    except (ValueError,KeyError,OSError,RuntimeError,csv.Error) as e:p.exit(2,f'ERROR: {e}\n')

if __name__=='__main__':main()
