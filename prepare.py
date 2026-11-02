"""Validate native nf-core samplesheets and build an immutable per-run plan.

Standard library only. This validates the toolkit's supported subset; nf-core
performs complete upstream schema validation at execution time.
"""
import argparse, csv, hashlib, json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REFERENCE_FILES = ('fasta','gtf','gff','blacklist','transcript_fasta','barcode_whitelist','additional_fasta')
REFERENCE_INDEXES = ('star_index','rsem_index','salmon_index','bowtie2_index','bismark_index')
WORKFLOWS = {'wgbs':'methylseq','rrbs':'methylseq','atacseq':'atacseq','chipseq':'chipseq','rnaseq':'rnaseq','scrnaseq':'scrnaseq'}
HEADERS = {
    'wgbs':['sample','fastq_1','fastq_2','genome'],
    'rrbs':['sample','fastq_1','fastq_2','genome'],
    'atacseq':['sample','fastq_1','fastq_2','replicate'],
    'chipseq':['sample','fastq_1','fastq_2','replicate','antibody','control','control_replicate'],
    'rnaseq':['sample','fastq_1','fastq_2','strandedness'],
    'scrnaseq':['sample','fastq_1','fastq_2'],
}

def require(ok, message):
    if not ok: raise ValueError(message)

def local_file(value, base):
    require(isinstance(value,str) and value, 'Expected a nonempty local file path')
    require('://' not in value, 'Stage remote files locally before running')
    require(not re.match(r'^[A-Za-z]:',value) or sys.platform=='win32', 'Use Linux paths under WSL, not Windows drive paths')
    path=Path(value).expanduser()
    path=(base/path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_file() and path.stat().st_size>0,f'Missing or empty file: {path}')
    require(not any(c.isspace() for c in str(path)),f'Upstream workflows require paths without whitespace: {path}. Stage inputs in a Linux directory without spaces.')
    return str(path)

def properties(schema):
    result=dict(schema.get('properties',{}))
    for key in ('$defs','definitions'):
        for section in schema.get(key,{}).values(): result.update(properties(section))
    return result

def local_index(value, base):
    require(isinstance(value,str) and value and '://' not in value, 'Stage reference indexes locally')
    path=Path(value).expanduser()
    path=(base/path).resolve() if not path.is_absolute() else path.resolve()
    require(not any(c.isspace() for c in str(path)),f'Index path contains whitespace: {path}')
    if path.is_file():
        require(path.name.endswith('.tar.gz'),'Index archive must end in .tar.gz')
        return local_file(value,base)
    require(path.is_dir() and any(p.is_file() and p.stat().st_size for p in path.rglob('*')),f'Missing or empty index: {path}')
    return str(path)

def index_inventory(path):
    """Hash prebuilt index contents so same-size edits cannot silently resume."""
    root=Path(path)
    files=sorted(p for p in root.rglob('*') if p.is_file()) if root.is_dir() else [root]
    inventory={}
    for p in files:
        h=hashlib.sha256()
        with p.open('rb') as src:
            for block in iter(lambda:src.read(8*1024**2),b''): h.update(block)
        key=str(p.relative_to(root)) if root.is_dir() else p.name
        inventory[key]={'bytes':p.stat().st_size,'sha256':h.hexdigest()}
    return inventory

def validate_samples(assay, sheet, params):
    with sheet.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f); header=reader.fieldnames
        require(header is not None and len(header)==len(set(header)), 'Missing or duplicate CSV headers')
        require(header[:len(HEADERS[assay])]==HEADERS[assay],f'Expected first columns: {HEADERS[assay]}')
        rows=list(reader)
    require(rows,'Samplesheet is empty')
    seen_files=set(); libraries={}; groups={}
    for i,row in enumerate(rows,2):
        require(None not in row and all(v is not None for v in row.values()),f'CSV row {i} has wrong column count')
        require(all(v==v.strip() for v in row.values()),f'Whitespace around values at row {i}')
        name=row['sample'];require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',name),f'Unsafe sample name at row {i}')
        if assay=='scrnaseq': require(row['fastq_2'],'10x GEX requires barcode/UMI R1 and cDNA R2')
        for field in ('fastq_1','fastq_2'):
            if field=='fastq_2' and not row[field]: continue
            path=local_file(row[field],sheet.parent)
            require(path.endswith(('.fq.gz','.fastq.gz')),f'Expected gzipped FASTQ: {path}')
            require(path not in seen_files,f'FASTQ assigned more than once: {path}')
            seen_files.add(path);row[field]=path
        if assay in ('wgbs','rrbs'): require(not row['genome'],'Leave sample genome column blank: one explicit reference per run')
        rep=row.get('replicate','1')
        require(rep.isdigit() and int(rep)>=1,f'Invalid replicate at row {i}')
        groups.setdefault(name,set()).add(int(rep))
        signature=(bool(row['fastq_2']),row.get('antibody',''),row.get('control',''),row.get('control_replicate',''),row.get('strandedness',''))
        key=(name,rep)
        require(key not in libraries or libraries[key]==signature,f'Inconsistent lane metadata for {key}')
        libraries[key]=signature
        if assay=='rnaseq': require(row['strandedness'] in ('auto','forward','reverse','unstranded'),f'Invalid RNA strandedness at row {i}')
    if assay in ('atacseq','chipseq'):
        for name,reps in groups.items(): require(reps==set(range(1,max(reps)+1)),f'Replicates for {name} must start at 1 and be contiguous')
        for row in rows:
            control=row.get('control',''); rep=row.get('control_replicate','')
            if assay=='chipseq':
                require(bool(row['antibody'])==bool(control),f'ChIP antibody and matched control required together: {row["sample"]}')
            if control:
                key=(control,rep)
                require(key in libraries,f'Unknown control library: {key}')
                require(control!=row['sample'], 'An IP cannot be its own control')
                require(not libraries[key][1] and not libraries[key][2],f'Control {key} must be an input/IgG library, not another IP')
                require(bool(row['fastq_2'])==libraries[key][0],f'Control layout differs: {key}')
            else: require(not rep,'control_replicate provided without control')
        if assay=='chipseq': require(any(r['antibody'] for r in rows),'No ChIP IP libraries')
        if assay=='atacseq': require(not any(r.get('control') for r in rows),'ATAC controls need a separately reviewed preset; unsupported in initial bulk preset')
    return header,rows

def prepare(args):
    require(args.assay in WORKFLOWS,'Unsupported preprocessing assay')
    require(1<=args.cpus<=256 and 2<=args.memory_gb<=4096,'Invalid resource limits')
    if Path('/proc/meminfo').exists():
        physical_kib=int(re.search(r'MemTotal:\s+(\d+)',Path('/proc/meminfo').read_text()).group(1))
        require(args.memory_gb <= physical_kib/1048576-2,'Memory limit leaves less than 2 GiB for Nextflow/OS; lower --memory-gb or use a larger Linux runtime')
    sheet=Path(args.samples).resolve(); paramfile=Path(args.params).resolve(); run=Path(args.run).resolve()
    require(not any(c.isspace() for c in str(run)), 'Run directory must have no whitespace; use a Linux workspace, e.g. ~/analyses/client/run1')
    params=json.loads(paramfile.read_text(encoding='utf-8-sig'));require(isinstance(params,dict),'Parameters must be a JSON object')
    pipeline=WORKFLOWS[args.assay]
    version=json.loads((ROOT/'versions.json').read_text())[pipeline]
    schema_path=ROOT/'schemas'/f'{pipeline}.json'
    require(hashlib.sha256(schema_path.read_bytes()).hexdigest()==version['schema_sha256'],'Vendored parameter schema checksum mismatch')
    schema=properties(json.loads(schema_path.read_text()))
    for key,value in params.items():
        require(key in schema,f'Unknown upstream parameter: {key}')
        spec=schema[key]; kind=spec.get('type')
        if kind=='boolean': require(type(value) is bool,f'{key} must be boolean')
        if kind=='string': require(isinstance(value,str),f'{key} must be a string')
        if kind in ('integer','number'): require(type(value) in (int,float) and (kind!='integer' or type(value) is int),f'{key} must be {kind}')
        if 'enum' in spec: require(value in spec['enum'],f'Invalid {key}: {value}')
    for key in ('input','outdir','genome','max_cpus','max_memory','max_time'):
        require(key not in params,f'{key} is managed by launcher or prohibited: use local references and CLI limits')
    require(params.get('fasta'),'Explicit local fasta required')
    require(params.get('igenomes_ignore') is True,'Set igenomes_ignore:true to avoid implicit reference downloads')
    for key in REFERENCE_FILES:
        if params.get(key): params[key]=local_file(params[key],paramfile.parent)
    for key in REFERENCE_INDEXES:
        if params.get(key): params[key]=local_index(params[key],paramfile.parent)
    if args.assay in ('wgbs','rrbs'):
        require(params.get('aligner') in ('bismark','bwameth'),'Choose bismark or bwameth explicitly; comparative benchmark required before a production default')
        require(type(params.get('non_directional')) is bool,'Declare non_directional from library protocol')
        require(not (params['aligner']=='bwameth' and params['non_directional']),'Non-directional bwa-meth preset is not supported')
        require(params.get('rrbs',False)==(args.assay=='rrbs'),'RRBS preset must set rrbs:true; WGBS must not')
        if args.assay=='rrbs':
            require(params['aligner']=='bismark','Initial conventional RRBS preset supports Bismark only')
            require(params.get('skip_deduplication') is True,'Conventional RRBS preset retains positional duplicates; set skip_deduplication:true')
        params.setdefault('save_reference',True)
    if args.assay in ('atacseq','chipseq','rnaseq'):
        require(params.get('gtf'),'Explicit assembly-matched GTF required')
        require(params.get('skip_deseq2_qc') is True,'Set skip_deseq2_qc:true; downstream study-specific modeling is separate')
    if args.assay in ('atacseq','chipseq'):
        require(params.get('blacklist'),'Assembly-matched blacklist required for this preset')
        require(type(params.get('macs_gsize')) in (int,float) and params['macs_gsize']>0,'Positive numeric macs_gsize required')
        require(type(params.get('narrow_peak')) is bool,'Choose narrow_peak explicitly for the assay/target')
    if args.assay=='scrnaseq':
        require(params.get('aligner')=='star','Initial 10x raw-read preset supports STARsolo; other engines need separate validation')
        require(params.get('protocol') in ('10XV2','10XV3','10XV4'),'Declare 10x chemistry explicitly')
        require(params.get('gtf') and params.get('barcode_whitelist'),'Explicit matched GTF and chemistry barcode whitelist required')
        require(params.get('star_feature') in ('Gene','GeneFull','Gene Velocyto'),'Declare Gene, GeneFull or Gene Velocyto')
        require(params.get('skip_cellbender') is True,'Initial 10x preset requires skip_cellbender:true; ambient correction needs a separately reviewed branch')
    params.setdefault('save_align_intermeds',True)
    params.setdefault('save_reference',True)
    header,rows=validate_samples(args.assay,sheet,params)
    out=run/'provenance'
    params.update(input=str(out/'samples.resolved.csv'),outdir=str(run/'results'))
    input_files={r[k] for r in rows for k in ('fastq_1','fastq_2') if r[k]}
    input_files.update(params[k] for k in REFERENCE_FILES if params.get(k))
    stats={v:dict(bytes=Path(v).stat().st_size,mtime_ns=Path(v).stat().st_mtime_ns) for v in sorted(input_files)}
    code_hash=hashlib.sha256(Path(__file__).read_bytes()+(ROOT/'run.sh').read_bytes()).hexdigest()
    identity=dict(assay=args.assay,workflow=pipeline,release=version['version'],engine_required='25.10.4',toolkit_sha256=code_hash,params=params,rows=rows,profile=args.profile,cpus=args.cpus,memory_gb=args.memory_gb,input_file_stats=stats)
    identity['reference_index_contents']={k:index_inventory(params[k]) for k in REFERENCE_INDEXES if params.get(k)}
    fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    if out.exists():
        require(args.resume,'Run already exists; use a new directory or --resume')
        require((out/'plan.sha256').read_text().strip()==fingerprint,'Resume rejected: inputs, parameters or resources changed')
    else:
        require(not run.exists() or not any(run.iterdir()),'Run directory is not empty')
        out.mkdir(parents=True)
    with (out/'samples.resolved.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,header);w.writeheader();w.writerows(rows)
    (out/'params.resolved.json').write_text(json.dumps(params,indent=2)+'\n')
    (out/'workflow.tsv').write_text(f'{pipeline}\t{version["version"]}\n')
    (out/'resources.config').write_text(f"process {{\n executor = 'local'\n resourceLimits = [cpus: {args.cpus}, memory: {args.memory_gb}.GB, time: 48.h]\n}}\nexecutor.queueSize = 1\nexecutor.cpus = {args.cpus}\nexecutor.memory = {args.memory_gb}.GB\n")
    (out/'plan.sha256').write_text(fingerprint+'\n')
    (out/'plan.json').write_text(json.dumps(identity,indent=2)+'\n')
    (out/'samples.original.csv').write_bytes(sheet.read_bytes())
    (out/'params.original.json').write_bytes(paramfile.read_bytes())
    print(f'Validated {len(rows)} FASTQ rows for {args.assay}; pinned nf-core/{pipeline} {version["version"]}')
    return params

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('assay');p.add_argument('samples');p.add_argument('params');p.add_argument('run')
    p.add_argument('--cpus',type=int,default=4);p.add_argument('--memory-gb',type=int,default=24)
    p.add_argument('--profile',choices=('docker','apptainer','singularity'),default='docker');p.add_argument('--resume',action='store_true')
    try: prepare(p.parse_args())
    except (ValueError,OSError,KeyError,csv.Error) as e: p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__': main()
