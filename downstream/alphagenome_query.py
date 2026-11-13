"""Optional, explicit AlphaGenome 2 variant scoring. Planning is offline.

No general-purpose LLM is called. Remote predictions are annotations, not
measurements, independent replication, p-values or clinical classifications.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

SDK_VERSION = '2.0.0'
SCORERS = {'ATAC', 'CAGE', 'CHIP_HISTONE', 'CHIP_TF', 'DNASE', 'METHYLATION',
           'RNA_SEQ', 'SC_ATAC', 'SC_GENE_EXPRESSION',
           'SC_GENE_EXPRESSION_DENOISED', 'SPLICE_SITES', 'SPLICE_SITE_USAGE',
           'SPLICE_JUNCTIONS', 'POLYADENYLATION', 'CONTACT_MAPS'}
LENGTHS = {2**14, 2**17, 2**19, 2**20}
# GRCh38 primary chromosomes; mitochondrial and alternate contigs are not routed.
CHROM_LENGTHS = dict(zip(['chr'+str(x) for x in range(1, 23)]+['chrX', 'chrY'],
    [248956422,242193529,198295559,190214555,181538259,170805979,159345973,
     145138636,138394717,133797422,135086622,133275309,114364328,107043718,
     101991189,90338345,83257441,80373285,58617616,64444167,46709983,
     50818468,156040895,57227415]))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(path)


def read_plan(config_path):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    allowed = {'variants','assembly','condition','cell_counts','cell_condition_name',
               'scorers','sequence_length','model_version','max_variants'}
    if not isinstance(config, dict) or set(config)-allowed:
        raise ValueError('Unknown configuration field; credentials do not belong in configuration files')
    if config.get('assembly') != 'GRCh38':
        raise ValueError('This adapter requires GRCh38; liftOver and allele normalization are upstream steps')
    length = config.get('sequence_length', 2**20)
    if type(length) is not int or length not in LENGTHS:
        raise ValueError('sequence_length must be 16384, 131072, 524288 or 1048576')
    checkpoint = config.get('model_version', 'ALL_FOLDS')
    if checkpoint not in {'ALL_FOLDS', 'FOLD_1'}:
        raise ValueError('model_version must be ALL_FOLDS or FOLD_1')
    scorers = config.get('scorers')
    if (not isinstance(scorers, list) or not scorers or
            any(not isinstance(s, str) or s not in SCORERS for s in scorers) or
            len(set(scorers)) != len(scorers)):
        raise ValueError('Supply a nonempty, unique list of supported scorer names')
    condition, cells = config.get('condition'), config.get('cell_counts')
    if bool(condition) == bool(cells):
        raise ValueError('Choose exactly one reference condition or cell_counts file')
    if condition and (not isinstance(condition, str) or
            not re.fullmatch(r'[A-Z][A-Z0-9]*:[0-9]+(?:_[A-Z][A-Z0-9]*:[0-9]+)*', condition)):
        raise ValueError('condition must be a provider-supported ontology CURIE or composite CURIE')
    if 'SC_GENE_EXPRESSION_DENOISED' in scorers and not cells:
        raise ValueError('SC_GENE_EXPRESSION_DENOISED requires raw single-cell counts')
    if cells and not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,63}', str(config.get('cell_condition_name',''))):
        raise ValueError('cell_condition_name must be a generic label, not a donor identifier')
    def resolve(name):
        p = Path(config[name]).expanduser()
        return (config_path.parent/p).resolve() if not p.is_absolute() else p.resolve()
    path = resolve('variants')
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t' if path.suffix.lower()=='.tsv' else ',')
        if reader.fieldnames != ['variant_id','chromosome','position','ref','alt']:
            raise ValueError('Variants need exactly: variant_id,chromosome,position,ref,alt; no donor data')
        variants = list(reader)
    cap = config.get('max_variants',100)
    if type(cap) is not int or not 1 <= cap <= 1000 or not 1 <= len(variants) <= cap:
        raise ValueError('Use 1..max_variants rows (default 100; hard limit 1000)')
    ids, alleles = set(), set()
    for row in variants:
        ident = row['variant_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', ident or '') or ident in ids:
            raise ValueError('variant_id must be unique and filename-safe')
        ids.add(ident)
        chrom, ref, alt = row['chromosome'], row['ref'], row['alt']
        try:
            pos = int(row['position'])
        except (TypeError, ValueError):
            raise ValueError('position must be a one-based integer') from None
        if (chrom not in CHROM_LENGTHS or pos < 1 or
                not re.fullmatch('[ACGT]+',ref or '') or not re.fullmatch('[ACGT]+',alt or '') or
                ref == alt or max(len(ref),len(alt)) > 1000):
            raise ValueError('Use normalized biallelic A/C/G/T alleles on primary GRCh38 chromosomes')
        key = (chrom,pos,ref,alt)
        if key in alleles:
            raise ValueError('Duplicate genomic variant')
        alleles.add(key)
        # SDK rounds the center up for odd-width unstranded intervals.
        start = (2*(pos-1)+len(ref)+1)//2-length//2
        end = start+length
        if start < 0 or end > CHROM_LENGTHS[chrom] or pos-1+len(ref) > end:
            raise ValueError('Prediction window crosses a chromosome boundary; select another supported length')
        row.update(position=pos, interval_start=start, interval_end=end)
    plan = dict(adapter='kilnbio-alphagenome-variant-scoring-v1', sdk_version=SDK_VERSION,
                model='AlphaGenome 2', model_version=checkpoint, assembly='GRCh38',
                sequence_length=length, scorers=scorers, condition=condition,
                variants=variants, input_sha256={'variants':sha(path)},
                code_sha256=sha(__file__), reference_alleles_checked_locally=False,
                external_service='gdmscience.googleapis.com:443',
                prediction_is_independent_validation=False,
                output_use='Provider API terms apply; not a foundation-model training export')
    if cells:
        cell_path = resolve('cell_counts')
        plan.update(cell_counts=str(cell_path), cell_condition_name=config['cell_condition_name'])
        plan['input_sha256']['cell_counts'] = sha(cell_path)
    canonical = json.dumps(plan,sort_keys=True,separators=(',',':')).encode()
    plan['plan_sha256'] = hashlib.sha256(canonical).hexdigest()
    return plan


def prepare_cells(path):
    """Validate without silent truncation, subsampling or metadata upload."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    a = ad.read_h5ad(path)
    if not 1 <= a.n_obs <= 256 or a.n_vars == 0:
        raise ValueError('Provide 1..256 explicitly selected cells with raw UMI counts in X')
    genes = [str(g).split('.')[0] for g in a.var_names]
    if len(set(genes)) != len(genes) or any(not re.fullmatch(r'ENSG\d{11}',g) for g in genes):
        raise ValueError('Use unique human Ensembl gene IDs; version-stripped IDs must remain unique')
    x = sparse.csr_matrix(a.X, copy=True)
    x.sum_duplicates()
    if (not np.isfinite(x.data).all() or (x.data < 0).any() or
            (x.data > np.iinfo(np.int32).max).any() or (x.data != np.floor(x.data)).any()):
        raise ValueError('X must contain finite, nonnegative raw integer UMI counts, not normalized values')
    x.eliminate_zeros()
    if (np.diff(x.indptr)==0).any():
        raise ValueError('Each selected cell must have positive counts')
    return ad.AnnData(X=x.astype(np.int32),
                      obs=pd.DataFrame(index=['cell_'+str(i) for i in range(a.n_obs)]),
                      var=pd.DataFrame(index=genes))


def execute(plan, out, allow_cell_upload=False, resume=False):
    # All validation precedes client creation. No client exists during planning.
    if plan.get('cell_counts') and not allow_cell_upload:
        raise ValueError('Expression conditioning transmits counts; explicitly add --allow-cell-upload')
    if importlib.metadata.version('alphagenome') != SDK_VERSION:
        raise ValueError('Install the pinned optional alphagenome=='+SDK_VERSION+' environment')
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers
    import pandas as pd
    cell_data = prepare_cells(plan['cell_counts']) if plan.get('cell_counts') else None
    api_key = os.environ.get('ALPHA_GENOME_API_KEY','')
    if not api_key:
        raise ValueError('Set ALPHA_GENOME_API_KEY locally; never put the key in a sample sheet or Git')
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    pp = out/'plan.json'
    if pp.exists():
        if not resume or json.loads(pp.read_text(encoding='utf-8')) != plan:
            raise ValueError('Output exists; --resume requires exactly the same inputs and adapter')
    elif any(out.iterdir()):
        raise ValueError('Use an empty output directory')
    save_json(pp,plan)
    receipt = dict(status='running', completed=[], failed=[], plan_sha256=plan['plan_sha256'],
                   started_utc=datetime.now(timezone.utc).isoformat(),
                   validation_scope='Remote predictions only; no biological validation claim')
    save_json(out/'receipt.json',receipt)
    try:
        client = None
        kwargs = ({'cell_conditioning':cell_data,'cell_condition_name':plan['cell_condition_name']}
                  if cell_data is not None else {'condition':plan['condition']})
        scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[s] for s in plan['scorers']]
        for v in plan['variants']:
            path = out/(v['variant_id']+'.scores.tsv')
            check = path.with_suffix('.sha256')
            if resume and path.exists() and check.exists() and check.read_text().strip()==sha(path):
                receipt['completed'].append(v['variant_id'])
                continue
            variant = genome.Variant(chromosome=v['chromosome'],position=v['position'],
                                     reference_bases=v['ref'],alternate_bases=v['alt'])
            interval = variant.reference_interval.resize(plan['sequence_length'])
            if (interval.start,interval.end) != (v['interval_start'],v['interval_end']):
                raise ValueError('SDK coordinate convention differs from the frozen plan')
            try:
                if client is None:
                    client = dna_client.create(api_key, model_version=dna_client.ModelVersion[plan['model_version']], timeout=30)
                scores = client.score_variant(interval=interval,variant=variant,
                                               variant_scorers=scorers,**kwargs)
                tidy = variant_scorers.tidy_scores([scores])
                if tidy.empty:
                    raise ValueError('No scores returned')
                tidy.insert(0,'kilnbio_variant_id',v['variant_id'])
                tidy.insert(1,'kilnbio_evidence_type','model_prediction')
                tmp = path.with_suffix('.tmp')
                tidy.to_csv(tmp,sep='\t',index=False)
                tmp.replace(path)
                check.write_text(sha(path)+'\n',encoding='ascii')
                receipt['completed'].append(v['variant_id'])
            except Exception as exc:
                # Avoid persisting credentials, request bodies or participant data in remote error text.
                receipt['failed'].append({'variant_id':v['variant_id'],'error_type':type(exc).__name__})
                break
            finally:
                save_json(out/'receipt.json',receipt)
        receipt['status'] = 'complete' if len(receipt['completed'])==len(plan['variants']) else 'incomplete'
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        if receipt['status']=='complete':
            frames=[pd.read_csv(out/(v['variant_id']+'.scores.tsv'),sep='\t') for v in plan['variants']]
            pd.concat(frames,ignore_index=True).to_csv(out/'scores.tsv',sep='\t',index=False)
        else:
            (out/'scores.tsv').unlink(missing_ok=True)
        save_json(out/'receipt.json',receipt)
        if receipt['status']!='complete':
            raise RuntimeError('Prediction run incomplete; see receipt.json. Valid completed queries can be resumed.')
    except Exception as exc:
        receipt['status']='incomplete'
        receipt['error_type']=type(exc).__name__
        (out/'scores.tsv').unlink(missing_ok=True)
        save_json(out/'receipt.json',receipt)
        raise RuntimeError('Prediction run incomplete; see receipt.json for completion and error type') from None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('config');p.add_argument('out')
    p.add_argument('--execute',action='store_true',help='Explicitly send the configured variants to AlphaGenome')
    p.add_argument('--allow-cell-upload',action='store_true',help='Also allow sending the selected raw expression counts')
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    try:
        plan=read_plan(args.config)
        if args.execute:
            execute(plan,args.out,args.allow_cell_upload,args.resume)
            print('Completed model annotations. These are not independent assay observations.')
        else:
            if args.resume or args.allow_cell_upload:
                raise ValueError('--resume and --allow-cell-upload require --execute')
            save_json(Path(args.out).with_suffix('.plan.json'),plan)
            print('Offline plan written; no API call made. Review the plan before adding --execute.')
    except Exception as exc:
        # Known local errors are useful; remote exceptions may contain request details.
        msg=str(exc) if isinstance(exc,(ValueError,FileNotFoundError,RuntimeError,KeyError)) else type(exc).__name__
        print('ERROR: '+msg,file=sys.stderr)
        return 2
    return 0


if __name__=='__main__':
    raise SystemExit(main())
