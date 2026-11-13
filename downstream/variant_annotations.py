"""Offline, allele-exact research annotation. Never assigns a clinical diagnosis."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from datetime import datetime, timezone


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def textfile(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8-sig')


def allele(chrom, pos, ref, alt):
    chrom = str(chrom).removeprefix('chr')
    chrom = 'MT' if chrom == 'M' else chrom
    if chrom not in {str(x) for x in range(1, 23)} | {'X', 'Y', 'MT'}:
        raise ValueError('Only primary human chromosomes are supported')
    if not re.fullmatch(r'[1-9][0-9]*', str(pos)):
        raise ValueError('Position must be a positive 1-based integer')
    if not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt) or ref == alt:
        raise ValueError('Use normalized biallelic A/C/G/T variants; no symbolic alleles')
    return chrom, int(pos), ref, alt


def queries(path):
    with textfile(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        if r.fieldnames != ['variant_id', 'chromosome', 'position', 'ref', 'alt']:
            raise ValueError('TSV needs exactly variant_id,chromosome,position,ref,alt')
        rows = list(r)
    if not rows:
        raise ValueError('No query variants')
    seen, keys = set(), {}
    for row in rows:
        ident = row['variant_id']
        if not ident or ident in seen or any(v is None for v in row.values()):
            raise ValueError('Missing fields or duplicate variant_id')
        k = allele(row['chromosome'], row['position'], row['ref'], row['alt'])
        if k in keys:
            raise ValueError('Duplicate allele; provide one row per unique variant, without donor identifiers')
        seen.add(ident)
        keys[k] = ident
    return rows, keys


def vcf_matches(path, wanted, assembly, kind):
    headers, csq = [], None
    with textfile(path) as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.rstrip())
                if line.startswith('##INFO=<ID=CSQ,'):
                    m = re.search(r'Format: ([^\"]+)', line)
                    if m:
                        csq = m.group(1).split('|')
                if line.startswith('##reference=') or line.startswith('##assembly='):
                    other = 'GRCh37' if assembly == 'GRCh38' else 'GRCh38'
                    if other in line:
                        raise ValueError('Resource header contradicts declared assembly')
                continue
            a = line.rstrip('\r\n').split('\t')
            if len(a) < 8:
                raise ValueError('Malformed VCF row')
            # A split record is required to attribute Number=A/R and CSQ correctly.
            if ',' in a[4]:
                if any((a[0].removeprefix('chr'), int(a[1]), a[3], alt) in wanted for alt in a[4].split(',')):
                    raise ValueError('Matched multiallelic VCF: normalize/split with bcftools first')
                continue
            try:
                key = allele(a[0], a[1], a[3], a[4])
            except ValueError:
                continue  # Structural, nonprimary and unspecified alleles are outside this interface.
            if key not in wanted:
                continue
            info = dict(x.split('=', 1) if '=' in x else (x, True) for x in a[7].split(';') if x != '.')
            record = {'vcf_id': a[2], 'filter': a[6], 'info': info}
            if kind == 'vep':
                if not csq or 'CSQ' not in info:
                    raise ValueError('Matched VEP record lacks CSQ annotation/header')
                annotations = []
                for c in info['CSQ'].split(','):
                    fields = c.split('|')
                    if len(fields) != len(csq):
                        raise ValueError('VEP CSQ field count differs from header')
                    annotations.append(dict(zip(csq, fields)))
                record['transcript_consequences'] = annotations
            yield wanted[key], record


def missense_matches(path, wanted, assembly):
    columns = None
    with textfile(path) as f:
        for line in f:
            if line.startswith('#CHROM\t'):
                columns = line.rstrip().lstrip('#').split('\t')
                required = {'CHROM','POS','REF','ALT','genome','uniprot_id','transcript_id','protein_variant','am_pathogenicity','am_class'}
                if not required.issubset(columns):
                    raise ValueError('AlphaMissense columns missing')
                continue
            if line.startswith('#') or not line.strip():
                continue
            if columns is None:
                raise ValueError('AlphaMissense #CHROM header missing')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != len(columns):
                raise ValueError('Malformed AlphaMissense row')
            row = dict(zip(columns, fields))
            if row['genome'] != {'GRCh38':'hg38', 'GRCh37':'hg19'}[assembly]:
                raise ValueError('AlphaMissense assembly mismatch')
            key = allele(row['CHROM'], row['POS'], row['REF'], row['ALT'])
            if key not in wanted:
                continue
            score = float(row['am_pathogenicity'])
            if not math.isfinite(score) or not 0 <= score <= 1 or len(row['REF']) != 1 or len(row['ALT']) != 1:
                raise ValueError('Invalid AlphaMissense missense-SNV score')
            yield wanted[key], row


def write_tsv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def annotate(config_path, output):
    cp = Path(config_path).resolve()
    cfg = json.loads(cp.read_text(encoding='utf-8'))
    if set(cfg) - {'assembly','variants','normalized_biallelic','resources'}:
        raise ValueError('Unknown annotation configuration')
    assembly = cfg.get('assembly')
    if assembly not in {'GRCh37','GRCh38'} or cfg.get('normalized_biallelic') is not True:
        raise ValueError('Declare GRCh37/GRCh38 and normalized_biallelic=true; no liftover or normalization is performed')
    def path(s):
        p = Path(s).expanduser()
        return (cp.parent/p).resolve()
    qp = path(cfg['variants'])
    rows, wanted = queries(qp)
    resources = cfg.get('resources', {})
    if not resources or set(resources) - {'clinvar','alphamissense','vep'}:
        raise ValueError('Select one or more of clinvar, alphamissense, vep')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    receipt = {'status':'running', 'research_use_only':True, 'assembly':assembly,
               'recorded_utc':datetime.now(timezone.utc).isoformat(), 'query_sha256':digest(qp),
               'query_variants':len(rows), 'resources':{},
               'absence_is_not_benign':True, 'clinical_classification_performed':False}
    def save():
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    save()
    try:
        for kind, spec in resources.items():
            if set(spec) != {'path','version','source_url','assembly'} or spec['assembly'] != assembly:
                raise ValueError('Each resource requires path, version, source_url and matching assembly')
            if not spec['version'] or not re.match(r'https://[^/]+/', spec['source_url']):
                raise ValueError('Provide resource version/date and HTTPS provenance URL')
            p = path(spec['path'])
            counts = dict.fromkeys(wanted.values(), 0)
            iterator = missense_matches(p, wanted, assembly) if kind == 'alphamissense' else vcf_matches(p, wanted, assembly, kind)
            with (out/(kind+'.jsonl')).open('w',encoding='utf-8') as f:
                for ident, annotation in iterator:
                    counts[ident] += 1
                    f.write(json.dumps({'variant_id':ident,'source':kind,'annotation':annotation},allow_nan=False)+'\n')
            for row in rows:
                n = counts[row['variant_id']]
                row[kind+'_records'] = n
                row[kind+'_status'] = 'not_found_in_resource' if n == 0 else 'matched_multiple_records' if n > 1 else 'matched'
            receipt['resources'][kind] = {**spec, 'path':str(p), 'sha256':digest(p),
                'matched_variants':sum(n > 0 for n in counts.values()), 'annotation_records':sum(counts.values())}
            save()
        write_tsv(out/'variants.tsv', rows, list(rows[0]))
        receipt['status'] = 'complete'
        save()
    except Exception as exc:
        receipt.update(status='failed', error=str(exc))
        save()
        raise
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config'); parser.add_argument('output')
    args = parser.parse_args()
    try:
        print(json.dumps(annotate(args.config, args.output), indent=2))
    except (OSError, ValueError, KeyError) as exc:
        sys.exit(str(exc))
