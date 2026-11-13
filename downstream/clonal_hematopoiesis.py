"""Research interfaces for blood clonal hematopoiesis (CH), distinct from ChIP-seq.

screen consumes biallelic somatic VCF or an annotated TSV, an explicit curated
driver table and a specimen roster. import audits existing provider calls.
Neither path establishes a clinical CHIP diagnosis. See docs/EXOME_AND_CH.md.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path


def read_table(path, required):
    with Path(path).open(encoding='utf-8-sig',newline='') as handle:
        delimiter = ',' if str(path).lower().endswith('.csv') else '\t'
        reader=csv.DictReader(handle,delimiter=delimiter)
        if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames) or not set(required)<=set(reader.fieldnames):
            raise ValueError('Required columns: '+','.join(required))
        rows=list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('Malformed table row')
    return rows


def number(value):
    if value in ('','NA','.','nan',None):
        return None
    result=float(value)
    if not math.isfinite(result):
        raise ValueError('Nonfinite numeric value')
    return result


def count(value):
    v=number(value)
    if v is not None and (v<0 or v!=int(v)):
        raise ValueError('Read depth must be a nonnegative integer')
    return None if v is None else int(v)


def allele_key(row):
    chrom=row['chrom']; pos=int(row['pos']); ref=row['ref']; alt=row['alt']
    if pos<1 or not chrom or ref==alt or any(c not in 'ACGT' for c in ref+alt) or not ref or not alt:
        raise ValueError('Need normalized biallelic A/C/G/T variants')
    return chrom,pos,ref,alt


def vcf_rows(path, af_field, assembly):
    opener=gzip.open if str(path).endswith(('.gz','.bgz')) else open
    header=None
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            if line.lower().startswith('##reference='):
                aliases={'grch38':'38','hg38':'38','grch37':'37','hg19':'37'}
                expected=aliases.get(assembly.lower())
                declared={aliases[x.lower()] for x in re.findall(r'GRCh38|hg38|GRCh37|hg19',line,re.I)}
                if expected and declared and declared!={expected}:
                    raise ValueError('VCF reference header contradicts the declared assembly')
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header=line.rstrip().split('\t')
                if len(header)<10 or len(set(header[9:]))!=len(header[9:]):
                    raise ValueError('VCF needs unique sample columns')
                continue
            if line.startswith('#') or not line.strip():
                continue
            if header is None:
                raise ValueError('VCF sample header is missing')
            fields=line.rstrip().split('\t')
            if len(fields)!=len(header) or ',' in fields[4]:
                raise ValueError('Split and normalize multiallelic records with bcftools before screening')
            info=dict(s.split('=',1) for s in fields[7].split(';') if '=' in s)
            for sample,text in zip(header[9:],fields[9:]):
                fmt=dict(zip(fields[8].split(':'),text.split(':')))
                def alternate(name):
                    values=fmt.get(name,'.').split(',')
                    return values[1] if len(values)==2 else ''
                # ADF/ADR count read strands. F1R2/F2R1 describe pair orientation.
                strand='read_strand' if 'ADF' in fmt and 'ADR' in fmt else 'pair_orientation'
                forward=alternate('ADF' if strand=='read_strand' else 'F1R2')
                reverse=alternate('ADR' if strand=='read_strand' else 'F2R1')
                yield dict(sample=sample,assembly=assembly,chrom=fields[0],pos=fields[1],ref=fields[3],alt=fields[4],
                           filter=fields[6],dp=fmt.get('DP',''),alt_count=alternate('AD'),vaf=fmt.get('AF',''),
                           alt_forward=forward,alt_reverse=reverse,orientation_kind=strand,
                           population_af=info.get(af_field,''))
    if header is None:
        raise ValueError('No VCF header')


def quality_reasons(row,args):
    reasons=[]
    if row.get('filter')!='PASS':
        reasons.append('not_PASS')
    values={k:count(row.get(k)) for k in ('dp','alt_count','alt_forward','alt_reverse')}
    for key,minimum in [('dp',args.min_depth),('alt_count',args.min_alt_reads),('alt_forward',1),('alt_reverse',1)]:
        if values[key] is None:
            reasons.append('missing_'+key)
        elif values[key]<minimum:
            reasons.append('low_'+key)
    dp,alt=values['dp'],values['alt_count']
    if dp is not None and alt is not None and alt>dp:
        raise ValueError('Alternate reads exceed total depth')
    if alt is not None and any(values[k] is not None and values[k]>alt for k in ('alt_forward','alt_reverse')):
        raise ValueError('Orientation support exceeds alternate read count')
    vaf=number(row.get('vaf'))
    if vaf is None and dp and alt is not None:
        vaf=alt/dp
    if vaf is None:
        reasons.append('missing_vaf')
    elif not 0<=vaf<=1:
        raise ValueError('VAF outside [0,1]')
    elif vaf<args.min_vaf:
        reasons.append('low_vaf')
    af=number(row.get('population_af'))
    if af is None:
        reasons.append('population_frequency_unknown')
    elif not 0<=af<=1:
        raise ValueError('Population AF outside [0,1]')
    elif af>args.max_population_af:
        reasons.append('population_frequency_high')
    return reasons,vaf


def screen(args, roster):
    drivers=read_table(args.drivers,['assembly','chrom','pos','ref','alt','gene','evidence'])
    driver_map={}
    for row in drivers:
        key=allele_key(row)
        if row['assembly']!=args.assembly or key in driver_map or not row['gene'] or not row['evidence']:
            raise ValueError('Driver list requires unique assembly-matched alleles, gene and public evidence')
        driver_map[key]=row
    source=vcf_rows(args.variants,args.population_af_field,args.assembly) if str(args.variants).endswith(('.vcf','.vcf.gz','.vcf.bgz')) else read_table(args.variants,['sample','assembly','chrom','pos','ref','alt','filter','dp','alt_count','vaf','alt_forward','alt_reverse','population_af','orientation_kind'])
    requested={r['sample'] for r in roster}; audit=[]; seen=set()
    for row in source:
        if row['sample'] not in requested:
            continue
        key=allele_key(row); record=(row['sample'],key)
        if record in seen or row['assembly']!=args.assembly:
            raise ValueError('Duplicate allele/sample or assembly mismatch')
        seen.add(record)
        if row.get('orientation_kind') not in ('read_strand','pair_orientation'):
            raise ValueError('Declare read_strand or pair_orientation support')
        reasons,vaf=quality_reasons(row,args)
        driver=driver_map.get(key)
        if driver is None:
            reasons.append('not_in_curated_driver_set')
        audit.append(dict(sample=row['sample'],chrom=key[0],pos=key[1],ref=key[2],alt=key[3],
                          gene=driver['gene'] if driver else '',vaf=vaf,orientation_kind=row['orientation_kind'],
                          evidence=driver['evidence'] if driver else '',candidate=not reasons,reasons=';'.join(reasons),
                          review='somatic_origin_and_clinical_context_required' if not reasons else ''))
    by_sample=defaultdict(list)
    for row in audit:
        if row['candidate']: by_sample[row['sample']].append(row)
    calls=[]
    for sample in roster:
        variants=by_sample[sample['sample']]
        calls.append(dict(sample=sample['sample'],donor=sample['donor'],specimen=sample['specimen'],
                          ch_status='candidate_driver_detected' if variants else 'no_qualifying_variant_in_supplied_records',
                          n_candidates=len(variants),max_vaf=max((v['vaf'] for v in variants),default=None),
                          genes=';'.join(sorted({v['gene'] for v in variants})),
                          assessment='screening_only_no_negative_call'))
    return audit,calls


def import_calls(args,roster):
    required=['Biobank_Subject_ID','hasCHIP','nCHIP','Gene','VAF','Batch_CHIP_call']
    rows=read_table(args.calls,required)
    source={}
    for row in rows:
        key=row['Biobank_Subject_ID']
        if key in source:
            raise ValueError('Repeated source subject')
        source[key]=row
    variants=read_table(args.variants,['Biobank_Subject_ID','varID','Gene','VAF','DP','AD_Alt','ADF_Alt','ADR_Alt'])
    selected={r['source_subject']:r for r in roster}
    audit=[]; seen=set()
    for row in variants:
        if row['Biobank_Subject_ID'] not in selected:
            continue
        target=selected[row['Biobank_Subject_ID']]
        key=(target['sample'],row['varID'])
        if key in seen:
            raise ValueError('Repeated source variant')
        seen.add(key)
        vaf=number(row['VAF'])
        if vaf is None or not 0<=vaf<=1:
            raise ValueError('Invalid provider VAF')
        audit.append(dict(sample=target['sample'],variant=row['varID'],gene=row['Gene'],vaf=vaf,dp=count(row['DP']),
                          alt_count=count(row['AD_Alt']),alt_forward=count(row['ADF_Alt']),alt_reverse=count(row['ADR_Alt'])))
        record=audit[-1]
        if record['dp'] is not None and record['alt_count'] is not None and record['alt_count']>record['dp']:
            raise ValueError('Provider alternate reads exceed depth')
        if record['alt_count'] is not None and any(record[k] is not None and record[k]>record['alt_count'] for k in ('alt_forward','alt_reverse')):
            raise ValueError('Provider strand support exceeds alternate reads')
    by_sample=defaultdict(list)
    for row in audit: by_sample[row['sample']].append(row)
    calls=[]
    for target in roster:
        source_row=source.get(target['source_subject'])
        records=by_sample[target['sample']]
        if source_row is None:
            if records: raise ValueError('Provider variants exist without a subject classification')
            status='unknown'; n=None; batch=''
        else:
            n=count(source_row['nCHIP']); flag=source_row['hasCHIP']; batch=source_row['Batch_CHIP_call']
            if flag not in ('0','1') or n is None or (flag=='1')!=(n>0) or n!=len(records):
                raise ValueError('Provider call, variant count or mutation count disagree')
            status='provider_reported_positive' if flag=='1' else 'provider_reported_negative'
        calls.append(dict(sample=target['sample'],donor=target['donor'],specimen=target['specimen'],
                          ch_status=status,n_candidates=n,max_vaf=max((r['vaf'] for r in records),default=None),
                          genes=';'.join(sorted({r['gene'] for r in records})),batch=batch,
                          assessment='Imported provider classification; no independent somatic re-calling'))
    return audit,calls


def run(args):
    required=['sample','donor','specimen'] + (['source_subject'] if args.action=='import' else ['tissue'])
    roster=read_table(args.samples,required)
    if not roster or any(not row[k] for row in roster for k in required):
        raise ValueError('Nonempty specimen roster required')
    for field in ('sample','donor'):
        if len({r[field] for r in roster})!=len(roster):
            raise ValueError('Choose one independent specimen per donor')
    if args.action=='import' and len({r['source_subject'] for r in roster})!=len(roster):
        raise ValueError('Source subject mapping is not unique')
    if args.action=='screen' and any(r['tissue']!='blood' for r in roster):
        raise ValueError('CH screening requires an explicit blood specimen')
    if Path(args.out).exists():
        raise ValueError('Choose a new output directory')
    audit,calls=import_calls(args,roster) if args.action=='import' else screen(args,roster)
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.ch-',dir=out.parent))
    def write(rows,filename,empty_header):
        with (stage/filename).open('w',encoding='utf-8',newline='') as handle:
            writer=csv.DictWriter(handle,list(rows[0]) if rows else empty_header,delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
    try:
        fields=(['sample','variant','gene','vaf','dp','alt_count','alt_forward','alt_reverse'] if args.action=='import' else
                ['sample','chrom','pos','ref','alt','gene','vaf','orientation_kind','evidence','candidate','reasons','review'])
        write(audit,'variants.tsv',fields)
        write(calls,'samples.tsv',['sample','donor','specimen','ch_status'])
        # Preserve provider negative/unknown labels. Candidate-only screens do not create reference controls.
        events=[]
        for row in calls:
            state={'provider_reported_positive':'altered','provider_reported_negative':'reference'}.get(row['ch_status'],'unknown')
            events.append(dict(event='CH_provider_status' if args.action=='import' else 'CH_candidate_screen',origin='somatic',
                               donor=row['donor'],specimen=row['specimen'],state=state,source_status=row['ch_status']))
        write(events,'genomic_events.tsv',[])
        paths=[Path(args.samples),Path(args.variants),Path(args.calls if args.action=='import' else args.drivers)]
        write([dict(input=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],'inputs.tsv',[])
        summary=dict(action=args.action,samples=len(calls),variant_records=len(audit),
                     settings={k:v for k,v in vars(args).items() if k not in ('samples','variants','drivers','calls','out')},
                     limitations=['No clinical diagnosis','No inferred negative from missing VCF records',
                                  'Blood tissue and somatic origin need verification','U2AF1 reference duplication needs a separate validated approach'])
        (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        if out.exists():raise ValueError('Output appeared during execution')
        stage.rename(out)
    finally:
        if stage.exists():shutil.rmtree(stage)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    for name in ('screen','import'):
        s=sub.add_parser(name)
        for opt in ('samples','variants','out'):s.add_argument('--'+opt,required=True)
        if name=='import':s.add_argument('--calls',required=True)
        else:
            s.add_argument('--drivers',required=True);s.add_argument('--assembly',required=True)
            s.add_argument('--population-af-field',default='gnomAD_AF')
            s.add_argument('--min-depth',type=int,default=20);s.add_argument('--min-alt-reads',type=int,default=3)
            s.add_argument('--min-vaf',type=float,default=.02);s.add_argument('--max-population-af',type=float,default=.001)
    args=p.parse_args()
    if args.action=='screen' and (args.min_depth<1 or args.min_alt_reads<1 or not 0<args.min_vaf<1 or not 0<=args.max_population_af<=1):
        p.error('Invalid screening thresholds')
    try:print(json.dumps(run(args)))
    except (OSError,ValueError,KeyError) as error:p.exit(2,str(error)+'\n')


if __name__=='__main__':main()
