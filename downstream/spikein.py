"""Derive relative exposure factors from matched experimental spike-in counts.

Counts must use the same read/fragment unit and mapping/filtering rules. These
factors do not identify technical variation without an experimental reference.
No differential tests or biological uncertainty estimates are produced here.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
from analysis import digest, need, rows

FIELDS = ['sample','normalization_group','target_reads','reference_reads',
          'target_amount','reference_amount','amount_unit','reference_material']


def factors(records, group, min_reference):
    selected=[r for r in records if r['normalization_group']==group]
    need(selected,'Normalization group absent')
    need(len({r['sample'] for r in selected})==len(selected),'Duplicate sample')
    need(all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',r['sample']) for r in selected),'Unsafe sample ID')
    need(len({r['amount_unit'] for r in selected})==1 and selected[0]['amount_unit'] in ('ng_chromatin','cells'),'Use one declared input unit: ng_chromatin or cells')
    need(len({r['reference_material'] for r in selected})==1 and selected[0]['reference_material'].strip(),'Use the same reference stock/assay within a normalization group')
    result=[]
    for row in selected:
        target,reference=[float(row[k]) for k in ('target_reads','reference_reads')]
        a,b=[float(row[k]) for k in ('target_amount','reference_amount')]
        need(all(math.isfinite(v) for v in (target,reference,a,b)),'Nonfinite spike-in measurement')
        need(target>=0 and target.is_integer() and reference>=min_reference and reference.is_integer(),'Nonnegative integer target counts and sufficient integer reference counts required')
        need(a>0 and b>0,'Positive pre-IP input amounts required')
        # C / exposure = (C / spike counts) * (spike amount / target amount).
        exposure=reference*a/b
        need(math.isfinite(exposure) and exposure>0,'Invalid reference-derived exposure')
        result.append(dict(sample=row['sample'],target_reads=int(target),reference_reads=int(reference),
                           exposure=exposure,target_per_reference_adjusted_input=target/exposure))
    center=math.exp(sum(math.log(r['exposure']) for r in result)/len(result))
    for r in result:r['size_factor']=r['exposure']/center
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('counts',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--group',required=True)
    p.add_argument('--min-reference',type=int,default=1000)
    args=p.parse_args()
    try:
        need(args.min_reference>=1,'Positive minimum reference count required')
        header,records=rows(args.counts)
        need(header==FIELDS,f'Expected columns: {FIELDS}')
        result=factors(records,args.group,args.min_reference)
        need(not args.output.exists(),'Choose a new output directory')
        args.output.mkdir(parents=True)
        with (args.output/'size_factors.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,['sample','size_factor']);w.writeheader()
            w.writerows({k:r[k] for k in w.fieldnames} for r in result)
        with (args.output/'diagnostics.tsv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,result[0].keys(),delimiter='\t');w.writeheader();w.writerows(result)
        record={'source':str(args.counts.resolve()),'source_sha256':digest(args.counts),'group':args.group,
                'min_reference':args.min_reference,'formula':'exposure = reference_reads * target_amount / reference_amount; size_factor = exposure / geometric_mean(exposure)',
                'scope':'Relative signal per declared input unit; descriptive, no biological variance or significance',
                'assumptions':['Comparable reference recovery and stock across the contrast',
                    'Competitive mapping and identical count units/filters for both species',
                    'Pre-IP target/reference amounts known; amount units do not imply per-cell measurements',
                    'Cross-reactive antibodies and unequal recovery can bias scaling; inspect controls'],
                'citations':['https://doi.org/10.1371/journal.pone.0166438','https://doi.org/10.1038/s41587-024-02377-y']}
        (args.output/'provenance.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        print(args.output/'size_factors.csv')
    except (ValueError,KeyError,OSError) as e:p.exit(2,f'ERROR: {e}\n')


if __name__=='__main__':main()
