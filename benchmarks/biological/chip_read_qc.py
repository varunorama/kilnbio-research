"""Reconcile independently reported read counts, lengths and QC flags before mapping."""
import argparse
import csv
import json
import zipfile
from pathlib import Path
from chip_spikein import R1,SAMPLES

def fastqc(path):
    with zipfile.ZipFile(path) as z:
        entry,=[n for n in z.namelist() if n.endswith('/fastqc_data.txt')]
        summary,=[n for n in z.namelist() if n.endswith('/summary.txt')]
        lines=z.read(entry).decode().splitlines();end=lines.index('>>END_MODULE')
        stats=dict(l.split('\t',1) for l in lines[3:end])
        flags={row.split('\t')[1]:row.split('\t')[0] for row in z.read(summary).decode().splitlines()}
    return stats,flags

def reconcile(r1,r2,out):
    rows=[];allflags=[]
    for campaign,root,records in [('R1',r1,R1),('R2',r2,SAMPLES)]:
        for sample,(mark,condition) in records.items():
            folder=root/sample
            if not (folder/'fastqc.complete').is_file():raise ValueError('Incomplete FastQC: '+sample)
            trim=json.loads((folder/(sample+'.fastq.gz_trimming_report.json')).read_text())
            raw,before=fastqc(folder/(sample+'_fastqc.zip'));clean,after=fastqc(folder/(sample+'_trimmed_fastqc.zip'))
            n=int(raw['Total Sequences']);nkeep=int(clean['Total Sequences']);proc=trim['read_processing']
            if (trim['mode']!='single-end' or n!=proc['total_reads'] or n!=5000000 or nkeep!=proc['reads_written']
                    or not 0<nkeep<=n):raise ValueError('FastQC and trimmer count mismatch: '+sample)
            rows.append(dict(sample=sample,campaign=campaign,mark=mark,condition=condition,raw_reads=n,trimmed_reads=nkeep,
                retained_percent=100*nkeep/n,raw_length=raw['Sequence length'],trimmed_length=clean['Sequence length'],
                raw_GC_percent=int(raw['%GC']),trimmed_GC_percent=int(clean['%GC']),
                trimmed_per_base_quality=after['Per base sequence quality'],trimmed_adapter_content=after['Adapter Content'],
                flags_after_trim='; '.join(k+':'+v for k,v in after.items() if v!='PASS')))
            for stage,flags in [('raw',before),('trimmed',after)]:
                allflags.extend(dict(sample=sample,stage=stage,module=k,status=v) for k,v in flags.items())
    out.mkdir(parents=True,exist_ok=False)
    for name,data in [('read_qc.tsv',rows),('fastqc_flags.tsv',allflags)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,data[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['campaign1','campaign2','out']:p.add_argument(k,type=Path)
    a=p.parse_args();reconcile(a.campaign1,a.campaign2,a.out)
