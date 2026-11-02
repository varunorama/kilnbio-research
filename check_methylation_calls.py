"""Numerical integrity checks for Bismark coverage output; not biological truth."""
import argparse, gzip, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('results');p.add_argument('output');a=p.parse_args()
report=[]
for path in sorted(Path(a.results).rglob('*.cov.gz')):
    n=0;errors=0
    with gzip.open(path,'rt') as f:
        for line in f:
            fields=line.rstrip().split('\t')
            try:
                assert len(fields)==6
                start,end=int(fields[1]),int(fields[2]);percent=float(fields[3]);m,u=int(fields[4]),int(fields[5])
                assert start>=0 and end>=start and m>=0 and u>=0 and m+u>0 and 0<=percent<=100
                assert abs(percent-100*m/(m+u))<=0.011
            except (AssertionError,ValueError): errors+=1
            n+=1
    report.append(dict(file=str(path),rows=n,invalid_rows=errors))
result=dict(files=report,status='passed' if report and all(x['rows']>0 and x['invalid_rows']==0 for x in report) else 'failed',scope='Coverage-table numerical integrity only, not methylation accuracy or conversion efficiency')
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
raise SystemExit(0 if result['status']=='passed' else 1)
