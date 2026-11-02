"""Score peak callers without using ENCODE ATAC or ChIP peak calls as truth."""
import argparse
from array import array
import csv
import gzip
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
NAMES = {'MACS2':'calls_peaks.narrowPeak', 'MACS3':'calls_peaks.narrowPeak',
         'HMMRATAC':'calls_accessible_regions.narrowPeak','Omnipeak':'calls.broadPeak'}


def intervals(path):
    opener = gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as stream:
        for line in stream:
            if line.startswith(('#','track','browser')) or not line.strip():
                continue
            fields = line.split()
            start,end = int(fields[1]),int(fields[2])
            if not 0<=start<end:
                raise ValueError(f'Invalid BED interval in {path}')
            yield fields[0],start,end


def grid(sizes, width):
    return {c:np.zeros(size//width,dtype=bool) for c,size in sizes.items()}


def interval_grid(path, sizes, width, mask=False):
    result=grid(sizes,width)
    lengths={c:[] for c in sizes}
    for c,s,e in intervals(path):
        if c not in result:
            continue
        if mask:
            lo,hi=s//width,(e+width-1)//width
        else:
            lo=(s-width//2+width-1)//width
            hi=(e-width//2+width-1)//width
        result[c][max(0,lo):max(0,hi)]=True
        lengths[c].append(e-s)
    return result,lengths


def tag_grid(path, sizes, width, atac):
    counts={c:array('I',[0])*(size//width) for c,size in sizes.items()}
    with path.open() as stream:
        for line in stream:
            f=line.split()
            s,e=int(f[1]),int(f[2])
            pos=(s+e)//2 if atac else (s+100 if f[5]=='+' else e-1-100)
            b=pos//width
            if 0<=b<len(counts[f[0]]):
                counts[f[0]][b]+=1
    return {c:np.frombuffer(data,dtype=np.uint32) for c,data in counts.items()}


def ratio(x,y):
    return float(x/y) if y else float('nan')


def intersect(a,b):
    return int(np.count_nonzero(a & b))


def score_experiment(base,row):
    atac=row['assay']=='ATAC'
    w=100 if atac else 1000
    sizes=dict((c,int(n)) for c,n in (s.split() for s in (base/'ref/panel.sizes').read_text().splitlines()))
    masks,_=interval_grid(base/'ref/mask.bed',sizes,w,mask=True)
    available={c:~m for c,m in masks.items()}
    exp=base/row['experiment']
    callers=['MACS2','MACS3','HMMRATAC'] if atac else ['MACS3','Omnipeak']
    chip_ctrl=[]
    if atac:
        refs=[interval_grid(base/'input'/(row[f'dnase_rep{r}']+'.bed.gz'),sizes,w)[0] for r in [1,2]]
        truth={c:refs[0][c] & refs[1][c] & available[c] for c in sizes}
    else:
        chip_ctrl=[tag_grid(exp/f'control.{part}.bed',sizes,w,False) for part in ['train','test']]
    signal={r:tag_grid(exp/f'rep{r}.bed',sizes,w,atac) for r in [1,2]}
    calls={}; lengths={}
    for r in [1,2]:
        for name in callers:
            fn='calls_peaks.broadPeak' if not atac and name=='MACS3' else NAMES[name]
            raw,lens=interval_grid(exp/f'r{r}-{name}'/fn,sizes,w)
            calls[r,name,'native']={c:raw[c] & available[c] for c in sizes}
            calls[r,name,'equal_bp']={}
            lengths[r,name]=lens
        for ci,c in enumerate(sizes):
            # Each bin has one common training score and one common tie key for all callers.
            smooth=np.convolve(signal[r][c]*available[c],np.ones(5),'same')
            if not atac:
                control=np.convolve(chip_ctrl[0][c]*available[c],np.ones(5),'same')
                # Within-chromosome ranking is invariant to a constant library-size factor.
                smooth=(smooth+1)/(control+1)
            tie=np.random.RandomState(20261101+10*r+ci).random_sample(len(smooth))
            k=min(np.count_nonzero(calls[r,n,'native'][c]) for n in callers)
            for n in callers:
                bins=np.flatnonzero(calls[r,n,'native'][c])
                selected=bins[np.lexsort((tie[bins],-smooth[bins]))[:k]]
                keep=np.zeros(len(smooth),dtype=bool)
                keep[selected]=True
                calls[r,n,'equal_bp'][c]=keep
    records=[]
    for r in [1,2]:
        other=3-r
        for n in callers:
            for mode in ['native','equal_bp']:
                for region in ['panel',*sizes]:
                    chroms=list(sizes) if region=='panel' else [region]
                    a,b=calls[r,n,mode],calls[other,n,mode]
                    predicted=sum(int(a[c].sum()) for c in chroms)
                    common=sum(intersect(a[c],b[c]) for c in chroms)
                    union=predicted+sum(int(b[c].sum()) for c in chroms)-common
                    allbins=sum(int(available[c].sum()) for c in chroms)
                    widths=[v for c in chroms for v in lengths[r,n][c]]
                    rec=dict(experiment=row['experiment'],cell=row['cell'],assay=row['assay'],rep=r,caller=n,view=mode,region=region,
                             bin_bp=w,selected_bins=predicted,eligible_bins=allbins,coverage=ratio(predicted,allbins),
                             replicate_shared_bins=common,replicate_union_bins=union,
                             replicate_jaccard=ratio(common,union),native_n_peaks=len(widths),
                             native_median_width_bp=float(np.median(widths)) if widths else 0,
                             native_p90_width_bp=float(np.percentile(widths,90)) if widths else 0,
                             dnase_reference_bins='',dnase_overlap_bins='',dnase_precision='',dnase_recall='',dnase_f1='',
                             heldout_treatment_hits='',heldout_treatment_total='',heldout_input_hits='',heldout_input_total='',
                             heldout_enrichment='')
                    if atac:
                        positive=sum(int(truth[c].sum()) for c in chroms)
                        tp=sum(intersect(a[c],truth[c]) for c in chroms)
                        rec.update(dnase_reference_bins=positive,dnase_overlap_bins=tp,dnase_precision=ratio(tp,predicted),
                                   dnase_recall=ratio(tp,positive),dnase_f1=ratio(2*tp,predicted+positive))
                        if mode=='equal_bp' and not predicted:
                            rec.update(dnase_precision='',dnase_recall='',dnase_f1='')
                    else:
                        fractions=[]
                        for label,track in [('treatment',signal[other]),('input',chip_ctrl[1])]:
                            hit=sum(int(track[c][a[c]].sum()) for c in chroms)
                            total=sum(int(track[c][available[c]].sum()) for c in chroms)
                            rec[f'heldout_{label}_hits']=hit
                            rec[f'heldout_{label}_total']=total
                            fractions.append(ratio(hit,total))
                        rec['heldout_enrichment']=ratio(*fractions)
                    records.append(rec)
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('experiment')
    p.add_argument('output',type=Path)
    args=p.parse_args()
    with (ROOT/'experiments.tsv').open() as f:
        row=next(r for r in csv.DictReader(f,delimiter='\t') if r['experiment']==args.experiment)
    data=score_experiment(args.directory,row)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,data[0].keys(),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)
