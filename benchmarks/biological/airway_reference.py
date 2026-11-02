"""Build and audit the full human + ERCC reference used by the airway benchmark.

Never concatenate GTFs with unconditional blank separators: STAR may tolerate
them while RSEM rejects them. Keep only coordinates on the supplied FASTA and
record every filtered contig/row. This recipe expects the official ERCC92 ZIP.
"""
import argparse, collections, gzip, hashlib, json, re, zipfile
from pathlib import Path

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def make(fasta,gtf,ercc,out):
    out.mkdir(parents=True,exist_ok=False)
    lengths={};stats=collections.Counter();missing=collections.Counter();spike_exons=set()
    with zipfile.ZipFile(ercc) as z:
        with (out/'genome.fa').open('wb') as dst:
            for label,src in [('human',gzip.open(fasta,'rb')),('ERCC',z.open('ERCC92.fa'))]:
                name=None
                with src:
                    for line in src:
                        line=line.strip()
                        if not line:continue
                        if line.startswith(b'>'):
                            name=line[1:].split()[0].decode('ascii')
                            if not name or name in lengths:raise ValueError('Empty or duplicate FASTA contig')
                            lengths[name]=0;stats[label+'_contigs']+=1
                        elif name is None:raise ValueError('FASTA data before header')
                        else:lengths[name]+=len(line)
                        dst.write(line+b'\n')
        if stats['ERCC_contigs']!=92 or min(lengths.values())<=0:raise ValueError('FASTA/ERCC completeness failure')
        with (out/'genome.gtf').open('wb') as dst:
            for label,src in [('human',gzip.open(gtf,'rb')),('ERCC',z.open('ERCC92.gtf'))]:
                with src:
                    for number,line in enumerate(src,1):
                        line=line.strip()
                        if not line or line.startswith(b'#'):
                            stats['blank_or_comment_rows_removed']+=1;continue
                        fields=line.split(b'\t')
                        if len(fields)!=9:raise ValueError(f'{label} malformed GTF row {number}')
                        chrom=fields[0].decode('ascii')
                        if chrom not in lengths:
                            if label=='ERCC':raise ValueError('ERCC annotation/FASTA mismatch')
                            missing[chrom]+=1;continue
                        if not 1<=int(fields[3])<=int(fields[4])<=lengths[chrom]:raise ValueError(f'Out-of-bounds GTF row {number}')
                        if fields[2]==b'exon':
                            if fields[6] not in (b'+',b'-') or not all(re.search(k+b' "[^"]+"',fields[8]) for k in (b'gene_id',b'transcript_id')):raise ValueError(f'Invalid exon {number}')
                            if label=='ERCC':spike_exons.add(chrom)
                        dst.write(line+b'\n');stats[label+'_gtf_rows_written']+=1
    if len(spike_exons)!=92:raise ValueError('Expected exon annotation on each ERCC contig')
    record={'input_sha256':{p.name:sha256(p) for p in [fasta,gtf,ercc]},'output_sha256':{n:sha256(out/n) for n in ['genome.fa','genome.gtf']},
            'stats':dict(stats),'unrepresented_contig_rows_removed':dict(missing),'contig_lengths':lengths}
    (out/'reference-audit.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'stats':stats,'unrepresented_rows':sum(missing.values())}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('human_fasta_gz','human_gtf_gz','ercc_zip','output'):p.add_argument(name,type=Path)
    a=p.parse_args();make(a.human_fasta_gz,a.human_gtf_gz,a.ercc_zip,a.output)
