"""Structural checks of published tiny-test outputs, not assay accuracy."""
import argparse,gzip,json,math,struct
from pathlib import Path
from downstream.analysis import check_rsem
p=argparse.ArgumentParser();p.add_argument('assay',choices=['atacseq','chipseq','rnaseq']);p.add_argument('results');p.add_argument('output')
p.add_argument('--rna-quantifier', choices=['salmon','rsem'], default='salmon')
a=p.parse_args()
files=[f for f in Path(a.results).rglob('*') if f.is_file()]
counts={'bam':0,'bigwig':0,'peak_files':0,'peak_rows':0,'salmon_quant_files':0,'salmon_transcripts':0,'star_count_files':0,'featurecount_tables':0};errors=[]
for f in files:
    try:
        if f.suffix=='.bam':
            with gzip.open(f,'rb') as h: assert h.read(4)==b'BAM\1'
            counts['bam']+=1
        if f.suffix.lower() in ('.bw','.bigwig'):
            magic=f.open('rb').read(4)
            assert len(magic)==4 and 0x888FFC26 in (struct.unpack('<I',magic)[0],struct.unpack('>I',magic)[0])
            counts['bigwig']+=1
        if f.name.endswith(('.narrowPeak','.broadPeak','.narrowPeak.gz','.broadPeak.gz')):
            opener=gzip.open if f.suffix=='.gz' else open
            with opener(f,'rt') as h:
                for line in h:
                    if not line.strip() or line.startswith(('track','browser','#')):continue
                    v=line.rstrip().split('\t');assert len(v)>=9
                    assert int(v[1])>=0 and int(v[2])>int(v[1])
                    counts['peak_rows']+=1
            counts['peak_files']+=1
        if f.name=='quant.sf' or f.name.endswith('.quant.sf'):
            with f.open() as h:
                header=h.readline().rstrip().split('\t')
                assert header==['Name','Length','EffectiveLength','TPM','NumReads']
                for line in h:
                    v=line.rstrip().split('\t');assert len(v)==5
                    assert all(math.isfinite(float(x)) and float(x)>=0 for x in v[1:])
                    counts['salmon_transcripts']+=1
            counts['salmon_quant_files']+=1
        if f.name.endswith('ReadsPerGene.out.tab'):
            with f.open() as h:
                for line in h:
                    v=line.rstrip().split('\t');assert len(v)==4 and all(int(x)>=0 for x in v[1:])
            counts['star_count_files']+=1
        if f.name.endswith('.featureCounts.txt'):
            with f.open() as h:
                lines=(line for line in h if not line.startswith('#'))
                header=next(lines).rstrip().split('\t');assert header[0]=='Geneid'
                for line in lines:
                    v=line.rstrip().split('\t');assert len(v)==len(header) and all(math.isfinite(float(x)) and float(x)>=0 for x in v[6:])
            counts['featurecount_tables']+=1
    except (AssertionError,ValueError,OSError,EOFError,StopIteration) as e:errors.append(dict(file=str(f),error=type(e).__name__))
multiqc=any(f.match('*multiqc*.html') and f.stat().st_size>1000 for f in files)
required=counts['bam']>0 and multiqc
if a.assay in ('atacseq','chipseq'):required=required and counts['bigwig']>0 and counts['peak_files']>0 and counts['peak_rows']>0
if a.assay=='rnaseq' and a.rna_quantifier=='salmon':
    required=required and counts['salmon_quant_files']>0 and counts['salmon_transcripts']>0
if a.assay=='rnaseq' and a.rna_quantifier=='rsem':
    gene_files=[str(f) for f in files if f.name.endswith(('.genes.results','.genes.results.gz'))]
    counts['rsem_gene_files']=len(gene_files)
    required=required and bool(gene_files)
    try:
        check_rsem(gene_files)
    except (ValueError, OSError, EOFError, TypeError) as exc:
        errors.append(dict(error=str(exc)))
result=dict(status='passed' if required and not errors else 'failed',assay=a.assay,counts=counts,multiqc=multiqc,errors=errors,scope='Basic file/count/coordinate integrity, not mapping accuracy, biological QC thresholds or a method benchmark')
Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
raise SystemExit(0 if result['status']=='passed' else 1)
