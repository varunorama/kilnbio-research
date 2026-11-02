"""Stream CpG counts into chr,pos,N,X without estimating counts from percentages."""
import argparse, csv, gzip, json, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from analysis import need, file_at, digest

class Reference:
    """Random access via samtools .fai; uncompressed FASTA only."""
    def __init__(self,path,index):
        need(path.suffix!='.gz','Use uncompressed indexed FASTA')
        self.handle=path.open('rb');self.index={}
        for line in index.read_text(encoding='utf-8').splitlines():
            f=line.split('\t');need(len(f)>=5,'Invalid FASTA index')
            need(f[0] not in self.index,'Duplicate reference contig')
            self.index[f[0]]=tuple(map(int,f[1:5]))
        need(self.index,'Empty FASTA index')
    def base(self,chrom,pos):
        need(chrom in self.index,f'Contig absent from reference: {chrom}')
        length,offset,bases,width=self.index[chrom]
        need(1<=pos<=length and bases>0 and width>=bases,'Position/index outside reference bounds')
        q=pos-1;self.handle.seek(offset+(q//bases)*width+q%bases)
        return self.handle.read(1).decode('ascii').upper()
    def cpg(self,chrom,pos,strand=None):
        b=self.base(chrom,pos)
        if b=='C':
            need(strand in (None,'+','.'),'Strand disagrees with reference cytosine')
            need(self.base(chrom,pos+1)=='G','Non-CpG context is unsupported by this CpG interface')
            return pos,'+'
        need(b=='G' and strand in (None,'-'),'Site must be a reference C or reverse-strand C at G')
        need(self.base(chrom,pos-1)=='C','Non-CpG reverse context')
        return pos-1,'-'

def integer(value):
    need(str(value).isdigit(),'Counts/positions must be nonnegative decimal integers')
    return int(value)

def convert(config,output):
    cfg=json.loads(config.read_text(encoding='utf-8-sig'));out=output.resolve()
    allowed={'format','input','reference','reference_index','assembly','modification','chemistry','pipeline_version','measurement_reviewed'}
    need(isinstance(cfg,dict) and not set(cfg)-allowed,'Unknown import configuration keys')
    fmt=cfg.get('format');need(fmt in ('bismark_coverage','bismark_cytosine','biomodal_cxreport','modkit_bedmethyl'),'Unsupported counts format')
    mod=cfg.get('modification');need(mod in ('modC','5mC','5hmC'),'Declare modification channel')
    need(cfg.get('assembly') and cfg.get('measurement_reviewed') is True,'Declare assembly and review source measurement semantics')
    if fmt.startswith('bismark'):
        need(mod=='modC' and cfg.get('chemistry')=='bisulfite','Conventional bisulfite measures combined modC; specialized chemistry needs a dedicated adapter')
    if fmt=='biomodal_cxreport':
        need(cfg.get('chemistry') in ('duet_modC','duet_evoC'),'Declare duet_modC or duet_evoC')
        need(cfg.get('pipeline_version') and (cfg['chemistry']=='duet_evoC' or mod=='modC'),'duet +modC does not distinguish 5mC/5hmC')
    if fmt=='modkit_bedmethyl':
        need(mod in ('5mC','5hmC'),'Import m/h channels separately; combined modC is not inferred from a single channel')
        need(cfg.get('chemistry')=='nanopore' and cfg.get('pipeline_version'),'Record modkit version and chemistry')
    need(not out.exists(),'Use a new output directory')
    source=file_at(cfg['input'],config.parent);fasta=file_at(cfg['reference'],config.parent)
    fai=file_at(cfg.get('reference_index',str(fasta)+'.fai'),config.parent)
    ref=Reference(fasta,fai)
    out.mkdir(parents=True);status=out/'status.json'
    processed=written=zero=ignored=0;closed=set();last_chrom=None;last_raw=-1;pending=None
    handle=None
    try:
        opener=gzip.open if source.suffix=='.gz' else open
        with opener(source,'rt',encoding='utf-8') as src,(out/'counts.tsv').open('w',encoding='utf-8',newline='') as dst:
            writer=csv.writer(dst,delimiter='\t',lineterminator='\n');writer.writerow(['chr','pos','N','X'])
            def flush(record):
                nonlocal written,zero
                if record is None:return
                chrom,pos,N,X,_=record
                if N==0:zero+=1;return
                writer.writerow([chrom,pos,N,X]);written+=1
            for number,line in enumerate(src,1):
                if not line.strip() or line.startswith('#'):continue
                f=line.rstrip('\r\n').split('\t');processed+=1
                if fmt in ('biomodal_cxreport','bismark_cytosine'):
                    need(len(f)==7,f'Expected seven cytosine-report columns at line {number}')
                    chrom=f[0];pos=integer(f[1]);strand=f[2]
                    need(strand in ('+','-') and f[5]=='CG','Expected stranded CpG cytosine report')
                    X=integer(f[3]);N=X+integer(f[4])
                elif fmt=='bismark_coverage':
                    need(len(f)==6,f'Expected six Bismark coverage columns at line {number}')
                    chrom=f[0];pos=integer(f[1]);need(integer(f[2])==pos,'Coverage sites must be single-base, one-based')
                    X=integer(f[4]);N=X+integer(f[5]);strand=None
                else:
                    need(len(f)==18,f'Expected modkit 18-column bedMethyl at line {number}')
                    code=f[3].split(',')[0]
                    if code!=('m' if mod=='5mC' else 'h'):ignored+=1;continue
                    chrom=f[0];pos=integer(f[1])+1;strand=f[5]
                    need(integer(f[2])==pos and strand=='.','Require modkit --cpg --combine-strands, one-base records')
                    N=integer(f[9]);X=integer(f[11])
                    need(N==X+integer(f[12])+integer(f[13]),'N_valid must include modified, canonical and OTHER modification calls')
                    if ',' in f[3]:need(f[3].split(',')[1]=='CG','Only CpG motif is supported')
                need(N>=X>=0 and pos>=1,'Invalid coverage or modified count')
                canonical,observed=ref.cpg(chrom,pos,strand)
                if chrom!=last_chrom:
                    need(chrom not in closed,'Contig blocks must not recur; sort input by contig/position')
                    if last_chrom is not None:closed.add(last_chrom)
                    last_chrom=chrom;last_raw=-1
                need(pos>=last_raw,'Input must be coordinate sorted within each contig')
                last_raw=pos
                key=(chrom,canonical)
                if pending is None or key!=tuple(pending[:2]):
                    flush(pending);pending=[chrom,canonical,0,0,set()]
                need(observed not in pending[4],'Duplicate strand/site record; resolve lanes upstream')
                pending[2]+=N;pending[3]+=X;pending[4].add(observed)
            flush(pending)
        need(written>0,'No covered CpG records for the requested modification')
        record={'complete':True,'config':cfg,'records_read':processed,'sites_written':written,
                'zero_coverage_omitted':zero,'other_modification_rows_ignored':ignored,
                'coordinates':'1-based-cytosine','context':'CpG','strand_policy':'sum forward and reverse observations at forward C',
                'missingness':'Omitted sites are missing coverage, never zero methylation',
                'source_sha256':digest(Path(__file__)),
                'input_sha256':{str(p):digest(p) for p in (config,source,fasta,fai)}}
        status.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        return record
    except Exception as e:
        status.write_text(json.dumps({'complete':False,'error':str(e)})+'\n',encoding='utf-8');raise
    finally:ref.handle.close()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('config');p.add_argument('output');a=p.parse_args()
    try:print(json.dumps(convert(Path(a.config).resolve(),Path(a.output)),indent=2))
    except (ValueError,KeyError,OSError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()
