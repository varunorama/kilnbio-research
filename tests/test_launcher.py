import csv, gzip, json, os, subprocess, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare import HEADERS
ROOT=Path(__file__).resolve().parents[1]

class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='omics_')
        self.base=Path(self.temp.name)
        for name in ('ref.fa','ref.gtf','blacklist.bed','whitelist.txt'): (self.base/name).write_text('fixture\n')
        for name in ('r1','r2','c1','c2'):
            with gzip.open(self.base/f'{name}.fq.gz','wt') as f: f.write('@fixture\nACGT\n+\nIIII\n')
    def tearDown(self): self.temp.cleanup()
    def inputs(self,assay,preset=None):
        params=json.loads((ROOT/'examples'/f'{preset or assay}.json').read_text())
        for key,name in [('fasta','ref.fa'),('gtf','ref.gtf'),('blacklist','blacklist.bed'),('barcode_whitelist','whitelist.txt')]:
            if key in params: params[key]=str(self.base/name)
        row=dict(sample='S1',fastq_1='r1.fq.gz',fastq_2='r2.fq.gz')
        if assay in ('wgbs','rrbs'): row['genome']=''
        if assay=='rnaseq': row['strandedness']='auto'
        if assay in ('atacseq','chipseq'): row['replicate']='1'
        if assay=='chipseq': row.update(antibody='H3K27ac',control='INPUT',control_replicate='1')
        rows=[row]
        if assay=='chipseq': rows.append(dict(sample='INPUT',fastq_1='c1.fq.gz',fastq_2='c2.fq.gz',replicate='1',antibody='',control='',control_replicate=''))
        return params,rows
    def run_case(self,assay,params=None,rows=None,extra=(),expect=0,runname='run'):
        default_p,default_r=self.inputs(assay);params=default_p if params is None else params;rows=default_r if rows is None else rows
        sheet=self.base/'samples.csv';param=self.base/'params.json'
        with sheet.open('w',newline='') as f:
            w=csv.DictWriter(f,HEADERS[assay]);w.writeheader();w.writerows(rows)
        param.write_text(json.dumps(params))
        result=subprocess.run(['bash',str(ROOT/'run.sh'),assay,str(sheet),str(param),str(self.base/runname),'--dry-run',*extra],capture_output=True,text=True)
        self.assertEqual(result.returncode,expect,result.stdout+result.stderr)
        return result
    def test_all_assays(self):
        for assay in HEADERS:
            with self.subTest(assay=assay):
                self.run_case(assay,runname=assay)
                resolved=json.loads((self.base/assay/'provenance/params.resolved.json').read_text())
                self.assertTrue(resolved['save_align_intermeds'])
    def test_bad_control(self):
        p,r=self.inputs('chipseq');r[0]['control']='MISSING';self.run_case('chipseq',p,r,expect=2)
    def test_broad_histone_preset(self):
        params,rows=self.inputs('chipseq','chipseq_broad')
        rows[0]['antibody']='H3K36me3'
        self.run_case('chipseq',params,rows)
        resolved=json.loads((self.base/'run/provenance/params.resolved.json').read_text())
        self.assertEqual((resolved['narrow_peak'],resolved['macs_fdr'],resolved['broad_cutoff']),(False,0.05,0.1))
    def test_rrbs_dedup_rejected(self):
        p,r=self.inputs('rrbs');p['skip_deduplication']=False;self.run_case('rrbs',p,r,expect=2)
    def test_mixed_lane_metadata(self):
        p,r=self.inputs('rnaseq');r.append(dict(r[0],fastq_1='c1.fq.gz',fastq_2='',strandedness='auto'));self.run_case('rnaseq',p,r,expect=2)
    def test_duplicate_fastq(self):
        p,r=self.inputs('wgbs');r.append(dict(r[0],sample='OTHER'));self.run_case('wgbs',p,r,expect=2)
    def test_unknown_param(self):
        p,r=self.inputs('wgbs');p['imaginary_parameter']=True;self.run_case('wgbs',p,r,expect=2)
    def test_sc_ambient_correction_not_implicit(self):
        p,r=self.inputs('scrnaseq');p['skip_cellbender']=False;self.run_case('scrnaseq',p,r,expect=2)
    def test_missing_reference(self):
        p,r=self.inputs('wgbs');p['fasta']='missing.fa';self.run_case('wgbs',p,r,expect=2)
    def test_resume_consistency(self):
        self.run_case('wgbs');self.run_case('wgbs',extra=('--resume',))
        self.run_case('wgbs',extra=('--resume','--cpus','8'),expect=2)
    def test_changed_file_resume_rejected(self):
        self.run_case('wgbs');(self.base/'ref.fa').write_text('changed reference\n')
        self.run_case('wgbs',extra=('--resume',),expect=2)
    def test_prebuilt_index_resume_identity(self):
        p,r=self.inputs('rnaseq','rnaseq_rsem')
        index=self.base/'star';index.mkdir();(index/'Genome').write_bytes(b'AAAA')
        p['star_index']='star'
        self.run_case('rnaseq',p,r)
        resolved=json.loads((self.base/'run/provenance/params.resolved.json').read_text())
        self.assertEqual(resolved['star_index'],str(index.resolve()))
        self.run_case('rnaseq',p,r,extra=('--resume',))
        old=(index/'Genome').stat()
        (index/'Genome').write_bytes(b'BBBB')
        os.utime(index/'Genome',ns=(old.st_atime_ns,old.st_mtime_ns))
        self.run_case('rnaseq',p,r,extra=('--resume',),expect=2)
    def test_missing_and_empty_prebuilt_index(self):
        p,r=self.inputs('rnaseq','rnaseq_rsem');p['rsem_index']='absent'
        self.run_case('rnaseq',p,r,expect=2)
        (self.base/'absent').mkdir();self.run_case('rnaseq',p,r,expect=2)
    def test_single_end(self):
        p,r=self.inputs('wgbs');r[0]['fastq_2']='';self.run_case('wgbs',p,r)
    def test_no_shell_execution_from_sample(self):
        p,r=self.inputs('wgbs');r[0]['sample']='S$(touch_BAD)';self.run_case('wgbs',p,r,expect=2)
    def test_mock_nextflow_dispatch(self):
        self.run_case('wgbs')
        bindir=self.base/'bin';bindir.mkdir()
        (bindir/'docker').write_text('#!/usr/bin/env bash\nexit 0\n')
        (bindir/'nextflow').write_text('#!/usr/bin/env bash\nif [[ $1 == -version ]]; then echo "nextflow version 25.10.4"; else printf "%s\\n" "$@" > "$MOCK_ARGS"; fi\n')
        for f in bindir.iterdir(): f.chmod(0o755)
        env=dict(os.environ,PATH=str(bindir)+os.pathsep+os.environ['PATH'],MOCK_ARGS=str(self.base/'args.txt'))
        cmd=['bash',str(ROOT/'run.sh'),'wgbs',str(self.base/'samples.csv'),str(self.base/'params.json'),str(self.base/'run'),'--resume']
        result=subprocess.run(cmd,env=env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        argv=(self.base/'args.txt').read_text().splitlines()
        self.assertEqual(argv[:4],['run','nf-core/methylseq','-r','4.2.0'])
        self.assertIn('-resume',argv)
        first_report=argv[argv.index('-with-report')+1]
        result=subprocess.run(cmd,env=env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        second_argv=(self.base/'args.txt').read_text().splitlines()
        self.assertNotEqual(first_report,second_argv[second_argv.index('-with-report')+1])
        (bindir/'nextflow').write_text('#!/usr/bin/env bash\nif [[ $1 == -version ]]; then echo "nextflow version 25.10.4"; else exit 17; fi\n')
        result=subprocess.run(cmd,env=env,capture_output=True,text=True)
        self.assertEqual(result.returncode,17,result.stderr)

if __name__=='__main__': unittest.main(verbosity=2)
