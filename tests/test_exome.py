"""Dependency-free preflight regressions for the Sarek adapter."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from downstream.exome_prepare import validate_params,validate_samples

class ExomePreflight(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        (self.base/'ref.fai').write_text('chr1\t100\t0\t50\t51\n')
        (self.base/'target.bed').write_text('chr1\t0\t50\n')
        self.p=dict(wes=True,aligner='bwa-mem2',intervals=str(self.base/'target.bed'),fasta_fai=str(self.base/'ref.fai'),tools='haplotypecaller',skip_tools='baserecalibrator')
    def tearDown(self):self.temp.cleanup()
    def test_known_bed_and_no_proprietary_aligner(self):
        validate_params('exome',self.p)
        self.p['aligner']='sentieon-bwamem'
        with self.assertRaises(ValueError):validate_params('exome',self.p)
    def test_unknown_ploidy_is_not_silently_diploid(self):
        (self.base/'target.bed').write_text('chrX\t0\t50\n')
        with self.assertRaisesRegex(ValueError,'autosomal'):validate_params('exome',self.p)
    def test_bqsr_resources_must_be_explicit(self):
        self.p.pop('skip_tools')
        with self.assertRaisesRegex(ValueError,'BQSR'):validate_params('exome',self.p)
    def test_somatic_resources_must_be_explicit(self):
        self.p['tools']='mutect2'
        with self.assertRaisesRegex(ValueError,'panel of normals'):validate_params('ch-somatic',self.p)
    def test_sample_lane_and_file_identity(self):
        h=['patient','sex','status','sample','lane','fastq_1','fastq_2']
        row=dict(zip(h,['p1','XX','0','s1','L1','a.fastq.gz','b.fastq.gz']))
        p=self.base/'samples.csv'
        def write(rows):
            with p.open('w',newline='') as f:
                w=csv.DictWriter(f,h);w.writeheader();w.writerows(rows)
        write([row]);validate_samples('exome',p,self.p,lambda v,b:v)
        write([row,row])
        with self.assertRaisesRegex(ValueError,'Duplicate lane'):validate_samples('exome',p,self.p,lambda v,b:v)
        write([dict(row,fastq_2='a.fastq.gz')])
        with self.assertRaises(ValueError):validate_samples('exome',p,self.p,lambda v,b:v)

if __name__=='__main__':unittest.main(verbosity=2)
