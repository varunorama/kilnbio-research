import importlib.util, json, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mods',ROOT/'modifications/import_counts.py')
mods=importlib.util.module_from_spec(spec);spec.loader.exec_module(mods)

class ModificationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        (self.base/'ref.fa').write_bytes(b'>chr1\nACGTCGAA\n')
        (self.base/'ref.fa.fai').write_text('chr1\t8\t6\t8\t9\n')
        self.cfg={'format':'biomodal_cxreport','input':'sites.tsv','reference':'ref.fa','assembly':'fixture','modification':'5hmC','chemistry':'duet_evoC','pipeline_version':'1.6.0','measurement_reviewed':True}
    def tearDown(self):self.temp.cleanup()
    def run_import(self,text):
        (self.base/'sites.tsv').write_text(text);p=self.base/'config.json';p.write_text(json.dumps(self.cfg))
        return mods.convert(p,self.base/'out')
    def test_reverse_strand_collapses_and_other_mods_stay_denominator(self):
        self.run_import('chr1\t2\t+\t2\t8\tCG\tCGT\nchr1\t3\t-\t3\t7\tCG\tCGT\n')
        self.assertEqual((self.base/'out/counts.tsv').read_text(),'chr\tpos\tN\tX\nchr1\t2\t20\t5\n')
    def test_modkit_other_modification_in_denominator(self):
        self.cfg.update(format='modkit_bedmethyl',chemistry='nanopore')
        self.run_import('chr1\t1\t2\th\t10\t.\t1\t2\t0\t10\t20\t2\t3\t5\t0\t0\t0\t0\n')
        self.assertIn('chr1\t2\t10\t2',(self.base/'out/counts.tsv').read_text())
    def test_five_base_cannot_be_called_hydroxymethylation(self):
        self.cfg['chemistry']='duet_modC'
        with self.assertRaisesRegex(ValueError,'does not distinguish'):self.run_import('')
    def test_non_cpg_fails(self):
        with self.assertRaises(ValueError):self.run_import('chr1\t4\t-\t2\t8\tCG\tCGT\n')
    def test_duplicate_sites_fail_and_mark_incomplete(self):
        with self.assertRaisesRegex(ValueError,'Duplicate'):self.run_import('chr1\t2\t+\t2\t8\tCG\tCGT\n'*2)
        self.assertFalse(json.loads((self.base/'out/status.json').read_text())['complete'])
    def test_zero_coverage_is_omitted(self):
        record=self.run_import('chr1\t2\t+\t0\t0\tCG\tCGT\nchr1\t5\t+\t2\t3\tCG\tCGA\n')
        self.assertEqual(record['zero_coverage_omitted'],1);self.assertEqual(record['sites_written'],1)

if __name__=='__main__':unittest.main(verbosity=2)
