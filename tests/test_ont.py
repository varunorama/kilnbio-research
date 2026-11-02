import importlib.util,json,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ont',ROOT/'longread/ont.py');ont=importlib.util.module_from_spec(spec);spec.loader.exec_module(ont)
class ONTPlanTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        for f in ('ref.fa','ref.fa.fai','calls.bam'):(self.base/f).write_text('fixture')
        self.cfg={'input_stage':'aligned_modbam','input':'calls.bam','reference':'ref.fa','assembly':'fixture','channels':['5mC','5hmC'],'models_reviewed':True,'duplicate_policy':'retain_pcr_free','threads':2,'filter_threshold':.8,'versions':{'dorado':'2.1.2','modkit':'0.6.4','samtools':'1.24'}}
    def tearDown(self):self.temp.cleanup()
    def plan(self):
        p=self.base/'config.json';p.write_text(json.dumps(self.cfg));return ont.plan(p,self.base/'out')
    def test_aligned_input_not_realigned_or_converted_to_fastq(self):
        _,out,record=self.plan();self.assertFalse(out.exists())
        argv=[step['argv'] for step in record['commands']]
        self.assertFalse(any(a[0]=='dorado' for a in argv));self.assertIn('--combine-strands',argv[-1]);self.assertIn('--filter-threshold',argv[-1])
    def test_uncertain_model_not_silently_selected(self):
        self.cfg['models_reviewed']=False
        with self.assertRaisesRegex(ValueError,'model compatibility'):self.plan()
    def test_total_modification_not_inferred(self):
        self.cfg['channels']=['modC']
        with self.assertRaisesRegex(ValueError,'channels'):self.plan()
if __name__=='__main__':unittest.main(verbosity=2)
