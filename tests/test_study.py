import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from study import execute, plan


class StudyTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        (self.p/'input.tsv').write_text('data\n')
        self.cfg={'schema_version':1,'research_use':True,'inputs':['input.tsv'],'configs':{},'steps':[
            {'id':'first','command':'variant-annotations','args':['config','{{run}}/first'],'depends_on':[],'outputs':['first/data.tsv']},
            {'id':'second','command':'variant-annotations','args':['config','{{run}}/second'],'depends_on':['first'],'outputs':['second/data.tsv']}]}
        self.cp=self.p/'study.json';self.out=self.p/'run'
    def tearDown(self):self.tmp.cleanup()
    def write(self):self.cp.write_text(json.dumps(self.cfg));return self.cp
    def fake(self,argv,**kwargs):
        p=Path(argv[-1]);p.mkdir();(p/'data.tsv').write_text('result\n')
        return subprocess.CompletedProcess(argv,0)
    def test_plan_does_not_execute(self):
        with patch('study.subprocess.run',side_effect=AssertionError('launched')):
            self.assertEqual(execute(self.write(),self.out,True)['status'],'planned')
    def test_dependency_order_and_resume_skips_verified_outputs(self):
        self.cfg['steps'].reverse()
        with patch('study.subprocess.run',side_effect=self.fake) as runner:
            r=execute(self.write(),self.out)
            self.assertEqual(list(r['steps']),['first','second'])
            self.assertEqual(runner.call_count,2)
        with patch('study.subprocess.run',side_effect=AssertionError('reran')):
            self.assertEqual(execute(self.cp,self.out,resume=True)['status'],'complete')
    def test_changed_input_refuses_resume(self):
        with patch('study.subprocess.run',side_effect=self.fake):execute(self.write(),self.out)
        (self.p/'input.tsv').write_text('changed\n')
        with self.assertRaisesRegex(ValueError,'changed'):execute(self.cp,self.out,resume=True)
    def test_changed_output_refuses_resume(self):
        with patch('study.subprocess.run',side_effect=self.fake):execute(self.write(),self.out)
        (self.out/'first/data.tsv').write_text('corrupted\n')
        with self.assertRaisesRegex(ValueError,'output changed'):execute(self.cp,self.out,resume=True)
    def test_failure_blocks_descendants_and_preserves_failed_attempt(self):
        def fail(argv,**kwargs):
            r=self.fake(argv,**kwargs);r.returncode=1;return r
        with patch('study.subprocess.run',side_effect=fail) as runner:
            with self.assertRaisesRegex(ValueError,'Step failed'):execute(self.write(),self.out)
            self.assertEqual(runner.call_count,1)
        self.assertFalse((self.out/'second').exists())
        with patch('study.subprocess.run',side_effect=self.fake):r=execute(self.cp,self.out,resume=True)
        self.assertEqual(r['steps']['first']['attempt'],2)
        self.assertEqual(len(list((self.out/'attempts/first').glob('*/data.tsv'))),1)
    def test_cycles_and_unknown_commands_rejected(self):
        self.cfg['steps'][0]['depends_on']=['second']
        with self.assertRaisesRegex(ValueError,'cycle'):plan(self.write(),self.out)
        self.cfg['steps'][0]['depends_on']=[];self.cfg['steps'][0]['command']='arbitrary-shell'
        with self.assertRaisesRegex(ValueError,'allowed'):plan(self.write(),self.out)
    def test_output_escape_rejected(self):
        self.cfg['steps'][0]['outputs']=['../outside.tsv']
        with self.assertRaisesRegex(ValueError,'under their step'):plan(self.write(),self.out)
    def test_remote_model_explicit_opt_in(self):
        self.cfg['steps'][0]['command']='alphagenome';self.cfg['steps'][0]['args'].append('--execute')
        with self.assertRaisesRegex(ValueError,'Remote scientific'):plan(self.write(),self.out)
        self.cfg['allow_remote_scientific_services']=True
        self.assertTrue(plan(self.write(),self.out)['allow_remote_scientific_services'])
    def test_template_resolution_and_fingerprint(self):
        (self.p/'template.json').write_text('{"input":"{{study}}/input.tsv","output":"{{run}}/first"}')
        self.cfg['configs']={'contrast':'template.json'}
        p=plan(self.write(),self.out)
        self.assertEqual(Path(p['configs']['contrast']['input']),self.p/'input.tsv')
        self.assertIn(str(self.p/'template.json'),p['inputs_sha256'])


if __name__=='__main__':unittest.main(verbosity=2)
