"""Offline planning and runtime-policy regressions; no prediction credentials."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from downstream.alphagenome_query import read_plan
from ci.check_runtime_policy import inspect,scan


class OfflinePlan(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.table=self.p/'v.csv'
        self.table.write_text('variant_id,chromosome,position,ref,alt\nv1,chr2,219243701,C,T\n',encoding='utf-8')
        self.c={'variants':'v.csv','assembly':'GRCh38','condition':'UBERON:0000178_CL:0001054','scorers':['RNA_SEQ','ATAC']}
    def tearDown(self):self.tmp.cleanup()
    def plan(self):
        p=self.p/'c.json';p.write_text(json.dumps(self.c),encoding='utf-8');return read_plan(p)
    def test_offline_coordinates_and_input_hash(self):
        p=self.plan();v=p['variants'][0]
        self.assertEqual((v['interval_start'],v['interval_end']),(218719413,219767989))
        self.assertFalse(p['prediction_is_independent_validation'])
        self.table.write_text(self.table.read_text()+'v2,chr2,219243702,A,G\n')
        self.assertNotEqual(p['plan_sha256'],self.plan()['plan_sha256'])
    def test_no_wrong_assembly_or_unknown_configuration(self):
        self.c['assembly']='GRCh37'
        with self.assertRaisesRegex(ValueError,'GRCh38'):self.plan()
        self.c['assembly']='GRCh38';self.c['api_key']='not-a-real-key'
        with self.assertRaisesRegex(ValueError,'credentials'):self.plan()
    def test_canonical_alleles_and_boundary(self):
        for row in ['v1,chr2,1,C,T','v1,chr2,219243701,C,<DEL>','v1,chr2,219243701,C,C','v1,2,219243701,C,T']:
            self.table.write_text('variant_id,chromosome,position,ref,alt\n'+row+'\n')
            with self.assertRaises(ValueError):self.plan()
    def test_duplicates_and_request_budget(self):
        self.table.write_text(self.table.read_text()+'v2,chr2,219243701,C,T\n')
        with self.assertRaisesRegex(ValueError,'Duplicate genomic'):self.plan()
        self.c['max_variants']=1
        with self.assertRaisesRegex(ValueError,'max_variants'):self.plan()
    def test_cell_only_scorer_requires_cells(self):
        self.c['scorers']=['SC_GENE_EXPRESSION_DENOISED']
        with self.assertRaisesRegex(ValueError,'requires raw'):self.plan()
    def test_runtime_guard_detects_known_integrations(self):
        p=self.p/'bad.py';p.write_text('import openai\n')
        self.assertTrue(inspect(p))
        p.write_text('import numpy\n');self.assertEqual(inspect(p),[])
        p.write_text('import subprocess\nsubprocess.run(["claude", "-p", "x"])\n')
        self.assertTrue(inspect(p))
    def test_shipped_runtime_has_no_known_llm_calls(self):
        self.assertEqual(scan(Path(__file__).resolve().parents[1]),[])


if __name__=='__main__':unittest.main(verbosity=2)
