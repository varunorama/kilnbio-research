"""Offline contract tests against the real pinned SDK; service calls are mocked."""
import json
from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest.mock import patch,create_autospec
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client,gene_counts_condition
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from downstream.alphagenome_query import read_plan,prepare_cells,execute,SCORERS
from alphagenome.models.v2 import variant_scorers


class SDKContract(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.root=Path(__file__).resolve().parents[1]
        self.plan=read_plan(self.root/'examples/alphagenome/reference.json')
    def tearDown(self):self.tmp.cleanup()
    def cells(self,x,genes=None):
        a=ad.AnnData(X=sparse.csr_matrix(x),obs=pd.DataFrame({'donor':['PRIVATE']*len(x)},index=[f'private_{i}' for i in range(len(x))]),
                     var=pd.DataFrame({'private_annotation':['secret','secret']},index=genes or ['ENSG00000000001.1','ENSG00000000002']))
        f=self.p/'cells.h5ad';a.write_h5ad(f);return f
    def test_sdk_coordinates_and_scorers(self):
        self.assertTrue(SCORERS <= set(variant_scorers.RECOMMENDED_VARIANT_SCORERS))
        for pos in [1000000,1000001]:
            for ref in ['A','AG','AGG']:
                v=genome.Variant(chromosome='chr1',position=pos,reference_bases=ref,alternate_bases='C')
                win=v.reference_interval.resize(1048576)
                self.assertEqual(win.start,(2*(pos-1)+len(ref)+1)//2-524288)
    def test_cells_remove_metadata_and_serialize(self):
        a=prepare_cells(self.cells([[0,2],[3,0]]))
        self.assertEqual(list(a.obs.columns),[]);self.assertEqual(list(a.var.columns),[])
        self.assertEqual(list(a.obs_names),['cell_0','cell_1'])
        p=gene_counts_condition.from_anndata(a)
        self.assertEqual(len(p.cells),2)
        self.assertEqual(list(p.cells[0].gene_counts),[2])
    def test_fractional_empty_and_duplicate_counts_rejected(self):
        for x,genes in [([[.2,2]],None),([[0,0]],None),([[1,2]],['ENSG00000000001.1','ENSG00000000001.2']),([[-1,2]],None)]:
            with self.assertRaises(ValueError):prepare_cells(self.cells(x,genes))
    def test_explicit_upload_gate_precedes_client(self):
        plan=dict(self.plan,cell_counts=str(self.cells([[1,2]])),cell_condition_name='cells')
        with patch.object(dna_client,'create') as create:
            with self.assertRaisesRegex(ValueError,'allow-cell-upload'):execute(plan,self.p/'out')
            create.assert_not_called()
    def test_completion_resume_and_changed_inputs(self):
        client=create_autospec(dna_client.DnaClient,instance=True)
        client.score_variant.return_value=[]
        table=pd.DataFrame({'raw_score':[.12],'gene_name':['GLB1L']})
        out=self.p/'out'
        with patch.dict(os.environ,{'ALPHA_GENOME_API_KEY':'offline-test-only'}),patch.object(dna_client,'create',return_value=client) as create,patch.object(variant_scorers,'tidy_scores',side_effect=lambda _:table.copy()):
            execute(self.plan,out)
            self.assertEqual(json.loads((out/'receipt.json').read_text())['status'],'complete')
            self.assertTrue((out/'scores.tsv').exists());self.assertEqual(client.score_variant.call_count,1)
            execute(self.plan,out,resume=True)
            self.assertEqual(client.score_variant.call_count,1);self.assertEqual(create.call_count,1)
            with self.assertRaisesRegex(ValueError,'same inputs'):execute(dict(self.plan,condition='CL:0000000'),out,resume=True)
    def test_failed_query_is_not_success_or_secret_leak(self):
        client=create_autospec(dna_client.DnaClient,instance=True)
        client.score_variant.side_effect=ValueError('fake-secret-from-remote')
        with patch.dict(os.environ,{'ALPHA_GENOME_API_KEY':'offline-test-only'}),patch.object(dna_client,'create',return_value=client):
            with self.assertRaises(RuntimeError):execute(self.plan,self.p/'out')
        text=(self.p/'out'/'receipt.json').read_text()
        self.assertEqual(json.loads(text)['status'],'incomplete')
        self.assertNotIn('fake-secret',text);self.assertFalse((self.p/'out'/'scores.tsv').exists())


if __name__=='__main__':unittest.main(verbosity=2)
