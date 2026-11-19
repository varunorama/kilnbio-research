import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'downstream'))
import multiomic_meta
import pathway_compare
from visual_reports import multiomic_report
from study import plan


class VisualizationTests(unittest.TestCase):
    def test_matrix_annotations_preserve_scale_label(self):
        import matplotlib.pyplot as plt
        from visual_reports import matrix
        fig,ax=plt.subplots()
        try:
            matrix(fig,ax,np.array([[1,.876],[.876,1]]),['A','B'],['A','B'],'Pearson correlation',1)
            self.assertEqual(fig.axes[-1].get_ylabel(),'Pearson correlation')
            self.assertIn('0.88',[t.get_text() for t in ax.texts])
        finally:plt.close(fig)

    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()

    def fixture(self,nlayers=3):
        layers=[]
        for i in range(nlayers):
            pd.DataFrame({'feature':['strong','opposite','no_signal','missing','<script>alert(1)</script>'],
                          'effect':[1,-1 if i==1 else 1,0,np.nan,1],
                          'p':[.0001,.001,1,np.nan,.2]}).to_csv(self.p/f'L{i}.tsv',sep='\t',index=False)
            layers.append(dict(layer=f'L{i}',results=f'L{i}.tsv',feature_map='identity',orientation='1'))
        pd.DataFrame(layers).to_csv(self.p/'manifest.tsv',sep='\t',index=False)
        return SimpleNamespace(manifest=str(self.p/'manifest.tsv'),out=str(self.p/'output'),min_layers=2,alpha=.05)

    def assert_report(self,out,panels):
        record=json.loads((out/'visualization.json').read_text())
        self.assertEqual(len(record['panels']),panels)
        for path,digest in record['artifact_sha256'].items():
            self.assertEqual(hashlib.sha256((out/path).read_bytes()).hexdigest(),digest)
        self.assertTrue((out/'visualizations.pdf').read_bytes().startswith(b'%PDF'))
        self.assertNotIn('src="https://',(out/'report.html').read_text())

    def test_reports_leave_statistical_outputs_unchanged(self):
        args=self.fixture()
        with patch('visual_reports.multiomic_report'):multiomic_meta.run(args)
        baseline=(self.p/'output/genes.tsv').read_bytes();nominations=(self.p/'output/nominations.tsv').read_bytes()
        args.out=str(self.p/'rendered');multiomic_meta.run(args)
        self.assertEqual(baseline,(self.p/'rendered/genes.tsv').read_bytes())
        self.assertEqual(nominations,(self.p/'rendered/nominations.tsv').read_bytes())
        self.assert_report(self.p/'rendered',5)
        elig=pd.read_csv(self.p/'rendered/plotted_data/rank_eligibility.tsv',sep='\t')
        self.assertEqual(list(elig.common_genes),[4,4,4])
        changes=pd.read_csv(self.p/'rendered/plotted_data/rank_sensitivity.tsv',sep='\t')
        self.assertTrue(np.isfinite(changes.rank_gain_with_assay).all())
        self.assertNotIn('<script>alert',(self.p/'rendered/report.html').read_text())
        conc=pd.read_csv(self.p/'rendered/plotted_data/concordance.tsv',sep='\t')
        self.assertTrue((conc[['same','opposite','neutral','not_both']].sum(axis=1)==5).all())

    def test_two_assays_do_not_fabricate_leave_one_out_gain(self):
        args=self.fixture(2);multiomic_meta.run(args)
        elig=pd.read_csv(self.p/'output/plotted_data/rank_eligibility.tsv',sep='\t')
        self.assertTrue((elig.status=='not_estimable_fewer_layers_than_r').all())
        self.assertTrue(pd.read_csv(self.p/'output/plotted_data/rank_sensitivity.tsv',sep='\t').empty)

    def test_missing_only_report_and_atomic_failure(self):
        args=self.fixture(2)
        for f in self.p.glob('L*.tsv'):
            pd.DataFrame({'feature':['G'],'effect':[np.nan],'p':[np.nan]}).to_csv(f,sep='\t',index=False)
        result=multiomic_meta.run(args);self.assertEqual(result['testable_genes'],0)
        args.out=str(self.p/'failed')
        with patch('visual_reports.multiomic_report',side_effect=ValueError('render failed')):
            with self.assertRaisesRegex(ValueError,'render failed'):multiomic_meta.run(args)
        self.assertFalse((self.p/'failed').exists());self.assertFalse(list(self.p.glob('.multiomic-*')))

    def pathway_fixture(self):
        spec=importlib.util.spec_from_file_location('pathway_demo',ROOT/'examples/pathways/make_demo.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module.create(self.p/'pathways')

    def test_pathway_example_and_study_plan(self):
        data=self.pathway_fixture();out=self.p/'comparison';summary=pathway_compare.run(data/'compare.json',out)
        self.assertEqual(summary['pathway_family'],5);self.assertEqual(summary['selected_pathways'],2)
        result=pd.read_csv(out/'pathways.tsv',sep='\t').set_index('pathway')
        self.assertAlmostEqual(result.loc['SYNTHETIC_RESPONSE','p_pc'],.0004)
        # Opposing NES can still support a NONdirectional pathway test.
        self.assertLess(result.loc['SYNTHETIC_DISCORDANT','q_pc'],.05)
        self.assertTrue(np.isnan(result.loc['SYNTHETIC_UNMEASURED','p_pc']))
        self.assert_report(out,4);plan(data/'study.json',self.p/'study')

    def test_pathway_rejects_selected_only_and_wrong_gene_support(self):
        data=self.pathway_fixture();path=data/'RNA.tsv';tab=pd.read_csv(path,sep='\t')
        tab.iloc[:2].to_csv(path,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'every GMT'):pathway_compare.run(data/'compare.json',self.p/'bad')
        tab.loc[0,'contributing_genes']='OUTSIDE';tab.to_csv(path,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'Contributing genes'):pathway_compare.run(data/'compare.json',self.p/'bad')

    def test_pc_fixed_family_and_null_behavior(self):
        p=np.array([[.001,.002,np.nan],[.001,np.nan,np.nan],[.2,.3,.4]])
        combined,q=pathway_compare.partial_conjunction(p,2)
        self.assertAlmostEqual(combined[0],.004);self.assertTrue(np.isnan(combined[1]));self.assertAlmostEqual(combined[2],.6)
        self.assertAlmostEqual(q[0],.004*3*(1+1/2+1/3))
        # At the boundary null, one non-null may be arbitrarily strong; two
        # identical uniform null p values impose perfect dependence.
        u=np.linspace(.000005,.999995,100000)
        combined,_=pathway_compare.partial_conjunction(np.column_stack([np.full(len(u),1e-10),u,u]),2)
        self.assertLessEqual(float((combined<=.05).mean()),.05)
        reordered,_=pathway_compare.partial_conjunction(p[::-1,::-1],2)
        np.testing.assert_allclose(reordered[::-1],pathway_compare.partial_conjunction(p,2)[0],equal_nan=True)


if __name__=='__main__':unittest.main()
