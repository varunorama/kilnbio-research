"""Statistical oracles, failure cases and workflow integration for public functions."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'downstream'))
import decision_methods as dm
import study


class DecisionTests(unittest.TestCase):
    def test_signed_rank_cancellation_and_ties(self):
        result = dm.rank_assays(['B','A','C'], {'rna':[4,-4,1], 'protein':[-4,4,1]}, 2, 'signed')
        np.testing.assert_allclose(result['arms']['rna']['scores'], [.625,-.625,.25])
        np.testing.assert_allclose(result['arms']['combined']['scores'], [0,0,.25])
        self.assertEqual(result['arms']['combined']['ranking'], ['C','A','B'])

    def test_unsigned_ignores_sign_and_preserves_input_ties(self):
        a = dm.rank_assays(['B','A','C'], {'rna':[4,-4,1], 'meth':[-4,4,1]}, 2, 'unsigned')
        self.assertEqual(a['arms']['combined']['nominations'], ['B','A'])
        np.testing.assert_allclose(a['arms']['combined']['scores'], [5/6,5/6,1/3])

    def test_rank_rejects_changing_universe_missing_and_bad_budget(self):
        for args in [(['A','A'], {'rna':[1,2]}, 1, 'signed'),
                     (['A','B'], {'rna':[1,None]}, 1, 'signed'),
                     (['A','B'], {'rna':[1]}, 1, 'signed'),
                     (['A','B'], {'rna':[1,2]}, 3, 'signed')]:
            with self.assertRaises(ValueError): dm.rank_assays(*args)

    def test_nonadditive_replacements(self):
        r = dm.compare_nominations({'A':['x','y','z'],'B':['y','z','x'],'C':['x','y','z']}, 1, ['A','B','C'])
        self.assertEqual([t['replacements'] for t in r['transitions']], [1,1])
        self.assertEqual(sum(t['rank_change']['x'] for t in r['transitions']), 0)
        with self.assertRaises(ValueError): dm.compare_nominations({'A':['x'],'B':['z']},1,['A','B'])

    def test_validation_budget_keeps_missing_and_wrong_direction(self):
        r = dm.validation_yield(['a','b','c'], [1,1,-1], [2,None,2], [.01,None,.001], [True,False,True])
        self.assertEqual((r['budget'],r['testable'],r['confirmed']), (3,2,1))
        self.assertAlmostEqual(r['yield_per_slot'], 1/3)
        with self.assertRaises(ValueError): dm.validation_yield(['a'],[1],[1],[.01],[False])

    def test_donor_partition_groups_and_order(self):
        r = dm.donor_partition(['D1','D2','D1'], 'locked:')
        self.assertEqual(r['fold'][0], r['fold'][2])
        rev = dm.donor_partition(['D2','D1'], 'locked:')
        self.assertEqual(rev['fold'], r['fold'][1:])
        with self.assertRaises(ValueError): dm.donor_partition(['D1'],'x',5,[0,1,2,3,4])

    def test_promoter_fixed_map_transform_and_missing(self):
        r = dm.promoter_summary(['p1','p2'], [[.2,.5],[.8,None]], {'g':['p1','p2']}, 'beta')
        self.assertAlmostEqual(r['values'][0][0], 0)
        self.assertFalse(r['testable'][0][1])
        missing = dm.promoter_summary(['p1','p2'], [[.2,.5],[.8,.5]], {'g':['p1','absent']}, 'beta')
        self.assertFalse(any(missing['testable'][0]))
        self.assertEqual(missing['missing_mapped_probes'], [['absent']])

    def test_promoter_median_and_leave_one_out(self):
        args = (['p1','p2','p3'], [[0],[1],[20]], {'g':['p1','p2','p3']}, 'M')
        self.assertEqual(dm.promoter_summary(*args,method='median')['values'], [[1]])
        self.assertEqual(dm.promoter_summary(*args,drop_probe='p3')['values'], [[.5]])
        two = dm.promoter_summary(['p1','p2'],[[0],[1]],{'g':['p1','p2']},'M',drop_probe='p1')
        self.assertFalse(two['testable'][0][0])
        with self.assertRaises(ValueError): dm.promoter_summary(['p'],[[1.2]],{'g':['p','q']},'beta')

    def test_hc3_against_independent_statsmodels(self):
        import statsmodels.api as sm
        rng = np.random.default_rng(418)
        x = np.column_stack([np.ones(45),rng.normal(size=(45,2))])
        y = x @ rng.normal(size=(3,4)) + rng.normal(size=(45,4))*(1+abs(x[:,1,None]))
        c = np.array([0,1,-.5])
        r = dm.fit_hc3(list('abcd'),[f'd{i}' for i in range(45)],y,x,c)
        for j in range(4):
            oracle = sm.OLS(y[:,j],x).fit(cov_type='HC3',use_t=True).t_test(c)
            np.testing.assert_allclose([r['effect'][j],r['se'][j],r['t'][j],r['p'][j]],
                                      [oracle.effect.item(),oracle.sd.item(),oracle.tvalue.item(),oracle.pvalue.item()],rtol=1e-9,atol=1e-12)

    def test_hc3_refuses_invalid_design_and_leverage(self):
        y = np.arange(16.)[:,None]; ids=[str(i) for i in range(16)]
        with self.assertRaises(ValueError): dm.fit_hc3(['g'],ids,y,np.ones((16,2)),[0,1])
        x=np.column_stack([np.ones(16),[1]+[0]*15])
        with self.assertRaises(ValueError): dm.fit_hc3(['g'],ids,y,x,[0,1])
        with self.assertRaises(ValueError): dm.fit_hc3(['g'],ids,y,np.ones((16,1)),[0])

    def test_locked_conjunction_matches_max_tail_and_by(self):
        from statsmodels.stats.multitest import multipletests
        r = dm.directional_conjunction(['a','b','c'],{'R':[4,-3,None],'M':[-5,2,1]},
              {'R':20,'M':15},{'R':[1,-1,1],'M':[-1,-1,1]},
              {'R':[True,True,False],'M':[True,True,True]},'BY')
        p = np.maximum([stats.t.sf(4,20),stats.t.sf(3,20),1],stats.t.sf([5,-2,1],15))
        np.testing.assert_allclose(r['conjunction_p'],p)
        np.testing.assert_allclose(r['q'],multipletests(p,method='fdr_by')[1])
        self.assertEqual(r['testable'],[True,True,False])
        with self.assertRaises(ValueError):
            dm.directional_conjunction(['a'],{'R':[4],'M':[4]},{'R':20,'M':20},
                                       {'R':[0],'M':[1]},{'R':[True],'M':[True]})

    def test_summary_sensitivity_has_controlled_denominator(self):
        r = dm.promoter_sensitivity(['a','b','c'],[3,2,1],[3,2,1],[1,2,3],1)
        self.assertEqual(r['comparison']['replacements'],0)
        self.assertEqual(r['mean']['nominations'],['a'])
        with self.assertRaises(ValueError): dm.promoter_sensitivity(['a','b'],[3,2],[3,2],[1,None],1)

    def test_probe_coherence_retains_missing_and_zero(self):
        r = dm.probe_coherence(['a','b','c','d'],[2,-1,0,None],1)
        self.assertEqual((r['same_direction'],r['opposite_direction'],r['zero_effect'],r['missing']),(1,1,1,1))
        self.assertFalse(r['all_declared_same_direction'])
        self.assertAlmostEqual(r['fraction_of_tested'],1/3)

    def test_element_weights_within_replicate_not_across(self):
        r = dm.element_methylation([[0,0,0],[10,10,10]],[[5,5,5],[10,10,10]])
        self.assertAlmostEqual(r['methylation'],.5)  # pooled fraction would be 2/3
        no = dm.element_methylation([[0,0,0],[10,10,1]],[[5,5,5],[10,10,1]])
        self.assertFalse(no['testable']); self.assertIsNone(no['methylation'])
        with self.assertRaises(ValueError): dm.element_methylation([[2,0,0],[0,0,0]],[[1,5,5],[5,5,5]])

    def test_cli_output_provenance_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result.json'
            cmd=[sys.executable,str(ROOT/'downstream/decision_methods.py'),'run',str(ROOT/'examples/decision_rank.json'),'--out',str(out)]
            proc=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(proc.returncode,0,proc.stderr)
            payload=json.loads(out.read_text())
            self.assertEqual(payload['result']['arms']['combined']['nominations'],['GENE_A','GENE_C'])
            self.assertEqual(len(payload['input_sha256']),64)
            before=out.read_bytes()
            self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,2)
            self.assertEqual(out.read_bytes(),before)

    def test_study_plan_accepts_function_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=study.plan(ROOT/'examples/decision_study.json',Path(tmp)/'run')
            self.assertEqual(p['steps'][0]['command'],'decision-methods')

    @unittest.skipUnless(os.name == 'posix' and shutil.which('bash'), 'Linux/WSL dispatcher integration')
    def test_real_study_execution_resume_and_input_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)
            shutil.copy(ROOT/'examples/decision_rank.json',base/'decision_rank.json')
            shutil.copy(ROOT/'examples/decision_study.json',base/'decision_study.json')
            cmd=[sys.executable,str(ROOT/'downstream/study.py'),str(base/'decision_study.json'),str(base/'run')]
            env=dict(os.environ,KILNBIO_PYTHON=sys.executable)
            first=subprocess.run(cmd,env=env,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stdout+first.stderr)
            out=base/'run/rank/result.json'; before=out.read_bytes()
            again=subprocess.run(cmd+['--resume'],env=env,capture_output=True,text=True)
            self.assertEqual(again.returncode,0,again.stdout+again.stderr)
            self.assertEqual(out.read_bytes(),before)
            spec=json.loads((base/'decision_rank.json').read_text());spec['arguments']['k']=1
            (base/'decision_rank.json').write_text(json.dumps(spec))
            changed=subprocess.run(cmd+['--resume'],env=env,capture_output=True,text=True)
            self.assertNotEqual(changed.returncode,0)
            self.assertEqual(out.read_bytes(),before)

    def test_catalog_covers_all_dispatched_functions(self):
        c=dm.catalog()
        self.assertEqual(len(c['methods']),10)
        for name,context in c['methods'].items():
            self.assertTrue(callable(getattr(dm,name)))
            self.assertTrue(all(context[k] for k in ('purpose','assumptions','interpretation')))


if __name__ == '__main__': unittest.main()
