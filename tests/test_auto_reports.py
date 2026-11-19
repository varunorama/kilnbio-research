import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'downstream'))
from auto_reports import output_for,render
from analysis_visualizations import Builder,read,methods
import numpy as np
import pandas as pd


class AutomaticReports(unittest.TestCase):
    def test_plans_and_catalogs_have_no_biological_report(self):
        for c,a in [('de',['cfg','out','--dry-run']),('clocks',['catalog','--out','x']),('assay-clocks',['features','--out','x']),('alphagenome',['cfg','out']),('vep-annotate',['cfg','out']),('erythroid-stage',['--out','x','--list-features'])]:
            self.assertIsNone(output_for(c,a),(c,a))
        self.assertEqual(output_for('genomic-associate',['a','b','c','d','specimen','out']),Path('out').resolve())
        self.assertEqual(output_for('capture-panel',['--out=panel']),Path('panel').resolve())

    def test_bounded_input_and_missing_values(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'x.tsv';pd.DataFrame({'a':[1,np.nan,3,4]}).to_csv(p,sep='\t',index=False)
            t=read(p,limit=2);self.assertTrue(t.attrs['truncated']);self.assertEqual(len(t),2);self.assertTrue(pd.isna(t.a.iloc[1]))

    def test_custom_agreement_threshold_is_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            b=Builder(Path(td),'method-agreement');t=pd.DataFrame({'feature':['A'],'one_q':[.02],'two_q':[.001],'one_effect':[1],'two_effect':[1]})
            methods(b,t,alpha=.01);b.finish()
            values=pd.read_csv(Path(td)/'plotted_data/method_calls.tsv',sep='\t')
            self.assertEqual(values.one.iloc[0],0);self.assertEqual(values.two.iloc[0],1)

    def test_matrix_loading_limits_columns(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'matrix.tsv';pd.DataFrame(np.ones((2,350)),columns=['feature',*[f'S{i}' for i in range(349)]]).to_csv(p,sep='\t',index=False)
            b=Builder(Path(td),'matrix');t=b.load_matrix(p);b.finish()
            self.assertEqual(t.shape,(2,301));self.assertEqual(b.sources[0]['columns_available'],350)

    def test_report_keeps_statistics_and_records_real_missingness(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'out';(root/'results').mkdir(parents=True)
            p=root/'results/all_features.tsv';pd.DataFrame({'feature':['A','B','C'],'effect':[1,-1,np.nan],'p':[.001,.9,np.nan],'q':[.003,.9,np.nan]}).to_csv(p,sep='\t',index=False)
            digest=hashlib.sha256(p.read_bytes()).hexdigest();report=render('de',[],root)
            self.assertEqual(digest,hashlib.sha256(p.read_bytes()).hexdigest())
            r=json.loads((report/'visualization.json').read_text());self.assertIn('all_features_volcano',{x['name'] for x in r['panels']})
            data=pd.read_csv(report/'plotted_data/all_features_volcano.tsv',sep='\t');self.assertTrue(pd.isna(data.loc[2,'display_p']))
            for f,sha in r['artifact_sha256'].items():self.assertEqual(hashlib.sha256((report/f).read_bytes()).hexdigest(),sha)

    def test_failed_report_is_not_published(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'analysis';root.mkdir();(root/'result.tsv').write_text('real analysis\n')
            with patch('report_profiles.build',side_effect=ValueError('invalid plotting input')):
                with self.assertRaisesRegex(RuntimeError,'outputs are preserved'):render('de',[],root)
            self.assertTrue((root/'result.tsv').exists());self.assertFalse((root/'visualization').exists());self.assertEqual(len(list(root.glob('.visualization.pending-*/report_failure.json'))),1)

    @unittest.skipIf(os.name=='nt','Bash launcher runs on Linux/WSL')
    def test_command_and_study_generate_and_hash_reports(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);cfg={'schema_version':1,'research_use':True,'inputs':[str(ROOT/'examples/assay_roster.tsv')],'configs':{},'steps':[{'id':'audit','command':'assay-audit','args':['reach',str(ROOT/'examples/assay_roster.tsv'),'--out','{{run}}/audit/result.json'],'depends_on':[],'outputs':['audit/result.json']}]}
            # Direct command test covers the public launch boundary independently of study schema.
            out=td/'audit.json';env=dict(os.environ,KILNBIO_PYTHON=sys.executable)
            subprocess.run(['bash',str(ROOT/'run.sh'),'assay-audit','reach',str(ROOT/'examples/assay_roster.tsv'),'--out',str(out)],env=env,check=True,capture_output=True,text=True)
            self.assertTrue(Path(str(out)+'.visualization/report.html').is_file())
            from study import execute
            cp=td/'study.json';cp.write_text(json.dumps(cfg))
            with patch.dict(os.environ,env):state=execute(cp,td/'study')
            hashed=state['steps']['audit']['outputs_sha256']
            self.assertIn('audit/result.json.visualization/report.html',hashed)
            report=td/'study/audit/result.json.visualization/report.html';report.write_text('changed')
            with patch.dict(os.environ,env):
                with self.assertRaisesRegex(ValueError,'output changed'):execute(cp,td/'study',resume=True)


if __name__=='__main__':unittest.main()
