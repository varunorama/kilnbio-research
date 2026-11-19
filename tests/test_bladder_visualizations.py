import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from auto_reports import render


class BladderVisualizations(unittest.TestCase):
    def test_comparisons_and_budget_are_preserved(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            pd.DataFrame([
                {'gene':g,'arm':arm,'rank':rank,'score':score}
                for arm,rows in {
                    'A':[('X',1,.8),('Y',2,-.4),('Z',3,.1)],
                    'B':[('X',2,.3),('Y',1,-.7),('Z',3,.1)],
                    'C':[('X',1,.6),('Y',2,-.3),('Z',3,.1)],
                    'P':[('X',1,.9),('Y',3,.1),('Z',2,.2)]}.items()
                for g,rank,score in rows]).to_csv(root/'rankings.tsv',sep='\t',index=False)
            pd.DataFrame([
                {'comparison':'A_to_B','gene':'X','status':'exited'},
                {'comparison':'A_to_B','gene':'Y','status':'entered'},
                {'comparison':'B_to_C','gene':'Y','status':'exited'},
                {'comparison':'B_to_C','gene':'X','status':'entered'},
                {'comparison':'A_to_C','gene':'X','status':'retained'}]).to_csv(root/'changes.tsv',sep='\t',index=False)
            (root/'report.json').write_text(json.dumps({'budget':1,'common_genes':3,'groups':{
                'RNA_all':{'patients':8,'positive':3,'no_qualifying_call':5},
                'matched':{'patients':6,'positive':2,'no_qualifying_call':4},
                'RNA_only':{'patients':2,'positive':1,'no_qualifying_call':1}}}))
            before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            out=render('bladder-example',[],root)
            self.assertEqual(before,{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in before})
            counts=pd.read_csv(out/'plotted_data/shortlist_status.tsv',sep='\t').set_index('comparison')
            self.assertEqual(counts.loc['A_to_B','entered'],1)
            self.assertEqual(counts.loc['A_to_C','retained'],1)
            returning=pd.read_csv(out/'plotted_data/bladder_returning_candidates.tsv',sep='\t')
            self.assertEqual(returning.gene.tolist(),['X'])
            ranks=pd.read_csv(out/'plotted_data/shortlist_rank_changes.tsv',sep='\t')
            self.assertEqual(set(ranks.loc[ranks.displayed,'gene']),{'X','Y'})
            evidence=pd.read_csv(out/'plotted_data/bladder_rank_evidence.tsv',sep='\t').set_index('gene')
            self.assertEqual(evidence.loc['X','C'],.6)
            self.assertEqual(evidence.loc['Y','P'],.1)
            receipt=json.loads((out/'visualization.json').read_text())
            self.assertIn('bladder_visualizations.py',receipt['report_source_sha256'])
            self.assertEqual(len(receipt['panels']),6)


if __name__=='__main__':unittest.main()
