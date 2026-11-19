import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'downstream'))
from circos_report import load_config, run
from multiomic_meta import run as integrate
spec=importlib.util.spec_from_file_location('circos_demo',ROOT/'examples/circos/make_demo.py')
demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)


class CircosTests(unittest.TestCase):
    def test_tracks_follow_supplied_assays_and_coordinates(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);cp=demo.create(root/'input');cfg=json.loads(cp.read_text())
            run(cp,root/'five')
            cfg['tracks']=cfg['tracks'][:2];cfg.pop('links');cp.write_text(json.dumps(cfg))
            run(cp,root/'two')
            for name,n in [('five',5),('two',2)]:
                receipt=json.loads((root/name/'visualization.json').read_text())
                self.assertEqual(len(receipt['details']['tracks']),n)
                t=pd.read_csv(root/name/'plotted_data/circos_tracks.tsv',sep='\t')
                self.assertEqual(t.track.nunique(),n)
                self.assertEqual(len(t),n*480)
                self.assertTrue(t.value.isna().any())
            self.assertEqual(json.loads((root/'two/visualization.json').read_text())['details']['links'],0)

    def test_reject_build_mismatch_and_invalid_coordinate(self):
        with tempfile.TemporaryDirectory() as td:
            cp=demo.create(Path(td)/'input');cfg=json.loads(cp.read_text())
            cfg['tracks'][0]['genome_build']='other-build';cp.write_text(json.dumps(cfg))
            with self.assertRaisesRegex(ValueError,'genome build'):load_config(cp)
            cfg['tracks'][0]['genome_build']=cfg['genome_build'];cp.write_text(json.dumps(cfg))
            t=pd.read_csv(cp.parent/'RNA.tsv',sep='\t');t.loc[0,'start']=-1;t.to_csv(cp.parent/'RNA.tsv',sep='\t',index=False)
            with self.assertRaisesRegex(ValueError,'out-of-bounds'):load_config(cp)

    def test_optional_panel_does_not_change_integration_statistics(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);cp=demo.create(root/'input');cfg=json.loads(cp.read_text());cfg['tracks']=cfg['tracks'][:2];cp.write_text(json.dumps(cfg))
            for layer in ['RNA','protein']:
                pd.DataFrame({'feature':['A','B','C'],'effect':[1,-1,1],'p':[.001,.02,.7]}).to_csv(root/(layer+'.tsv'),sep='\t',index=False)
            pd.DataFrame([dict(layer=n,results=n+'.tsv',feature_map='identity',orientation=1) for n in ['RNA','protein']]).to_csv(root/'manifest.tsv',sep='\t',index=False)
            for name,config in [('plain',None),('circular',str(cp))]:
                integrate(SimpleNamespace(manifest=str(root/'manifest.tsv'),out=str(root/name),min_layers=2,alpha=.05,circos_config=config))
            self.assertEqual((root/'plain/genes.tsv').read_bytes(),(root/'circular/genes.tsv').read_bytes())
            self.assertEqual((root/'plain/nominations.tsv').read_bytes(),(root/'circular/nominations.tsv').read_bytes())
            self.assertTrue((root/'circular/figures/circos_tracks.svg').is_file())


if __name__=='__main__':unittest.main()
