import sys, unittest, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"downstream"))
from method_agreement import agreement,run

class Agreement(unittest.TestCase):
    def test_same_fraction_different_evidence(self):
        a=agreement({"A":(1,.01),"B":(2,.02),"C":(None,None)},.05)
        b=agreement({"A":(1,.01),"B":(2,.02),"C":(-3,.001)},.05)
        self.assertEqual(a[5],b[5]);self.assertAlmostEqual(a[5],2/3)
        self.assertEqual((a[1],a[7],b[1],b[7]),(2,"partial_up",3,"opposite_calls"))
    def test_different_chemistry_fails(self):
        with tempfile.TemporaryDirectory() as d:
            m=Path(d)/"m.tsv";x=Path(d)/"x.tsv"
            m.write_text("method\tstudy\tcontrast\tfeature_family\tfeature_unit\tassembly\tchannel\nA\ts\tc\tf\tCpG\thg38\tmodC\nB\ts\tc\tf\tCpG\thg38\t5hmC\n")
            x.write_text("feature\tA_effect\tA_q\tB_effect\tB_q\nx\t1\t0.001\t1\t0.01\n")
            with self.assertRaisesRegex(ValueError,"channel"):run(x,m,.05)
