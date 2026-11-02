import csv
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"downstream"))
from graft_readouts import run


class GraftReadoutTests(unittest.TestCase):
    def test_serial_pools_missing_recipients_and_technical_replicates(self):
        base=ROOT/"examples"/"graft_readouts"
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"result"
            audit=run(base/"lots.tsv",base/"grafts.tsv",base/"readouts.tsv",
                      ("hsc_identity",15,"culture"),("hCD45_of_live",140,"marrow"),"secondary",out)
            self.assertEqual((audit["lots_with_both"],audit["distinct_donors_with_both"]),(2,1))
            self.assertEqual(audit["multi_lot_recipients_excluded"],1)
            with (out/"lots.tsv").open() as f:
                a=next(csv.DictReader(f,delimiter="\t"))
            self.assertEqual([float(a[k]) for k in ("predictor","enrolled_recipients","measured_recipients","median_outcome")],[.5,2,1,15])

    def test_incompatible_units_fail_before_publication(self):
        base=ROOT/"examples"/"graft_readouts"
        with tempfile.TemporaryDirectory() as tmp:
            inp=Path(tmp)/"input.tsv"
            inp.write_text((base/"readouts.tsv").read_text().replace("percent","fraction",1))
            out=Path(tmp)/"result"
            with self.assertRaises(ValueError):
                run(base/"lots.tsv",base/"grafts.tsv",inp,("hsc_identity",15,"culture"),("hCD45_of_live",140,"marrow"),"primary",out)
            self.assertFalse(out.exists())
