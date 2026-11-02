import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"downstream"))
from assay_audit import pairs,policies
from candidate_evidence import summarize

class Audits(unittest.TestCase):
    def test_cross_tissue_ambiguous_and_unknown(self):
        roster=[{"donor":x} for x in "ABCD"]
        def row(d,s,t,a,i):return dict(donor=d,specimen=s,tissue=t,assay=a,measurement=str(i))
        r=pairs(roster,[row("A","a","tum","r",1),row("A","a","tum","m",2),
                        row("B","b","tum","r",3),row("B","b","tum","m",4),row("B","b","tum","r",5),
                        row("C","c","normal","r",6),row("C","cc","tum","m",7),
                        row("D","","tum","m",8),row("D","","tum","r",9)],"r","m","tum")
        self.assertEqual([r["summary"][k] for k in ("both_any_tissue","both_requested_tissue","both_same_specimen","unique_one_to_one_donors","ambiguous_donors")],[4,3,2,1,1])
        self.assertEqual([v["donor"] for v in r["unique_pairs"]],["A"])
    def test_coverage_changes_apparent_gain(self):
        roster=[dict(donor=str(i),base_ok=b,add_ok=a) for i,(b,a) in enumerate([("1","1"),("1","0"),("0","0"),("1","1")])]
        sc=[dict(donor=str(i),observed=y,baseline=b,augmented=a) for i,(y,b,a) in enumerate([("0","3","1"),("7","7",""),("6","",""),("","3","1")])]
        r=policies(roster,sc)
        self.assertEqual((r["both_scored"]["gain"],r["all_scored"]["gain"],r["labeled_abstentions"],r["missing_outcomes"]),(2,1,1,1))
        with self.assertRaises(ValueError):policies(roster,sc[:-1])
        with self.assertRaises(ValueError):policies(roster,sc[:2]+[dict(sc[2],baseline="0")]+sc[3:])
    def test_direction_missingness_and_dependent_endpoints(self):
        e=[dict(assay=a,cohort="same") for a in ("r","m","p")]
        n=[dict(gene="G",assay=a,direction="-1") for a in ("r","m","p")]
        o=[dict(gene="G",**x,testable=t,effect=v,q=q) for x,(t,v,q) in zip(e,[("1","-3",".01"),("1","2",".01"),("0","","")])]
        r=summarize(n,e,o)["candidates"][0]
        self.assertEqual(r["confirmed_assays"],["r"])
        self.assertEqual(r["confirmed_cohort_labels"],["same"])
        self.assertEqual(len(r["endpoints"]),3)
        with self.assertRaises(ValueError):summarize(n,e,o[:-1])
        with self.assertRaises(ValueError):summarize(n,e,o+[o[0]])
