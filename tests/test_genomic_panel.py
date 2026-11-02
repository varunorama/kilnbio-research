import sys,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"downstream"))
from genomic_panel import extract

class GenomicPanel(unittest.TestCase):
    def test_absent_not_reference_or_low_quality(self):
        samples=[dict(sample="S",donor="P",specimen="blood")]
        panel=[dict(event="v"+str(pos),assembly="GRCh38",chrom="chr1",pos=str(pos),ref="A",alt="G") for pos in range(1,6)]
        header="##fileformat=VCFv4.3\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        v=["chr1\t%d\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t%s"%(i,g) for i,g in [(1,"0/0:20:40"),(2,"1|0:18:25"),(3,"./1:30:30"),(4,"1/1:6:60")]]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"tiny.vcf";p.write_text(header+"\n".join(v)+"\n")
            r=extract(p,panel,samples,"GRCh38",10,20)
            self.assertEqual([x[4] for x in r],["reference","altered","unknown","unknown","unknown"])
            self.assertEqual(r[1][5],1)
            self.assertEqual(r[-1][6],"no_exact_record")
            p.write_text(header.replace("##reference=GRCh38","##reference=hg19")+"\n".join(v)+"\n")
            with self.assertRaisesRegex(ValueError,"contradicts"):extract(p,panel,samples,"GRCh38",10,20)
