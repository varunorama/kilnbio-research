import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "downstream" / "capture_panel.py"


class CapturePanel(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "sizes").write_text("chr1\t20\n")
        (self.root / "target").write_text("chr1\t1\t5\nchr1\t3\t6\nchr1\t6\t7\n")
        self.args = [sys.executable, str(SCRIPT), "--targets", "target", "--name", "test",
                     "--sizes", "sizes", "--assembly", "GRCh38", "--data-assembly", "GRCh38", "--out", "out"]

    def invoke(self, *args, success=True):
        result = subprocess.run(self.args + list(args), cwd=self.root, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            return json.loads((self.root / "out/summary.json").read_text())
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / "out").exists())
        return result.stderr

    def test_union_boundaries_cpg_and_baits(self):
        (self.root / "baits").write_text("chr1\t0\t8\n")
        (self.root / "cpg").write_text("chr\tpos\tN\tX\nchr1\t1\t2\t1\nchr1\t6\t2\t0\nchr1\t7\t1\t1\n")
        info = self.invoke("--measurements", "cpg", "--format", "cpg", "--baits", "baits")
        self.assertEqual((info["target_rows"],info["target_union_bases"],info["baits"]["union_bases"]), (3,6,8))
        self.assertEqual((info["features"]["fully_on_target"],info["features"]["tested"]), (1,3))

    def test_vcf_full_span_pass_only(self):
        (self.root / "calls").write_text("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\t.\tA\tT\t.\tPASS\t.\nchr1\t7\t.\tAT\tA\t.\tPASS\t.\n"
            "chr1\t3\t.\tA\t<DEL>\t.\tPASS\t.\nchr1\t3\t.\tA\tC\t.\tLowQual\t.\n")
        info = self.invoke("--measurements", "calls", "--format", "vcf")
        self.assertEqual((info["features"]["fully_on_target"],info["features"]["tested"]), (1,2))

    def test_wrong_build_duplicate_cpg_and_unknown_contig(self):
        (self.root / "x").write_text("##reference=hg19\n")
        self.assertIn("contradicts", self.invoke("--measurements","x","--format","vcf",success=False))
        (self.root / "x").write_text("chr\tpos\tN\tX\nchr1\t2\t2\t1\nchr1\t2\t2\t1\n")
        self.assertIn("unique", self.invoke("--measurements","x","--format","cpg",success=False))
        (self.root / "x").write_text("1\t0\t1\n")
        self.assertIn("Unknown contig", self.invoke("--measurements","x",success=False))
        self.assertIn("original vendor",self.invoke("--preset","twist-exome-2-hg38",success=False))

    @unittest.skipUnless(shutil.which("samtools"), "requires samtools")
    def test_zero_coverage_in_denominator(self):
        sam = self.root / "tiny.sam"
        sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\n"
                       "read1\t0\tchr1\t2\t60\t3M\t*\t0\t0\tAGC\tIII\n")
        subprocess.run(["samtools","view","-b","-o",str(self.root/"tiny.bam"),str(sam)],check=True)
        info = self.invoke("--bam","tiny.bam","--min-depth","1")["alignment_depth"]
        self.assertEqual((info["target_bases"],info["covered_bases"],info["zero_depth_bases"]), (6,3,3))
        self.assertEqual(info["fraction_at_min_depth"],0.5)


if __name__ == "__main__":
    unittest.main()
