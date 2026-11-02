import tempfile
import unittest
import math
from pathlib import Path
from score import interval_grid, tag_grid, score_experiment


class GridTest(unittest.TestCase):
    def test_midpoints_half_open_and_mask_overlaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed=Path(tmp)/'regions.bed'
            bed.write_text('chr1\t50\t150\nchr1\t275\t279\n')
            self.assertEqual(interval_grid(bed,{'chr1':300},100)[0]['chr1'].tolist(),[True,False,False])
            self.assertEqual(interval_grid(bed,{'chr1':300},100,mask=True)[0]['chr1'].tolist(),[True,True,True])

    def test_chromatin_fragment_centers_and_stranded_chip(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed=Path(tmp)/'tags.bed'
            bed.write_text('chr1\t0\t20\t.\t60\t+\nchr1\t180\t200\t.\t60\t-\nchr1\t280\t299\t.\t60\t+\n')
            self.assertEqual(tag_grid(bed,{'chr1':300},100,True)['chr1'].tolist(),[1,1,1])
            self.assertEqual(tag_grid(bed,{'chr1':300},100,False)['chr1'].tolist(),[1,1,0])

    def test_heldout_enrichment_and_equal_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            ref = base/'ref'
            ref.mkdir()
            (ref/'panel.sizes').write_text('chr1\t10000\n')
            (ref/'mask.bed').touch()
            exp = base/'toy'
            exp.mkdir()

            def tags(bins):
                return ''.join(f'chr1\t{1000*b}\t{1000*b+20}\t.\t60\t+\n' for b in bins)

            for rep in (1, 2):
                (exp/f'rep{rep}.bed').write_text(tags([0, 2, 2, 2, 8]))
                for caller, filename, interval in [('MACS3', 'calls_peaks.broadPeak', (0,3000)),
                                                   ('Omnipeak', 'calls.broadPeak', (2000,4000))]:
                    folder = exp/f'r{rep}-{caller}'
                    folder.mkdir()
                    (folder/filename).write_text(f'chr1\t{interval[0]}\t{interval[1]}\n')
            (exp/'control.test.bed').write_text(tags(range(10)))
            (exp/'control.train.bed').write_text(tags([2,8,8,8,8]))
            rows = score_experiment(base, dict(assay='toyMark', experiment='toy', cell='toyCell'))
            panel = {(r['caller'], r['view']): r for r in rows if r['rep']==1 and r['region']=='panel'}
            # MACS3 recovers 4/5 treatment and 3/10 withheld control tags.
            self.assertAlmostEqual(panel['MACS3','native']['heldout_enrichment'], 8/3)
            self.assertAlmostEqual(panel['Omnipeak','native']['heldout_enrichment'], 3)
            self.assertEqual(panel['MACS3','native']['replicate_jaccard'], 1)
            for caller in ('MACS3', 'Omnipeak'):
                self.assertEqual(panel[caller,'equal_bp']['selected_bins'], 2)
            (exp/'r1-MACS3'/'calls_peaks.broadPeak').write_text('chr1\t4\t75\n')
            sparse = score_experiment(base, dict(assay='toyMark', experiment='toy', cell='toyCell'))
            for entry in sparse:
                if entry['rep']==1 and entry['view']=='equal_bp':
                    self.assertEqual(entry['selected_bins'], 0)
                    self.assertTrue(math.isnan(entry['heldout_enrichment']))


if __name__=='__main__':
    unittest.main()
