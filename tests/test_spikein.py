import importlib.util, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'downstream'))
from spikein import factors

class SpikeInTests(unittest.TestCase):
    def records(self):
        return [dict(sample=s,normalization_group='mark',target_reads=10000,reference_reads=n,
                     target_amount=30000,reference_amount=750,amount_unit='ng_chromatin',reference_material='stock1')
                for s,n in [('control',1000),('treated',4000)]]
    def test_global_loss_and_amount_correction(self):
        r=self.records();out=factors(r,'mark',1000)
        self.assertAlmostEqual(out[1]['size_factor']/out[0]['size_factor'],4)
        self.assertAlmostEqual(out[1]['target_per_reference_adjusted_input']/out[0]['target_per_reference_adjusted_input'],.25)
        # Four-fold extra material, rather than occupancy, explains extra spike reads.
        r[1]['reference_amount']=3000
        out=factors(r,'mark',1000)
        self.assertAlmostEqual(out[1]['size_factor']/out[0]['size_factor'],1)
    def test_inadequate_reference_and_mixed_units_fail(self):
        r=self.records();r[1]['reference_reads']=0
        with self.assertRaisesRegex(ValueError,'sufficient'):factors(r,'mark',1000)
        r=self.records();r[1]['amount_unit']='cells'
        with self.assertRaisesRegex(ValueError,'unit'):factors(r,'mark',1000)
    def test_mixed_stock_fails(self):
        r=self.records();r[1]['reference_material']='different'
        with self.assertRaisesRegex(ValueError,'stock'):factors(r,'mark',1000)

if __name__=='__main__':unittest.main()
