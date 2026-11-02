"""Offline rejection tests for biological design and assay semantics."""
import csv, importlib.util, json, os, random, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('analysis',ROOT/'downstream/analysis.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)

class DesignTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        self.samples=[{'sample':f'S{i}','biological_id':f'D{i}','condition':'A' if i<3 else 'B','batch':'X' if i%2 else 'Y'} for i in range(6)]
        self.cfg={'samples':'samples.csv','numerator':'B','denominator':'A','assay':'rnaseq','engine':'deseq2','count_source':'raw_counts','counts':'counts.tsv'}
        (self.base/'counts.tsv').write_text('feature_id\t'+'\t'.join(f'S{i}' for i in range(6))+'\nG1\t'+ '\t'.join(['10']*6)+'\n')
    def tearDown(self):self.temp.cleanup()
    def prepare(self):
        with (self.base/'samples.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,self.samples[0].keys());w.writeheader();w.writerows(self.samples)
        (self.base/'config.json').write_text(json.dumps(self.cfg))
        return analysis.prepare('de',self.base/'config.json',self.base/'out')
    def test_valid_design_dry_run_has_no_writes(self):
        cfg,plan,out=self.prepare();self.assertFalse(out.exists());self.assertEqual(len(plan['design_matrix']),6)
    def test_batch_confounding_is_rejected(self):
        for s in self.samples:s['batch']=s['condition']
        self.cfg['covariates']=[{'name':'batch','type':'factor'}]
        with self.assertRaisesRegex(ValueError,'rank deficient'):self.prepare()
    def test_technical_replicates_cannot_inflate_n(self):
        self.samples[1]['biological_id']='D0'
        with self.assertRaisesRegex(ValueError,'technical replicates'):self.prepare()
    def test_repeated_donors_require_block(self):
        for i in range(3):self.samples[i+3]['biological_id']=self.samples[i]['biological_id']
        with self.assertRaisesRegex(ValueError,'Repeated donors'):self.prepare()
        self.cfg['covariates']=[{'name':'biological_id','type':'factor'}];self.prepare()
    def test_tpm_fraction_not_rounded_to_counts(self):
        p=self.base/'counts.tsv';p.write_text(p.read_text().replace('10','10.5',1))
        with self.assertRaisesRegex(ValueError,'integers'):self.prepare()
    def test_rsem_reordered_gene_rows_are_rejected_before_tximport(self):
        header='gene_id\ttranscript_id(s)\teffective_length\texpected_count\tTPM\n'
        genes=['G1\tT1\t404\t12.5\t100\n','G2\tT2\t560\t7.5\t95\n']
        p=self.base/'A.genes.results';q=self.base/'B.genes.results'
        p.write_text(header+''.join(genes));q.write_text(header+''.join(reversed(genes)))
        with self.assertRaisesRegex(ValueError,'order'):analysis.check_rsem([str(p),str(q)])
        q.write_text(p.read_text());analysis.check_rsem([str(p),str(q)])
    def test_unknown_parameter_fails(self):
        self.cfg['fdrr']=.01
        with self.assertRaisesRegex(ValueError,'Unknown configuration'):self.prepare()
    def test_rsem_collects_zero_length_genes_across_libraries(self):
        header = 'gene_id\ttranscript_id(s)\teffective_length\texpected_count\tTPM\n'
        first, second = self.base/'first.results', self.base/'second.results'
        first.write_text(header + 'G1\tT1\t0\t0.93\t0\nG2\tT2\t40\t15\t90\n')
        second.write_text(header + 'G1\tT1\t20\t12\t80\nG2\tT2\t0\t0\t0\n')
        self.assertEqual(analysis.check_rsem([str(first), str(second)]), ['G1', 'G2'])
    @unittest.skipUnless(os.environ.get('KILNBIO_RSCRIPT'), 'requires the KilnBio Bioconductor environment')
    def test_rsem_fit_excludes_zero_offset_even_if_abundant_elsewhere(self):
        rng = random.Random(218)
        names = [s['sample'] for s in self.samples]
        (self.base/'quant.csv').write_text('sample,genes_results\n' + ''.join(f'{s},{s}.results\n' for s in names))
        for i, sample in enumerate(names):
            with (self.base/f'{sample}.results').open('w', newline='') as out:
                writer = csv.writer(out, delimiter='\t')
                writer.writerow(['gene_id','transcript_id(s)','length','effective_length','expected_count','TPM','FPKM'])
                for gene in range(160):
                    count = max(1, round(rng.gammavariate(6, (30 + gene % 40) / 6)))
                    invalid = i == 0 and gene == 0
                    writer.writerow([f'G{gene}',f'T{gene}',500,0 if invalid else 400,
                                     0.9 if invalid else count,0 if invalid else count,0 if invalid else count])
        del self.cfg['counts']
        self.cfg.update(count_source='rsem_tximport',quant_manifest='quant.csv',min_count=10,min_samples=2)
        self.prepare()
        subprocess.run([sys.executable,str(ROOT/'downstream/analysis.py'),'de',
                        str(self.base/'config.json'),str(self.base/'out')],check=True)
        results = self.base/'out/results'
        with (results/'feature_filter.tsv').open() as table:
            flags = {r['feature_id']:r for r in csv.DictReader(table,delimiter='\t')}
        self.assertEqual(flags['G0']['valid_effective_length'], 'FALSE')
        self.assertEqual(flags['G0']['included'], 'FALSE')
        with (results/'all_features.tsv').open() as table:
            fitted = {r['feature_id'] for r in csv.DictReader(table,delimiter='\t')}
        self.assertNotIn('G0', fitted)
        self.assertIn('G1', fitted)
    def test_chromatin_requires_reviewed_normalization(self):
        self.cfg.update(assay='atacseq',feature_universe='label_blind_consensus')
        with self.assertRaisesRegex(ValueError,'normalization'):self.prepare()
    def test_no_independent_replicates_fails(self):
        self.samples=self.samples[:1]+self.samples[3:4]
        with self.assertRaisesRegex(ValueError,'two independent'):self.prepare()
    def test_external_factors_match_samples_and_are_positive(self):
        self.cfg.update(external_size_factors='factors.csv',normalization_reviewed=True)
        p=self.base/'factors.csv'
        p.write_text('sample,size_factor\n'+''.join(f'S{i},1\n' for i in reversed(range(6))))
        self.prepare()
        p.write_text(p.read_text().replace('S0,1','S0,0'))
        with self.assertRaisesRegex(ValueError,'positive'):self.prepare()
        p.write_text(p.read_text().replace('S0,0','MISSING,1'))
        with self.assertRaisesRegex(ValueError,'exactly once'):self.prepare()
    def test_external_factors_do_not_override_tximport(self):
        self.cfg.update(count_source='rsem_tximport',external_size_factors='factors.csv',normalization_reviewed=True)
        with self.assertRaisesRegex(ValueError,'length offsets'):self.prepare()

if __name__=='__main__':unittest.main(verbosity=2)
