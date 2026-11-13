"""Scientific edge cases; run explicitly in the NumPy/pandas/SciPy environment."""
import argparse
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from downstream import olink_import as olink, multiomic_meta as meta
from downstream import clonal_hematopoiesis as ch, exome_prepare as exome


class Platforms(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='kilnbio_platform_')
        self.base=Path(self.temp.name)
    def tearDown(self): self.temp.cleanup()
    def write(self,name,rows):
        p=self.base/name
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
        return str(p)
    def olink_args(self,**changes):
        rows=[dict(SampleID=s,OlinkID=f,NPX=v,ExtNPX=v+4,PlateID='P1',WellID='A1' if s=='s1' else 'B2',SampleQC='PASS',AssayQC='PASS')
              for s,k in [('s1',0),('s2',1),('outsider',100)] for f,v in [('P1',1.5+k),('P2',3.5+k)]]
        args=dict(export=self.write('npx.tsv',rows),samples=self.write('roster.tsv',[dict(sample='s1',group='0'),dict(sample='s2',group='1')]),out=str(self.base/'out'),renormalize='vendor',lod_file=None,min_detect=None)
        args.update(changes); return argparse.Namespace(**args)
    def test_npx_preserved_no_foreign_samples(self):
        a=self.olink_args(); summary=olink.run(a)
        x=pd.read_csv(Path(a.out)/'npx.tsv',sep='\t',index_col=0)
        np.testing.assert_array_equal(x.values,[[1.5,2.5],[3.5,4.5]])
        self.assertEqual(list(x.columns),['s1','s2']); self.assertEqual(summary['lod_source'],'unavailable')
    def test_missing_lod_cannot_pass_detection_filter(self):
        a=self.olink_args(min_detect=.5); s=olink.run(a); self.assertEqual(s['eligible_assays'],0)
    def test_study_centering_transforms_lod_in_same_units(self):
        a=self.olink_args(renormalize='study')
        a.lod_file=self.write('lod.tsv',[dict(OlinkID='P1',LOD_ExtNPX=6),dict(OlinkID='P2',LOD_ExtNPX=8)])
        olink.run(a); x=pd.read_csv(Path(a.out)/'measurement_qc.tsv',sep='\t')
        np.testing.assert_allclose(x.LOD_out,0)
        self.assertEqual(list(x.above_LOD),[False,False,True,True])
    def test_fail_is_not_measured_and_no_silent_duplicate(self):
        a=self.olink_args(); d=pd.read_csv(a.export,sep='\t'); d.loc[0,'AssayQC']='FAIL'; d.to_csv(a.export,sep='\t',index=False)
        self.assertEqual(olink.run(a)['eligible_assays'],1)
        a.out=str(self.base/'bad'); pd.concat([d,d.iloc[:1]]).to_csv(a.export,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'Duplicate'): olink.run(a)
    def test_legacy_controls_are_not_participants(self):
        a=self.olink_args(); d=pd.read_csv(a.export,sep='\t'); d['Sample_Type']='CONTROL'; d.to_csv(a.export,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'no biological'): olink.run(a)
    def test_distinct_panels_retain_all_plate_identities(self):
        a=self.olink_args(); d=pd.read_csv(a.export,sep='\t'); d.loc[d.OlinkID=='P2','PlateID']='P2'; d.to_csv(a.export,sep='\t',index=False)
        olink.run(a); s=pd.read_csv(Path(a.out)/'samples.tsv',sep='\t')
        self.assertTrue((s.olink_PlateID=='P1|P2').all())
    def test_exact_confounding_distinguished_from_association(self):
        v,full=olink.rank_association(pd.Series(['A','A','B','B']),pd.Series(['1','1','2','2']))
        self.assertAlmostEqual(v,1); self.assertTrue(full)
        _,full=olink.rank_association(pd.Series(['A','A','B','B','B']),pd.Series(['1','2','1','2','2']))
        self.assertFalse(full)
    def test_pc_concordant_discordant_and_missing(self):
        up=np.array([[.001,.002],[.001,.999],[.001,np.nan]])
        p,d=meta.directional_pc(up,1-up,2)
        np.testing.assert_allclose(p[:2],[.004,1]); self.assertTrue(np.isnan(p[2]))
        self.assertEqual(d[0],1)
    def test_partial_null_not_just_global_null(self):
        # One arbitrarily strong layer and two correlated null layers: r=2 null is true.
        rng=np.random.default_rng(20261112); n=100000
        for rho in [0,.8,1]:
            shared=rng.normal(size=n)
            z=np.column_stack([np.full(n,12),shared,rho*shared+np.sqrt(1-rho*rho)*rng.normal(size=n)])
            p,_=meta.directional_pc(norm.sf(z),norm.cdf(z),2)
            self.assertLessEqual(np.mean(p<=.05),.052)
    def test_feature_direction_is_not_chosen_before_bonferroni(self):
        self.write('results.tsv',[dict(feature='a',effect=1,p=.01),dict(feature='b',effect=-1,p=.02)])
        self.write('map.tsv',[dict(feature='a',gene='G'),dict(feature='b',gene='G')])
        t,_=meta.load_layer(dict(results='results.tsv',feature_map='map.tsv',orientation='1',layer='a'),self.base)
        self.assertAlmostEqual(t.loc['G','up'],.01); self.assertAlmostEqual(t.loc['G','down'],.02)
    def test_by_matches_closed_form(self):
        np.testing.assert_allclose(meta.adjust([.01,.04,.9]),[.055,.11,1])
    def test_tiny_acat_not_zero(self):
        p=meta.acat(np.array([[1e-100,.4]]),np.array([1,1]))[0]
        self.assertGreater(p,0); self.assertAlmostEqual(p/1e-100,2)
    def ch_args(self):
        return argparse.Namespace(min_depth=20,min_alt_reads=3,min_vaf=.02,max_population_af=.001,assembly='GRCh38',population_af_field='gnomAD_AF')
    def test_ch_quality_boundaries_and_missing_frequency(self):
        a=self.ch_args(); row=dict(filter='PASS',dp='150',alt_count='3',alt_forward='1',alt_reverse='2',vaf='.02',population_af='.001')
        self.assertEqual(ch.quality_reasons(row,a)[0],[])
        row['population_af']=''; self.assertIn('population_frequency_unknown',ch.quality_reasons(row,a)[0])
        row['alt_count']='151'
        with self.assertRaisesRegex(ValueError,'exceed'): ch.quality_reasons(row,a)
    def test_ch_absence_is_not_negative(self):
        a=self.ch_args(); a.drivers=self.write('drivers.tsv',[dict(assembly='GRCh38',chrom='1',pos=10,ref='A',alt='G',gene='FAKE',evidence='synthetic test only')])
        a.variants=self.write('variants.tsv',[dict(sample='s1',assembly='GRCh38',chrom='1',pos=10,ref='A',alt='G',filter='PASS',dp=150,alt_count=3,alt_forward=1,alt_reverse=2,vaf=.02,population_af=0,orientation_kind='read_strand')])
        _,calls=ch.screen(a,[dict(sample=s,donor=s,specimen=s) for s in ['s1','s2']])
        self.assertEqual(calls[0]['ch_status'],'candidate_driver_detected')
        self.assertEqual(calls[1]['ch_status'],'no_qualifying_variant_in_supplied_records')
    def test_vcf_orientation_and_multiallelic_rejection(self):
        p=self.base/'x.vcf'; header='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
        p.write_text(header+'1\t10\t.\tA\tG\t50\tPASS\tgnomAD_AF=0\tDP:AD:F1R2:F2R1\t100:97,3:50,1:47,2\n')
        row=list(ch.vcf_rows(p,'gnomAD_AF','GRCh38'))[0]
        self.assertEqual(row['orientation_kind'],'pair_orientation'); self.assertEqual(ch.quality_reasons(row,self.ch_args())[1],.03)
        p.write_text(p.read_text().replace('##fileformat=VCFv4.2','##fileformat=VCFv4.2\n##reference=hg19'))
        with self.assertRaisesRegex(ValueError,'contradicts'): list(ch.vcf_rows(p,'gnomAD_AF','GRCh38'))
        p.write_text(p.read_text().replace('##reference=hg19\n',''))
        p.write_text(p.read_text().replace('\tG\t','\tG,T\t'))
        with self.assertRaisesRegex(ValueError,'multiallelic'): list(ch.vcf_rows(p,'gnomAD_AF','GRCh38'))
    def test_exome_rejects_assembly_mismatch_and_missing_resources(self):
        (self.base/'ref.fai').write_text('chr1\t100\t0\t50\t51\n')
        (self.base/'targets.bed').write_text('chr1\t0\t50\n')
        p=dict(wes=True,aligner='bwa-mem2',intervals=str(self.base/'targets.bed'),fasta_fai=str(self.base/'ref.fai'),tools='haplotypecaller',skip_tools='baserecalibrator')
        exome.validate_params('exome',p); self.assertTrue(p['normalize_vcfs'])
        (self.base/'targets.bed').write_text('1\t0\t50\n')
        with self.assertRaisesRegex(ValueError,'contig'):exome.validate_params('exome',p)
        (self.base/'targets.bed').write_text('chr1\t0\t50\n');p['tools']='mutect2'
        with self.assertRaisesRegex(ValueError,'panel of normals'):exome.validate_params('ch-somatic',p)
    def test_provider_import_preserves_negative_and_unknown(self):
        calls=[dict(Biobank_Subject_ID='private1',hasCHIP='1',nCHIP='1',Gene='SIM',VAF='.03',Batch_CHIP_call='batch1'),dict(Biobank_Subject_ID='private2',hasCHIP='0',nCHIP='0',Gene='',VAF='',Batch_CHIP_call='batch1')]
        variants=[dict(Biobank_Subject_ID='private1',varID='synthetic',Gene='SIM',VAF='.03',DP='100',AD_Alt='3',ADF_Alt='1',ADR_Alt='2')]
        args=argparse.Namespace(calls=self.write('calls.tsv',calls),variants=self.write('variants.tsv',variants))
        roster=[dict(sample=f's{i}',donor=f'd{i}',specimen=f'b{i}',source_subject=f'private{i}') for i in [1,2,3]]
        _,result=ch.import_calls(args,roster)
        self.assertEqual([r['ch_status'] for r in result],['provider_reported_positive','provider_reported_negative','unknown'])
        self.assertFalse(any('private' in str(r) for r in result))
        calls[0]['nCHIP']='2';self.write('calls.tsv',calls)
        with self.assertRaisesRegex(ValueError,'disagree'): ch.import_calls(args,roster)

if __name__=='__main__': unittest.main(verbosity=2)
