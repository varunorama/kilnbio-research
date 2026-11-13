"""Check allele attribution, missingness, mapping and explicit execution boundaries."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'downstream'))
from variant_annotations import annotate, queries
from alphafold_lookup import lookup, map_protein
from vep_annotate import run as run_vep


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name)
        (self.p/'q.tsv').write_text('variant_id\tchromosome\tposition\tref\talt\na\tchr1\t100\tA\tC\nb\t1\t101\tT\tG\n')
        self.c={'assembly':'GRCh38','normalized_biallelic':True,'variants':'q.tsv','resources':{}}
    def tearDown(self):self.temp.cleanup()
    def source(self,kind,name,text):
        (self.p/name).write_text(text)
        self.c['resources'][kind]={'path':name,'version':'fixture-v1','source_url':'https://example.org/fixture','assembly':'GRCh38'}
    def run_annotation(self):
        (self.p/'c.json').write_text(json.dumps(self.c))
        return annotate(self.p/'c.json',self.p/'out')
    def test_clinvar_retains_conflict_and_separate_somatic_assertion(self):
        self.source('clinvar','cv.vcf','##reference=GRCh38\n1\t100\tCV1\tA\tC\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Benign(1)|Pathogenic(1);CLNREVSTAT=criteria_provided,_conflicting_classifications;ONC=Oncogenic;SCI=Tier_II\n')
        r=self.run_annotation();self.assertEqual(r['resources']['clinvar']['matched_variants'],1)
        data=json.loads((self.p/'out/clinvar.jsonl').read_text())['annotation']['info']
        self.assertEqual(data['ONC'],'Oncogenic');self.assertIn('Pathogenic(1)',data['CLNSIGCONF'])
        with (self.p/'out/variants.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(rows[1]['clinvar_status'],'not_found_in_resource')
    def test_assembly_contradiction_fails_with_receipt(self):
        self.source('clinvar','cv.vcf','##reference=GRCh37\n')
        with self.assertRaisesRegex(ValueError,'assembly'):self.run_annotation()
        self.assertEqual(json.loads((self.p/'out/receipt.json').read_text())['status'],'failed')
    def test_matched_multiallelic_record_not_misattributed(self):
        self.source('clinvar','cv.vcf','1\t100\t.\tA\tC,G\t.\t.\tCLNSIG=Benign,Pathogenic\n')
        with self.assertRaisesRegex(ValueError,'multiallelic'):self.run_annotation()
    def test_missense_all_transcripts_retained(self):
        self.source('alphamissense','am.tsv','#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n1\t100\tA\tC\thg38\tP04637\tT1\tR175H\t0.8\tlikely_pathogenic\n1\t100\tA\tC\thg38\tP04637\tT2\tR175H\t0.8\tlikely_pathogenic\n')
        r=self.run_annotation();self.assertEqual(r['resources']['alphamissense']['annotation_records'],2)
        self.assertEqual(r['resources']['alphamissense']['matched_variants'],1)
    def test_invalid_missense_score_fails(self):
        self.source('alphamissense','am.tsv','#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n1\t100\tA\tC\thg38\tP04637\tT1\tR175H\tNaN\tunknown\n')
        with self.assertRaisesRegex(ValueError,'score'):self.run_annotation()
    def test_vep_preserves_transcript_consequences(self):
        self.source('vep','v.vcf','##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence|Feature">\n1\t100\t.\tA\tC\t.\tPASS\tCSQ=C|missense_variant|T1,C|intron_variant|T2\n')
        self.run_annotation()
        x=json.loads((self.p/'out/vep.jsonl').read_text())['annotation']['transcript_consequences']
        self.assertEqual([r['Feature'] for r in x],['T1','T2'])
    def test_duplicate_normalized_alleles_and_unsafe_alleles_rejected(self):
        for text in ['b\t1\t100\tA\tC','b\t1\t101\tT\t<DEL>','b\t1\t0\tT\tG']:
            (self.p/'q.tsv').write_text('variant_id\tchromosome\tposition\tref\talt\na\tchr1\t100\tA\tC\n'+text+'\n')
            with self.assertRaises(ValueError):queries(self.p/'q.tsv')
    def test_no_silent_liftover(self):
        self.c['normalized_biallelic']=False
        with self.assertRaisesRegex(ValueError,'normalization'):self.run_annotation()
    def test_fragment_mapping_and_mismatch(self):
        model={'sequence':'ACD','sequenceStart':10,'sequenceEnd':12}
        self.assertEqual(map_protein(model,'C11H'),('reference_residue_matches',2))
        self.assertEqual(map_protein(model,'R11H')[0],'reference_residue_mismatch')
        self.assertEqual(map_protein(model,'A1C')[0],'outside_model_fragment')
    def test_alphafold_offline_plan_never_contacts_api(self):
        q=self.p/'proteins.tsv';q.write_text('uniprot_id\tprotein_variant\nP04637\tR175H\n')
        with patch('urllib.request.urlopen',side_effect=AssertionError('network')):
            self.assertEqual(lookup(q,self.p/'out')['status'],'planned')
    def test_alphafold_cache_current_schema_retains_model_version(self):
        q=self.p/'proteins.tsv';q.write_text('uniprot_id\tprotein_variant\nP04637\tR2H\n')
        (self.p/'P04637.json').write_text(json.dumps([{'uniprotAccession':'P04637','modelEntityId':'AF-P04637-F1','latestVersion':6,'sequence':'MRA','sequenceStart':1,'sequenceEnd':3}]))
        r=lookup(q,self.p/'out',cache=self.p)
        self.assertEqual(r['output_rows'],1)
        with (self.p/'out/structures.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(row['model_version'],'6');self.assertEqual(row['mapping_status'],'reference_residue_matches')
    def test_alphafold_wrong_accession_not_silently_mapped(self):
        q=self.p/'proteins.tsv';q.write_text('uniprot_id\tprotein_variant\nP04637\tR2H\n')
        (self.p/'P04637.json').write_text(json.dumps([{'uniprotAccession':'P04637-2'}]))
        lookup(q,self.p/'out',cache=self.p)
        self.assertIn('no_matching_model',(self.p/'out/structures.tsv').read_text())
    def vep_config(self):
        (self.p/'genome.fa').write_text('>1\nACGT\n')
        (self.p/'v.vcf').write_text('##fileformat=VCFv4.2\n')
        (self.p/'cache/homo_sapiens/116_GRCh38').mkdir(parents=True)
        (self.p/'cache/homo_sapiens/116_GRCh38/info.txt').write_text('fixture only\n')
        c={'input_vcf':'v.vcf','fasta':'genome.fa','cache_dir':'cache','assembly':'GRCh38','forks':1}
        p=self.p/'vep.json';p.write_text(json.dumps(c));return p
    def test_vep_plan_is_offline_and_keeps_all_transcripts(self):
        with patch('subprocess.run',side_effect=AssertionError('execute')):
            p=run_vep(self.vep_config(),self.p/'out')
        self.assertIn('--check_ref',p['argv']);self.assertNotIn('--pick',p['argv']);self.assertIn('--offline',p['argv'])
    def test_vep_rejects_wrong_installed_version(self):
        from subprocess import CompletedProcess
        with patch('subprocess.run',return_value=CompletedProcess([],0,'ensembl-vep 115.2','')):
            with self.assertRaisesRegex(ValueError,'116.3'):run_vep(self.vep_config(),self.p/'out',True)
        self.assertEqual(json.loads((self.p/'out/receipt.json').read_text())['status'],'failed')


if __name__=='__main__':unittest.main(verbosity=2)
