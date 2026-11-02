import importlib.util
import unittest
from pathlib import Path

spec=importlib.util.spec_from_file_location('ont_platform',Path(__file__).resolve().parents[1]/'longread/ont.py')
ont=importlib.util.module_from_spec(spec);spec.loader.exec_module(ont)

class PlatformTests(unittest.TestCase):
    def test_wrong_platform_rejected(self):
        with self.assertRaisesRegex(ValueError,'different sequencing platform'):
            ont.check_platform('@RG\tID:HG002\tPL:PACBIO\n')
        with self.assertRaisesRegex(ValueError,'different sequencing platform'):
            ont.check_platform('@RG\tID:a\tPL:ONT\n@RG\tID:b\tPL:ILLUMINA\n')

    def test_unknown_platform_is_not_confirmation(self):
        self.assertFalse(ont.check_platform('@HD\tSO:coordinate\n')['ONT_declared'])
        self.assertTrue(ont.check_platform('@RG\tID:a\tPL:ONT\n')['ONT_declared'])

if __name__=='__main__':unittest.main()
