"""Regression checks for profile selection and failed-install activation."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("kilnbio_installer", ROOT / "environments/install.py")
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class InstallerTests(unittest.TestCase):
    def test_minimal_default(self):
        self.assertEqual(installer.profiles_requested(None), ["core"])

    def test_singlecell_includes_pseudobulk_statistics(self):
        self.assertEqual(installer.profiles_requested(["singlecell"]),
                         ["analysis", "bladder", "core", "singlecell"])

    def test_all_includes_each_managed_profile(self):
        self.assertEqual(installer.profiles_requested(["all"]), sorted(installer.PROFILES))

    def test_failed_smoke_preserves_existing_command(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory)
            command = prefix / "bin/kilnbio"
            command.parent.mkdir()
            command.write_text("previous working installation")
            tools = prefix / "tools" / installer.MAMBA_VERSION
            tools.mkdir(parents=True)
            mamba = tools / "micromamba"
            mamba.write_text("fixture")
            (prefix / "environments").mkdir()
            digest = installer.digest
            with patch.object(installer.platform, "system", return_value="Linux"), \
                 patch.object(installer.platform, "machine", return_value="x86_64"), \
                 patch.object(installer.shutil, "which", return_value="fixture"), \
                 patch.object(installer, "digest", side_effect=lambda p: installer.MAMBA_SHA256 if p == mamba else digest(p)), \
                 patch.object(installer, "execute"), \
                 patch.object(installer, "smoke", side_effect=RuntimeError("bad result")):
                with self.assertRaisesRegex(RuntimeError, "bad result"):
                    installer.install(SimpleNamespace(prefix=prefix, container="docker"), ["core"])
            self.assertEqual(command.read_text(), "previous working installation")
            self.assertFalse((prefix / "installation.json").exists())
            self.assertFalse((prefix / ".setup-lock").exists())


if __name__ == "__main__":
    unittest.main()
