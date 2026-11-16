"""Build the uploadable Claude skill from the project skill and its license."""
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "assistants" / "kilnbio-skill.zip"


def build():
    with ZipFile(OUTPUT, "w") as archive:
        for name, source in (("SKILL.md", ROOT / ".claude/skills/kilnbio/SKILL.md"),
                             ("LICENSE", ROOT / "LICENSE")):
            entry = ZipInfo(f"kilnbio/{name}")
            entry.compress_type = ZIP_DEFLATED
            entry.create_system = 3
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, source.read_bytes().replace(b"\r\n", b"\n"))


if __name__ == "__main__":
    build()
