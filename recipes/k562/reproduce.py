"""Run the K562 worked example from public deposited assay files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen


HERE = Path(__file__).resolve().parent


def checksum(path, algorithm):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def acquire(cache, sources):
    cache.mkdir(parents=True, exist_ok=True)
    for row in sources:
        p=cache/row["file"]
        alg="sha256" if "sha256" in row else "md5"
        if not p.exists():
            with TemporaryDirectory(dir=cache) as tmp:
                candidate=Path(tmp)/"download"
                with urlopen(Request(row["url"],headers={"User-Agent":"KilnBio-research-reproduction"}),timeout=240) as resp, candidate.open("wb") as out:
                    shutil.copyfileobj(resp,out,1024*1024)
                if checksum(candidate,alg)!=row[alg]:raise ValueError("Source checksum mismatch for "+p.name)
                candidate.replace(p)
        if checksum(p,alg)!=row[alg]:raise ValueError("Source checksum mismatch for "+p.name)


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cache",type=Path,default=Path("public-sources"))
    cli.add_argument("--out",type=Path,required=True)
    a=cli.parse_args()
    if a.out.exists():cli.error("Choose a new directory for --out")
    src=json.loads((HERE/"public_sources.json").read_text())
    cache=a.cache.resolve()
    out=a.out.resolve()
    acquire(cache,src)
    out.mkdir(parents=True)
    shutil.copy2(HERE/"compare.py",out/"compare.py")
    shutil.copy2(HERE/"sources.json",out/"sources.json")
    config=out/"benchmark";config.mkdir()
    shutil.copy2(HERE/"discovery_files.json",config/"sources.json")
    for row in src:
        if row["file"].startswith("enhancer-"):
            shutil.copy2(cache/row["file"],config/row["file"])
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",PYTHONNOUSERSITE="1")
    cmd=[sys.executable,"-I",str(out/"compare.py"),"--benchmark",str(config)]
    cmd.extend(["--source-dir",str(cache)])
    subprocess.run(cmd,check=True,env=env)


if __name__=="__main__":
    main()
