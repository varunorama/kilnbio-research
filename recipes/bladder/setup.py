"""Install only two pinned R packages into a new, isolated library."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request

CRAN = "https://cran.r-project.org"
BIOC = ["https://bioconductor.posit.co/packages/3.23/bioc","https://bioconductor.org/packages/3.23/bioc"]
SOURCES = {
    "statmod_1.5.2.tar.gz": ("5f83d82043c5e26f1048fc0139580da1b7e40becf95fb4cb26ed3afcd46a581d",
                             [CRAN+"/src/contrib", CRAN+"/src/contrib/Archive/statmod"]),
    "limma_3.68.5.tar.gz": ("40b1d2059a76ad98db4bf80139760acf8e8bae92a8265aaaf16ac199831533f4",
                            [b+s for s in ["/src/contrib","/src/contrib/Archive/limma"] for b in BIOC]),
    "statmod_1.5.2.zip": ("29ba741681c6f9b1c011de09e84ecebb75711a85344ef67fd7ad0a6393c9d891",
                         [CRAN+"/bin/windows/contrib/4.6"]),
    "limma_3.68.5.zip": ("8feaaf7819b8d9fc88375f8dc07bff5d7160249260644e9f283a1c22f345f4e9",
                        [b+"/bin/windows/contrib/4.6" for b in BIOC]),
}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rscript",default="Rscript")
    p.add_argument("--library",required=True,type=Path)
    args=p.parse_args()
    lib=args.library.resolve()
    if lib.exists():
        p.error("Choose a new library directory; existing R packages will not be overwritten")
    script=str(Path(__file__).with_name("setup.R"))
    subprocess.run([args.rscript,"--vanilla",script,"check"],check=True)
    win=sys.platform=="win32"
    with tempfile.TemporaryDirectory() as tmp:
        files=[]
        for name,(digest,bases) in SOURCES.items():
            if name.endswith(".zip")!=win:
                continue
            for base in bases:
                try:
                    with urllib.request.urlopen(base+"/"+name,timeout=90) as f:
                        data=f.read()
                    if hashlib.sha256(data).hexdigest()!=digest:
                        raise ValueError("Package checksum differs: "+name)
                    break
                except OSError as err:
                    if base==bases[-1]:
                        raise RuntimeError("Could not download pinned "+name+"; see README.md") from err
            path=Path(tmp)/name
            path.write_bytes(data)
            files.append(str(path))
        lib.mkdir(parents=True)
        subprocess.run([args.rscript,"--vanilla",script,str(lib),"win.binary" if win else "source",*files],check=True)


if __name__=="__main__":
    main()
