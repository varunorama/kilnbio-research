"""Stage the pinned public OmniAge PC reference archive in an isolated installation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from importlib.util import find_spec
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from clocks import PC


URL = "https://zenodo.org/api/records/18287372/files/PCClocks.tar.gz/content"
SIZE = 2934797397
MD5 = "b4840fd9d0b10b86a78c0b402966f9cf"
CHUNK = 16 * 1024 * 1024


def fetch_piece(cache, start):
    end = min(start + CHUNK, SIZE) - 1
    piece = cache / f"pc-{start:010d}.chunk"
    if piece.exists() and piece.stat().st_size == end - start + 1:
        return piece
    for attempt in range(5):
        try:
            req = Request(URL, headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"})
            with urlopen(req, timeout=70) as response, tempfile.NamedTemporaryFile(dir=cache, delete=False) as out:
                tmp = Path(out.name)
                try:
                    if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{SIZE}":
                        raise ValueError("Unexpected partial-content response")
                    shutil.copyfileobj(response, out)
                except BaseException:
                    out.close()
                    tmp.unlink(missing_ok=True)
                    raise
            if tmp.stat().st_size != end - start + 1:
                tmp.unlink()
                raise OSError("Incomplete chunk")
            tmp.replace(piece)
            return piece
        except (HTTPError, URLError, TimeoutError, OSError) as err:
            if attempt == 4:
                raise
            wait = int(err.headers.get("Retry-After", 0)) if isinstance(err, HTTPError) and str(err.headers.get("Retry-After", "")).isdigit() else 0
            time.sleep(max(wait, 2 ** attempt))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cache", required=True, type=Path, help="dedicated resumable PC weights cache")
    cli.add_argument("--archive", type=Path, help="use an already downloaded archive")
    cli.add_argument("--workers", type=int, choices=range(1, 9), default=4)
    a = cli.parse_args()
    a.cache.mkdir(parents=True, exist_ok=True)
    archive = a.archive or a.cache / "PCClocks.tar.gz"
    if not archive.exists():
        if a.archive:
            raise FileNotFoundError(archive)
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            parts = list(pool.map(lambda start: fetch_piece(a.cache, start), range(0, SIZE, CHUNK)))
        with tempfile.TemporaryDirectory(dir=a.cache) as tmp:
            staging = Path(tmp) / archive.name
            with staging.open("wb") as out:
                for part in parts:
                    with part.open("rb") as stream:
                        shutil.copyfileobj(stream, out)
            with staging.open("rb") as stream:
                if staging.stat().st_size != SIZE or hashlib.file_digest(stream, "md5").hexdigest() != MD5:
                    for part in parts:
                        part.unlink(missing_ok=True)
                    raise ValueError("PC archive failed verification. Retry to refetch its chunks.")
            staging.replace(archive)
        for part in parts:
            part.unlink(missing_ok=True)
    else:
        with archive.open("rb") as stream:
            if archive.stat().st_size != SIZE or hashlib.file_digest(stream, "md5").hexdigest() != MD5:
                raise ValueError("PC weights do not match the Zenodo record")
    root = Path(find_spec("omniage").origin).parent / "data"
    target = root / "PCClocks"
    if target.exists():
        raise FileExistsError(f"Use a clean OmniAge installation. Already exists: {target}")
    with tempfile.TemporaryDirectory(dir=root) as tmp, tarfile.open(archive, mode="r|gz") as bundle:
        stage = Path(tmp) / "PCClocks"
        for entry in bundle:
            p = PurePosixPath(entry.name)
            if not entry.isfile() or len(p.parts) not in (2, 3) or p.parts[0] != "PCClocks":
                continue
            if p.name == "PC_Impute_Means.parquet" and len(p.parts) == 2:
                dest = stage / p.name
            elif len(p.parts) == 3 and p.parts[1] in PC and p.name in {
                    "center.parquet", "rotation.parquet", "model.parquet", "intercept.txt", "surrogate_weights.parquet", "final_model.parquet"}:
                dest = stage / p.parts[1] / p.name
            else:
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            with bundle.extractfile(entry) as source, dest.open("xb") as stream:
                shutil.copyfileobj(source, stream)
        for name in PC:
            required = ("surrogate_weights.parquet", "final_model.parquet") if name == "PCGrimAge1" else ("model.parquet", "intercept.txt")
            for filename in ("center.parquet", "rotation.parquet", *required):
                if not (stage / name / filename).is_file():
                    raise ValueError(f"Missing {name}/{filename}")
        if not (stage / "PC_Impute_Means.parquet").is_file():
            raise ValueError("Missing shared probe reference")
        stage.rename(target)
    print("Prepared six PC-clock resources at", target)


if __name__ == "__main__":
    main()
