"""Stage cBioPortal files from a pinned Datahub revision, or verify a local cache."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import time
import urllib.error
import urllib.request


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def get(url, body=None, headers=None):
    fields = {"User-Agent": "KilnBio-public-research"}
    if body:
        fields.update({"Content-Type": "application/vnd.git-lfs+json", "Accept": "application/vnd.git-lfs+json"})
    fields.update(headers or {})
    return urllib.request.urlopen(urllib.request.Request(url, data=body, headers=fields), timeout=120)


def stage(cache, offline=False):
    spec = json.loads(Path(__file__).with_name("sources.json").read_text())
    base = "https://raw.githubusercontent.com/cBioPortal/datahub/{}/public/{}/".format(spec["revision"], spec["study"])
    for name, (digest, size) in spec["files"].items():
        dest = cache / name
        if dest.exists():
            if dest.stat().st_size != size or sha(dest) != digest:
                raise ValueError(f"Checksum failed: {dest}; remove this cache file and retry online.")
            continue
        if offline:
            raise FileNotFoundError(f"Offline cache missing {name}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        for attempt in range(3):
            temp = None
            try:
                with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as out:
                    temp = Path(out.name)
                    with get(base + name) as initial:
                        prefix = initial.read(256)
                        if prefix.startswith(b"version https://git-lfs.github.com/spec/v1"):
                            pointer = (prefix + initial.read()).decode().splitlines()
                            if f"oid sha256:{digest}" not in pointer or f"size {size}" not in pointer:
                                raise ValueError(f"Unexpected upstream LFS pointer: {name}")
                            payload = {"operation": "download", "transfers": ["basic"], "objects": [{"oid": digest, "size": size}]}
                            with get(spec["lfs_batch"], json.dumps(payload).encode()) as reply:
                                obj = json.load(reply)["objects"][0]
                            action = obj["actions"]["download"]
                            with get(action["href"], headers=action.get("header")) as stream:
                                shutil.copyfileobj(stream, out, 1024 * 1024)
                        else:
                            out.write(prefix)
                            shutil.copyfileobj(initial, out, 1024 * 1024)
                if temp.stat().st_size != size or sha(temp) != digest:
                    raise ValueError(f"Downloaded file failed size/SHA-256 check: {name}")
                temp.replace(dest)
                print(f"Verified {name} ({size:,} bytes)", flush=True)
                break
            except (OSError, EOFError):
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
            finally:
                if temp:
                    temp.unlink(missing_ok=True)
    return spec
