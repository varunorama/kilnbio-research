"""Complete-input CPU scoring for four pinned human non-DNAm aging models."""
import argparse
import csv
import hashlib
import shutil
import tempfile
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen


# Model revision, SHA256 of its PyAging weight file, input contract, reported endpoint.
MODELS = {
    "pasta": ("d4aba6bbc39be2554c7e9f3cbdba95ff7b21f058","e317da12eadf4e6f9ce008c4f4e09164a72c718d2fa5f498a4af04f9232be596",
              "human-gene-expression","relative_age_shift_years"),
    "ocampoatac2": ("9b44d9708dd2443698e0119dba8e69a6b720009b","f29d303fb6bec7db50ab00808f266f670016cfdfe7f159a92eec52e67a6b3a15",
                    "pbmc-grch38-tn5-counts","pbmc_chronological_age_years"),
    "camilloh3k27ac": ("16d6afb27970dc74069f281221bdec905dfe3da6","35d9df32d6ddd23aade05e0560415b07044c12bf637d2a0e82fb3484f2b75c54",
                       "h3k27ac-gene-enrichment","histone_chronological_age_years"),
    "pac": ("3fdf87ce233f7d0ea8308fd0b95620bd0dca54d6","62374a890ca1d130b8fac59b58fe142dba22d881286b63d3542f347dc2c430a7",
            "plasma-olink-explore3072-npx-and-age","mortality_equivalent_age_years"),
}


def verified(f, name):
    if hashlib.file_digest(f,"sha256").hexdigest() != MODELS[name][1]:
        raise ValueError(f"Weight checksum mismatch for {name}")
    f.seek(0)


def fetch(cache,name):
    cache.mkdir(parents=True,exist_ok=True)
    path = cache/(name+".pt")
    if path.exists():
        with path.open("rb") as f:
            verified(f,name)
        return
    url = f"https://huggingface.co/pyaging/{name}/resolve/{MODELS[name][0]}/{name}.pt"
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        dest = Path(tmp)/path.name
        with urlopen(url,timeout=120) as src, dest.open("wb") as dst:
            shutil.copyfileobj(src,dst)
        with dest.open("rb") as f:
            verified(f,name)
        dest.replace(path)


def load(cache,name):
    import torch
    if version("pyaging") != "0.5.7":
        raise ValueError("Use the pinned pyaging==0.5.7 environment")
    with (cache/(name+".pt")).open("rb") as f:
        verified(f,name)  # Only deserialize the pinned public artifact.
        model = torch.load(f,map_location="cpu",weights_only=False).to(torch.float64).eval()
    if model.metadata["clock_name"] != name or len(set(model.features)) != len(model.features):
        raise ValueError("Unexpected model identity or features")
    return model


def score(model,name,matrix,output):
    import numpy as np
    import torch
    with matrix.open(encoding="utf-8-sig",newline="") as f:
        reader = csv.reader(f,delimiter="\t")
        header = next(reader)
        if not header or header[0] != "sample" or len(set(header)) != len(header) or "" in header:
            raise ValueError("Expected sample followed by unique, nonempty feature names")
        lookup = {key:i for i,key in enumerate(header)}
        missing = set(model.features) - set(header)
        if missing:
            raise ValueError(f"{len(missing)} required features absent; first: {sorted(missing)[:4]}")
        indices = [lookup[key] for key in model.features]
        w = csv.writer(output,delimiter="\t",lineterminator="\n")
        w.writerow(["sample",MODELS[name][3]])
        seen = set()
        for row in reader:
            if len(row) != len(header) or not row[0].strip() or row[0] in seen:
                raise ValueError("Unexpected row length, duplicate or empty sample ID")
            seen.add(row[0])
            x = np.array([float(row[i]) for i in indices], dtype=np.float64)
            if not np.isfinite(x).all():
                raise ValueError("Every required measurement must be finite")
            if name == "ocampoatac2" and ((x < 0).any() or (x != np.floor(x)).any() or x.sum()==0):
                raise ValueError("Supply nonnegative integer Tn5 counts with positive total")
            if name == "pac" and x[model.features.index("age")] < 0:
                raise ValueError("Chronological age must be nonnegative in years")
            if name in {"pasta","camilloh3k27ac"} and (x.min() == x.max() or (x < 0).any()):
                raise ValueError("Expected variable nonnegative expression or enrichment")
            # No missing values. Pasta's cohort fill has no effect, and these four
            # frozen models have no cohort-level postprocessing.
            with torch.inference_mode():
                value = model(torch.from_numpy(x[None,:])).item()
            if not np.isfinite(value):
                raise ValueError(f"Nonfinite prediction for {row[0]}")
            w.writerow([row[0],value])
        if not seen:
            raise ValueError("No sample rows")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command",required=True)
    for cmd in ("fetch","features","score"):
        c = sub.add_parser(cmd)
        c.add_argument("--model",choices=MODELS,required=True)
        c.add_argument("--cache",type=Path,required=True)
        if cmd != "fetch":
            c.add_argument("--out",type=Path,required=True)
        if cmd == "score":
            c.add_argument("--input-kind",required=True,choices=[v[2] for v in MODELS.values()])
            c.add_argument("--matrix",type=Path,required=True)
    a = p.parse_args()
    if a.command == "fetch":
        return fetch(a.cache,a.model)
    if a.command == "score" and a.input_kind != MODELS[a.model][2]:
        p.error("Input kind disagrees with model contract")
    if a.out.exists():
        p.error("Choose a new output file")
    model = load(a.cache,a.model)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.out.parent) as tmp:
        path = Path(tmp)/"result"
        with path.open("w",newline="",encoding="utf-8") as f:
            if a.command == "features":
                f.write("feature\n"+"\n".join(model.features)+"\n")
            else:
                score(model,a.model,a.matrix,f)
        path.replace(a.out)


if __name__ == "__main__":
    main()
