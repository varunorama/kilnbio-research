"""Inventory installed clock engines and run the supported local DNAm projections."""
import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.request import urlopen


GADD = "https://raw.githubusercontent.com/DanniGadd/EpiScores-for-protein-levels/a5130fab3895a0d95f0dcc8826aa9fb5e8c0fa86/Weights_files/"
T2D = "https://raw.githubusercontent.com/marioni-group/episcores-diabetes-prediction/22eebea219b23cc0d80a1333326b11b08fdbe91c/kora_validation/"
WEIGHTS = {
    "KORA_weights_finalised_120121.csv": (GADD, "d41753e26d0618e60cbf38264ea8991c405e8357b163bb7e3c37d471aa0927b5"),
    "LBC_weights_finalised_120121.csv": (GADD, "782e5165cd2faa339aa53bee2e17ab596ee9314ecdf5f7c479437bb084cc5f4f"),
    "directEpiScoreCoefficients.rds": (T2D, "c8de0320e724668a45b2e6da0bdba1132e19462ed8edef9321769584f6d146f0"),
    "dModelCoefficients.rds": (T2D, "9d48be5ff9512da3eff9ef93bdfdb7c3164000b12fab10d4fd6d108098f376ae"),
}
LOCAL = ("Hannum", "Horvathv1", "PhenoAge", "DunedinPACE")
PC = ("PCHorvath2013", "PCHorvath2018", "PCHannum", "PCPhenoAge", "PCDNAmTL", "PCGrimAge1")
OMNI = ("Hannum", *PC)


def local_linear_models():
    # Load Biolearn's Torch runtime before the other native libraries on Windows.
    import torch
    from biolearn.model import model_definitions
    from biolearn.util import get_data_file
    names = set()
    for name, definition in model_definitions.items():
        spec = definition["model"]
        filename = spec.get("file", "")
        if (spec["type"] == "LinearMethylationModel" and filename
                and Path(filename).name == filename and Path(get_data_file(filename)).is_file()):
            names.add(name)
    return names


def resources(cache, names, download=False):
    cache.mkdir(parents=True, exist_ok=True)
    for name in names:
        base, digest = WEIGHTS[name]
        dest = cache / name
        if not dest.is_file() and download:
            with tempfile.TemporaryDirectory(dir=cache) as staging:
                staged = Path(staging) / name
                with urlopen(base + name, timeout=60) as response, staged.open("wb") as stream:
                    shutil.copyfileobj(response, stream)
                if hashlib.sha256(staged.read_bytes()).hexdigest() != digest:
                    raise ValueError(f"Unexpected upstream digest for {name}")
                staged.replace(dest)
        if not dest.exists() or hashlib.sha256(dest.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Missing or unverified {name}. See fetch-marioni.")


def inventory():
    available = local_linear_models()
    from biolearn.model import model_definitions
    import omniage
    return {
        "biolearn_version": version("biolearn"),
        "biolearn": {name: {"species": m.get("species"), "tissue": m.get("tissue"),
                           "output": m.get("output"), "paper": m.get("source"),
                           "type": m["model"]["type"], "kilnbio_local_adapter": name in available}
                     for name, m in model_definitions.items()},
        "omniage_version": version("omniage"),
        "omniage_installed_groups": omniage.list_available_clocks(check_installed=True, verbose=False),
        "omniage_paper": "https://doi.org/10.1038/s41467-026-76038-w",
        "marioni_protein_panel": "109 DNAm proxies, cohort-standardized",
        "EpiDiff": "DNAm differentiation index, 0 to 1, preprint; R package 0.1.0",
        "marioni_diabetes": "145-probe score, optional clinical-plus-DNAm ten-year prediction",
    }


def read_betas(path):
    import pandas as pd
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream, delimiter="\t"))
    if not header or header[0] != "probe" or len(header) < 2 or len(header) != len(set(header)) or "" in header:
        raise ValueError("Expected unique sample names following the probe column")
    x = pd.read_csv(path, sep="\t", index_col=0)
    if not x.index.is_unique or x.index.hasnans:
        raise ValueError("Expected unique probe IDs")
    x.columns = header[1:]
    return x


def complete(x, needed, name):
    import numpy as np
    import pandas as pd
    missing = set(needed).difference(x.index)
    if missing:
        raise ValueError(f"{name} missing {len(missing)} of {len(needed)} required scoring/reference probes")
    x = x.loc[sorted(needed)].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(x.to_numpy()).all() or x.min().min() < 0 or x.max().max() > 1:
        raise ValueError("Supply measured, finite beta values in [0,1] at all required probes")
    return x


def score_biolearn(name, path):
    if name not in local_linear_models():
        raise ValueError("Expected a catalogued local linear methylation model")
    return biolearn_from_matrix(name, read_betas(path))


def biolearn_from_matrix(name, x):
    import numpy as np
    import pandas as pd
    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery
    from biolearn.util import get_data_file
    if name not in local_linear_models():
        raise ValueError("Expected a catalogued local linear methylation model")
    model = ModelGallery().get(name, imputation_method="none")
    needed = set(model.methylation_sites())
    if name == "DunedinPACE":
        reference = pd.read_csv(get_data_file("DunedinPACE_Gold_Means.csv"), index_col=0)
        needed.update(reference.index)
    x = complete(x, needed, name)
    scores = model.predict(GeoData(pd.DataFrame(index=x.columns), dnam=x))
    if not scores.index.equals(x.columns) or scores.shape[1] != 1 or not np.isfinite(scores.to_numpy()).all():
        raise ValueError("Unexpected upstream result or sample alignment")
    scores.columns = [name]
    return scores


def score_omni(name, path, clinical):
    import numpy as np
    import pandas as pd
    import omniage.models as models
    if name in PC:
        assets = Path(models.__file__).resolve().parent.parent / "data" / "PCClocks" / name
        if not (assets / "rotation.parquet").is_file():
            raise ValueError(f"Provision the public OmniAge PC loadings first (see README). Missing {assets}")
    model = getattr(models, name)()
    x = read_betas(path)
    x = complete(x, model.rotation.index if name in PC else model.weights.index, name)
    kwargs = {}
    if name == "PCGrimAge1":
        if clinical is None:
            raise ValueError("PCGrimAge1 requires --clinical with sample, age and sex")
        meta = pd.read_csv(clinical, sep="\t", dtype={"sample": str}).set_index("sample", verify_integrity=True)
        if set(meta.index) != set(x.columns):
            raise ValueError("Clinical samples must match the beta matrix")
        meta = meta.loc[x.columns]
        if not meta.sex.isin(["Female", "Male"]).all() or not np.isfinite(meta.age).all() or (meta.age < 0).any():
            raise ValueError("Expected nonnegative ages and upstream Female/Male sex coding")
        kwargs = {"ages": meta.age, "sex": meta.sex}
    scores = model.predict(x, verbose=False, **kwargs)
    if isinstance(scores, pd.Series):
        scores = scores.to_frame(name=name)
    if not scores.index.equals(x.columns) or scores.empty or not np.isfinite(scores.to_numpy()).all():
        raise ValueError("Unexpected upstream score or sample alignment")
    return scores


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    sub = cli.add_subparsers(dest="command", required=True)
    sub.add_parser("catalog").add_argument("--out", type=Path, required=True)
    sub.add_parser("fetch-marioni").add_argument("--cache", type=Path, required=True)
    score = sub.add_parser("score")
    score.add_argument("--model", required=True, help="A local Biolearn model, a documented omni- selector, or a Marioni panel")
    score.add_argument("--beta", type=Path, required=True)
    score.add_argument("--out", type=Path, required=True)
    score.add_argument("--cache", type=Path)
    score.add_argument("--clinical", type=Path, help="clinical covariates for diabetes risk or PCGrimAge1")
    score.add_argument("--rscript", default="Rscript")
    score.add_argument("--epidiff-upstream-mean", action="store_true",
                       help="Allow the published EpiDiff missing-row fill up to the example's 85 absent CpGs")
    a = cli.parse_args()
    if a.command == "fetch-marioni":
        resources(a.cache, WEIGHTS, download=True)
        return
    if a.out.exists():
        cli.error("Choose a new output file")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.out.parent) as tmp:
        staged = Path(tmp) / "result"
        if a.command == "catalog":
            staged.write_text(json.dumps(inventory(), indent=2, default=str) + "\n", encoding="utf-8")
        elif a.epidiff_upstream_mean and a.model != "EpiDiff":
            cli.error("--epidiff-upstream-mean requires --model EpiDiff")
        elif a.clinical and a.model not in ("marioni-diabetes", "omni-PCGrimAge1"):
            cli.error("--clinical applies only to marioni-diabetes or omni-PCGrimAge1")
        elif a.model == "EpiDiff":
            policy = "upstream-mean" if a.epidiff_upstream_mean else "complete"
            subprocess.run([a.rscript, "--vanilla", Path(__file__).with_name("epidiff.R"),
                            a.beta, staged, policy], check=True)
        elif a.model.startswith("omni-"):
            if a.model[5:] not in OMNI:
                cli.error("Unsupported OmniAge selector. See aging/README.md")
            score_omni(a.model[5:], a.beta, a.clinical).to_csv(staged, sep="\t", index_label="sample")
        elif a.model.startswith("marioni-"):
            if a.model not in ("marioni-proteins", "marioni-diabetes"):
                cli.error("Expected marioni-proteins or marioni-diabetes")
            if a.cache is None:
                cli.error("Supply the verified --cache from fetch-marioni")
            protein = a.model == "marioni-proteins"
            files = [n for n in WEIGHTS if n.endswith(".csv") == protein]
            resources(a.cache, files)
            script = Path(__file__).with_name("marioni.R")
            args = [a.rscript, "--vanilla", script, "proteins" if protein else "diabetes", a.cache, a.beta, staged]
            subprocess.run(args + ([a.clinical] if a.clinical else []), check=True)
        else:
            score_biolearn(a.model, a.beta).to_csv(staged, sep="\t", index_label="sample")
        staged.rename(a.out)


if __name__ == "__main__":
    main()
