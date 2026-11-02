"""Run the public matched-cohort, aging-repeatability and noise examples."""
import argparse
import csv
import hashlib
import html
import math
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile
from urllib.request import urlopen

import openpyxl

from paired_predictions import compare


HERE = Path(__file__).resolve().parent
MOFA = "https://raw.githubusercontent.com/bioFAM/MOFAdata/5c1e71bdd77188cdac00e8ce974e0a5768972b1f/data/"
FIG3 = ("https://media.springernature.com/original/springer-static/esm/"
        "art%3A10.1038%2Fs43587-022-00248-2/MediaObjects/43587_2022_248_MOESM5_ESM.xlsx")
SOURCES = {
    "CLL_data.RData": (MOFA + "CLL_data.RData", "a85c19b2a8160006fa42bf836d48745e3959567a8723eaabe73ee6efd904893a"),
    "CLL_covariates.RData": (MOFA + "CLL_covariates.RData", "3abfd16dbca5bb0afc2a500edddb8e36f2123a529e81438f06a734b15dbc1cd7"),
    "aging_fig3.xlsx": (FIG3, "c753754d60eee9d23b45c07811768cc5d6c419f9ce3efd77541f7968b34ec975"),
}
ATAC = "https://raw.githubusercontent.com/SunScript0/ATAC-clock/1b1b73881a435baeade20606823f26ef9e4a0f08/clocks/parallel/"
for model, folder, digest in [
    ("atac", "2023-02-04_15-10_atac_vs_rna", "4998791d00db2c986ba474226a705276ccca99f7754e033102d4459a91846464"),
    ("rna", "2023-02-04_15-10_rna_vs_atac", "65175571b727255bf8004ad8036227654d4a34ff2ff793428bf28f1036e3056a"),
    ("combined", "2023-02-27_12-09_multiomic", "0b7482b1a1540ed3751be1ba534fc284fa5746e117bac13c22ab59cbe884be73"),
]:
    SOURCES[f"age_{model}.tsv"] = (ATAC + folder + "/preds.tsv", digest)


def get_inputs(cache, offline):
    cache.mkdir(parents=True, exist_ok=True)
    for name, (url, digest) in SOURCES.items():
        dest = cache / name
        if not dest.exists():
            if offline:
                raise FileNotFoundError(dest)
            with tempfile.TemporaryDirectory(dir=cache) as staging:
                staged = Path(staging) / name
                with urlopen(url, timeout=90) as stream, staged.open("wb") as out:
                    shutil.copyfileobj(stream, out)
                verify(staged, digest)
                staged.replace(dest)
        verify(dest, digest)


def verify(path, digest):
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Unexpected SHA-256 for {path}")


def aging(cache, out):
    book = openpyxl.load_workbook(cache / "aging_fig3.xlsx", read_only=True, data_only=True)
    try:
        rows = iter(book["3a-f"].values)
        if tuple(next(rows)) != ("sample", 1, 2, "diff", "abs_diff", "Clock"):
            raise ValueError("Unknown Figure 3a-f source layout")
        clocks = {}
        for donor, first, second, reported_diff, reported_abs, clock in rows:
            if clock not in {"PhenoAge", "PCPhenoAge"}:
                continue
            paired = clocks.setdefault(clock, {})
            diff = second - first
            if donor in paired or not math.isclose(diff, -reported_diff, abs_tol=1e-8) or not math.isclose(abs(diff), reported_abs, abs_tol=1e-8):
                raise ValueError("Unexpected replicate pair in source data")
            paired[donor] = diff
    finally:
        book.close()
    if len(clocks["PhenoAge"]) != 36 or clocks["PhenoAge"].keys() != clocks["PCPhenoAge"].keys():
        raise ValueError("Expected identical sets of 36 paired specimens")
    with (out / "aging_repeatability.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["clock", "pairs", "median_absolute_years", "maximum_absolute_years",
                         "mean_signed_years_rep2_minus_rep1", "lower_agreement_years", "upper_agreement_years"])
        for clock, pairs in clocks.items():
            differences = list(pairs.values())
            absolute = list(map(abs, differences))
            bias = statistics.mean(differences)
            width = 1.96 * statistics.stdev(differences)
            writer.writerow([clock, len(pairs), statistics.median(absolute), max(absolute),
                             bias, bias - width, bias + width])


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def matched_age_predictions(cache, out):
    records = []
    for model in ("atac", "rna", "combined"):
        with (cache / f"age_{model}.tsv").open() as stream:
            for donor, age, prediction, fold in csv.reader(stream, delimiter="\t"):
                records.append(dict(donor=donor, fold=fold, observed=age, model=model, predicted=prediction))
    summary = compare(records, "atac")
    with (out / "matched_predictions.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    with (out / "atac_rna_aging.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["assay", "donors", "MAE_years", "MAE_reduction_vs_ATAC", "paired_lower", "paired_upper"])
        writer.writerows(summary)


def table(rows, columns):
    def render(value):
        try:
            return f"{float(value):.4f}" if "." in value else html.escape(value)
        except ValueError:
            return html.escape(value)
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = "".join("<tr>" + "".join(f"<td>{render(row[key])}</td>" for key, _ in columns) + "</tr>" for row in rows)
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def report(out):
    cohort = read_tsv(out / "cohort.tsv")
    counts = {k: sum(row[k] == "TRUE" for row in cohort) for k in ("RNA", "DNAm", "outcome", "eligible")}
    metrics = table(read_tsv(out / "performance.tsv"),
                    [("model", "Inputs"), ("n", "Test donors"), ("MAE", "Mean absolute error"), ("RMSE", "Root mean squared error")])
    differences = table(read_tsv(out / "paired_differences.tsv"),
                        [("reference", "Comparator"), ("improvement_MAE", "MAE reduction with both assays"),
                         ("lower", "Lower 95% bound"), ("upper", "Upper 95% bound")])
    clocks = table(read_tsv(out / "aging_repeatability.tsv"),
                   [("clock", "Clock"), ("pairs", "Pairs"), ("median_absolute_years", "Median absolute gap (years)"),
                    ("maximum_absolute_years", "Maximum absolute gap (years)")])
    age_model = table(read_tsv(out / "atac_rna_aging.tsv"), [("assay", "Assay"), ("donors", "Matched donors"),
                      ("MAE_years", "Mean absolute error (years)"), ("MAE_reduction_vs_ATAC", "Improvement over ATAC (years)")])
    noise = table(read_tsv(out / "noise_challenge.tsv"), [("assay", "Assay"), ("method", "Adjustment"),
                  ("technical_RMSE", "Distance to clean matrix (RMSE)"),
                  ("effect_RMSE", "Biological contrast error (RMSE)")])
    (out / "report.html").write_text(f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width"><title>KilnBio | Public modeling examples</title>
<style>body{{font:17px/1.55 system-ui;margin:3em auto;max-width:65em;padding:0 1em;color:#17303f}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;margin:1em 0}}th,td{{text-align:left;padding:.42em .8em;border-bottom:1px solid #ccd8dc}}th{{background:#eaf3f2}}a{{color:#175c6a}}</style>
<main><h1>KilnBio public modeling examples</h1>
<h2>Does methylation add predictive information to RNA?</h2>
<p><strong>KilnBio model fits.</strong> The MOFA team supplied processed measurements.
KilnBio selected the endpoint and trained and evaluated the four models below.</p>
<p>The CLL cohort contains {len(cohort)} donors. {counts['RNA']} have RNA, {counts['DNAm']} have methylation,
and {counts['outcome']} have the selected ex vivo drug response (D_001_3). {counts['eligible']} have all required data.
All models include the source gender covariate and use identical patients. Smaller errors are better.</p>{metrics}
<p>Positive differences favor measuring both assays. Intervals are percentile bootstraps of paired test-donor
errors conditional on this development fit. They omit training-set and between-study variation.</p>{differences}
<p>The authors selected variable molecular features across their cohort before publishing these matrices.
This is an internal demonstration on processed data. Plate and lab annotations are unavailable in these inputs.
<a href="https://pubmed.ncbi.nlm.nih.gov/29925568/">Original MOFA study</a></p>
<h2>How far apart can two aging scores be on the same blood specimen?</h2>
<p><strong>Reproduction.</strong> The clock authors published the paired scores and the repeatability finding.
KilnBio recalculated the paired gaps from their Figure 3 source table.</p>{clocks}
<p>These are remeasurements of the same 36 specimens on the 450K platform. The gaps describe repeatability under
the reported protocol. A treatment effect or an individual's aging rate requires a different design.
<a href="https://www.nature.com/articles/s43587-022-00248-2">Higgins-Chen et al., Figure 3</a></p>
<h2>Does adding RNA improve a blood ATAC age clock?</h2>
<p><strong>Published predictions reanalyzed by KilnBio.</strong> The source study fitted all three clocks.
KilnBio computed paired error differences and their conditional intervals.</p>{age_model}
<p>These errors were recomputed from the authors' published predictions after matching every donor, age, and fold.
This reproduces a comparison within the original study. It does not refit the clocks or establish performance in an independent cohort.
<a href="https://doi.org/10.1007/s11357-023-00986-0">ATAC-clock study</a></p>
<h2>Can we reduce a known artifact with one run per donor?</h2>
<p><strong>KilnBio simulation.</strong> The group and technical effects used for scoring were introduced by KilnBio.</p>
<p>The underlying RNA and methylation values come from 80 matched CLL donors. We planted the biological contrasts,
batch offsets, drift and control features. The clean matrix is known in this challenge only.
Lower scores are better. The first column of errors measures recovery of values.
The second measures distortion of the group difference. See <a href="noise_challenge.tsv">the full results</a> for assay-specific contrasts.</p>{noise}
<p>When no technical anchor exists, or biological group and plate coincide, the script declines adjustment.
This exercise tests known additive effects. Natural technical variation in unreplicated clinical samples remains unresolved.</p>
</main></html>""", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="new results directory")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    cache, dest = args.cache.resolve(), args.out.resolve()
    if dest.exists() or cache == dest:
        parser.error("Choose a new results directory distinct from the cache")
    get_inputs(cache, args.offline)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        out = Path(tmp) / "results"
        out.mkdir()
        subprocess.run([args.rscript, "--vanilla", HERE / "cll.R", cache, out], check=True)
        subprocess.run([args.rscript, "--vanilla", HERE / "noise_demo.R", HERE / "adjust.R", cache, out], check=True)
        aging(cache, out)
        matched_age_predictions(cache, out)
        report(out)
        with (out / "sources.tsv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["file", "public_url", "sha256"])
            writer.writerows((name, *entry) for name, entry in SOURCES.items())
        (out / "python_version.txt").write_text(f"Python {platform.python_version()}\nopenpyxl {openpyxl.__version__}\n")
        out.rename(dest)
    print(dest / "report.html")


if __name__ == "__main__":
    main()
