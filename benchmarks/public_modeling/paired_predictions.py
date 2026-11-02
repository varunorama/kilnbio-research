"""Compare frozen continuous-outcome predictions on identical held-out donors."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import tempfile


def compare(rows, baseline, seed=17043):
    design = {}
    errors = {}
    for row in rows:
        donor, fold, model = (row[key] for key in ("donor", "fold", "model"))
        outcome, predicted = (float(row[key]) for key in ("observed", "predicted"))
        if not all((donor, fold, model)) or not all(map(math.isfinite, (outcome, predicted))):
            raise ValueError("Each row needs a donor, held-out fold, model and finite measurements")
        if design.setdefault(donor, (fold, outcome)) != (fold, outcome):
            raise ValueError(f"Inconsistent outcome or fold for donor {donor}")
        block = errors.setdefault(model, {})
        if donor in block:
            raise ValueError(f"Repeated donor {donor} in {model}")
        block[donor] = abs(predicted - outcome)
    if baseline not in errors or len(errors) < 2 or len(design) < 2:
        raise ValueError("Supply the named baseline and a comparator for at least two donors")
    donors = sorted(design)
    if any(block.keys() != design.keys() for block in errors.values()):
        raise ValueError("Every model must have the same held-out donors")
    rng = random.Random(seed)
    summary = []
    for name, block in errors.items():
        delta = [errors[baseline][donor] - block[donor] for donor in donors]
        boot = sorted(statistics.mean(rng.choices(delta, k=len(delta))) for _ in range(2000))
        summary.append((name, len(donors), statistics.mean(block.values()), statistics.mean(delta),
                        boot[49], boot[1949]))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("predictions", type=Path, help="TSV: donor, fold, observed, model, predicted")
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=17043)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Choose a new output directory")
    with args.predictions.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not {"donor","fold","observed","model","predicted"} <= set(fields):
            parser.error("Missing or duplicated TSV column")
        summary = compare(reader, args.baseline, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out.parent) as tmp:
        staged = Path(tmp) / "results"
        staged.mkdir()
        with (staged / "paired_mae.tsv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["model","donors","MAE","reduction_vs_baseline","paired_lower","paired_upper"])
            writer.writerows(summary)
        info = dict(baseline=args.baseline, seed=args.seed, resamples=2000,
                    input_sha256=hashlib.sha256(args.predictions.read_bytes()).hexdigest(),
                    interval="95% percentile interval over independent donors, conditional on frozen predictions")
        (staged / "design.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
        staged.rename(args.out)


if __name__ == "__main__":
    main()
