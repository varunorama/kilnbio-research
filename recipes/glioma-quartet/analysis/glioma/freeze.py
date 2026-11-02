"""Build the three panels using only the pinned discovery summaries."""
import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "public_modeling"))
from glioma_evidence import FILES, LAYERS, read_table


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--discovery-cache", required=True, type=Path)
    cli.add_argument("--out", required=True, type=Path)
    a = cli.parse_args()
    p, fc, joint = [read_table(a.discovery_cache, n, True) for n in FILES]
    if not p.index.is_unique or not p.index.equals(fc.index) or not p.index.equals(joint.index):
        raise ValueError("Mismatched discovery genes")
    q = np.column_stack([false_discovery_control(p[col]) for col in LAYERS])
    only = (false_discovery_control(joint.DPM) <= .05) & (q.min(axis=1) > .05)
    if (len(p), only.sum()) != (9209, 31):
        raise ValueError("Unexpected discovery universe or joint-only set")
    panels = {"joint_only": joint.loc[only, "DPM"], "RNA_only": p.rna, "protein_only": p.protein}
    with a.out.open("x", encoding="utf-8", newline="\n") as out:
        out.write("panel\trank\tgene\trna_sign\tprotein_sign\tpromoter_methylation_sign\tselection_p\n")
        for panel, metric in panels.items():
            for rank, gene in enumerate(sorted(metric.index, key=lambda g: (metric[g], g))[:31], 1):
                signs = np.sign(fc.loc[gene, LAYERS]).astype(int)
                out.write("\t".join(map(str, (panel, rank, gene, *signs, format(metric[gene], ".15g")))) + "\n")


if __name__ == "__main__":
    main()
