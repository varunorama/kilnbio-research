"""Test whether the deployed ChIP caller version changes the native interval boundaries."""
import argparse
import csv
import json
import subprocess
from pathlib import Path

from run import MACS3, ROOT

DEPLOYED = 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'


def boundaries(path):
    with path.open() as file:
        return {tuple(line.split()[:3]) for line in file if not line.startswith(('#','track','browser'))}


def check(base, output):
    with (ROOT/'experiments.tsv').open() as file:
        rows = list(csv.DictReader(file, delimiter='\t'))
    results = []
    for row in rows:
        if row['assay'] == 'ATAC':
            continue
        exp = row['experiment']
        for rep in (1,2):
            current = base/exp/f'r{rep}-MACS3'
            deployed = base/exp/f'r{rep}-MACS3-3.0.1'
            deployed.mkdir(exist_ok=True)
            status = json.loads((current/'status.json').read_text())
            if status['exit'] or MACS3 not in status['command']:
                raise ValueError(f'Invalid original command in {current}')
            cmd = [DEPLOYED if a == MACS3 else a.replace(f'/b/{exp}/{current.name}', f'/b/{exp}/{deployed.name}')
                   for a in status['command']]
            with (deployed/'tool.log').open('w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
            name = 'calls_peaks.broadPeak'
            old, new = boundaries(deployed/name), boundaries(current/name)
            results.append(dict(experiment=exp, assay=row['assay'], cell=row['cell'], rep=rep,
                                macs_3_0_1=len(old), macs_3_1_0=len(new),
                                identical_boundaries=len(old & new), only_3_0_1=len(old-new), only_3_1_0=len(new-old)))
    with output.open('w', newline='') as file:
        writer = csv.DictWriter(file, results[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    check(args.directory.resolve(), args.output)
