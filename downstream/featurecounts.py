"""Count fragments in annotated exons grouped by gene_id using Subread 2.1.1."""
import argparse
import csv
import json
import os
import re
import subprocess
from pathlib import Path

from analysis import digest, file_at, need, rows

IMAGE = 'quay.io/biocontainers/subread@sha256:a0f7ea960dd0f55240b330b5c4515388a493e266a7b427d6a59a39164bc14865'
STRAND = {'unstranded': '0', 'forward': '1', 'reverse': '2'}


def plan(config, out):
    cfg = json.loads(config.read_text(encoding='utf-8-sig'))
    need(set(cfg) <= {'gtf', 'libraries', 'layout', 'runtime', 'threads'}, 'Unknown featureCounts option')
    layout = cfg['layout']
    need(layout in ('paired', 'single'), 'Declare paired or single layout')
    runtime = cfg.get('runtime', 'docker')
    need(runtime in ('docker', 'local'), 'Choose docker or local runtime')
    threads = cfg.get('threads', 1)
    need(type(threads) is int and 1 <= threads <= 64, 'Use 1..64 threads')
    gtf = file_at(cfg['gtf'], config.parent)
    sheet = file_at(cfg['libraries'], config.parent)
    header, libs = rows(sheet)
    need(header == ['sample', 'bam', 'strandedness'], 'Expected sample,bam,strandedness')
    ids = [row['sample'] for row in libs]
    need(len(set(ids)) == len(ids) and all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', i) for i in ids), 'Unique safe sample IDs required')
    need(all(row['strandedness'] in STRAND for row in libs), 'Declare resolved strandedness per library')
    bams = [file_at(row['bam'], sheet.parent) for row in libs]
    need(len(set(bams)) == len(bams) and all(p.suffix == '.bam' for p in bams), 'Supply a distinct genomic BAM per sample')
    need(not out.exists(), 'Choose a new output directory')
    prefix = []
    annotation, destination = str(gtf), str(out / 'featureCounts.tsv')
    inputs = list(map(str, bams))
    if runtime == 'docker':
        annotation = '/input/' + gtf.name
        destination = '/output/featureCounts.tsv'
        inputs = [f'/bam/{name}.bam' for name in ids]
        prefix = ['docker', 'run', '--rm', '--network', 'none', '--user', f'{os.getuid()}:{os.getgid()}']
        for src, target in [(gtf, annotation), *zip(bams, inputs)]:
            prefix.extend(['-v', f'{src}:{target}:ro'])
        prefix.extend(['-v', f'{out}:/output', IMAGE])
    args = ['featureCounts', '-t', 'exon', '-g', 'gene_id', '-T', str(threads),
            '-s', ','.join(STRAND[row['strandedness']] for row in libs)]
    if layout == 'paired':
        args.extend(['-p', '--countReadPairs'])
    command = prefix + args + ['-a', annotation, '-o', destination, *inputs]
    return ids, inputs, command, prefix + ['featureCounts', '-v'], [config, sheet, gtf, *bams]


def publish_matrix(out, ids, inputs):
    assigned = [0] * len(ids)
    seen = set()
    with (out / 'featureCounts.tsv').open() as source, (out / 'counts.tmp').open('w', newline='') as target:
        table = csv.reader((l for l in source if not l.startswith('#')), delimiter='\t')
        expected = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', *inputs]
        need(next(table, []) == expected, 'Unexpected featureCounts columns or BAM order')
        writer = csv.writer(target, delimiter='\t', lineterminator='\n')
        writer.writerow(['feature_id', *ids])
        for row in table:
            need(len(row) == len(expected) and row[0] and row[0] not in seen, 'Malformed or duplicate gene row')
            values = [int(x) for x in row[6:]]
            need(all(v >= 0 for v in values), 'Invalid gene counts')
            seen.add(row[0])
            assigned = [n + v for n, v in zip(assigned, values)]
            writer.writerow([row[0], *values])
    header, stats = rows(out / 'featureCounts.tsv.summary', '\t')
    need(header == ['Status', *inputs], 'Assignment summary does not match BAM order')
    totals = next((row for row in stats if row['Status'] == 'Assigned'), {})
    need(assigned == [int(totals[p]) for p in inputs], 'Gene totals disagree with assignment summary')
    need(all(assigned), 'A library has zero assigned fragments/reads; check genomic reference, annotation and strand')
    (out / 'counts.tmp').replace(out / 'counts.tsv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    try:
        ids, inputs, command, version_cmd, files = plan(args.config.resolve(), out)
        record = {'samples': ids, 'command': command, 'inputs': list(map(str, files))}
        if args.dry_run:
            print(json.dumps(record, indent=2))
            return
        out.mkdir(parents=True, exist_ok=False)
        record['input_sha256'] = {str(p): digest(p) for p in files}
        (out / 'plan.json').write_text(json.dumps(record, indent=2) + '\n')
        version = subprocess.run(version_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
        (out / 'version.txt').write_text(version)
        need(re.search(r'featureCounts v2\.1\.1\b', version), 'Use Subread featureCounts 2.1.1')
        with (out / 'featureCounts.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        publish_matrix(out, ids, inputs)
        print(out / 'counts.tsv')
    except (KeyError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(2, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
