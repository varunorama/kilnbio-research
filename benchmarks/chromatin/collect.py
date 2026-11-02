"""Score completed runs and export compact evidence that can be shared without the alignments."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean
import numpy as np

from score import NAMES, ROOT, intervals, score_experiment


def write_tsv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: 'NA' if value == '' or isinstance(value, float) and not math.isfinite(value)
                             else value for key, value in row.items()})


def check_scores(rows, callers):
    expected = 2 * len(callers) * 2 * 4   # replicates, callers, views, regions
    if len(rows) != expected:
        raise ValueError('Incomplete scores')
    for rep in (1, 2):
        for chrom in ('panel', 'chr1', 'chr19', 'chr21'):
            budgets = [r['selected_bins'] for r in rows if r['rep'] == rep
                       and r['region'] == chrom and r['view'] == 'equal_bp']
            if len(budgets) != len(callers) or len(set(budgets)) != 1:
                raise ValueError(f'Coverage mismatch for replicate {rep}, {chrom}')


def contrasts(rows, methods, atac):
    metric = 'dnase_f1' if atac else 'heldout_enrichment'
    matched = [r for r in rows if r['view'] == 'equal_bp']
    baseline = 'MACS2' if atac else 'MACS3'

    def panel_mean(method, field):
        values = [r[field] for r in matched if r['caller'] == method and r['region'] == 'panel']
        return mean(values) if len(values) == 2 and all(isinstance(v, (int,float)) for v in values) else math.nan

    result = []
    for method in methods:
        if method == baseline:
            continue
        base_value = panel_mean(baseline, metric)
        value = panel_mean(method, metric)
        jaccard_change = panel_mean(method, 'replicate_jaccard') - panel_mean(baseline, 'replicate_jaccard')
        gain = value - base_value if atac else (value/base_value - 1 if base_value > 0 else math.nan)
        min_gain = 0.02 if atac else 0.10
        if math.isfinite(gain) and math.isfinite(jaccard_change):
            passed = 'yes' if gain >= min_gain and jaccard_change >= -0.02 else 'no'
        else:
            passed = 'unscorable'
        result.append(dict(experiment=rows[0]['experiment'], cell=rows[0]['cell'], assay=rows[0]['assay'],
                           baseline=baseline, candidate=method, metric=metric, baseline_mean=base_value,
                           candidate_mean=value, gain_type='absolute' if atac else 'fractional',
                           gain=gain, jaccard_delta=jaccard_change, passes_cell_numeric_margins=passed,
                           zero_budget_chromosome_replicates=sum(r['selected_bins']==0 for r in matched
                               if r['caller']==baseline and r['region']!='panel')))
    return result


def check_native_atac(rows, audit):
    fields = ['selected_bins','dnase_reference_bins','dnase_overlap_bins']
    expected = {(r['experiment'], str(r['rep']), r['caller']): tuple(r[f] for f in fields)
                for r in rows if r['assay']=='ATAC' and r['view']=='native' and r['region']=='panel'}
    with audit.open() as file:
        observed = list(csv.DictReader(file, delimiter='\t'))
    counts = {(r['experiment'], r['rep'], r['caller']): tuple(int(r[f]) for f in fields) for r in observed}
    if len(counts) != len(observed) or counts != expected:
        raise ValueError('Native ATAC counts disagree with the independent BEDTools calculation')


def collect(base, output, native_audit=None):
    output.mkdir(exist_ok=True, parents=True)
    scores, runs, libraries, comparisons, fragments, agreement = [], [], [], [], [], []
    with (ROOT/'experiments.tsv').open() as file:
        experiments = list(csv.DictReader(file, delimiter='\t'))
    for row in experiments:
        exp = row['experiment']
        atac = row['assay'] == 'ATAC'
        methods = ['MACS2', 'MACS3', 'HMMRATAC'] if atac else ['MACS3', 'Omnipeak']
        counts = json.loads((base/exp/'depths.json').read_text())
        for acc, n in counts['eligible'].items():
            role = next((r for r in ('rep1', 'rep2') if row[r] == acc), 'control')
            libraries.append(dict(experiment=exp, file=acc, role=role, eligible=n))
        for rep in (1, 2):
            if atac:
                boundaries = [set(intervals(base/exp/f'r{rep}-{caller}'/NAMES[caller])) for caller in ('MACS2','MACS3')]
                old, new = boundaries
                agreement.append(dict(experiment=exp, cell=row['cell'], rep=rep, macs2=len(old), macs3=len(new),
                                      identical_boundaries=len(old & new), only_macs2=len(old-new), only_macs3=len(new-old)))
                histogram = np.zeros(2001, dtype=np.int64)
                with (base/exp/f'rep{rep}.bed').open() as bed:
                    for line in bed:
                        _, start, end = line.split()
                        histogram[int(end)-int(start)] += 1
                total = int(histogram.sum())
                if total != counts['balanced_treatment']:
                    raise ValueError('ATAC fragment histogram does not match the calling depth')
                middle = np.searchsorted(histogram.cumsum(), [(total+1)//2, (total+2)//2])
                fragments.append(dict(experiment=exp, cell=row['cell'], rep=rep, fragments=total,
                                      median_bp=middle.mean(), fraction_1_100=histogram[:101].sum()/total,
                                      fraction_150_250=histogram[150:251].sum()/total,
                                      fraction_300_500=histogram[300:501].sum()/total))
            for caller in methods:
                folder = base/exp/f'r{rep}-{caller}'
                state = json.loads((folder/'status.json').read_text())
                if state['exit'] != 0:
                    raise ValueError(f'Unsuccessful caller in {folder}')
                name = 'calls_peaks.broadPeak' if not atac and caller == 'MACS3' else NAMES[caller]
                peakfile = folder/name
                if not peakfile.is_file():
                    raise FileNotFoundError(peakfile)
                log = (folder/'tool.log').read_text()
                omni_quality = ('yes' if 'Low quality of data detected' in log else 'no') if caller=='Omnipeak' else 'NA'
                hmm_em_cap = ('yes' if 'Reached maximum number (20) of iterations' in log else 'no') if caller=='HMMRATAC' else 'NA'
                runs.append(dict(experiment=exp, assay=row['assay'], cell=row['cell'], rep=rep,
                                 caller=caller, balanced_treatment=counts['balanced_treatment'],
                                 control_training=counts.get('control_training',''),
                                 control_evaluation=counts.get('control_evaluation',''),
                                 seconds=state['seconds'], peaks_sha256=hashlib.sha256(peakfile.read_bytes()).hexdigest(),
                                 empty_peak_file='yes' if peakfile.stat().st_size==0 else 'no',
                                 omnipeak_low_quality_warning=omni_quality, hmm_fragment_em_iteration_cap=hmm_em_cap))
        metrics = score_experiment(base, row)
        check_scores(metrics, methods)
        scores.extend(metrics)
        comparisons.extend(contrasts(metrics, methods, atac))
    if native_audit:
        check_native_atac(scores, native_audit)
    write_tsv(output/'metrics.tsv', scores)
    write_tsv(output/'runs.tsv', runs)
    write_tsv(output/'library_counts.tsv', libraries)
    write_tsv(output/'atac_fragment_lengths.tsv', fragments)
    write_tsv(output/'atac_macs_agreement.tsv', agreement)
    write_tsv(output/'contrasts.tsv', comparisons)
    (output/'references.json').write_bytes((base/'ref/references.json').read_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--native-audit', type=Path)
    args = parser.parse_args()
    collect(args.directory, args.output, args.native_audit)
