"""Reusable research comparisons from KilnBio recipes; no automatic biological claims.

Arrays are explicitly aligned by supplied feature/sample identifiers. See
docs/METHODOLOGY.md for theory, assumptions, provenance and use cases.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy import stats

API_VERSION = 1


def _ids(values, name, unique=True):
    out = list(values)
    if not out or any(not isinstance(x, str) or not x or x.strip() != x or any(c in x for c in '\r\n\t') for x in out):
        raise ValueError(f'{name} requires nonempty string identifiers without outer whitespace')
    if unique and len(set(out)) != len(out):
        raise ValueError(f'{name} contains duplicate identifiers')
    return out


def _integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)


def _array(values, name, ndim=1, finite=True):
    x = np.asarray(values, dtype=float)
    if x.ndim != ndim or any(s == 0 for s in x.shape):
        raise ValueError(f'{name} must be a nonempty {ndim}-dimensional array')
    if np.isinf(x).any() or (finite and not np.isfinite(x).all()):
        raise ValueError(f'{name} contains invalid or missing values')
    return x


def _flags(values, size, name):
    raw = np.asarray(values)
    if raw.shape != (size,) or any(type(v) not in (bool, np.bool_) for v in raw):
        raise ValueError(f'{name} requires {size} explicit booleans')
    return raw.astype(bool)


def _adjust(p, method):
    if method not in ('BH', 'BY'):
        raise ValueError('adjustment must be BH or BY')
    order = np.argsort(p, kind='stable')
    count = len(p)
    factor = np.sum(1 / np.arange(1, count + 1)) if method == 'BY' else 1.
    values = p[order] * count * factor / np.arange(1, count + 1)
    q = np.empty(count)
    q[order] = np.minimum(1., np.minimum.accumulate(values[::-1])[::-1])
    return q


def rank_assays(features, statistics, k, mode):
    """Rank a fixed common universe using equal-weight signed or unsigned ranks.

    `statistics` maps assay names to finite statistics aligned to `features`.
    signed uses rank(abs(t))/(G+1) * sign(t), absolute-score ordering and
    UTF-8 identifier tie breaking. unsigned uses rank(abs(t))/G and stable
    input-order ties. Directions are descriptive, not inferred probabilities.
    """
    features = _ids(features, 'features')
    if not isinstance(statistics, dict) or len(statistics) < 1:
        raise ValueError('statistics must map at least one assay to a vector')
    _ids(statistics, 'assays')
    if 'combined' in statistics:
        raise ValueError('combined is reserved for the equal-weight result')
    k = _integer(k, 'k')
    if k > len(features) or mode not in ('signed', 'unsigned'):
        raise ValueError('Require k <= number of features and mode signed or unsigned')
    scores = {}
    for assay, values in statistics.items():
        x = _array(values, assay)
        if len(x) != len(features):
            raise ValueError('Every assay must use the same ordered feature family')
        s = stats.rankdata(abs(x), method='average') / (len(x) + (mode == 'signed'))
        scores[assay] = s * np.sign(x) if mode == 'signed' else s
    scores['combined'] = np.mean(list(scores.values()), axis=0)
    result = {}
    for assay, score in scores.items():
        order = (sorted(range(len(features)), key=lambda j: (-abs(score[j]), features[j].encode('utf-8')))
                 if mode == 'signed' else np.argsort(-score, kind='stable').tolist())
        result[assay] = {'scores': score.tolist(), 'ranking': [features[j] for j in order],
                         'nominations': [features[j] for j in order[:k]]}
        if mode == 'signed':
            result[assay]['directions'] = np.sign(score).astype(int).tolist()
    return {'features': features, 'mode': mode, 'k': k, 'arms': result,
            'interpretation': 'Descriptive priorities on a fixed feature family; no P value or pooled effect.'}


def compare_nominations(rankings, k, path):
    """Compare consecutive rankings on one feature universe without causal claims."""
    path = _ids(path, 'path')
    if len(path) < 2 or not isinstance(rankings, dict) or any(a not in rankings for a in path):
        raise ValueError('path needs at least two named rankings')
    orders = {a: _ids(rankings[a], a) for a in path}
    universe = set(orders[path[0]])
    if any(set(v) != universe for v in orders.values()):
        raise ValueError('Ranking comparisons require exactly the same feature family')
    k = _integer(k, 'k')
    if k > len(universe):
        raise ValueError('k exceeds feature family')
    transitions = []
    for a, b in zip(path[:-1], path[1:]):
        old, new = set(orders[a][:k]), set(orders[b][:k])
        ra, rb = ({g: i + 1 for i, g in enumerate(orders[x])} for x in (a, b))
        transitions.append({'from': a, 'to': b, 'retained': sorted(old & new),
                            'entered': [g for g in orders[b][:k] if g not in old],
                            'exited': [g for g in orders[a][:k] if g not in new],
                            'replacements': len(new - old),
                            'rank_change': {g: rb[g] - ra[g] for g in sorted(universe)}})
    return {'k': k, 'transitions': transitions,
            'interpretation': 'Changes conditional on the supplied rankings; replacement counts are not additive.'}


def validation_yield(nominations, directions, effects, adjusted_p, testable, alpha=0.05):
    """Count directional confirmations without removing untestable follow-up slots.

    Effects and adjusted P values must be computed externally in the declared
    validation family, independently of selection. Use null/NaN for untestable rows.
    """
    ids = _ids(nominations, 'nominations')
    if not 0 < alpha < 1:
        raise ValueError('alpha must be between zero and one')
    d, e, q = (_array(x, name, finite=(name == 'directions')) for x, name in
               [(directions, 'directions'), (effects, 'effects'), (adjusted_p, 'adjusted_p')])
    ok = _flags(testable, len(ids), 'testable')
    if any(len(x) != len(ids) for x in (d, e, q)) or not np.isin(d, [-1, 1]).all():
        raise ValueError('Aligned values and directions -1 or +1 are required')
    if (not np.isfinite(e[ok]).all() or not np.isfinite(q[ok]).all()
            or ((q[ok] < 0) | (q[ok] > 1)).any()
            or np.isfinite(e[~ok]).any() or np.isfinite(q[~ok]).any()):
        raise ValueError('Testable rows need finite effects and valid q; untestable rows must be missing')
    hits = ok & (d * e > 0) & (q <= alpha)
    return {'budget': len(ids), 'testable': int(ok.sum()), 'confirmed': int(hits.sum()),
            'yield_per_slot': float(hits.mean()), 'alpha': alpha,
            'status': {g: ('unmeasured' if not ok[i] else 'confirmed' if hits[i] else 'not_confirmed')
                       for i, g in enumerate(ids)}}


def donor_partition(donors, salt, folds=5, validation_folds=(0, 1)):
    """Hash whole donors; duplicate donor entries intentionally get the same fold.

    Uses first 16 hexadecimal digits of SHA256(salt + donor) modulo folds.
    This does not stratify covariates or guarantee adequate model designs.
    """
    donors = _ids(donors, 'donors', unique=False)
    folds = _integer(folds, 'folds', 2)
    if not isinstance(salt, str) or not salt:
        raise ValueError('Provide a nonempty prespecified salt')
    vf = [_integer(x, 'validation fold', 0) for x in validation_folds]
    if len(set(vf)) != len(vf) or not vf or len(vf) >= folds or any(x >= folds for x in vf):
        raise ValueError('Specify a nonempty proper subset of validation folds')
    ids = [int(hashlib.sha256((salt + d).encode()).hexdigest()[:16], 16) % folds for d in donors]
    return {'donors': donors, 'fold': ids, 'phase': ['validation' if f in vf else 'discovery' for f in ids],
            'salt': salt, 'folds': folds, 'validation_folds': vf}


def promoter_summary(probe_ids, values, feature_map, scale, method='mean', min_probes=2, drop_probe=None):
    """Summarize a fixed probe map in M-value space, optionally omitting one probe.

    values has probe rows and sample columns. Any missing required probe/value
    makes that gene/sample untestable. beta clips at 1e-6; scale M uses supplied
    M-values. Mapping is supplied by the caller, not inferred from gene names.
    """
    probes = _ids(probe_ids, 'probe_ids')
    x = _array(values, 'values', ndim=2, finite=False)
    if len(probes) != x.shape[0] or scale not in ('beta', 'M') or method not in ('mean', 'median'):
        raise ValueError('Require matching probe rows, scale beta/M and method mean/median')
    _integer(min_probes, 'min_probes', 2)
    if not isinstance(feature_map, dict):
        raise ValueError('feature_map must map genes to fixed probe identifiers')
    genes = _ids(feature_map, 'genes')
    mapping = {g: _ids(feature_map[g], f'probes for {g}') for g in genes}
    if drop_probe is not None and drop_probe not in {p for ps in mapping.values() for p in ps}:
        raise ValueError('drop_probe is not in the declared map')
    if scale == 'beta':
        finite = x[np.isfinite(x)]
        if ((finite < 0) | (finite > 1)).any():
            raise ValueError('Beta values must lie in [0,1] or be missing')
        beta = np.clip(x, 1e-6, 1 - 1e-6)
        x = np.log2(beta / (1 - beta))
    lookup = {p: i for i, p in enumerate(probes)}
    summaries, counts, missing = [], [], []
    for gene in genes:
        required = [p for p in mapping[gene] if p != drop_probe]
        absent = [p for p in required if p not in lookup]
        counts.append(len(required)); missing.append(absent)
        if len(required) < min_probes or absent:
            summaries.append(np.full(x.shape[1], np.nan))
        else:
            # np.mean/median deliberately propagate missing values; no changing map.
            summaries.append((np.mean if method == 'mean' else np.median)(x[[lookup[p] for p in required]], axis=0))
    summary = np.asarray(summaries)
    return {'genes': genes, 'values': summary.tolist(), 'testable': np.isfinite(summary).tolist(),
            'required_probe_counts': counts, 'missing_mapped_probes': missing,
            'method': method, 'scale': 'M', 'drop_probe': drop_probe}


def fit_hc3(features, sample_ids, response, design, contrast, min_residual_df=10):
    """OLS contrast with HC3 covariance and a residual-df t reference.

    response is samples x features; design is samples x terms and must explicitly
    include an intercept when intended. Rows must be independent biological units.
    All covariates, filtering and the contrast are caller-specified before testing.
    """
    features, samples = _ids(features, 'features'), _ids(sample_ids, 'sample_ids')
    y, x = _array(response, 'response', 2), _array(design, 'design', 2)
    c = _array(contrast, 'contrast')
    minimum = _integer(min_residual_df, 'min_residual_df', 1)
    if y.shape != (len(samples), len(features)) or x.shape[0] != len(samples) or len(c) != x.shape[1]:
        raise ValueError('Response, design, sample identifiers and contrast are misaligned')
    if not np.any(c) or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Contrast must be nonzero and design full rank')
    df = len(samples) - x.shape[1]
    if df < minimum:
        raise ValueError('Insufficient residual degrees of freedom')
    pinv = np.linalg.pinv(x)
    coef = pinv @ y
    residual = y - x @ coef
    leverage = (x * pinv.T).sum(axis=1)
    if np.any(leverage >= 1 - 1e-8):
        raise ValueError('HC3 is undefined at unit leverage; revise the design/eligibility')
    effect = c @ coef
    se = np.sqrt(((c @ pinv) / (1 - leverage)) ** 2 @ residual ** 2)
    valid = np.isfinite(se) & (se > 1e-12)
    t = np.divide(effect, se, out=np.full(len(features), np.nan), where=valid)
    p = 2 * stats.t.sf(abs(t), df)
    return {'features': features, 'sample_ids': samples, 'effect': effect.tolist(), 'se': se.tolist(),
            't': t.tolist(), 'p': p.tolist(), 'testable': valid.tolist(), 'df': df,
            'max_leverage': float(leverage.max()),
            'interpretation': 'HC3 protects against heteroskedasticity, not confounding or dependent rows.'}


def directional_conjunction(features, t_statistics, degrees_freedom, directions, eligibility, adjustment='BH'):
    """Validate all supplied assay-specific directions, frozen before validation.

    Every assay maps to aligned t statistics, directions and eligibility flags;
    degrees_freedom maps assays to positive residual df. Ineligible tests receive
    P=1 and stay in the family. The family must be the full frozen nomination union.
    """
    features = _ids(features, 'features')
    if not isinstance(t_statistics, dict) or len(t_statistics) < 2:
        raise ValueError('Require at least two assays')
    assays = _ids(t_statistics, 'assays')
    if any(not isinstance(m, dict) or set(m) != set(assays) for m in (degrees_freedom, directions, eligibility)):
        raise ValueError('All validation inputs must name the same assays')
    tails, valid_all = {}, np.ones(len(features), dtype=bool)
    for assay in assays:
        tt, dd = _array(t_statistics[assay], assay, finite=False), _array(directions[assay], 'directions')
        ok = _flags(eligibility[assay], len(features), 'eligibility')
        df = degrees_freedom[assay]
        if not np.isscalar(df) or isinstance(df, (bool, str)) or not np.isfinite(df) or df <= 0:
            raise ValueError('Each assay needs positive finite residual df')
        if len(tt) != len(features) or len(dd) != len(features) or not np.isin(dd, [-1, 1]).all() or not np.isfinite(tt[ok]).all():
            raise ValueError('Statistics/directions must align; eligible t must be finite and directions +/-1')
        p = np.ones(len(features)); p[ok] = stats.t.sf(tt[ok] * dd[ok], df)
        tails[assay] = p
        valid_all &= ok
    joint = np.max(list(tails.values()), axis=0)
    return {'features': features, 'directional_p': {a: p.tolist() for a, p in tails.items()},
            'conjunction_p': joint.tolist(), 'q': _adjust(joint, adjustment).tolist(),
            'testable': valid_all.tolist(), 'adjustment': adjustment,
            'interpretation': 'Support in every frozen assay direction; no mediation or causal claim.'}


def promoter_sensitivity(features, rna_t, mean_t, median_t, k=25):
    """Compare only promoter-summary changes on identical genes and RNA statistics."""
    a = rank_assays(features, {'RNA': rna_t, 'methylation': mean_t}, k, 'unsigned')
    b = rank_assays(features, {'RNA': rna_t, 'methylation': median_t}, k, 'unsigned')
    ar, br = a['arms']['combined'], b['arms']['combined']
    change = compare_nominations({'mean': ar['ranking'], 'median': br['ranking']}, k, ['mean', 'median'])
    # Correlation is undefined for a constant score; expose missing, never zero.
    sa, sb = np.asarray(ar['scores']), np.asarray(br['scores'])
    rho = float(stats.spearmanr(sa, sb).statistic) if np.ptp(sa) and np.ptp(sb) else float('nan')
    return {'k': k, 'mean': ar, 'median': br, 'comparison': change['transitions'][0],
            'score_spearman': rho, 'interpretation': 'Decision sensitivity, not independent validation.'}


def probe_coherence(probe_ids, effects, summary_direction):
    """Describe probe effect signs; does not fit effects or test mediation."""
    probes = _ids(probe_ids, 'probe_ids')
    e = _array(effects, 'effects', finite=False)
    if len(e) != len(probes) or summary_direction not in (-1, 1) or isinstance(summary_direction, bool):
        raise ValueError('Align probes/effects and specify summary_direction +/-1')
    tested = np.isfinite(e); concordant = tested & (e * summary_direction > 0)
    return {'probes': probes, 'tested': int(tested.sum()), 'missing': int((~tested).sum()),
            'same_direction': int(concordant.sum()), 'opposite_direction': int((tested & (e * summary_direction < 0)).sum()),
            'zero_effect': int((tested & (e == 0)).sum()),
            'fraction_of_tested': float(concordant.sum() / tested.sum()) if tested.any() else float('nan'),
            'all_declared_same_direction': bool(concordant.all())}


def element_methylation(methylated, depth, min_depth=5, min_cpgs=3):
    """Average replicate depth-weighted fractions over common covered CpG dyads.

    Inputs have replicate rows and the same pre-aligned CpG dyad columns. Methylated
    read equivalents may be fractional for deposited percentage x depth inputs.
    """
    m, d = _array(methylated, 'methylated', 2, False), _array(depth, 'depth', 2, False)
    _integer(min_depth, 'min_depth'); _integer(min_cpgs, 'min_cpgs')
    if m.shape != d.shape or m.shape[0] < 2:
        raise ValueError('Require at least two replicates with aligned CpG columns')
    present = np.isfinite(m) & np.isfinite(d)
    if (np.any(m[np.isfinite(m)] < 0) or np.any(d[np.isfinite(d)] < 0)
            or np.any(m[present] > d[present])):
        raise ValueError('Require 0 <= methylated <= depth')
    keep = np.all(present & (d >= min_depth), axis=0)
    count = int(keep.sum())
    fractions = (m[:, keep].sum(axis=1) / d[:, keep].sum(axis=1)) if count >= min_cpgs else None
    return {'testable': fractions is not None, 'common_cpgs': count, 'eligible_columns': np.flatnonzero(keep).tolist(),
            'replicate_fractions': fractions.tolist() if fractions is not None else None,
            'methylation': float(fractions.mean()) if fractions is not None else None}


CATALOG = {
    'rank_assays': ('Prioritize a fixed common feature family', 'Comparable contrasts; prespecified eligible features; equal assay weights', 'Signed bladder or unsigned GTEx priorities; scores are not P values'),
    'compare_nominations': ('Separate patient-set and added-assay changes', 'Same feature universe and fixed nomination budget across the requested path', 'Entered/exited targets and rank changes; descriptive rather than causal'),
    'validation_yield': ('Evaluate frozen follow-up nominations', 'Independent validation; frozen directions; q values from a declared full family', 'Confirmed count per original budget; untestable slots remain'),
    'donor_partition': ('Keep all observations of each donor in one partition', 'Prespecify salt/folds before outcomes; verify balance and adequate designs', 'Deterministic allocation; does not guarantee independent related donors'),
    'promoter_summary': ('Build mean/median or leave-one-probe-out promoter measurements', 'Fixed assembly-compatible annotation; matched assay scale; no outcome-driven mapping', 'M-value summaries and explicit untestable states; no transcriptional mechanism claim'),
    'fit_hc3': ('Estimate an adjusted contrast in continuous measurements', 'Independent rows; full-rank prespecified design; adequate residual df; no unit leverage', 'OLS effects with HC3 SE and t reference; not a count model or confounding remedy'),
    'directional_conjunction': ('Test both or all frozen assay directions in held-out data', 'Directions fixed independently; full nomination union; valid t reference', 'Maximum directional P followed by BH/BY; all-assay association rather than mediation'),
    'promoter_sensitivity': ('Measure nomination changes from mean versus median summaries', 'Same donors, RNA statistics, gene family and budget; summary is the changed input', 'Replacements and score agreement; no independent replication'),
    'probe_coherence': ('Inspect whether probe signs agree with a promoter direction', 'Effects from the same prespecified contrast and scale; report missing probes', 'Descriptive signs; a high fraction does not establish significance'),
    'element_methylation': ('Summarize shared covered CpG dyads within an element', 'Identical dyad order, compatible chemistry, assembly and replicate units', 'Equal replicate average of depth-weighted fractions; unmeasured is not zero'),
}


def catalog():
    """Machine-readable chooser; selecting a method remains a study-design decision."""
    return {'api_version': API_VERSION, 'methods': {name: {'purpose': v[0], 'assumptions': v[1], 'interpretation': v[2]}
                                                   for name, v in CATALOG.items()},
            'existing_interfaces': {'multiomic-meta': 'Directional partial conjunction and BY evidence combination',
                                    'candidate-evidence': 'Endpoint-complete validation evidence records',
                                    'assay-audit': 'Specimen matching, reach and fallback-policy comparisons',
                                    'method-agreement': 'Agreement and conflicts among planned methods',
                                    'adjust_assay / project_control_adjustment': 'Reusable R functions in recipes/glioma-quartet/analysis/quartet/adjust.R'},
            'status': 'New reusable interfaces extracted from research recipes; archived results retain their original code.'}


def _clean(value):
    if isinstance(value, dict): return {k: _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_clean(v) for v in value]
    if isinstance(value, (np.integer, np.floating, np.bool_)): value = value.item()
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('catalog', help='Show available methods, assumptions and interpretation')
    run = sub.add_parser('run', help='Run one declared method from JSON; use Python API for large arrays')
    run.add_argument('specification', type=Path)
    run.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.action == 'catalog':
            print(json.dumps(catalog(), indent=2)); return
        if args.out.exists(): raise ValueError('Choose a new output file')
        raw = args.specification.read_bytes()
        spec = json.loads(raw)
        if (not isinstance(spec, dict) or set(spec) != {'schema_version', 'method', 'arguments'}
                or type(spec['schema_version']) is not int or spec['schema_version'] != 1):
            raise ValueError('Require schema_version=1, method and arguments')
        name = spec['method']
        if name not in CATALOG or not isinstance(spec['arguments'], dict): raise ValueError('Unknown method or invalid arguments')
        result = globals()[name](**spec['arguments'])
        payload = {'api_version': API_VERSION, 'method': name, 'method_context': catalog()['methods'][name],
                   'input_sha256': hashlib.sha256(raw).hexdigest(),
                   'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'result': result}
        rendered = json.dumps(_clean(payload), indent=2, allow_nan=False) + '\n'
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as handle: handle.write(rendered)
    except (ValueError, TypeError, KeyError, OSError) as error:
        parser.exit(2, f'{parser.prog}: {error}\n')


if __name__ == '__main__':
    main()
