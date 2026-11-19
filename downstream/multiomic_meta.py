"""Directional gene-level partial conjunction across molecular assays.

Primary test: Bonferroni feature aggregation within each direction, Bonferroni
partial conjunction across layers, two-direction correction, then BY across genes.
This tolerates arbitrary dependence given valid directional input p-values. It
tests at least r associated layers, not independent replication. Existing methods
are cited in docs/MULTIOMIC_AND_OLINK.md. ACAT is an exploratory global-null result.
"""
import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def adjust(p, method='BY'):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    mask = np.isfinite(p)
    n = int(mask.sum())
    if n:
        order = np.argsort(p[mask], kind='stable')
        factor = np.sum(1 / np.arange(1, n + 1)) if method == 'BY' else 1
        ranked = p[mask][order] * n * factor / np.arange(1, n + 1)
        values = np.empty(n)
        values[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
        q[mask] = values
    return q


def directional_pc(up, down, r):
    """Union over the two directions of a Bonferroni partial conjunction test."""
    m = np.isfinite(up).sum(axis=1)
    factor = np.maximum(1, m - r + 1)
    pu = np.sort(np.where(np.isfinite(up), up, np.inf), axis=1)[:, r - 1] * factor
    pdn = np.sort(np.where(np.isfinite(down), down, np.inf), axis=1)[:, r - 1] * factor
    p = np.minimum(1, 2 * np.minimum(pu, pdn))
    p[m < r] = np.nan
    direction = np.where(pu <= pdn, 1, -1).astype(float)
    direction[m < r] = np.nan
    return p, direction


def acat(p, weights):
    """Stable Cauchy combination, an asymptotic global-null diagnostic."""
    x = np.asarray(p, float)
    w = np.where(np.isfinite(x), weights, 0)
    bounded = np.clip(x, 1e-300, 1-1e-15)
    transformed = np.tan(np.pi * (.5-bounded))
    tiny = bounded < 1e-12
    transformed[tiny] = 1 / (np.pi * bounded[tiny])
    total = w.sum(axis=1)
    stat = np.nansum(w * transformed, axis=1) / np.where(total > 0, total, 1)
    result = .5 - np.arctan(stat) / np.pi
    large = stat > 1e12
    result[large] = 1 / (np.pi * stat[large])
    result[total == 0] = np.nan
    return np.clip(result, 0, 1)


def load_layer(row, base):
    path = base/row['results']
    table = pd.read_csv(path, sep='\t', dtype={'feature':str})
    if not {'feature','effect','p'} <= set(table) or table.feature.isna().any() or table.feature.duplicated().any():
        raise ValueError('Each result needs unique feature, effect and two-sided p')
    table['p'] = pd.to_numeric(table.p, errors='raise')
    table['effect'] = pd.to_numeric(table.effect, errors='raise')
    observed = table.p.notna() & table.effect.notna()
    if ((table.loc[observed,'p'] <= 0) | (table.loc[observed,'p'] > 1)).any() or np.isinf(table[['p','effect']]).any().any():
        raise ValueError('p must be in (0,1], effect finite, or missing')
    if ((table.effect == 0) & (table.p < 1)).any():
        raise ValueError('Zero effect with p<1 cannot specify a direction')
    inputs = [path]
    if row['feature_map'] == 'identity':
        mapping = pd.DataFrame({'feature':table.feature,'gene':table.feature})
    else:
        inputs.append(base/row['feature_map'])
        mapping = pd.read_csv(inputs[-1], sep='\t', dtype=str)
        if not {'feature','gene'} <= set(mapping) or mapping[['feature','gene']].isna().any().any():
            raise ValueError('Feature map needs nonempty feature and gene')
        mapping = mapping[['feature','gene']].drop_duplicates()
    joined = mapping.merge(table, on='feature', how='left', validate='many_to_one')
    orient = int(row['orientation'])
    sign = np.sign(joined.effect) * orient
    # A signed two-sided p defines directional tails for symmetric feature tests.
    joined['up'] = np.where(sign >= 0, joined.p/2, 1-joined.p/2)
    joined['down'] = np.where(sign <= 0, joined.p/2, 1-joined.p/2)
    joined.loc[joined.effect.isna(), ['up','down']] = np.nan
    rows = []
    for gene, group in joined.groupby('gene', sort=True):
        k = len(group)
        complete = group[['up','down']].notna().all().all()
        up = min(1, k * group.up.min()) if complete else np.nan
        down = min(1, k * group.down.min()) if complete else np.nan
        rows.append(dict(gene=gene, up=up, down=down, p=min(1, 2*min(up,down)) if complete else np.nan, features=k))
    if not rows:
        raise ValueError('No mapped genes in layer ' + row['layer'])
    return pd.DataFrame(rows).set_index('gene'), inputs


def run(args):
    out = Path(args.out)
    if out.exists():
        raise ValueError('Choose a new output directory')
    manifest = pd.read_csv(args.manifest, sep='\t', dtype=str, keep_default_na=False)
    required = {'layer','results','feature_map','orientation'}
    if not required <= set(manifest) or len(manifest) < 2 or manifest.layer.duplicated().any():
        raise ValueError('Need at least two unique layers with results, feature_map and orientation')
    if not manifest.layer.str.fullmatch('[A-Za-z][A-Za-z0-9_]*').all() or not set(manifest.orientation) <= {'1','-1'}:
        raise ValueError('Use simple layer identifiers and orientation 1 or -1')
    if not 2 <= args.min_layers <= len(manifest) or not 0 < args.alpha < 1:
        raise ValueError('Invalid alpha or number of layers')
    base = Path(args.manifest).resolve().parent
    loaded = [load_layer(row, base) for _,row in manifest.iterrows()]
    genes = sorted(set().union(*(set(t.index) for t,_ in loaded)))
    up = np.column_stack([t.up.reindex(genes) for t,_ in loaded])
    down = np.column_stack([t.down.reindex(genes) for t,_ in loaded])
    p_layers = np.column_stack([t.p.reindex(genes) for t,_ in loaded])
    p, direction = directional_pc(up,down,args.min_layers)
    # All declared genes stay in the multiplicity family, including untestable ones.
    q = adjust(np.where(np.isfinite(p),p,1))
    q[~np.isfinite(p)] = np.nan
    weights = pd.to_numeric(manifest.get('weight',pd.Series(1,index=manifest.index)).replace('',1)).to_numpy(float)
    if not (np.isfinite(weights) & (weights > 0)).all():
        raise ValueError('Weights must be positive and finite')
    p_acat = acat(p_layers,weights)
    within = np.column_stack([adjust(np.where(np.isfinite(p_layers[:,i]),p_layers[:,i],1)) for i in range(len(manifest))])
    result = pd.DataFrame(dict(gene=genes,layers_measured=np.isfinite(p_layers).sum(axis=1),p_pc=p,q_pc=q,
                               pc_direction=direction,p_acat_exploratory=p_acat))
    for i,name in enumerate(manifest.layer):
        result['p_'+name] = p_layers[:,i]
        result['p_up_'+name] = up[:,i]
        result['p_down_'+name] = down[:,i]
        result['q_within_'+name] = np.where(np.isfinite(p_layers[:,i]),within[:,i],np.nan)
    result['class'] = np.where(q<=args.alpha,'convergent',np.where(np.isfinite(p),'not_selected','untestable'))
    nominations=[]
    for k in np.flatnonzero(q<=args.alpha):
        for i,row in manifest.iterrows():
            favored = 1 if up[k,i] < down[k,i] else -1 if down[k,i] < up[k,i] else 0
            if favored == direction[k]:
                nominations.append(dict(gene=genes[k],assay=row['layer'],direction=int(direction[k])*int(row['orientation'])))
    summary = dict(primary_test='directional_bonferroni_partial_conjunction',multiplicity='Benjamini-Yekutieli',
                   feature_aggregation='directional_Bonferroni',family_genes=len(genes),testable_genes=int(np.isfinite(p).sum()),
                   min_layers=args.min_layers,alpha=args.alpha,convergent_genes=int((q<=args.alpha).sum()),
                   assumptions='Valid symmetric directional feature tests; prespecified feature maps, orientation and assay selection.',
                   interpretation='At least r associated layers in the selected direction. Same-donor agreement requires external validation.',
                   exploratory_test='ACAT global null; no convergence claim or discovery selection from ACAT')
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.multiomic-',dir=out.parent))
    try:
        result.sort_values(['p_pc','gene']).to_csv(stage/'genes.tsv',sep='\t',index=False,na_rep='NA')
        pd.DataFrame(nominations,columns=['gene','assay','direction']).to_csv(stage/'nominations.tsv',sep='\t',index=False)
        (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        paths=[Path(args.manifest),*(p for _,inputs in loaded for p in inputs)]
        pd.DataFrame([dict(input=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]).to_csv(stage/'inputs.tsv',sep='\t',index=False)
        from visual_reports import multiomic_report
        multiomic_report(stage, result, manifest, summary, circos_config=getattr(args,'circos_config',None))
        if out.exists():
            raise ValueError('Output appeared during execution')
        stage.rename(out)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest')
    parser.add_argument('--out',required=True)
    parser.add_argument('--min-layers',type=int,default=2)
    parser.add_argument('--alpha',type=float,default=.05)
    parser.add_argument('--circos-config',help='Optional coordinate/track manifest; appends a Circos-style panel without changing tests')
    args=parser.parse_args()
    try:
        print(json.dumps(run(args)))
    except (ValueError,OSError,KeyError) as error:
        parser.exit(2,str(error)+'\n')


if __name__=='__main__':
    main()
