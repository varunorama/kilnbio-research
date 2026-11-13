"""Import delivered Olink NPX without changing its vendor normalization by default.

Accepts CSV/TSV or parquet long tables. Select a study roster explicitly. Optional
study centering is a single-plate sensitivity analysis. No instrument conversion
or negative-control LOD formula is inferred. See docs/MULTIOMIC_AND_OLINK.md.
"""
import argparse
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_long(path):
    path = Path(path)
    if path.suffix.lower() == '.parquet':
        table = pd.read_parquet(path)
    else:
        table = pd.read_csv(path, sep=None, engine='python', dtype=str, keep_default_na=False)
    required = {'SampleID','OlinkID','NPX'}
    if not required <= set(table.columns) or table.empty:
        raise ValueError('Export needs SampleID, OlinkID and delivered NPX')
    # Documented older OlinkAnalyze exports use these spellings.
    for old, new in [('Sample_Type','SampleType'), ('Assay_Warning','AssayQC')]:
        if old in table and new not in table:
            table[new] = table[old]
    for name in required - {'NPX'}:
        if table[name].isna().any() or (table[name].astype(str).str.strip() == '').any():
            raise ValueError('Missing identifier: ' + name)
        table[name] = table[name].astype(str)
    for name in ('NPX','ExtNPX','LOD'):
        if name in table:
            raw = table[name].replace({'': np.nan, 'NA': np.nan, 'NaN': np.nan})
            table[name] = pd.to_numeric(raw, errors='raise')
            if np.isinf(table[name]).any():
                raise ValueError('Infinite values in ' + name)
    return table


def rank_association(x, group):
    """Report association and exact confounding separately for categorical designs."""
    from scipy.stats import chi2_contingency
    tab = pd.crosstab(x, group)
    if min(tab.shape) < 2:
        return None, False
    chi2 = chi2_contingency(tab, correction=False)[0]
    v = float(np.sqrt(chi2 / (tab.values.sum() * (min(tab.shape) - 1))))
    design = np.column_stack([np.ones(len(x)), pd.get_dummies(x, drop_first=True, dtype=float)])
    nuisance = np.column_stack([np.ones(len(group)), pd.get_dummies(group, drop_first=True, dtype=float)])
    confounded = np.linalg.matrix_rank(np.column_stack([nuisance, design])) == np.linalg.matrix_rank(nuisance)
    return v, bool(confounded)


def run(args):
    out = Path(args.out)
    if out.exists():
        raise ValueError('Choose a new output directory')
    samples = pd.read_csv(args.samples, sep='\t', dtype=str, keep_default_na=False)
    if samples.empty or samples.columns[0] != 'sample' or samples['sample'].duplicated().any():
        raise ValueError('Roster needs a unique first column sample')
    if any(not s or s != s.strip() for s in samples['sample']):
        raise ValueError('Empty or padded sample identifier')
    d = read_long(args.export)
    selected = d[d.SampleID.isin(samples['sample'])].copy()
    if 'SampleType' in selected:
        selected = selected[selected.SampleType.str.upper().isin(['SAMPLE','BIOLOGICAL_SAMPLE'])]
    if 'AssayType' in selected:
        selected = selected[selected.AssayType.str.lower() == 'assay']
    if set(selected.SampleID) != set(samples['sample']):
        raise ValueError('Some roster samples have no biological assay rows')
    if selected.duplicated(['SampleID','OlinkID']).any():
        raise ValueError('Duplicate sample/assay measurements; resolve repeated plates or runs explicitly')
    for field in ('PlateID','WellID'):
        if field not in selected:
            selected[field] = 'unreported'
    for field in ('SampleQC','AssayQC'):
        if field not in selected:
            selected[field] = 'UNREPORTED'
        selected[field] = selected[field].fillna('UNREPORTED').astype(str).str.upper()
        if not set(selected[field]) <= {'PASS','WARN','FAIL','WARNING','UNREPORTED',''}:
            raise ValueError('Unrecognized ' + field + ' value; preserve the vendor QC vocabulary')
    if 'QC_Warning' in selected and (selected['SampleQC'] == 'UNREPORTED').all():
        # Target exports supply a warning column rather than separate QC flags.
        selected['SampleQC'] = selected['QC_Warning'].astype(str).str.upper()
        selected.loc[~selected.SampleQC.isin(['PASS','WARN','WARNING','FAIL']), 'SampleQC'] = 'UNREPORTED'
    failed = set(selected.loc[selected.SampleQC == 'FAIL', 'SampleID'])
    kept = [s for s in samples['sample'] if s not in failed]
    if not kept:
        raise ValueError('All selected samples failed vendor QC')
    selected = selected[selected.SampleID.isin(kept)].copy()
    selected['NPX_out'] = selected['NPX']
    selected['LOD_out'] = selected['LOD'] if 'LOD' in selected else np.nan
    lod_source = 'vendor_NPX' if 'LOD' in selected else 'unavailable'
    if args.lod_file:
        lf = pd.read_csv(args.lod_file, sep='\t', dtype={'OlinkID':str})
        fields = [c for c in ('LOD_NPX','LOD_ExtNPX') if c in lf]
        if 'OlinkID' not in lf or len(fields) != 1 or lf.OlinkID.duplicated().any():
            raise ValueError('LOD table needs unique OlinkID and exactly one of LOD_NPX, LOD_ExtNPX')
        values = pd.to_numeric(lf.set_index('OlinkID')[fields[0]], errors='raise')
        if not np.isfinite(values).all():
            raise ValueError('LOD file must contain finite values')
        selected['LOD_out'] = selected.OlinkID.map(values)
        if fields[0] == 'LOD_ExtNPX':
            if 'ExtNPX' not in selected:
                raise ValueError('ExtNPX-scale LOD needs delivered ExtNPX')
            selected['LOD_out'] -= selected.ExtNPX - selected.NPX
        lod_source = 'explicit_' + fields[0]
    if args.renormalize == 'study':
        if selected.PlateID.nunique() != 1 or 'unreported' in set(selected.PlateID):
            raise ValueError('Study centering requires one declared plate')
        if 'ExtNPX' not in selected or not selected.ExtNPX.notna().all():
            raise ValueError('Study centering requires complete ExtNPX')
        offset = selected.ExtNPX - selected.NPX
        spread = offset.groupby(selected.OlinkID).agg(lambda x: x.max() - x.min())
        if (spread > 1e-5).any():
            raise ValueError('Delivered scale does not have a constant ExtNPX offset per assay')
        centered = selected.ExtNPX - selected.OlinkID.map(selected.groupby('OlinkID').ExtNPX.median())
        selected['LOD_out'] += centered - selected.NPX
        selected['NPX_out'] = centered
    selected.loc[selected.AssayQC == 'FAIL', 'NPX_out'] = np.nan
    above = (selected.NPX_out > selected.LOD_out).where(selected.LOD_out.notna() & selected.NPX_out.notna())
    selected['above_LOD'] = above
    wide = selected.pivot(index='OlinkID', columns='SampleID', values='NPX_out').reindex(columns=kept)
    annotations = [c for c in ('UniProt','Assay') if c in selected]
    for col in annotations:
        if selected.groupby('OlinkID')[col].nunique().max() > 1:
            raise ValueError('Inconsistent annotation for an OlinkID: ' + col)
    features = selected[['OlinkID',*annotations]].drop_duplicates('OlinkID').set_index('OlinkID').reindex(wide.index)
    features['complete'] = wide.notna().all(axis=1)
    features['LOD_observations'] = selected.groupby('OlinkID').above_LOD.count().reindex(wide.index)
    features['fraction_above_LOD'] = selected.groupby('OlinkID').above_LOD.mean().reindex(wide.index)
    features['eligibility'] = np.where(features.complete, 'retained', 'missing_or_failed_measurements')
    if args.min_detect is not None:
        unknown = features.LOD_observations < len(kept)
        features.loc[unknown, 'eligibility'] = 'LOD_unavailable'
        features.loc[~unknown & (features.fraction_above_LOD < args.min_detect), 'eligibility'] = 'below_detection_threshold'
    eligible = wide.loc[features.eligibility == 'retained']
    samp = samples.set_index('sample').loc[kept].copy()
    # A sample may legitimately span panels/plates with distinct OlinkIDs.
    # Retain all plate/well identities; duplicate sample/OlinkID remains forbidden.
    meta = selected.groupby('SampleID')[['PlateID','WellID','SampleQC']].agg(
        lambda x: '|'.join(sorted(set(x.astype(str))))).loc[kept]
    for key in ('PlateID','WellID','SampleQC'):
        samp['olink_' + key] = meta[key].values
    one_well = meta.WellID.where(~meta.WellID.str.contains('|',regex=False))
    samp['plate_row'] = one_well.str.extract(r'^([A-Za-z]+)', expand=False).fillna('multiple_or_unreported').values
    samp['plate_column'] = one_well.str.extract(r'([0-9]+)$', expand=False).fillna('multiple_or_unreported').values
    alias = []
    for field in samples.columns[1:]:
        values = samp[field]
        if values.nunique() > max(5, len(kept) // 2):
            continue
        for nuisance in ('olink_PlateID','plate_row','plate_column'):
            v, confounded = rank_association(values, samp[nuisance])
            alias.append(dict(covariate=field, versus=nuisance, cramers_v=v,
                              fully_confounded=confounded, associated=v is not None and v >= .5))
    qc = dict(samples_listed=len(samples), samples_written=len(kept), samples_failed=len(failed),
              assays=len(wide), eligible_assays=len(eligible), normalization=args.renormalize,
              lod_source=lod_source, min_detect=args.min_detect,
              qc_unreported_rows=int(((selected.SampleQC == 'UNREPORTED') | (selected.AssayQC == 'UNREPORTED')).sum()),
              scope='Delivered NPX import. Scientific comparisons require an appropriate design and comparable products/projects.')
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix='.olink-', dir=out.parent))
    try:
        for table, name, index in [(wide.rename_axis('feature'),'npx.tsv',True),
                                   (eligible.rename_axis('feature'),'eligible_npx.tsv',True),
                                   (features,'features.tsv',True),(samp,'samples.tsv',True),
                                   (pd.DataFrame(alias,columns=['covariate','versus','cramers_v','fully_confounded','associated']),'aliasing.tsv',False),
                                   (selected[['SampleID','OlinkID','PlateID','WellID','LOD_out','above_LOD','SampleQC','AssayQC']],'measurement_qc.tsv',False)]:
            table.to_csv(tmp/name, sep='\t', index=index, na_rep='NA')
        (tmp/'qc.json').write_text(json.dumps(qc, indent=2, allow_nan=False)+'\n')
        paths = [args.export,args.samples] + ([args.lod_file] if args.lod_file else [])
        pd.DataFrame([dict(input=Path(p).name,sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in paths]).to_csv(tmp/'inputs.tsv',sep='\t',index=False)
        if out.exists():
            raise ValueError('Output appeared during the run')
        tmp.rename(out)
    finally:
        if tmp.exists():
            shutil.rmtree(tmp)
    return qc


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('export')
    p.add_argument('--samples', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--renormalize', choices=['vendor','study'], default='vendor')
    p.add_argument('--lod-file')
    p.add_argument('--min-detect', type=float, help='Optional fraction in [0,1]; missing LOD makes a feature ineligible')
    args = p.parse_args()
    if args.min_detect is not None and not 0 <= args.min_detect <= 1:
        p.error('min-detect must lie in [0,1]')
    try:
        print(json.dumps(run(args)))
    except (OSError, ValueError, ImportError) as error:
        p.exit(2, str(error)+'\n')


if __name__ == '__main__':
    main()
