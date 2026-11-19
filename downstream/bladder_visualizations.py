"""Display the frozen bladder recipe's populations and ranking decisions."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from visual_reports import canvas, matrix
from plot_design import BLUE, ORANGE, TEAL, INK, MUTED, MISSING

ARM_NAMES = {'A': 'RNA / all eligible', 'B': 'RNA / matched',
             'C': 'RNA + protein / matched', 'P': 'Protein / matched'}
COMPARISONS = {'A_to_B': 'A → B\nRestrict the cohort',
               'B_to_C': 'B → C\nAdd protein evidence',
               'A_to_C': 'A → C\nOverall change'}


def bladder_report(b, root):
    path = root / 'report.json'
    obj = json.loads(path.read_text()) if path.is_file() else {}
    budget = obj.get('budget')
    if not isinstance(budget, int) or isinstance(budget, bool) or budget < 1:
        b.skip('shortlist_rank_changes', 'A positive integer shortlist budget is required')
        return
    b.report.title = 'KilnBio bladder cancer report: cohort and assay choices'
    b.report.notes.extend([
        'A = RNA in all eligible patients; B = RNA in the matched RNA–RPPA cohort; '
        'P = protein in that same cohort; C = the mean of the matched signed RNA and protein rank scores.',
        'The frozen recipe compares FGFR3 protein-altering somatic SNV/indel calls with no qualifying call. '
        'FGFR3 itself is excluded from candidate ranking. A negative grouping label does not mean mutation-free.',
        f'The recorded shortlist budget is {budget}. A → B isolates cohort restriction; '
        'B → C adds protein evidence on the same patients. All arms use the recorded common gene universe.'
    ])

    groups = obj.get('groups', {})
    population = []
    for key, label in [('RNA_all', 'A / All eligible RNA'), ('matched', 'B, C, P / Matched'), ('RNA_only', 'Excluded by matching')]:
        row = groups.get(key, {})
        if all(k in row for k in ['patients', 'positive', 'no_qualifying_call']):
            if row['positive'] + row['no_qualifying_call'] != row['patients']:
                raise ValueError('Population group counts do not sum to their recorded denominator')
            population.append(dict(population=label, **row))
    if population:
        data = pd.DataFrame(population)
        fig, ax = canvas(8, 11)
        y = np.arange(len(data))
        ax.barh(y, data.no_qualifying_call, color=BLUE, height=.53, label='No qualifying FGFR3 call')
        ax.barh(y, data.positive, left=data.no_qualifying_call, color=ORANGE, height=.53, label='Qualifying FGFR3 call')
        for i, row in data.iterrows():
            ax.text(row.no_qualifying_call/2, i, str(row.no_qualifying_call), ha='center', va='center', color='white')
            ax.text(row.no_qualifying_call+row.positive/2, i, str(row.positive), ha='center', va='center', color='white')
            ax.text(row.patients+data.patients.max()*.025, i, f'n = {row.patients:,}', va='center', color=INK)
        ax.set_yticks(y, data.population); ax.invert_yaxis()
        ax.set_xlim(0, data.patients.max()*1.23); ax.set_xlabel('Patients')
        ax.set_title('Who remains when both assays are required?')
        ax.legend(loc='upper left', bbox_to_anchor=(0, -.15), ncol=2)
        b.add('bladder_cohort_coverage', fig,
              'Recorded population counts by FGFR3 grouping. The matched and excluded populations partition '
              'the all-RNA population; the three rows are not independent cohorts.', data)
    else:
        b.skip('bladder_cohort_coverage', 'No recorded population counts')

    changes = b.load(root/'changes.tsv')
    if changes is not None and {'comparison', 'status'} <= set(changes):
        counts = pd.crosstab(changes.comparison, changes.status).reindex(
            index=[x for x in COMPARISONS if x in set(changes.comparison)],
            columns=['retained', 'entered', 'exited'], fill_value=0).fillna(0).astype(int)
        if len(counts):
            fig, ax = canvas(9, 11); y = np.arange(len(counts))
            for i, (status, color) in enumerate([('retained', BLUE), ('entered', TEAL), ('exited', ORANGE)]):
                ax.barh(y+(i-1)*.23, counts[status], height=.19, color=color, label=status.capitalize())
                for j, value in enumerate(counts[status]):
                    ax.text(value+.15, j+(i-1)*.23, str(value), va='center', fontsize=10)
            ax.set_yticks(y, [COMPARISONS[x] for x in counts.index]); ax.invert_yaxis()
            ax.set_xlim(0, max(budget, counts.to_numpy().max())+2)
            ax.set_xlabel('Genes'); ax.set_title(f'What changes within a {budget}-gene follow-up budget?')
            ax.legend(ncol=3, loc='upper left', bbox_to_anchor=(0, -.15))
            b.add('shortlist_status', fig,
                  'Gene counts are separated by comparison. Retained + entered describes the new shortlist; '
                  'retained + exited describes the previous shortlist. Counts are not pooled across comparisons.',
                  counts.rename_axis('comparison').reset_index())
    elif changes is not None:
        b.categories('shortlist_status', changes.get('status'), 'Shortlist status / comparison labels unavailable',
                     'Input records have no comparison identifier; counts cannot be attributed to a cohort or assay change.')

    t = b.load(root/'rankings.tsv')
    if t is None or not {'arm', 'gene', 'rank'} <= set(t):
        b.skip('shortlist_rank_changes', 'No gene ranks with arm labels')
        return
    t['rank'] = pd.to_numeric(t['rank'], errors='raise')
    if t.duplicated(['arm', 'gene']).any() or not np.isfinite(t['rank']).all() or (t['rank'] < 1).any():
        raise ValueError('Ranks must be finite and positive, with unique gene–arm pairs')
    ranks = t.pivot(index='gene', columns='arm', values='rank')
    arms = [arm for arm in ['A', 'B', 'C'] if arm in ranks]
    eligible = ranks[arms].le(budget).any(axis=1)
    order = ranks.loc[eligible, arms].min(axis=1).rename('best_rank').reset_index().sort_values(['best_rank', 'gene'])
    genes = order.gene.head(40).tolist()
    shown = ranks.loc[genes, arms]
    if not len(shown):
        b.skip('shortlist_rank_changes', 'No candidates within the recorded budget')
        return
    labels = [arm+'\n'+ARM_NAMES[arm].replace(' / ', '\n') for arm in arms]
    data = ranks.reset_index()
    data['displayed'] = data.gene.isin(genes)
    fig, ax = canvas(len(shown), 9)
    color_values = np.where(shown.isna(), np.nan, shown.le(budget).astype(float))
    cmap = ListedColormap(['#F0F3F6', BLUE]); cmap.set_bad(MISSING)
    ax.imshow(np.ma.masked_invalid(color_values), aspect='auto', cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(arms)), labels); ax.xaxis.tick_top()
    ax.set_yticks(range(len(genes)), genes)
    for y in range(len(genes)):
        for x in range(len(arms)):
            rank = shown.iloc[y,x]
            text = 'NA' if pd.isna(rank) else f'{rank:.0f}'
            ax.text(x,y,text,ha='center',va='center',color='white' if rank<=budget else MUTED)
    ax.set_xticks(np.arange(-.5,len(arms),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(genes),1),minor=True)
    ax.grid(which='minor',color='white',lw=2); ax.tick_params(which='minor',bottom=False,left=False)
    for spine in ax.spines.values():spine.set_visible(False)
    fig.suptitle(f'Shortlist membership and exact rank\nBlue: rank ≤ {budget} · Gray: outside the budget',x=.02,ha='left',fontsize=15,fontweight='bold')
    b.add('shortlist_rank_changes',fig,
          f'Genes nominated in any A/B/C arm; numbers are their original ranks among the recorded common universe. '
          f'{len(genes)} of {len(order)} nominated genes displayed, sorted by best A/B/C rank then gene. '
          'The TSV retains all ranked genes and records which rows were displayed.',data)

    if all(arm in ranks for arm in ['A','B','C']):
        returning = ranks.index[ranks.A.le(budget) & ranks.B.gt(budget) & ranks.C.le(budget)].sort_values().tolist()
        if returning:
            fig, ax = canvas(8, 10)
            for gene in returning[:8]:
                values = ranks.loc[gene,['A','B','C']].to_numpy(float)
                ax.plot(range(3), values, marker='o', lw=2, label=gene)
                for x, value in enumerate(values):ax.annotate(f'{gene}: {value:g}',(x,value),xytext=(0,10),textcoords='offset points',ha='center',fontsize=10)
            upper = max(budget+2, ranks.loc[returning[:8],['A','B','C']].max().max()+3)
            ax.axhspan(.5,budget+.5,color='#E9F2F8'); ax.axhline(budget+.5,color=MUTED,ls='--')
            ax.set_ylim(upper,.5); ax.set_xlim(-.25,2.25)
            ax.set_xticks(range(3),[a+'\n'+ARM_NAMES[a].replace(' / ','\n') for a in ['A','B','C']])
            ax.set_ylabel('Rank / smaller is higher'); ax.set_title('Candidates that leave and return to the shortlist')
            ax.text(.98,.03,f'Shaded area: {budget}-gene budget',transform=ax.transAxes,ha='right',color=MUTED)
            b.add('bladder_returning_candidates',fig,
                  'Display rule: rank within budget in A and C, but outside budget in B. Up to eight such genes '
                  'are drawn in alphabetical order; the TSV includes every qualifying gene. This shows the '
                  'recorded ranking behavior, not independent biological confirmation.',ranks.loc[returning].reset_index())

    if 'score' in t:
        scores=t.pivot(index='gene',columns='arm',values='score')
        columns=[x for x in ['A','B','P','C'] if x in scores]
        display=scores.reindex(index=genes,columns=columns)
        fig,ax=canvas(len(genes),11)
        matrix(fig,ax,display,genes,[a+'\n'+ARM_NAMES[a].replace(' / ','\n') for a in columns],
               'Signed within-arm rank score',limit=1)
        ax.set_title('Direction and rank evidence across RNA and protein')
        b.add('bladder_rank_evidence',fig,
              'A/B/P are signed rank scores from their original moderated t statistics; C is the recorded mean '
              'of B and P. These dimensionless scores preserve the recipe, and are not expression levels, '
              'fold changes or adjusted p-values. Display uses the same nominated genes as the rank panel.',
              display.reset_index())
