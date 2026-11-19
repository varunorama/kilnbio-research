"""Shared presentation layer. No statistical inference or data selection."""
import re

import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

INK = '#203447'
MUTED = '#637587'
BLUE = '#286FA3'
ORANGE = '#CE713B'
TEAL = '#14857E'
MISSING = '#E0E6EC'
PALETTE = [BLUE, ORANGE, TEAL, '#8B63A9', '#B5A34C', '#63819A', '#AB6475', '#586A48',
           '#7C9DCC', '#DBA16C', '#59A99F', '#B897CB', '#C9B878', '#8CA4B5',
           '#CF97A3', '#8D9E70', '#3B4E76', '#8B542E', '#405E5A', '#674775']
CMAP = LinearSegmentedColormap.from_list('kilnbio_evidence', [BLUE, '#F9FAFC', ORANGE])
CMAP.set_bad(MISSING)


def configure():
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 11,
        'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': '#C4CED7',
        'axes.titlecolor': INK, 'axes.titlesize': 14, 'axes.titleweight': 'bold',
        'axes.titlelocation': 'left', 'axes.titlepad': 16, 'axes.labelpad': 9,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.prop_cycle': mpl.cycler(color=PALETTE),
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'xtick.labelsize': 10, 'ytick.labelsize': 10,
        'xtick.major.size': 0, 'ytick.major.size': 0,
        'xtick.major.pad': 7, 'ytick.major.pad': 7,
        'legend.frameon': False, 'legend.fontsize': 9,
        'legend.title_fontsize': 9, 'legend.labelspacing': .7,
        'grid.color': '#E8EDF2', 'grid.linewidth': .7,
        'axes.axisbelow': True, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'savefig.dpi': 200, 'pdf.fonttype': 42,
        'svg.fonttype': 'none', 'lines.linewidth': 1.5,
        'boxplot.patchartist': True, 'boxplot.boxprops.color': BLUE,
        'boxplot.medianprops.color': ORANGE,
        'figure.constrained_layout.w_pad': .12, 'figure.constrained_layout.h_pad': .12,
        'figure.constrained_layout.wspace': .08, 'figure.constrained_layout.hspace': .12,
    })


def readable(value):
    text = str(value).replace('_', ' ')
    return re.sub(r'\b(qc|pca|umap|npx|qq|ma|vep)\b', lambda m:m[0].upper(), text, flags=re.I)


def finish_figure(fig):
    """Apply consistent spacing and table treatment to all typed reports."""
    for ax in fig.axes:
        if ax.get_label() == '<colorbar>':
            ax.tick_params(labelsize=9)
            continue
        if ax.axison and not ax.images:
            ax.grid(axis='y', alpha=.7)
        for title in (ax.title, ax._left_title, ax._right_title):
            if title.get_text():
                title.set_text(readable(title.get_text()))
        for table in ax.tables:
            for (row, col), cell in table.get_celld().items():
                cell.set_edgecolor('white'); cell.set_linewidth(1.2)
                cell.PAD = .08
                if row == 0:
                    cell.set_facecolor(INK); cell.set_text_props(color='white', weight='bold')
                else:
                    cell.set_facecolor('#F0F4F8' if row % 2 else '#FAFBFD')
                    cell.set_text_props(color=INK)
    if fig._supxlabel is None:
        fig.supxlabel('KILNBIO  /  ANALYSIS REPORT', fontsize=8, color=MUTED, ha='right', x=.985)


def panel_guide(name):
    """Editorial interpretation, independent of a run's observed findings."""
    rules = [
        (['circos_tracks'], 'Multi-omics in shared genomic context',
         'A common genomic view makes it easier to inspect where different assays have measurements and which regions are explicitly linked.',
         'Read rings from outside to inside using their separate scales. Dots mark supplied interval midpoints; curves represent declared links. Gray crosses are missing values, while blank positions have no supplied record. Co-location does not establish an interaction.'),
        (['bladder_cohort_coverage'], 'Which patients remain in the comparison?',
         'Adding an assay can change the study population before it changes the evidence used to prioritize genes.',
         'Compare all eligible RNA patients with the matched subset and the excluded patients. Colored segments retain the original FGFR3 grouping; matched and excluded patients partition the first row.'),
        (['shortlist_status'], 'Separate cohort effects from added assay evidence',
         'A different shortlist can reflect either a different cohort or additional molecular evidence.',
         'Read each comparison separately. A to B restricts the cohort; B to C adds protein on those same patients. Retained plus entered equals the new shortlist, while retained plus exited equals the previous shortlist.'),
        (['shortlist_rank_changes'], 'Which candidates remain within the follow-up budget?',
         'A fixed experimental budget makes it useful to see exactly which candidates each analysis would nominate.',
         'Each row is a gene and each column is a declared analysis arm. Numbers give exact ranks; blue cells fall within the recorded budget and gray cells fall outside it.'),
        (['bladder_returning_candidates'], 'Why a candidate can drop out and return',
         'Looking only at the first and last shortlist can hide changes caused by restricting the cohort.',
         'Follow a gene from all-RNA to matched-RNA to joint evidence. The shaded area marks the same follow-up budget at every stage. Returning reflects the recorded ranking, not a new validation experiment.'),
        (['bladder_rank_evidence'], 'How RNA and protein contribute to prioritization',
         'Separate assay scores reveal the direction and relative strength of the evidence behind a combined rank.',
         'Blue and copper indicate negative and positive signed rank scores. C averages the matched RNA and protein scores. The colors show ranking evidence, not expression magnitude or statistical significance.'),
        (['evidence_heatmap'], 'Evidence across assays',
         'A candidate can look convincing in one assay while being unmeasured or discordant in another.',
         'Read each row across assays, then check measured-layer coverage and joint adjusted evidence. Color shows the declared direction; a missing cell is not a negative result.'),
        (['pathway_evidence'], 'Pathway evidence across assays',
         'Different molecular layers can implicate the same process through different genes.',
         'Compare pathway support across assays, then inspect coverage and contributing genes. Signed and unsigned upstream tests have different interpretations.'),
        (['combined_ranking','pathway_ranking'], 'Ranked evidence',
         'A prioritized list needs both an evidence scale and an explanation of how much information supports each entry.',
         'Use the axis, threshold and size legend together. Prioritization helps choose follow-up candidates; it does not make the candidates experimentally validated.'),
        (['volcano'], 'Effect size and statistical evidence',
         'A small p-value and a large biological effect answer different questions. Seeing both helps avoid selecting candidates on only one criterion.',
         'Horizontal position is the reported effect; height is raw-p evidence. Colored points meet the stated adjusted-p display threshold. Labels identify a small set for inspection, not a new selection rule.'),
        (['_ma'], 'Abundance and effect size',
         'Low-abundance features often have different precision from well-measured features.',
         'Check whether the largest effects concentrate among weakly measured features before choosing them for follow-up.'),
        (['_qq'], 'P-value distribution check',
         'Unexpected p-value behavior can reveal a mismatch between the model and the data, as well as real associations.',
         'Compare points with the reference line. Departures deserve investigation alongside design and QC; this chart alone cannot identify their cause.'),
        (['manhattan','methylation_locus','methylation_sample_locus','panel_intervals'], 'Where the evidence occurs',
         'Genomic context helps distinguish an isolated observation from a regional pattern and connects results to measurable loci.',
         'Inspect coordinates and local coverage together. Gaps mean no displayed measurement; they do not establish absence of an effect.'),
        (['pca','correlation'], 'Sample relationships',
         'Unexpected sample groupings or unusually dissimilar samples can change the interpretation of a downstream comparison.',
         'Compare patterns with sample condition and batch metadata. Investigate outliers rather than removing them solely because of this plot.'),
        (['metadata_association'], 'Relationship between study variables',
         'When condition and technical variables overlap, their effects can be difficult or impossible to separate statistically.',
         'Inspect the recorded association between metadata variables before fitting a model. This is a design diagnostic; it does not remove confounding.'),
        (['spikein'], 'External-reference normalization',
         'An external reference can help assess global changes that ordinary library-size normalization might obscure.',
         'Compare target counts, reference counts and recorded factors. Interpretation depends on how the reference material was added to the experiment.'),
        (['oncoprint','variants_alleles','CH_'], 'Genomic states and candidate evidence',
         'Variant-level evidence helps connect inherited or acquired genomic differences to other molecular measurements.',
         'Check observed state, allele fraction, depth and review flags. Unknown calls stay distinct from reference calls, and research candidates still require appropriate confirmation.'),
        (['missing','coverage','reach','denominator','eligibility','measurement_patterns','states'], 'What was actually measurable?',
         'Apparent improvements can arise when an analysis quietly changes which samples or features it includes.',
         'Check the measured denominator and missing categories before comparing results. Coverage and eligibility are different from positive biological evidence.'),
        (['method_calls','method_intersections','method_effects'], 'Sensitivity to analysis method',
         'Different methods can produce different candidate lists from the same measurements.',
         'Look for stable directions and calls, then inspect disagreements. Method agreement is a robustness description, not an extra significance test.'),
        (['rank','shortlist'], 'How prioritization changes',
         'An additional assay or a different analysis choice is useful when its effect on the follow-up list is visible.',
         'Inspect which candidates move and compare their eligibility. Rank movement is meaningful only within the declared comparison set.'),
        (['concordance'], 'Agreement between molecular layers',
         'Combining assays is more informative when agreement, disagreement and missing evidence remain distinguishable.',
         'Compare direction counts with the total declared feature set. Matching directions across assays are not automatically independent confirmations.'),
        (['cell_expression_dotplot'], 'Expression and detection by cluster',
         'A cluster average can hide whether most cells or just a small fraction express a gene.',
         'Read color as mean expression and dot area as the fraction of cells with a nonzero measurement. These displayed genes are not automatically tested cell-type markers.'),
        (['umap','cell_composition'], 'Single-cell structure and composition',
         'Cell populations and sample composition provide context for interpreting donor-level differences.',
         'Check sample mixing and cluster composition. Embedding distance and cluster labels are exploratory; donors remain the replication unit for comparative inference.'),
        (['heldout','residual','intervals','confirmation'], 'Validation and uncertainty',
         'A promising training result needs to be assessed against measured outcomes and the uncertainty in the comparison.',
         'Inspect the reference line or interval, sample denominator and endpoint definition. These panels display the supplied validation outputs.'),
        (['leading_edge','shared_gene_support'], 'Genes behind pathway support',
         'Two assays can support the same pathway using entirely different genes.',
         'Check which genes contribute and whether support is shared across layers before interpreting pathway agreement as a common mechanism.'),
        (['pathway_overlap'], 'Redundancy between pathway annotations',
         'Closely overlapping gene sets can make one biological pattern appear to be many separate findings.',
         'Edges indicate shared annotation genes at the declared overlap threshold, not regulatory interactions. Use this view to organize related results.'),
        (['enrichment_curve'], 'Where enrichment is concentrated',
         'The running score shows which portion of the full gene ranking drives a gene-set result.',
         'Read the score peak alongside hit positions and the signed ranking. Enrichment direction refers to that ranking, not proven pathway activation.'),
        (['clock','score_'], 'Biomarker outputs and coverage',
         'Different biomarkers measure different constructs and may require different input features.',
         'Check model units and input coverage before comparing scores. Age, pace and other scores should not be interpreted on a common numeric scale.'),
        (['regulatory','prediction_output','structure','protein_positions','annotation','vep_'], 'Annotation and scientific-model context',
         'Variant and model annotations provide extra context for research prioritization when their coverage and score definitions are visible.',
         'Check matched records, model identity and score units. Different model outputs keep separate scales; annotation is distinct from experimental confirmation.'),
        (['feature_values','feature_heatmap'], 'Patterns in measured features',
         'Individual measurements reveal variation that is hidden by a summary statistic.',
         'Inspect variation within groups as well as between them. Display-selected features and row-standardized colors do not introduce new statistical tests.'),
        (['metrics','receipt','scope'], 'Analysis scope and recorded quantities',
         'Reliable interpretation depends on knowing what ran and which denominators and units were recorded.',
         'Use these recorded quantities to trace the analysis. Values with different units remain separate rather than sharing an arbitrary numeric scale.'),
    ]
    for tokens,title,why,read in rules:
        if any(name.endswith(token) if token in ['_ma','_qq'] else token in name for token in tokens):return dict(title=title,why=why,read=read)
    return dict(title=readable(name).capitalize(),
                why='Measurement quality and distribution determine which comparisons the data can support.',
                read='Inspect the units, counts and unusual values alongside the source QC flags. This display summarizes the recorded measurements without changing the analysis.')


configure()
