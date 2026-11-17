# KilnBio methods and interpretation

KilnBio exposes reusable methods for deciding which measurements and analysis choices change a research conclusion. This guide explains which function to use, its inputs, the statistical or logical basis, and the limits of its output. The functions require an explicit study design; the software does not choose a biological claim from the data.

The `decision-methods` Python API and command are new development interfaces extracted from the manuscript's analysis procedures. They are available in the source checkout containing this guide, not in the unchanged `v0.2.0-research.1`–`.3` archives. Those archives and the manuscript companion resources remain the provenance for the published analyses. New reusable interfaces have numerical and integration checks; they are not an additional biological replication or a demonstrated improvement in novice performance.

## Choose a function

| Function | Why use it | Required assumptions | Interpret the result |
|---|---|---|---|
| `rank_assays` | Prioritize a fixed common feature family | Comparable contrasts; prespecified eligible features; equal assay weights | Signed bladder or unsigned GTEx priorities; scores are not P values |
| `compare_nominations` | Separate patient-set and added-assay changes | Same feature universe and fixed nomination budget across the requested path | Entered/exited targets and rank changes; descriptive rather than causal |
| `validation_yield` | Evaluate frozen follow-up nominations | Independent validation; frozen directions; q values from a declared full family | Confirmed count per original budget; untestable slots remain |
| `donor_partition` | Keep all observations of each donor in one partition | Prespecify salt/folds before outcomes; verify balance and adequate designs | Deterministic allocation; does not guarantee independent related donors |
| `promoter_summary` | Build mean/median or leave-one-probe-out promoter measurements | Fixed assembly-compatible annotation; matched assay scale; no outcome-driven mapping | M-value summaries and explicit untestable states; no transcriptional mechanism claim |
| `fit_hc3` | Estimate an adjusted contrast in continuous measurements | Independent rows; full-rank prespecified design; adequate residual df; no unit leverage | OLS effects with HC3 SE and t reference; not a count model or confounding remedy |
| `directional_conjunction` | Test both or all frozen assay directions in held-out data | Directions fixed independently; full nomination union; valid t reference | Maximum directional P followed by BH/BY; all-assay association rather than mediation |
| `promoter_sensitivity` | Measure nomination changes from mean versus median summaries | Same donors, RNA statistics, gene family and budget; summary is the changed input | Replacements and score agreement; no independent replication |
| `probe_coherence` | Inspect whether probe signs agree with a promoter direction | Effects from the same prespecified contrast and scale; report missing probes | Descriptive signs; a high fraction does not establish significance |
| `element_methylation` | Summarize shared covered CpG dyads within an element | Identical dyad order, compatible chemistry, assembly and replicate units | Equal replicate average of depth-weighted fractions; unmeasured is not zero |

Existing reusable interfaces complement these ten functions:

| Interface | Why use it | Basis and limits |
|---|---|---|
| `multiomic-meta` | Require evidence in at least a declared number of assay layers | Directional Bonferroni partial conjunction, two-direction correction, BY adjustment; valid tests and fixed maps are essential. This is not pooled-effect meta-analysis. |
| `candidate-evidence` | Maintain complete validation records across genes, assays and cohorts | Frozen direction plus adjusted evidence; unmeasured endpoints remain explicit. Cohort labels do not prove independence. |
| `assay-audit pair` | Establish a unique specimen-level join | One donor/tissue per specimen and one measurement per assay in one qualifying specimen. Identifier checks do not authenticate material physically. |
| `assay-audit reach` and `policy` | Separate coverage from prediction performance | Availability counts and paired error differences on identical donors; these functions evaluate predictions rather than train them. |
| `method-agreement` | Inspect analysis-method sensitivity | Same-direction significant calls over planned methods, with conflicts/missing tests explicit; no independence assumption between methods is introduced. |
| R `adjust_assay` and `project_control_adjustment` | Fit and apply an exploratory control-feature projection | Protected biology plus nuisance control PC; apply training centers/loadings/coefficients to new samples. Controls and design separation govern applicability. |

## Run a method

Install the existing scientific Python dependencies in a dedicated environment. Use the source checkout's license and access terms, which are unchanged by this documentation.

```bash
python -m pip install -r environments/multiomic-requirements.txt
export KILNBIO_PYTHON="$(command -v python)"
bash run.sh decision-methods catalog
bash run.sh decision-methods run examples/decision_rank.json --out /tmp/kilnbio-rank-example.json
bash run.sh study examples/decision_study.json /tmp/kilnbio-decision-study
```

Use new output locations. The example is synthetic and has no biological interpretation. After installation, `kilnbio decision-methods ...` uses the same dispatcher. The study example declares its input, step, dependency list and expected output; its result is `rank/result.json` inside the study run directory. The older `workflows/analysis.nf` mode list is a separate interface and has not been extended for these functions.

On Windows, individual functions can be run directly with Python:

```powershell
python downstream/decision_methods.py catalog
python downstream/decision_methods.py run examples/decision_rank.json --out result.json
```

Each JSON request has `schema_version: 1`, a `method` from the catalog, and an `arguments` object using that function's parameter names. Each result records the function's assumptions, the input SHA256 and the implementation SHA256. Missing numerical outputs become JSON `null`. The command never overwrites an existing output file. JSON is useful for small examples and declared workflow steps; use the Python functions with NumPy arrays for large matrices.

```python
from downstream.decision_methods import rank_assays, compare_nominations

# All arrays use precisely this ordered, prespecified gene family.
features = ["GENE_A", "GENE_B", "GENE_C", "GENE_D"]
ranked = rank_assays(features,
    {"RNA": [4, -3, 2, 1], "protein": [3, 4, 1, -2]}, k=2, mode="signed")
print(ranked["arms"]["combined"]["nominations"])

# Supply full rankings for A = all RNA, B = matched RNA, C = matched joint.
changes = compare_nominations(
    {"A": ["GENE_B", "GENE_A", "GENE_C", "GENE_D"],
     "B": ranked["arms"]["RNA"]["ranking"],
     "C": ranked["arms"]["combined"]["ranking"]}, k=2, path=["A", "B", "C"])
```

Functions return dictionaries with explicit identifiers and named outputs. Array arguments must already be aligned; the functions cannot infer that a covariate row or CpG column represents the intended biological unit. Only donor identifiers may repeat in `donor_partition`; its purpose is to keep those repeated observations together.

## Signatures and input contracts

```python
rank_assays(features, statistics, k, mode)
compare_nominations(rankings, k, path)
validation_yield(nominations, directions, effects, adjusted_p, testable, alpha=0.05)
donor_partition(donors, salt, folds=5, validation_folds=(0, 1))
promoter_summary(probe_ids, values, feature_map, scale,
                 method="mean", min_probes=2, drop_probe=None)
fit_hc3(features, sample_ids, response, design, contrast, min_residual_df=10)
directional_conjunction(features, t_statistics, degrees_freedom,
                        directions, eligibility, adjustment="BH")
promoter_sensitivity(features, rna_t, mean_t, median_t, k=25)
probe_coherence(probe_ids, effects, summary_direction)
element_methylation(methylated, depth, min_depth=5, min_cpgs=3)
```

- `rank_assays`: `statistics` maps assay names to finite vectors. Input filtering is completed before the call. Every assay has equal weight; `combined` is a reserved output name. Signed mode reproduces the bladder rule; unsigned mode reproduces the GTEx rule. Other assay sets are a generalization and need study-specific evaluation.
- `compare_nominations`: each named ranking is a complete, unique permutation of the same feature set. `path` explicitly orders the comparisons. This function reports set/rank changes; to inspect sign reversals, compare the aligned `directions` returned by signed `rank_assays`.
- `validation_yield`: all vectors align to the frozen nomination list. `testable` contains booleans. Untestable effects and adjusted P values must be `None`/NaN; valid rows need finite values and directions ±1. The caller supplies adjusted P values from the complete prespecified validation family; the function never recalculates them on a favorable subset.
- `donor_partition`: a prespecified salt makes assignment deterministic. Repeating a donor repeats its fold. Check tissue/group representation and model estimability after splitting. The GTEx salt/folds are given below; a different salt is a sensitivity analysis rather than an independent cohort.
- `promoter_summary`: `values` has probe rows and sample columns; `feature_map` maps each gene to its fixed probe-ID list. `scale` is `beta` or `M`. Means/medians propagate missing required measurements. Omitting a probe may leave a gene untestable if fewer than two remain. The function does not fetch annotations or decide promoter boundaries.
- `fit_hc3`: `response` has sample rows and feature columns, `design` has the same sample rows and one column per covariate, and `contrast` has one entry per design column. Include an intercept explicitly. It returns effect, SE, t, unadjusted two-sided P, residual df, leverage audit and testability. Repeated measurements of the same donor are not independent rows merely because sample IDs differ; use an appropriate repeated-measure model instead.
- `directional_conjunction`: dictionaries use identical assay keys and feature order. `degrees_freedom` supplies positive scalar residual df per assay; `directions` contains frozen ±1 values; `eligibility` contains booleans. Ineligible assay tests become P=1 and remain in the family. Input t statistics are not arbitrary z scores. BH needs its usual dependence conditions; BY is available for arbitrary dependence among valid tests.
- `promoter_sensitivity`: RNA statistics and features are supplied once, so the API holds them fixed. Mean and median methylation statistics must come from the same donors, covariates and fit procedure. The software cannot establish that fact from the vectors alone.
- `probe_coherence`: finite probe effects count as tested; missing and zero effects have separate counts. `fraction_of_tested` excludes missing values, whereas `all_declared_same_direction` is false if any declared probe is missing, zero or opposed.
- `element_methylation`: arrays have replicate rows and matching CpG-dyad columns. Values are compatible methylated-read equivalents and total depth, not array beta values or modification probabilities. At least two replicates are required. It returns eligible columns, common-CpG count, per-replicate fractions and the mean, or an explicit untestable result.

## Equations and rationale

The numbering below matches the expanded manuscript. Each procedure retains its original scope: the glioma discovery statistics were supplied by the published directional P-value merging study; the new partial-conjunction implementation did not generate those nominations. Fixed/random-effects pooling and heterogeneity estimation across studies are not implemented by `multiomic-meta`.

### Validation at a fixed follow-up budget

For a frozen list Sₐ(k) of k nominations from strategy a, let dᵍ be the nominated direction, Δᵍ the validation effect and qᵍ its adjusted P value in the declared validation family. Define an endpoint-specific confirmation indicator and validation yield as

Equation 1:

$$
C_g=\mathbf{1}\{\mathrm{testable}_g,\ d_g\Delta_g>0,\ q_g\leq\alpha\},\qquad Y_a(k)=\sum_{g\in S_a(k)}C_g.
$$

Untestable nominations contribute zero confirmations and remain in the denominator k. Thus Yₐ(k)/k describes the return on a fixed follow-up budget, rather than the success rate after dropping unmeasurable targets. Contrasts between strategies use the same validation cohort and multiplicity family. The test and family are specific to each example: the glioma comparison uses the directional hypotheses defined above, whereas the lung comparison uses its externally measurable gene family. Resampling keeps the nominated lists fixed and preserves pairing between strategies. These yields evaluate held-out evidence and are not probabilities of causal or clinical validity.

The validation_yield function implements this fixed-budget accounting from frozen directions, validation effects, adjusted P values and explicit testability flags. The existing candidate-evidence command maintains the richer gene-by-assay-by-cohort record and requires every nominated endpoint, including unmeasured ones. This prevents a favorable yield being created by silently excluding unsuccessful or unavailable follow-up measurements.

### Covered CpG summaries for regulatory elements

For an eligible element e and replicate r, let cₑᵣⱼ be methylated-read equivalents and nₑᵣⱼ be total depth at common CpG dyad j. Released percentages multiplied by depth supply the methylated-read equivalents. With Jₑ denoting dyads passing the depth criterion in both replicates, the element summary is

Equation 2:

$$
M_e=\frac{1}{2}\sum_{r=1}^{2}\frac{\sum_{j\in J_e}c_{erj}}{\sum_{j\in J_e}n_{erj}},\qquad |J_e|\geq3.
$$

This averages replicate fractions after depth weighting within each replicate. Model comparisons distinguish the original eligible link set from the subset with a measurable methylation summary. Within that common subset, each model selects 25 links using its own held-out prediction scores, and yield is the number of released screen-positive labels among those links. The count of screen-positive links excluded by methylation coverage is reported separately. It cannot be interpreted as a loss caused by methylation biology. Standardization and logistic-regression fitting use training pairs outside the test chromosome; paired target-gene resampling quantifies uncertainty conditional on the fitted predictions.

The reusable element_methylation function implements the covered-dyad summary. Users supply CpG dyads already aligned across compatible replicates and declare the coverage thresholds. Depth weighting estimates each replicate’s measured methylation fraction; equal weighting across replicate fractions preserves the replicate as the averaging unit. Assembly, chemistry and dyad alignment must be established before calling the function.

### Signed ranks and specimen-selection comparisons

For the bladder recipe, A denotes RNA in all eligible tumors, B denotes RNA in matched tumors, and P denotes protein in those same matched tumors. All scores use the intersection U of genes eligible in A, B and P, with G = |U|. For gene g and arm a, the moderated statistic tᵍₐ is converted to an ascending average rank of its absolute value; the sign is then restored. The combined arm C is

Equation 3:

$$
s_{ga}=\operatorname{sgn}(t_{ga})\frac{\operatorname{rank}_{U}(|t_{ga}|)}{G+1},\qquad s_{gC}=\frac{s_{gB}+s_{gP}}{2}.
$$

Each arm orders genes by decreasing absolute score and uses gene identifiers to break exact score ties. The first 15 genes form its nomination set. The sign of the combined score supplies its direction; discordant assay signs can cancel. These scores are descriptive priorities, not P values, pooled effect estimates or probabilities of validation. The denominator G + 1 is the implemented bladder convention; it does not affect ordering within the shared family.

The A-to-B comparison measures the observed change after restricting the patient set; B-to-C measures the observed change after adding protein on that fixed set. If rᵍₐ is the final rank, their rank differences satisfy

Equation 4:

$$
r_{gC}-r_{gA}=(r_{gB}-r_{gA})+(r_{gC}-r_{gB}).
$$

This identity organizes the comparison and does not identify causal effects of selection or measurement. For a transition a to b, entering genes are Sᵦ(k) \ Sₐ(k), exiting genes are Sₐ(k) \ Sᵦ(k), and direction changes among retained genes are recorded separately. Counts of replacements need not add across transitions because a gene can leave and re-enter. The common feature family prevents a changed gene denominator from being mistaken for a changed ranking rule.

The rank_assays function exposes signed and unsigned modes with the denominators and tie rules specified here and in the GTEx section. Ranking places assay statistics with different measurement units on a common ordinal scale, while equal weighting makes the contribution of each supplied assay explicit. It sacrifices effect-size calibration and cannot correct an invalid assay model. The compare_nominations function applies the same budget to successive rankings and returns entries, exits and per-gene rank changes; users supply which transition represents specimen restriction or added evidence.

### Control projection and separate agreement objectives

The control projection can be written as Y = XB + zγᵀ + E, where rows of Y are training measurements, X contains the intercept and reference-material indicators, and z is the centered first principal-component score of standardized control features. Least squares estimates B and γ. For a new measurement y with control subset C, the fitted transformation is

Equation 5:

$$
z_{\mathrm{new}}=\left(\frac{y_C-\mu_C}{\sigma_C}\right)^{\!T}v-\bar z,\qquad y^{*}=y-z_{\mathrm{new}}\widehat\gamma.
$$

The training control means μ, standard deviations σ, loading vector v, score center and coefficients are frozen. No held-out material label is required. The implementation declines adjustment when the protected and nuisance terms are rank deficient, fewer than three residual degrees of freedom remain, or the condition number of the standardized non-intercept design exceeds 10. These are operational safeguards for this implementation, not evidence that every remaining adjustment preserves biology.

Let D(u,v) be the root mean squared difference across the scored non-control features, and let μₘ,₋ₗ denote the equal-laboratory reference for material m excluding held-out laboratory l. For each evaluation batch and replicate, the two objectives are

Equation 6:

$$
R_{\mathrm{tech}}=\frac{1}{3}\sum_m D(y_m^{*},\mu_{m,-l}^{*}),
$$

Equation 7:

$$
R_{\mathrm{bio}}=\frac{1}{3}\sum_{m<n}D(y_m^{*}-y_n^{*},\mu_{m,-l}^{\mathrm{raw}}-\mu_{n,-l}^{\mathrm{raw}}).
$$

The three materials generate three pairwise contrasts. Technical agreement uses references on the same processing scale; biological contrast agreement uses the unadjusted reference differences. Their separate reporting makes a reduction in laboratory discrepancy distinguishable from preservation of material differences. Subtracting the same concurrent D6 vector from both materials cancels algebraically in a within-batch contrast.

The existing R functions adjust_assay and project_control_adjustment provide fitting and application, respectively. The rationale is that control-feature covariance can supply a measured nuisance direction, while including material identity in training protects the specified biological contrast during coefficient estimation (Risso et al. 2014). Applicability depends on informative controls and separable design terms. Their use here remains a custom exploratory projection; it does not establish a universally transportable batch correction.

### Assay availability and prediction policy

Write bᵢ and aᵢ for eligibility of donor i for the baseline and additional assay. With baseline prediction fᵦ and two-assay prediction fᵦₐ, the frozen policy is

Equation 8:

$$
\widehat y_i^{\pi}=\begin{cases}f_{ba}(x_i),&b_i=1,\ a_i=1,\\f_b(x_i),&b_i=1,\ a_i=0,\\\mathrm{abstain},&b_i=0.\end{cases}
$$

The complete two-assay reach is the number with both eligibility indicators equal to one; policy reach is the number with baseline eligibility. For the same set L of scored donors with observed outcomes, the paired error gain is

Equation 9:

$$
\Delta_{\mathrm{MAE}}=\frac{1}{|L|}\sum_{i\in L}\left(|y_i-f_b(x_i)|-|y_i-\widehat y_i^{\pi}|\right).
$$

Positive gain favors the policy. Missing outcomes and abstentions are counted separately. The fallback donors contribute zero paired gain because both policies use the same baseline prediction. The reported resampling interval is conditional on the frozen fits and predictions; it does not include uncertainty from refitting a new development cohort.

The existing assay-audit reach and policy operations implement these comparisons. The policy function evaluates supplied frozen predictions; it does not train a replacement model or impute a missing assay. Comparing errors on the same scored donors isolates the performance of the stated policy on that population, while separate reach counts show who can receive a prediction.

The pairing rule first requires every nonmissing specimen identifier to map to one donor and one tissue. Within a donor and requested tissue, a pair is accepted only if exactly one specimen contains both requested assays and that specimen contains exactly one measurement per assay. Multiple qualifying specimens or files remain ambiguous, while missing specimen identifiers remain unknown. This rule audits the supplied identifiers; it does not independently authenticate the physical origin of a specimen.

The assay-audit pair operation implements this eligibility rule. Requiring a unique specimen-level join prevents shared donor identity alone from being interpreted as evidence that two assays measured the same material.

### Agreement among analysis methods

For a feature g, let M be the number of prespecified methods and Nᵍ⁺ and Nᵍ⁻ the numbers making positive and negative calls, respectively, at their within-method adjusted threshold. The reported agreement score is

Equation 10:

$$
A_g=\frac{\max(N_g^{+},N_g^{-})}{M}.
$$

Missing tests remain in the planned-method denominator. A second descriptive fraction divides by the number actually tested, when nonzero. Any positive and negative calls together produce an opposite-calls flag, including a tie; the score alone does not resolve that conflict. Agreement among methods fitted to the same observations is not independent replication, a combined P value or a calibrated confidence measure.

The method-agreement command exposes this calculation through the agreement function. Its use is to identify method-sensitive conclusions and conflicting directions under a common study, contrast, feature family, assembly and measurement channel. Because the methods share observations, their agreement is a robustness description rather than an increase in the number of independent biological replicates.

### Promoter summaries and matched GTEx validation

The GTEx ranking differs from the signed bladder rule. For assay a in a source tissue, it uses the unsigned percentile uᵍₐ of the absolute age statistic within the shared eligible gene family:

Equation 11:

$$
u_{ga}=\frac{\operatorname{rank}_{U}(|t_{ga}|)}{G},\qquad u_{g,\mathrm{combined}}=\frac{u_{g,\mathrm{RNA}}+u_{g,\mathrm{meth}}}{2}.
$$

Average ranks handle tied statistics; stable sorting of the fixed gene order resolves equal nomination scores. The two assay directions are retained separately. Thus nomination does not require RNA and methylation to change in the same direction. For CpG beta value β, the M-value is log₂[β/(1 − β)] after clipping β to [10⁻⁶, 1 − 10⁻⁶]. The gene summary takes the mean across its fixed annotated CpGs; the controlled alternative takes their median. The original RNA model used log₂ counts per million with a 0.5-count pseudocount and the corresponding library-size correction; TMM and voom were evaluated as sensitivities, as described below.

The donor_partition function assigns the first 16 hexadecimal digits of SHA256(salt + donor identifier), modulo the declared number of folds. The original GTEx salt was gtex-20261104: with five folds, of which 0 and 1 were assigned to validation. All occurrences of a donor therefore share a partition. Hashing provides a reproducible allocation independent of row order; it does not ensure covariate balance or independence between related donors.

The promoter_summary function implements a fixed feature-to-gene map on the beta or M-value scale, with mean, median and optional omission of a named probe. Missing required probes or values produce an untestable gene–sample summary rather than changing the mapped set. Mean summaries represent the average measured promoter signal; medians assess sensitivity to an extreme probe. Neither summary implies that all promoter CpGs have equivalent regulatory effects.

The fit_hc3 function estimates an explicitly supplied linear contrast in continuous measurements. For gene g, design X, residuals eᵢᵍ, leverages hᵢ and A = (XᵀX)⁻¹, its covariance estimate is

Equation 12:

$$
\widehat V_g=A X^{T}\operatorname{diag}\left\{\frac{e_{ig}^{2}}{(1-h_i)^2}\right\}X A.
$$

A contrast c has estimated variance cᵀV̂ᵍc, and the implementation uses a t reference with n − rank(X) residual degrees of freedom. The HC3 leverage correction addresses heteroskedastic residual variance (MacKinnon and White 1985). It does not remove unmeasured confounding, accommodate dependent donor rows, or make ordinary least squares appropriate for untransformed counts. The function requires a full-rank design, no unit-leverage observation and adequate residual degrees of freedom; zero standard-error contrasts are flagged untestable.

For frozen discovery direction dᵍₐ and validation t statistic with residual degrees of freedom ν, directional validation uses the upper-tail probability below. The joint test requires support in both assay-specific directions:

Equation 13:

$$
p_{ga}^{\mathrm{dir}}=1-F_{t_\nu}(d_{ga}t_{ga}^{\mathrm{val}}),\qquad p_g^{\mathrm{both}}=\max(p_{g,\mathrm{RNA}}^{\mathrm{dir}},p_{g,\mathrm{meth}}^{\mathrm{dir}}).
$$

The maximum is a conjunction test; it does not test whether one assay mediates the other. No extra search over direction occurs in validation because directions were frozen in discovery. A gene failing validation RNA expression eligibility receives RNA P = 1 and remains in the nomination-union adjustment family. Benjamini–Yekutieli adjustment supplies the additional dependence sensitivity.

At budget k = 25, the controlled summary comparison counts changed nominations as

Equation 14:

$$
D_k=k-|S_{\mathrm{mean}}(k)\cap S_{\mathrm{median}}(k)|.
$$

This is a sensitivity of the decision to a promoter summary, not a count of newly validated targets. Probe-level effect signs and leave-one-probe-out estimates distinguish an aggregate association from one dependent on a single CpG. Nomination overlap, directional coherence and adjusted validation evidence remain separate outputs; no composite robustness score was fitted.

The directional_conjunction function evaluates the frozen directions across the full nominated family; its maximum P value implements the requirement that every supplied assay support its own expected direction. The promoter_sensitivity function compares the mean and median rankings with the RNA statistics and gene family held fixed. The probe_coherence function reports positive, opposite, zero and missing probe effects relative to the summary direction. For a leave-one-probe-out check, promoter_summary is called with one omitted probe and the same design is refitted using fit_hc3. Summaries with fewer than two remaining probes are untestable. Repeated donor partitions or omitted probes reuse observations, so these operations assess sensitivity rather than create new independent validation cohorts.

### Directional evidence across assay layers

For a feature f mapped to gene g in layer l, let pᶠˡ be a valid two-sided P value from a symmetric test and eᶠˡ its effect multiplied by the prespecified layer orientation. The positive directional tail is pᶠˡ/2 when eᶠˡ ≥ 0 and 1 − pᶠˡ/2 otherwise; the negative tail reverses the inequality. A zero effect is accepted only with P = 1, giving both tails 0.5. If Kᵍˡ features map to the gene in that layer, directional feature aggregation is

Equation 15:

$$
P_{gl}^{d}=\min\left(1,K_{gl}\min_{f\mapsto g}p_{fl}^{d}\right),\qquad d\in\{+,-\}.
$$

Every mapped feature must have a finite effect and P value for that gene–layer result to be testable. An incomplete mapped feature set is not silently reduced. For a gene with mᵍ testable layers and a required minimum r ≥ 2, let Pᵍ₍ᵣ₎ᵈ be the rth smallest directional layer P value. The implemented Bonferroni partial-conjunction test and correction for choosing either direction are

Equation 16:

$$
P_g^{\mathrm{PC}}=\min\left\{1,\ 2(m_g-r+1)\min\left(P_{g,(r)}^{+},P_{g,(r)}^{-}\right)\right\},\qquad m_g\geq r.
$$

Genes with fewer than r testable layers are untestable. For multiplicity adjustment, all G declared genes remain in the family, with untestable results temporarily assigned P = 1. For ordered P values, Benjamini–Yekutieli adjustment is

Equation 17:

$$
q_{(i)}=\min\left\{1,\min_{j\geq i}\frac{G H_G P_{(j)}^{\mathrm{PC}}}{j}\right\},\qquad H_G=\sum_{h=1}^{G}\frac{1}{h}.
$$

Untestable q values are then restored to missing. A testable gene is labeled convergent when q ≤ α. This supports at least r associated layers in the selected orientation under the stated assumptions, including prespecified assay availability and valid input tests. It does not imply that every layer is individually significant or that the layers represent independent cohorts. Bonferroni and Benjamini–Yekutieli address dependence among valid tests; they do not repair invalid input models or outcome-driven feature selection (Benjamini and Heller 2008; Benjamini and Yekutieli 2001).

The module also returns an exploratory Cauchy-combination global-null P value using positive, prespecified layer weights (Liu and Xie 2020). That output is not used to label convergence or to select the reported glioma nominations. The term multiomic meta-analysis here denotes combination of evidence across assay layers; the module does not estimate a fixed-effect or random-effects pooled association across studies, or a between-study heterogeneity parameter.

The existing multiomic-meta command exposes the feature aggregation, directional_pc and adjustment operations as one validated input-to-output procedure. Partial conjunction is useful when the claim requires at least r assay layers instead of any one significant layer. Its order-statistic factor bounds the partial-conjunction null using established Bonferroni reasoning; the factor of two allows selection between the two orientations. Users must choose the required layers, feature mappings and orientations before inspecting the joint result.

### Operational evidence rules

Operational checks keep annotation evidence distinct from observed molecular effects. Missing QC or callability is retained as unknown rather than converted to a passing or negative state; conflicting identifiers, incompatible units and invalid model designs stop the affected analysis. A computational prediction or database assertion does not add a numerical validation bonus to the nomination score. The repository methodology guide links these rules and the equations above to their implementations and distinguishes implemented interfaces from empirically evaluated comparisons.

### Exploratory Cauchy combination

For the optional global-null diagnostic, each finite layer result is made nondirectional as $P_{gl} = \min(1,2\min(P_{gl}^{+},P_{gl}^{-}))$. With positive prespecified weights normalized over finite layers,

$$
T_g=\frac{\sum_l w_l\tan\{\pi(1/2-P_{gl})\}}{\sum_l w_l},\qquad
p_g^{\mathrm{ACAT}}=1/2-\arctan(T_g)/\pi.
$$

The implementation clips numerical endpoints and uses tail approximations for extremely small P values or large statistics. It returns missing when no layer is observed. The Cauchy tail approximation is established methodology (Liu and Xie 2020), and this diagnostic tests a global null rather than the partial-conjunction claim. Its weights do not alter the primary partial-conjunction selection. An assay's inclusion or missingness cannot be selected after inspecting favorable evidence without invalidating the intended inference.

## Reusing the existing R correction functions

```r
source("recipes/glioma-quartet/analysis/quartet/adjust.R")
# Training matrix: rows = independent measurements; columns = fixed features.
# biology is the protected training contrast; controls are prespecified/selected
# using training data only. Each input must satisfy adjust_assay's design checks.
fit <- adjust_assay(training_matrix, biology, controls=control_ids)
if (!is.null(fit$control_projection)) {
    corrected_test <- project_control_adjustment(fit, test_matrix)
}
```

Inspect `fit$status`; a declined fit is an analysis outcome, not permission to remove safeguards. The frozen projection is available only for a successful control-only nuisance fit. Batch covariates and categorical/continuous biology have additional overlap checks in `adjust_assay`. The reported Quartet recipe selects controls using training material means and residual variability and excludes those controls from scoring. A new study must define controls and biological contrasts for its own design.

## Connecting the promoter functions

Use `donor_partition` before model selection, then fit each tissue/partition on matched independent donors. Build the gene map and common eligibility family before ranking. Call `promoter_summary` with mean and median on the same mapped CpGs; transpose its gene-by-sample output when supplying sample-by-feature responses to `fit_hc3`. Drop untestable genes according to a rule fixed before comparing summaries and record those exclusions.

Fit RNA and each promoter summary with the same donor design. Pass their discovery t statistics to `rank_assays` or `promoter_sensitivity`. Freeze nomination lists and assay-specific directions before fitting validation effects. Use `directional_conjunction` on the full nomination union. Retain untestable validation hypotheses in that family. For each mapped probe, call `promoter_summary(..., drop_probe=probe_id)`, refit the same design, and record whether the direction and validation evidence survive. Use `probe_coherence` to summarize individual-probe signs separately.

These functions expose reusable calculations; they do not automatically acquire GTEx data, normalize RNA counts, construct a valid design, recover missing covariates, or establish consent/access. The source GTEx analysis and its detailed sensitivity protocol remain in manuscript Online Resources 6 and 7. No participant-level GTEx data are included here.

## Source map and validation

| Method | Implementation or original recipe |
|---|---|
| New Python functions and JSON interface | [decision_methods.py](../downstream/decision_methods.py) |
| Study orchestration | [study.py](../downstream/study.py), [run.sh](../run.sh), [study example](../examples/decision_study.json) |
| Signed ranking and specimen restriction | [bladder rank.R](../recipes/bladder/rank.R), [bladder protocol](../recipes/bladder/PROTOCOL.md) |
| Regulatory element measurement and fixed-budget comparison | [K562 compare.py](../recipes/k562/compare.py) |
| Control fitting and projection | [adjust.R](../recipes/glioma-quartet/analysis/quartet/adjust.R), [held-laboratory recipe](../recipes/glioma-quartet/analysis/quartet/quartet_hidden_replicates.R) |
| Multi-assay partial conjunction and ACAT | [multiomic_meta.py](../downstream/multiomic_meta.py), [shared-model guide](MULTIOMIC_AND_OLINK.md) |
| Specimen joins, reach, fallback predictions | [assay_audit.py](../downstream/assay_audit.py) |
| Endpoint-complete validation records | [candidate_evidence.py](../downstream/candidate_evidence.py) |
| Method support/conflict reporting | [method_agreement.py](../downstream/method_agreement.py) |
| Numerical and boundary checks | [test_decision_methods.py](../tests/test_decision_methods.py) |

The initial extraction was checked against all 129 genes in the archived bladder results: combined ranks matched exactly and the maximum score difference was below 5 × 10⁻¹⁶. HC3 effects, SEs, t statistics and P values were checked against statsmodels using an independently fitted contrast. BY-adjusted conjunction values were checked against statsmodels multiplicity correction. Tests also cover signed cancellation, tie handling, full nomination budgets, untestable promoters, unit leverage, invalid designs, shared-donor folds, missing probe effects, and the distinction between equal replicate averaging and pooled-read averaging. These are computational checks, not new biological validation.

```bash
python -m pip install -r environments/decision-tests-requirements.txt
python -m unittest discover -s tests -p test_decision_methods.py -v
```

The evidence-combination statistics, limma, HC3, control-factor rationale and validation metrics have established antecedents. KilnBio contributes explicit procedures, interfaces, assumptions and records connecting those components to research decisions. Writing an equation does not establish priority for a statistical method. Different rankings are not automatically better rankings; independent validation and the intended use determine the value of a change.

## Statistical sources

- Benjamini and Hochberg (1995), [false discovery rate control](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
- Benjamini and Yekutieli (2001), [FDR under dependence](https://doi.org/10.1214/aos/1013699998).
- Benjamini and Heller (2008), [partial conjunction hypotheses](https://doi.org/10.1111/j.1541-0420.2007.00984.x).
- Liu and Xie (2020), [Cauchy combination](https://doi.org/10.1080/01621459.2018.1554485).
- MacKinnon and White (1985), [heteroskedasticity-consistent covariance estimation](https://doi.org/10.1016/0304-4076(85)90158-7).
- Risso and colleagues (2014), [control-feature factor analysis](https://doi.org/10.1038/nbt.2931).
- Ritchie and colleagues (2015), [limma](https://doi.org/10.1093/nar/gkv007).
- Slobodyanyuk and colleagues (2024), [directional P-value merging](https://www.nature.com/articles/s41467-024-49986-4).
