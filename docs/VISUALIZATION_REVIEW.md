# Visualization review across KilnBio analyses

Reviewed November 18, 2026. This is a broad, structured methods and documentation
review across KilnBio's supported assays and adjacent extensions, with targeted
2024–2026 searches and foundational methods retained. It is not a claim to have
enumerated every visualization ever published or a preregistered systematic
review. PubMed/PMC, primary publisher articles, Bioconductor and official tool
documentation were used. Exact recorded search batches are in
[the search record](../evidence/visualization_search_2026-11-18.json).

## How to read this inventory

**Automatic** means the named KilnBio command now writes the figure. **Existing**
means a figure was already in a KilnBio command or named recipe. **Upstream**
means it is provided by a configured pinned processing workflow, conditional on
its enabled steps. **Extension** means researched and useful, but no tested
automatic KilnBio adapter is claimed here. **New analysis** means fitting another
model or supplying another data type is necessary before the plot can be valid.

The [report guide](VISUALIZATIONS.md) contains working commands, outputs and
interpretation. Installing KilnBio does not install every tool mentioned below.
Use the underlying method citation as well as KilnBio's version when publishing.

## 1. Sequencing and sample quality

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Read quality by base, GC and sequence length; adapter content | Are reads suitable for the intended assay? FASTQ QC metrics | FastQC/MultiQC [1] | Upstream |
| Mapping, duplication, retained reads and library complexity | Which libraries lose useful reads? Alignment metrics; optional complexity model | MultiQC; pinned workflow output docs below | Upstream, configured steps |
| Sample count/depth distributions and sample correlation heatmap | Which libraries differ in depth or global pattern? Normalized/transformed matrix | `de`, MultiQC | Existing `de` QC.pdf |
| PCA with biology, donor and batch annotations | Are major axes associated with design or technical factors? Transformed matrix and metadata | `de`; ComplexHeatmap [2] for annotations | Existing PCA; expanded metadata panels are Extension |
| Missing-measurement heatmap and exclusive coverage patterns | Who/what is excluded when an assay is required? Explicit measurement roster | `multiomic-meta`; UpSet design [3] | Automatic gene coverage patterns; donor display is Extension |
| Before/after normalization distributions | Did scaling remove gross library differences? Raw and normalized values with method labels | MultiQC; assay-specific models | Existing selected QC; general paired comparison is Extension |
| Technical replicate scatter, CV and Bland–Altman plots | How reproducible is measurement on the same specimen? Technical replicates | Assay-specific replicate models | Extension; biological replicates are not substitutes |
| Variance-partition and residual-diagnostic panels | How much variation is associated with donor, batch or covariates? Fitted variance model | Model-derived plots | New analysis; PCA separation alone cannot identify technical noise |

## 2. RNA-seq and differential expression

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Volcano plot | Effect size versus uncertainty; full tested gene set | DESeq2/edgeR/limma outputs | Existing `de` QC.pdf |
| MA/MD plot | Does differential signal vary with abundance? Mean abundance and model effect | DESeq2/edgeR/limma | Extension |
| Dispersion/mean-variance and voom trend | Does the count model fit the data? Fitted model object | Assay model diagnostics | Extension |
| Top-gene heatmap with sample annotations | Which samples drive the signature? Transformed counts, fixed gene choice | ComplexHeatmap [2] | Extension; label row scaling explicitly |
| Per-gene box/violin plus individual or paired points | Are effects consistent across biological samples? Normalized measurements and donor roster | Standard model-compatible plots | Extension |
| Method overlap matrix / UpSet and effect-effect scatter | Which results survive DESeq2, edgeR and limma choices? Full comparable result families | UpSet [3]; method-agreement tables | Extension; overlap is not an FDR procedure |
| Splice-junction sashimi and transcript-usage plots | Does regulation affect isoforms rather than total gene abundance? Junction/transcript evidence | Genome browser/Gviz [4] | New analysis for differential splicing; not inferred from gene DE |

## 3. Methylation arrays, WGBS, RRBS and capture methylation

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| M-bias/read-position methylation | Does conversion/calling vary along reads? Bismark/MethylDackel output | Pinned methylseq docs | Upstream, caller-dependent |
| CpG depth and methylation-fraction distributions | How much of the methylome is callable? Modified/total counts | methylKit [5]; processing reports | Upstream summaries; general downstream display is Extension |
| Conversion/control and CpG/CHG/CHH QC | Are chemistry and controls behaving appropriately? Assay-specific controls | Pinned methylseq docs | Upstream when supplied; non-CpG signal is not a universal conversion control |
| CpG/DMR effect versus evidence; Manhattan and QQ plots | Where are differential sites/regions, and how calibrated are tests? Complete tested table | methylKit [5]; DSS/dmrseq outputs | Extension |
| DMR locus tracks with per-sample fraction and depth | Is the region supported by several covered CpGs and donors? Coordinates/counts/annotation | Gviz [4], pyGenomeTracks [6] | Extension |
| Promoter/enhancer/gene-body annotation summaries | Where do DMRs fall relative to a frozen annotation? Genome-matched intervals | Genomic annotation tools | Extension; use an appropriate measured-region background |
| CpG/DMR-by-sample heatmaps | Is change consistent across samples? M/beta values or covered count estimates | ComplexHeatmap [2] | Extension; avoid equating missing coverage with zero methylation |
| Methylation pathway evidence | Which functions are implicated after accounting for measurement bias? Valid upstream tests | GOmeth/GOregion [7] for arrays; reviewed sequencing methods | Automatic comparison/report from imported valid pathway tests; upstream bias-aware adapter remains Extension |

GOmeth/GOregion explicitly address methylation-array probe/gene biases and
multi-gene annotations [7]. This does not automatically validate their use for
all WGBS/RRBS/capture experiments. A generic nearest-gene mapping followed by
unqualified overrepresentation analysis is not the default recommendation.

## 4. ATAC-seq and ChIP-seq

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Fragment lengths / nucleosomal periodicity | Does ATAC library structure look plausible? Paired-end fragments | Pinned atacseq docs | Upstream |
| TSS enrichment and signal heatmap/metaprofile | Is accessible or bound signal concentrated at relevant regions? BigWig/BAM plus fixed regions | deepTools [8], ChIPseeker [9] | Upstream configured QC; reusable arbitrary-region report is Extension |
| FRiP, enrichment, complexity and replicate QC | Which libraries support peak analysis? Assay-specific metrics | MultiQC [1]; pinned ATAC/ChIP workflows | Upstream, enabled steps |
| Peak width, annotation and TSS-distance distributions | Do narrow/broad calls suit the assay and mark? BED/narrowPeak/broadPeak | ChIPseeker [9] | Extension |
| Replicate/caller overlaps and concordance | How sensitive is the peak set to callers? Harmonized genomic intervals | Intervene [10], UpSet [3] | Existing comparison tables; automatic plots are Extension |
| Differential accessibility/binding volcano, MA and sample heatmaps | Which regions change under the specified design? Consensus count matrix and fit | `de`; upstream workflow | Existing volcano/PCA where count-based DE is used; additional views are Extension |
| Stacked locus tracks across ATAC, ChIP, methylation and RNA | Do measurements support the same regulatory hypothesis? Same genome build and declared scaling | pyGenomeTracks [6], Gviz [4] | Extension |
| Motif enrichment and motif activity heatmaps | Are TF motifs enriched or varying? Sequences or motif-deviation model | Signac [11] and reviewed motif tools | New analysis; accessible motif is not proof of TF binding |
| Footprints and enhancer–gene links | Where is protection or regulatory association supported? Bias-corrected cuts / link evidence | Signac [11]; dedicated models | New analysis; distance alone is not a validated enhancer–gene link |
| Regulatory-region pathway plots | Which processes are linked to regions? Appropriate regulatory domains and region background | GREAT [12], ChIPseeker [9] | Automatic imported pathway comparison; execution of GREAT is Extension |

## 5. Protein abundance, Olink and metabolomics

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| NPX distributions, plate/control effects, detection rates | Are Olink measurements comparable and measurable? Delivered NPX/QC/LOD metadata | OlinkAnalyze official docs [13] | Existing import/audit tables; plots are Extension |
| Protein/peptide missingness and CV plots | Does identification or quantification vary by sample/batch? Quantification and QC metadata | MSstats [14] | Extension; raw MS processing is not implied |
| Protein volcano, sample PCA and abundance heatmaps | Which proteins change? Fitted abundance model and matrix | limma/MSstats [14] | General pathway/evidence reports Automatic; protein-specific QC adapter Extension |
| RNA–protein effect scatter and discordance panels | Is transcription reflected in protein change? Matched eligible features and compatible contrasts | Paired result tables | Automatic directional concordance; effect-size scatter Extension |
| Peptide-to-protein or phosphosite-to-protein views | Which measurements support the protein result? Peptide/site mapping and parent abundance | MSstats family; reviewed PTM workflows | New analysis / Extension |
| Chromatograms, mass-error and peak-identification views | Is the measured molecule correctly identified? Raw MS and feature identification data | MetaboAnalyst 6.0 [15] | New analysis; not part of current KilnBio raw-processing support |
| Metabolite volcano, PCA, heatmap and pathway overlays | Are biochemical changes consistent with other assays? Identified metabolites and reviewed models | MetaboAnalyst [15], PaintOmics [16] | Extension; metabolite IDs and pathway topology require explicit mapping |

## 6. Single-cell RNA, chromatin and multimodal data

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| UMI/gene/mitochondrial QC violin and scatter | Which cells pass prespecified QC? Cell metrics | Scanpy [17] | Existing QC tables; additional QC plots Extension |
| UMAP colored by cluster and sample | What structure is visible? Fitted neighborhood graph/embedding | Scanpy [17] | Existing `scrna` automatic UMAP |
| Marker dotplots, stacked violin and cell-by-gene heatmaps | What supports cell annotation? Expression, marker model and cell labels | Scanpy [17], ComplexHeatmap [2] | Extension |
| Donor-by-cell-type composition and pseudobulk contrasts | Is a difference compositional or within cell type? Donor IDs and cell labels | Existing pseudobulk route | Extension figures; cells cannot replace independent donors |
| Paired RNA/ATAC embeddings and coverage tracks | How do modalities align at cell or population level? Truly paired/mapped modalities | Signac [11], MOFA+ [18] | New analysis beyond current 10x RNA route |
| Trajectories and RNA velocity streamlines | Is a temporal/dynamical hypothesis supported? Fitted trajectory; spliced/unspliced data for velocity | scVelo [19] | New analysis; an attractive UMAP is not evidence of lineage |
| Ligand–receptor dotplots and communication networks | Which cell populations have compatible signaling evidence? Appropriate expression/spatial model | LIANA+ (2024) [20] | New analysis; inferred communication remains a model output |
| Spatial overlays and tissue-linked multimodal tracks | Where do molecular patterns occur? Spatial coordinates and image registration | SGS (2025) [21] | New analysis; spatial support is not currently claimed |

## 7. Genomics, exomes, CH candidates, long reads and modified bases

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Capture coverage and callable-target summaries | Which targets could the experiment evaluate? Coverage/BED/callability | Existing capture/exome QC interfaces | Existing tables; coverage dashboards Extension |
| VAF/depth/quality distributions | Which variants have convincing read evidence? Variant calls and supporting metrics | maftools [22] / caller QC | Extension; low VAF alone does not establish CH |
| Oncoprint, mutation burden and co-occurrence maps | What is the cohort mutation landscape? Somatic MAF plus callable denominator | maftools [22] | Extension; germline and somatic interpretations differ |
| Protein-position lollipop and domain tracks | Where do variants affect protein sequence? Transcript-specific annotations | maftools [22] | Extension |
| CNV segments, structural-variant arcs and read evidence | Which genomic alterations are supported? Suitable caller output and alignments | Genome browser/Gviz [4] | New analysis for unsupported callers; annotation is not variant calling |
| Variant–RNA–chromatin–protein evidence panels | What functional context surrounds an alteration? Reviewed links and harmonized coordinates | Gviz [4], pyGenomeTracks [6], pathway reports | Extension locus view; shared pathway evidence Automatic |
| Nanopore length/quality/yield plots | Are long reads suitable? FASTQ/BAM/run metadata | NanoPack [23] | Extension visualization adapter |
| 5mC/5hmC fraction/depth and haplotype tracks | Which modification is measured, and with what support? Chemistry-aware count/BAM tags | Existing modification import; locus tools | Extension; do not label joint modC as resolved 5mC/5hmC |
| Biomodal chemistry/control and genomic coverage views | Do duet +modC/evoC outputs support the requested quantity? Provider outputs and chemistry labels | Existing import/QC schemas | Extension figures; different chemistries remain distinct |

## 8. Integrated evidence and validation

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Gene-by-assay signed evidence heatmap | Which molecular layers agree or conflict? Reviewed gene-level directional tests | `multiomic-meta`, heatmap principles [2] | Automatic |
| Joint ranking with measured-layer count | Which candidates have adjusted joint evidence? Partial conjunction + BY | `multiomic-meta` [24,25] | Automatic |
| Pairwise concordance and coverage patterns | How much evidence is shared versus absent? Complete declared family | `multiomic-meta`; set principles [3] | Automatic |
| Leave-one-assay-out rank sensitivity | Which rankings depend on an assay? Fixed r and common testable family | `multiomic-meta` | Automatic when estimable |
| Matched-versus-available cohort and eligibility diagrams | Is an apparent gain caused by changing patients/features? Explicit rosters | Existing assay-audit/compare-assays | Existing tables; general automatic diagram Extension |
| Shared and assay-specific latent factor loadings/variance | What structure is shared or modality-specific? Fitted joint model | MOFA+ [18], mixOmics [26] | New analysis; no causal technical/biological decomposition from a plot |
| Selection stability across resampling or methods | Which targets are reproducible? Proper repeated fitting | Model-specific resampling | Extension; correlated callers are not independent votes |
| Held-out yield curves and uncertainty intervals | Does integration improve follow-up performance? Independent prespecified endpoints | Existing glioma/Quartet recipe plots | Existing recipe; general report adapter Extension |
| ROC/precision–recall, calibration and error by subgroup | Does a biomarker predict reliably on new subjects? Held-out predictions and labels | Model-validation diagnostics | New analysis/Extension; train/test leakage invalidates the interpretation |
| Forest and funnel plots for study meta-analysis | Are compatible effects consistent across cohorts? Effect estimates, SEs and a study-level model | Conventional fixed/random-effects analysis | New analysis; not the current `multiomic-meta` test |

## 9. Aging clocks and reprogramming biomarkers

| Visualization | Question / required data | Tool or source | KilnBio status |
|---|---|---|---|
| Model-feature coverage and missing-input heatmap | Can this assay score the requested clock? Model CpGs/features and measured roster | Existing clock/coverage interfaces | Existing tables; automatic figures Extension |
| Predicted versus chronological age with identity line | Is a clock calibrated in this population? Scores and known ages | Model validation | Extension; report tissue/range and held-out error |
| Age-acceleration residuals and covariate plots | Which differences remain after the declared adjustment? Explicit fitted model | Model validation | New analysis / Extension |
| Clock-to-clock agreement and subgroup errors | Are clocks measuring consistent quantities? Compatible scores and cohorts | Comparative validation | Extension; units differ for pace versus age |
| Paired reprogramming trajectories | Does a score change within donors over treatment/time? Repeated specimens and model | Existing EpiDiff/lineage interfaces | Extension; score changes do not by themselves establish rejuvenation |

## 10. 16S and microbiome: a separate future route

**16S sequencing is not implemented in KilnBio.** RNA-seq or incidental microbial
QC is not a substitute for an amplicon workflow. This is included in the review
because it was requested, not listed as current platform support.

| Visualization | Required inputs / scientific purpose | Candidate tools | Status |
|---|---|---|---|
| Read filtering/denoising retention and quality profiles | Amplicon reads, primers and DADA2/denoising metrics | QIIME 2 [27] | New pipeline |
| ASV abundance/prevalence and taxonomy composition | Valid feature table and reference taxonomy | phyloseq [28], QIIME 2 | New pipeline |
| Alpha diversity and sequencing-depth sensitivity | Sample metadata and declared diversity estimator | QIIME 2/phyloseq | New pipeline; diversity estimation and differential abundance need different decisions |
| Beta-diversity PCoA/NMDS with dispersion checks | Appropriate distances/transforms and design | QIIME 2/phyloseq | New pipeline; apparent separation is not a group test |
| Differential taxon effect/error and composition-aware heatmaps | Valid compositional model, batch and repeated-measure handling | Reviewed microbiome differential method | New analysis; do not feed relative abundances to count DE unchanged |
| Host–microbiome pathway or association panels | Matched donors, host assays and microbial features | Reviewed cross-domain models | New analysis; inferred 16S functional potential is not measured metagenomic function |

## Pathway and ontology analyses across modalities

Start with an assay-specific model, then map results into a shared, versioned
pathway vocabulary. Preserve the genes or regions actually tested in each assay.
Use full rankings for rank-based enrichment and explicit backgrounds for
overrepresentation. Correct within the full declared testing family. Overlapping
GO terms/pathways are correlated and should not be counted as independent
biological confirmations [29].

| Visualization family | What it reveals | Current implementation or researched tool |
|---|---|---|
| NES/enrichment dot or lollipop plot | Strength, direction where meaningful, and uncertainty | Automatic `enrich`; clusterProfiler/enrichplot [30] |
| Running enrichment curve with gene positions | Where a gene set occurs in the ranking | Automatic fgsea `plotEnrichment`, not a re-created approximation [31] |
| Leading-edge heatmap | Which genes drive ranked enrichment | Automatic `enrich` |
| Pathway-by-assay heatmap | Shared, conflicting and absent pathway evidence | Automatic `pathway-compare` |
| Pathway coverage heatmap | Whether each assay could test the pathway | Automatic `pathway-compare` |
| Shared versus assay-specific contributors | Whether pathway agreement comes from the same genes | Automatic `pathway-compare`; full contributor TSV |
| Enrichment network / similarity graph | Redundant or related pathway annotations | Automatic GMT-Jaccard network; EnrichmentMap [29] offers more interaction |
| Gene–concept bipartite network and UpSet | Gene membership and exclusive pathway intersections | clusterProfiler/enrichplot [30], UpSet [3]; Extension |
| Ontology DAG and semantic clustering | Parent/child relationships and redundant terms | GO-aware packages; Extension; GMT alone does not contain ontology edges |
| Molecular overlays on pathway diagrams | Location of RNA/protein/metabolite effects in curated topology | Pathview [32], PaintOmics 4 [16]; Extension |
| Sample-by-pathway activity heatmap | How pathway-level scores vary between specimens | ReactomeGSA [33], decoupleR [34], PathIntegrate [35]; New analysis |
| Direction-aware integrated pathway map | Evidence respecting prespecified inter-assay direction relationships | ActivePathways/DPM [36]; distinct from current nondirectional pathway comparison |

This area already has substantial prior art: multiGSEA [37], ReactomeGSA [33],
PaintOmics [16], ActivePathways/DPM [36] and PathIntegrate [35]. KilnBio's new
report layer contributes automatic, traceable outputs tied to its own validated
interfaces. Neither combined pathway analysis nor heatmaps are claimed as new
statistical inventions.

For ontology ORA, a selected-only gene list with the entire genome as background
can overstate evidence when measurement is restricted. CpGs, peaks, transcripts,
proteins and metabolites have different mapping multiplicities and selection
biases. A common pathway name helps interpretation; it does not make their raw
values exchangeable. The optional joint test in `pathway-compare` combines valid
raw pathway p-values with a dependence-safe partial-conjunction rule, and leaves
score direction and shared-gene support visible. See [the statistical
specification](VISUALIZATIONS.md#optional-joint-pathway-evidence).

## What recent work changes about the priorities

The 2024 PathIntegrate and directional integration papers reinforce two distinct
needs: sample-level pathway modeling and direction-aware cross-assay evidence
[35,36]. They are not interchangeable with plotting one enrichment list. A 2025
benchmark evaluates multimodal integration across seven tasks, including batch
correction and feature selection [38]; consequently, a good-looking embedding
alone is insufficient validation. New 2025 genome-browser resources emphasize
linked genomic views [21,39], supporting the priority of a future locus report
that juxtaposes variants, methylation, chromatin, RNA and protein annotations.
A 2026 signaling-state method, PathwayEmbed [40], is a candidate for a separately
validated activity module, not grounds to relabel enrichment as activity.

The next implementation priorities from this review are: (1) a shared donor/QC
dashboard and method-agreement report; (2) genome-build-checked locus views;
(3) CpG/DMR and Olink QC reports; (4) richer cell-marker/pseudobulk displays;
(5) independently tested activity, trajectory, latent-factor and microbiome
modules. These are extension opportunities, not completed capabilities.

## Primary references and official documentation

1. Ewels et al. MultiQC. Bioinformatics (2016). [PMID 27312411](https://pubmed.ncbi.nlm.nih.gov/27312411/), DOI 10.1093/bioinformatics/btw354.
2. Gu et al. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics (2016). [PMID 27207943](https://pubmed.ncbi.nlm.nih.gov/27207943/), DOI 10.1093/bioinformatics/btw313.
3. Lex et al. UpSet. IEEE TVCG (2014). [PMID 26356912](https://pubmed.ncbi.nlm.nih.gov/26356912/), DOI 10.1109/TVCG.2014.2346248.
4. Hahne and Ivanek. Visualizing Genomic Data Using Gviz and Bioconductor (2016). [PMID 27008022](https://pubmed.ncbi.nlm.nih.gov/27008022/), DOI 10.1007/978-1-4939-3578-9_16.
5. Akalin et al. methylKit. Genome Biology (2012). [PMID 23034086](https://pubmed.ncbi.nlm.nih.gov/23034086/), DOI 10.1186/gb-2012-13-10-r87.
6. Lopez-Delisle et al. pyGenomeTracks. Bioinformatics (2021; online 2020). [PMID 32745185](https://pubmed.ncbi.nlm.nih.gov/32745185/), DOI 10.1093/bioinformatics/btaa692.
7. Maksimovic et al. Gene set enrichment analysis for genome-wide DNA methylation data. Genome Biology (2021). [PMID 34103055](https://pubmed.ncbi.nlm.nih.gov/34103055/), DOI 10.1186/s13059-021-02388-x. [Current missMethyl documentation](https://bioconductor.org/packages/release/bioc/html/missMethyl.html).
8. Ramirez et al. deepTools2. NAR (2016). [PMID 27079975](https://pubmed.ncbi.nlm.nih.gov/27079975/), DOI 10.1093/nar/gkw257.
9. Yu et al. ChIPseeker. Bioinformatics (2015). [PMID 25765347](https://pubmed.ncbi.nlm.nih.gov/25765347/), DOI 10.1093/bioinformatics/btv145.
10. Khan and Mathelier. Intervene (2017). [PMID 28569135](https://pubmed.ncbi.nlm.nih.gov/28569135/).
11. Stuart et al. Single-cell chromatin state analysis with Signac. Nature Methods (2021). [Primary article](https://doi.org/10.1038/s41592-021-01282-5).
12. McLean et al. GREAT. Nature Biotechnology (2010). [PMID 20436461](https://pubmed.ncbi.nlm.nih.gov/20436461/), DOI 10.1038/nbt.1630.
13. [OlinkAnalyze official source and documentation](https://github.com/Olink-Proteomics/OlinkRPackage). Package documentation, not a peer-reviewed accuracy claim.
14. Kohler et al. MSstats version 4.0. Journal of Proteome Research (2023). [PMID 37018319](https://pubmed.ncbi.nlm.nih.gov/37018319/), DOI 10.1021/acs.jproteome.2c00834.
15. Pang et al. MetaboAnalyst 6.0. NAR (2024). [Primary article](https://doi.org/10.1093/nar/gkae253).
16. Liu et al. PaintOmics 4. NAR (2022). [PMID 35609982](https://pubmed.ncbi.nlm.nih.gov/35609982/), DOI 10.1093/nar/gkac352.
17. Wolf et al. SCANPY. Genome Biology (2018). [PMID 29409532](https://pubmed.ncbi.nlm.nih.gov/29409532/), DOI 10.1186/s13059-017-1382-0.
18. Argelaguet et al. MOFA+. Genome Biology (2020). [PMID 32393329](https://pubmed.ncbi.nlm.nih.gov/32393329/), DOI 10.1186/s13059-020-02015-1.
19. Bergen et al. Generalizing RNA velocity to transient cell states through dynamical modeling. Nature Biotechnology (2020). [PMID 32747759](https://pubmed.ncbi.nlm.nih.gov/32747759/).
20. Dimitrov et al. LIANA+. Nature Cell Biology (2024). [PMID 39223377](https://pubmed.ncbi.nlm.nih.gov/39223377/), DOI 10.1038/s41556-024-01469-w.
21. Empowering integrative and collaborative exploration of single-cell and spatial multimodal data with SGS genome browser. Cell Genomics (2025). [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12143324/), DOI 10.1016/j.xgen.2025.100848.
22. Mayakonda et al. Maftools. Genome Research (2018). [PMID 30341162](https://pubmed.ncbi.nlm.nih.gov/30341162/).
23. De Coster et al. NanoPack. Bioinformatics (2018). [PMID 29547981](https://pubmed.ncbi.nlm.nih.gov/29547981/), DOI 10.1093/bioinformatics/bty149.
24. Benjamini and Heller. Screening for partial conjunction hypotheses. Biometrics (2008). [PMID 18261164](https://pubmed.ncbi.nlm.nih.gov/18261164/), DOI 10.1111/j.1541-0420.2007.00984.x.
25. Benjamini and Yekutieli. FDR under dependency. Annals of Statistics (2001). [DOI 10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998).
26. Rohart et al. mixOmics. PLoS Computational Biology (2017). [PMID 29099853](https://pubmed.ncbi.nlm.nih.gov/29099853/), DOI 10.1371/journal.pcbi.1005752.
27. Bolyen et al. QIIME 2. Nature Biotechnology (2019). [PMID 31341288](https://pubmed.ncbi.nlm.nih.gov/31341288/), DOI 10.1038/s41587-019-0209-9.
28. McMurdie and Holmes. phyloseq. PLoS ONE (2013). [PMID 23630581](https://pubmed.ncbi.nlm.nih.gov/23630581/), DOI 10.1371/journal.pone.0061217.
29. Reimand et al. Pathway enrichment analysis and visualization of omics data. Nature Protocols (2019). [PMID 30664679](https://pubmed.ncbi.nlm.nih.gov/30664679/), DOI 10.1038/s41596-018-0103-9.
30. Wu et al. clusterProfiler 4.0. The Innovation (2021). [PMID 34557778](https://pubmed.ncbi.nlm.nih.gov/34557778/), DOI 10.1016/j.xinn.2021.100141. [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html).
31. [fgsea official tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html). Implementation reference; the tutorial identifies the algorithm's preprint separately.
32. Luo and Brouwer. Pathview. Bioinformatics (2013). [PMID 23740750](https://pubmed.ncbi.nlm.nih.gov/23740750/), DOI 10.1093/bioinformatics/btt285.
33. Griss et al. ReactomeGSA. Molecular & Cellular Proteomics (2020). [PMID 32907876](https://pubmed.ncbi.nlm.nih.gov/32907876/), DOI 10.1074/mcp.TIR120.002155.
34. Badia-i-Mompel et al. decoupleR. Bioinformatics Advances (2022). [PMID 36699385](https://pubmed.ncbi.nlm.nih.gov/36699385/), DOI 10.1093/bioadv/vbac016.
35. Wieder et al. PathIntegrate. PLoS Computational Biology (2024). [PMID 38527092](https://pubmed.ncbi.nlm.nih.gov/38527092/), DOI 10.1371/journal.pcbi.1011814.
36. Slobodyanyuk et al. Directional integration and pathway enrichment analysis for multi-omics data. Nature Communications (2024). [PMID 38971800](https://pubmed.ncbi.nlm.nih.gov/38971800/), DOI 10.1038/s41467-024-49986-4.
37. Canzler and Hackermueller. multiGSEA. BMC Bioinformatics (2020). [PMID 33287694](https://pubmed.ncbi.nlm.nih.gov/33287694/), DOI 10.1186/s12859-020-03910-x.
38. Liu et al. Multitask benchmarking of single-cell multimodal omics integration methods. Nature Methods (2025). [PMID 41083898](https://pubmed.ncbi.nlm.nih.gov/41083898/), DOI 10.1038/s41592-025-02856-3.
39. Takahashi et al. Japan Omics Browser. BMC Genomics (2025). [PMID 40335902](https://pubmed.ncbi.nlm.nih.gov/40335902/), DOI 10.1186/s12864-025-11639-1.
40. PathwayEmbed: a computational tool to quantify intracellular signaling transduction states from transcriptomic data. Bioinformatics (2026). [PMID 42209442](https://pubmed.ncbi.nlm.nih.gov/42209442/), DOI 10.1093/bioinformatics/btag346.

Pinned processing-output documentation checked for this review:
[RNA-seq 3.27.0](https://nf-co.re/rnaseq/3.27.0/docs/output),
[ATAC-seq 2.1.2](https://nf-co.re/atacseq/2.1.2/docs/output),
[ChIP-seq 2.1.0](https://nf-co.re/chipseq/2.1.0/docs/output),
[methylseq 4.2.0](https://nf-co.re/methylseq/4.2.0/docs/output).
Tool-specific output depends on enabled workflow options; these links do not
assert that every optional module ran in every KilnBio study.
