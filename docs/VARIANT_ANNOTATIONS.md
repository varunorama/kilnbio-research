# Connect genomic variants to molecular context

These optional interfaces connect called alleles to public annotations and protein structure records for **research prioritization**. They do not assign a clinical diagnosis, establish a disease mechanism or run a general-purpose LLM. See [research-use scope](RESEARCH_USE.md).

| Interface | Input and output | What it actually does |
| --- | --- | --- |
| `variant-annotations` | Unique normalized allele TSV plus local resources → complete query roster and source-specific JSONL | Exact assembly/chromosome/position/REF/ALT joins to ClinVar VCF, AlphaMissense published scores and VEP-annotated VCF |
| `vep-annotate` | VCF, assembly, reference FASTA and local cache → plan or annotated VCF | Runs installed Ensembl VEP 116.3 with cache 116, offline, with reference checking; keeps all transcript consequences |
| `alphafold` | UniProt accessions and optional one-letter missense changes → structure metadata and mapping status | Looks up AlphaFold Database records or reads saved responses; records model version and checks the reference residue against the returned sequence |
| `alphagenome` | Variants plus condition → optional remote predictions | Separate [AlphaGenome 2 adapter](ALPHAGENOME.md); never a substitute for measured molecular evidence |

GATK HaplotypeCaller and Mutect2, and selectable DeepVariant/Strelka routes, belong to [variant discovery](EXOME_AND_CH.md), before annotation. AlphaFold Database lookup is **not local AlphaFold2 inference**, mutation-induced structure prediction or protein stability estimation. AlphaMissense uses released scores; KilnBio does not train or execute a new AlphaMissense model.

## Run local annotation

Use Linux/WSL2 with Python 3.10+. The annotation and structure interfaces use the Python standard library. Install VEP and its assembly-matched cache separately when running the local VEP step. Paths inside JSON resolve relative to that JSON file. Every output directory must be new.

```bash
kilnbio vep-annotate vep.json vep-plan
kilnbio vep-annotate vep.json vep-run --execute
kilnbio variant-annotations annotations.json annotated-research
kilnbio alphafold proteins.tsv structure-plan
kilnbio alphafold proteins.tsv structures --execute
# Alternatively, replay saved accession.json responses without contacting the API:
kilnbio alphafold proteins.tsv structures-replayed --cache saved-responses
```

Templates: [annotation configuration](../examples/variant-annotations.json), [VEP configuration](../examples/vep-annotate.json), [alleles](../examples/annotation-variants.tsv), [proteins](../examples/annotation-proteins.tsv). Template resource paths are placeholders; obtain public resources yourself and record their actual versions. Neither example silently downloads a large reference.

Before joining alleles, normalize and split both variant calls and resources against the **same reference** where needed, for example `bcftools norm -f genome.fa -m -any`. The interface does not perform liftover, infer a genome build or left-normalize an indel. It accepts GRCh37/GRCh38 primary human chromosomes and biallelic A/C/G/T alleles; structural and symbolic variants require a separate interpretation route. A matching multiallelic VCF row is rejected rather than guessing allele attribution. Query TSVs contain unique variants, not donor identifiers.

`variants.tsv` retains every query allele, including those absent from a resource. `not_found_in_resource` is **not a benign classification**. The JSONL files preserve all matched source records, transcript consequences and original ClinVar INFO fields, including review status, conflicts and separate somatic/germline assertions. No transcript is silently selected and no source classifications are merged into a clinical verdict. The receipt records assembly, source URL/version, file SHA-256, query denominator and match counts. Failed runs remain marked failed.

ClinVar and AlphaMissense files are streamed with bounded resource memory; scanning a full unindexed release can take time. A prefiltered local resource is acceptable, but its version must identify the subset, and missing matches then mean absence from that subset. Preserve the original source file and its license header. The AlphaMissense file retrieved during development carried a CC BY-NC-SA notice; verify the terms attached to your exact download rather than assuming that code and data licenses are identical. KilnBio does not redistribute these databases.

## Structure mapping and evidence limits

The structure lookup sends **only explicitly listed UniProt accessions** to the public AlphaFold Database API when `--execute` is set. Without that flag it writes an offline plan. Cached replay copies the response and hashes it; the original retrieval date should be retained alongside your cache. Accessions, isoforms, fragments, sequence ranges, model versions and URLs remain explicit. A missing model or a reference-residue mismatch remains an unsuccessful mapping. An offset within a fragment is a sequence position, not an assertion about PDB residue numbering. Mean pLDDT describes model confidence, not variant pathogenicity.

Map genomic variants to a documented transcript/protein/isoform before supplying protein substitutions. VEP, AlphaMissense and UniProt identifiers are not automatically interchangeable. KilnBio deliberately requires the protein mapping to be reviewed rather than guessing among incompatible isoforms.

The new interfaces have allele-matching, conflict/missingness, transcript-retention, sequence-mapping and failure-path tests. Small public resource subsets and a live public structure lookup provide software interoperability checks. The standalone local VEP route has plan/contract tests, **not a completed full-cache annotation benchmark**. These checks do not measure variant pathogenicity accuracy, whole-exome sensitivity or clinical performance. The live AlphaGenome prediction service remains untested here.

## Sources

- McLaren et al. The Ensembl Variant Effect Predictor. *Genome Biology* (2016). [Paper](https://doi.org/10.1186/s13059-016-0974-4); [VEP 116.3 source](https://github.com/Ensembl/ensembl-vep/tree/release/116.3).
- [NCBI ClinVar download documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/) and [data-use limitations](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).
- Cheng et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* (2023). [Paper](https://doi.org/10.1126/science.adg7492); [official implementation and data links](https://github.com/google-deepmind/alphamissense).
- Jumper et al. Highly accurate protein structure prediction with AlphaFold. *Nature* (2021). [Paper](https://doi.org/10.1038/s41586-021-03819-2).
- Bertoni et al. AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage. *Nucleic Acids Research* (2026 issue). [Database paper](https://doi.org/10.1093/nar/gkaf1226). This is the resource accessed by the lookup command; it does not run the structure-prediction algorithm.
- [AlphaFold Database API](https://alphafold.ebi.ac.uk/api-doc), [release notes](https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes) and [current response-field migration](https://www.ebi.ac.uk/pdbe/news/breaking-changes-afdb-predictions-api). The adapter uses `modelEntityId`, `sequence`, `sequenceStart` and `sequenceEnd`.

The contribution is the reusable, source-traceable connection to multi-omic research. These established tools and databases are credited to their developers; their use alone is not a new biological discovery.
