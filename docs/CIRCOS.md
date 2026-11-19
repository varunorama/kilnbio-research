# Circular genomic views across omics

KilnBio generates an offline **Circos-style** report from explicitly mapped assay
measurements. The number and order of rings follow the tracks in your configuration:
RNA + protein produces two rings; adding methylation, ATAC and ChIP produces five.
The renderer uses Matplotlib, already installed by every KilnBio setup profile.

```bash
python examples/circos/make_demo.py DEMO_INPUT
kilnbio circos DEMO_INPUT/config.json NEW_REPORT
```

Open `NEW_REPORT/report.html`. PDF, SVG, PNG and the exact plotted values are
generated automatically. The example uses fictional coordinates and values.

To include the circular view in the integrated-evidence report:

```bash
kilnbio multiomic-meta manifest.tsv --min-layers 2 \
  --circos-config mapped_tracks.json --out NEW_INTEGRATION
```

Track names must match assay layers in that integration manifest. The optional
view does not change the gene tests, multiplicity correction or nominations.
The `circos` command can also be a step in `kilnbio study`; list its configuration
and every referenced table in the study's `inputs` so resume detects changes.

## Required inputs

Use one declared genome build and **0-based, half-open** intervals throughout.
Paths resolve relative to the configuration file. There are no implicit genome
downloads, gene-symbol lookups, coordinate conversions or liftover operations.

```json
{
  "schema_version": 1,
  "genome_build": "GRCh38",
  "coordinate_system": "0-based-half-open",
  "chrom_sizes": "chrom_sizes.tsv",
  "tracks": [
    {
      "name": "RNA",
      "file": "rna_mapped.tsv",
      "genome_build": "GRCh38",
      "value_label": "log2FC",
      "value_min": -4,
      "value_max": 4
    }
  ],
  "links": {"file": "reviewed_links.tsv", "genome_build": "GRCh38"}
}
```

| Input | Required TSV columns | Meaning |
|---|---|---|
| Chromosome sizes | `chrom`, `length` | Positive integer lengths; row order defines chromosome order. |
| Each assay track | `chrom`, `start`, `end`, `value` | Mapped measurement and its original units. `NA` is preserved; extra annotation columns are retained. |
| Optional links | `chrom1`, `start1`, `end1`, `chrom2`, `start2`, `end2`, `relation` | Explicit interval pairs and a description of their relationship. Omit `links` when none are supplied. |

RNA and protein values need reviewed gene/protein-to-genome mappings. Methylation,
ATAC and ChIP intervals can use their native coordinates after checking the build
and coordinate convention. A protein-to-gene mapping supplies a display location;
it does not locate the protein molecule on the genome. Link labels should identify
the actual source, such as a measured contact or a prespecified feature map.

## Interpretation and safeguards

Rings show interval midpoints. Chromosome arcs are proportional to the supplied
lengths; each track has its own declared color scale and units. Gray crosses are
supplied missing values. Blank positions have no input record. Curves show only
provided links; neither proximity nor a similar color establishes regulation.
The plot is a descriptive view of your completed analyses.

Build disagreements, unknown chromosomes, out-of-bounds intervals and values
outside the declared scale stop rendering. No values are silently clipped or
imputed. The current display limits are 1–8 tracks, 64 chromosomes, 50,000 intervals
per track and 200 links with at most eight relation types. Larger inputs need an
explicit, recorded selection or aggregation before plotting. Crowded points may
overlap; TSV files preserve every supplied row. Input hashes and renderer version
are recorded in `visualization.json`.

The circular visualization concept is described by [Krzywinski et al., *Genome
Research* (2009), PMID 19541911](https://pubmed.ncbi.nlm.nih.gov/19541911/),
[doi:10.1101/gr.092759.109](https://doi.org/10.1101/gr.092759.109).
KilnBio's renderer does not invoke the original Circos executable.
