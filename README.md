# Configurable sgRNA Design Pipeline

`sgrna_design_pipeline_full.R` is an R/Bioconductor pipeline for designing and documenting CRISPR guide RNAs. It can design guides from a gene, genomic interval, FASTA file, or raw sequence, and can also re-process an existing guide-set RDS file.

The pipeline is configuration-driven. Most decisions are made in `default_config()` near the top of the script, or in a separate configuration file. The HTML report and tabular outputs are generated from the same configuration and guide table.

## Main capabilities

- SpCas9 and other nuclease definitions available from `crisprBase`.
- Custom PAMs, spacer lengths, PAM side, weights, and target type.
- BSgenome, FASTA, or raw-sequence input.
- Gene, genomic-region, and custom-sequence targets.
- CDS, transcript, or exon feature selection for gene targets.
- Mature, spliced-transcript coordinates for guide cleavage sites.
- Transcript/exon plot with selectable transcripts and guide markers.
- CRISPR knockout, CRISPRi, CRISPRa, base-editing, and prime-editing annotations.
- On-target scores, sequence features, optional off-target alignment/scoring, and VCF/SNP annotation.
- GC, homopolymer, strand, score, sequence-exclusion, TSS-window, and custom filtering.
- Individual guides and optional two-guide deletion pairs.
- CSV, RDS, HTML, BED, paired-guide CSV, and transcript-map CSV outputs.

## Requirements

The pipeline requires R plus Bioconductor packages used by `crisprDesign`:

- `crisprBase`
- `crisprDesign`
- `jsonlite`

Depending on the configuration, it also uses a BSgenome package, `GenomicFeatures` and a compatible TxDb package, an organism annotation package, and `Biostrings`/`Rsamtools` for FASTA input.

The script does not install missing packages unless `runtime$install_missing = TRUE`. Installing packages may require an internet connection and appropriate Bioconductor repositories.

## Running the pipeline

The simplest workflow is:

```bash
Rscript --vanilla sgrna_design_pipeline_full.R
```

This runs the values in `default_config()`.

For a separate configuration file, create `my_config.R` containing a list named `config`:

```r
config <- list(
 project_name = "MYGENE_sgRNA_design",
 target = list(
 type = "gene",
 name = "MYGENE",
 gene_id = "1234"
 ),
 output = list(
 directory = "/home/goro/R/output",
 prefix = "MYGENE_sgRNA_design"
 )
)
```

Run it with:

```bash
Rscript --vanilla sgrna_design_pipeline_full.R my_config.R
```

The supplied configuration is recursively merged into the defaults, so only values that need changing have to be included. An RDS configuration is also supported if it contains a list:

```bash
Rscript --vanilla sgrna_design_pipeline_full.R my_config.rds
```

Using `Rscript` is recommended because it gives a clear command-line run and writes the report to the configured output directory.

## Configuration reference

All configurable sections are at the top of the script inside `default_config()`.

### `mode` and project naming

```r
mode = "design" # "design" or "rds"
project_name = "TP53_sgRNA_design"
```

- `design` generates guides with `crisprDesign`.
- `rds` reads an existing guide object from `input$rds_file` and applies downstream annotation, filtering, pairing, and reporting.
- `project_name` is used in the report title and summary.
- `output$prefix` controls output filenames. Change both `project_name` and `output$prefix` when changing projects so an old report name is not reused.

### `runtime`

```r
runtime = list(
 install_missing = FALSE,
 quiet = FALSE
)
```

Set `install_missing = TRUE` only when automatic Bioconductor/CRAN installation is desired. Otherwise install packages manually and keep it `FALSE` for reproducibility.

### `genome`

#### BSgenome input

```r
genome = list(
 source = "bsgenome",
 package = "BSgenome.Hsapiens.UCSC.hg38",
 object = "BSgenome.Hsapiens.UCSC.hg38"
)
```

`package` is the installed package name and `object` is the exported genome object. The genome must match the TxDb/annotation build used for transcript coordinates.

#### FASTA input

```r
genome = list(
 source = "fasta",
 fasta = "/path/to/genome.fa"
)
```

FASTA targets require `target$seqnames`, `target$start`, and `target$end`. FASTA sequence names must match the target chromosome names.

#### Raw sequence input

```r
genome = list(
 source = "sequence",
 sequence = "ACGT..."
)
```

Alternatively put the sequence in `target$sequence`. Raw-sequence mode is useful for amplicons, plasmids, synthetic constructs, or regions without genomic annotation. Set `target$sequence_start` if sequence coordinates should be offset to a genomic or custom coordinate system.

### `target`

#### Gene target

```r
target = list(
 type = "gene",
 name = "TP53",
 gene_id = "7157",
 gene_column = "gene_id",
 feature_type = "cds"
)
```

Gene mode requires a BSgenome source and a compatible TxDb package in `annotation$txdb_package`.

`gene_id` must use the identifier system represented by the TxDb grouping. With the default human known-gene annotation, it is an Entrez identifier. `feature_type` can be `cds`, `transcripts`, or `exons`.

The `name` value is a display label. It does not identify a gene by itself; change `gene_id` as well.

#### Genomic-region target

```r
target = list(
 type = "region",
 name = "TP53_promoter",
 seqnames = "chr17",
 start = 7670000,
 end = 7675000
)
```

`start` and `end` are 1-based inclusive genomic coordinates. The pipeline validates that `start <= end`.

#### Custom-sequence target

```r
target = list(
 type = "sequence",
 name = "my_amplicon",
 sequence = "ACGTACGT...",
 sequence_start = 1L
)
```

Use this for any DNA sequence that can be represented as a character string. Transcript and exon annotation is not automatically available unless a compatible genomic model is also supplied.

### `annotation`

```r
annotation = list(
 enabled = TRUE,
 txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene",
 txdb_object = "TxDb.Hsapiens.UCSC.hg38.knownGene",
 organism_package = "org.Hs.eg.db",
 organism_keytype = "ENTREZID"
)
```

Set `enabled = FALSE` when transcript annotation is not wanted or is unavailable. Genomic guide design still works, but transcript positions, transcript plots, exon labels, and mature-transcript columns are unavailable.

The TxDb and BSgenome must describe the same genome assembly. For another species, replace the BSgenome, TxDb, and organism package with matching packages.

### `transcript`

```r
transcript = list(
 coordinate_system = "mature_transcript",
 selected_id = NULL,
 selector_in_report = TRUE
)
```

The supported coordinate system is `mature_transcript`, meaning the spliced transcript/cDNA sequence counted from 1.

- `selected_id = NULL` selects the first annotated transcript initially.
- Set `selected_id` to a transcript name such as `ENST00000269305.9` to choose it initially.
- `selector_in_report = TRUE` adds a transcript dropdown to the HTML report.

The report contains two related views: the upper genomic transcript/exon plot, with guide cut markers overlaid on the selected transcript row; and the lower plot, showing guide positions along the mature transcript.

An exonic cut receives a 1-based mature-transcript position. An intronic cut remains `NA` for mature-transcript position because introns are not part of the mature transcript, but it is labelled as `intron` in the transcript map.

### `cell_context`

```r
cell_context = list(
 name = "HEK293T",
 species = "human",
 strain = NULL,
 vcf = NULL,
 metadata = list()
)
```

`name`, `species`, `strain`, and `metadata` are recorded as context in the configuration/report. They do not by themselves change guide sequences.

If `vcf` is supplied, the pipeline attempts to add SNP annotations with `crisprDesign::addSNPAnnotation`. A cell-line-specific VCF is therefore the way to model sequence variation in a particular culture.

### `nuclease`

#### Built-in nuclease

```r
nuclease = list(
 name = "SpCas9",
 pams = NULL,
 weights = NULL,
 pam_side = "3prime",
 spacer_length = 20L,
 target_type = "DNA",
 metadata = list()
)
```

With `pams = NULL`, `name` must be one of the nuclease definitions exposed by `crisprBase::getAvailableCrisprNucleases()`.

#### Custom nuclease/PAM

```r
nuclease = list(
 name = "CustomCas",
 pams = c("NGG", "NAG"),
 weights = c(1, 0.5),
 pam_side = "3prime",
 spacer_length = 20L,
 target_type = "DNA",
 metadata = list(description = "custom PAM model")
)
```

The exact PAM syntax and weight interpretation follow `crisprBase::CrisprNuclease`. Check the installed package documentation for unusual PAM alphabets, RNA targets, or non-standard nuclease definitions.

### `design`

```r
design = list(
 canonical = TRUE,
 both_strands = TRUE,
 strict_overlap = TRUE,
 remove_ambiguities = TRUE,
 remove_duplicates = TRUE
)
```

- `both_strands = TRUE` searches both strands.
- `canonical` and `strict_overlap` control `crisprDesign` interval handling.
- `remove_ambiguities` removes ambiguous spacer sequences.
- `remove_duplicates` removes duplicate guide candidates.

### `scoring`

```r
scoring = list(
 add_on_target = TRUE,
 methods = c("ruleset1", "crisprater"),
 add_sequence_features = TRUE,
 add_composite = TRUE
)
```

The default methods are fast, reproducible SpCas9-compatible on-target scores. `add_sequence_features` adds sequence-derived features. The pipeline normalises numeric score columns and calculates `composite_score` as their row-wise mean when possible.

Set `add_on_target = FALSE` only when using an input guide object that already contains suitable scores or when scoring is intentionally omitted. If every score is missing, score plots and score filters cannot rank guides meaningfully.

### `off_target`

```r
off_target = list(
 enabled = FALSE,
 bowtie_index = NULL,
 mismatches = 3L,
 all_alignments = TRUE
)
```

When enabled, the pipeline attempts to call `addSpacerAlignments` and `addOffTargetScores`. `bowtie_index` must point to a compatible Bowtie index. If no index is supplied, the pipeline warns and skips off-target search.

Off-target analysis is genome-build-specific. Use an index made from the same genome assembly and chromosome naming scheme as the design genome.

### `modality`

```r
modality = list(
 type = "CRISPRko",
 tss_window = c(NA_real_, NA_real_)
)
```

Supported values are `CRISPRko`, `CRISPRi`, `CRISPRa`, `base_edit`, and `prime_edit`.

For CRISPRi/a, set a finite signed TSS window, for example:

```r
modality = list(type = "CRISPRi", tss_window = c(-100, 50))
```

The distance is strand-aware: positive-strand transcripts use `cut_site - TSS`, while negative-strand transcripts use `TSS - cut_site`.

### `editing`

```r
editing = list(
 window = 4:8,
 from = "C",
 to = "T",
 pbs_length = 13L,
 rtt_length = 20L,
 pbs_sequence = NULL,
 rtt_sequence = NULL
)
```

For base editing, `window` is 1-based protospacer positions. The output column `editing_candidates` lists matching positions such as `C4>T;C7>T`.

For prime editing, the output records the configured PBS and RTT lengths. The current pipeline annotates candidate guide sequences; it does not infer a complete pegRNA edit template from a desired variant.

### `cloning`

```r
cloning = list(
 enabled = FALSE,
 five_prime_adapter = "",
 three_prime_adapter = "",
 scaffold = "",
 restriction_sites = character()
)
```

When enabled, `cloning_oligo` is built as `five_prime_adapter + protospacer + three_prime_adapter + scaffold`. Strings in `restriction_sites` are searched literally and reported in `restriction_site_hit`. Verify strand orientation, overhangs, scaffold conventions, and enzyme compatibility before ordering oligos.

### `filters`

```r
filters = list(
 min_composite_score = NA_real_,
 min_gc = NA_real_,
 max_gc = NA_real_,
 max_homopolymer = 4L,
 strands = c("+", "-"),
 exclude_sequences = character(),
 custom = NULL
)
```

- `NA_real_` disables a numeric threshold.
- `min_gc` and `max_gc` are percentages, e.g. `30` and `70`.
- `max_homopolymer = 4` rejects runs longer than four identical bases.
- `strands` selects allowed guide strands.
- `exclude_sequences` removes exact protospacer matches, case-insensitively.
- `custom` can be a function returning one logical value per guide:

```r
filters = list(
 custom = function(tbl) {
 tbl$gc_percent >= 35 & tbl$gc_percent <= 65 &
 !grepl("TTTT", tbl$protospacer, fixed = TRUE)
 }
)
```

For CRISPRi/a, a finite `modality$tss_window` is applied as an additional filter.

### `paired_guides`

```r
paired_guides = list(
 enabled = TRUE,
 min_distance = 20L,
 max_distance = 1000L,
 opposite_strands_only = FALSE
)
```

Normal design creates individual guides. Pair design combines two filtered individual guides whose genomic cleavage sites are within the configured distance range.

- `min_distance` and `max_distance` are genomic cut-site distances.
- `opposite_strands_only = TRUE` keeps only pairs on opposite strands.
- Pair `mean_score` is the mean of the two guide composite scores.
- `deletion_bp` is the absolute distance between the two cut sites.

Pairing does not replace the normal guide output. It creates an additional `_pairs.csv` file and a paired-guide section in the HTML report. If no two guides satisfy the settings, the pairs file contains only a CSV header; this is expected and means the pair constraints or upstream filters need adjustment.

### `output`

```r
output = list(
 directory = "/home/goro/R/output",
 prefix = "TP53_sgRNA_design",
 write_rds = TRUE,
 write_csv = TRUE,
 write_html = TRUE,
 write_bed = TRUE,
 write_pairs = TRUE,
 write_transcript_map = TRUE,
 top_guides = 20L
)
```

The prefix determines every output filename. For `TP53_sgRNA_design`, the pipeline writes:

| File | Contents |
|---|---|
| `TP53_sgRNA_design.csv` | One row per filtered guide |
| `TP53_sgRNA_design.rds` | R data frame of the filtered guides |
| `TP53_sgRNA_design.html` | Interactive report |
| `TP53_sgRNA_design.bed` | BED-like guide intervals for genome browsers |
| `TP53_sgRNA_design_pairs.csv` | Two-guide deletion candidates |
| `TP53_sgRNA_design_transcript_map.csv` | All guide/transcript mappings |

`top_guides` controls how many guides and pairs are shown in the HTML tables. CSV/RDS outputs contain the full filtered result.

## Example configurations

### Change the gene

Change all project-specific identity fields and rerun the script:

```r
config <- list(
 project_name = "ITGA2B_sgRNA_design",
 target = list(
 type = "gene",
 name = "ITGA2B",
 gene_id = "3674",
 feature_type = "cds"
 ),
 output = list(
 directory = "/home/goro/R/output",
 prefix = "ITGA2B_sgRNA_design"
 )
)
```

Changing only `target$name` changes the display label but does not change the gene searched. Changing only `project_name` changes the report title but does not change the output filename. `target$gene_id` identifies the gene for gene-mode design.

### Design a genomic interval

```r
config <- list(
 project_name = "chr17_interval",
 target = list(
 type = "region",
 name = "chr17_interval",
 seqnames = "chr17",
 start = 7670000,
 end = 7675000
 ),
 annotation = list(enabled = FALSE),
 output = list(prefix = "chr17_interval")
)
```

### Design a custom sequence

```r
config <- list(
 project_name = "amplicon",
 genome = list(source = "sequence", sequence = "ACGTACGTACGT..."),
 target = list(type = "sequence", name = "amplicon"),
 annotation = list(enabled = FALSE),
 output = list(prefix = "amplicon")
)
```

### Select a transcript initially

```r
config <- list(
 transcript = list(
 selected_id = "ENST00000269305.9",
 selector_in_report = TRUE
 )
)
```

### Configure a deletion pair search

```r
config <- list(
 paired_guides = list(
 enabled = TRUE,
 min_distance = 100L,
 max_distance = 500L,
 opposite_strands_only = TRUE
 )
)
```

### Configure base editing

```r
config <- list(
 modality = list(type = "base_edit"),
 editing = list(window = 4:8, from = "C", to = "T")
)
```

## Understanding the outputs

### Individual guide table

The main CSV contains normalized fields including `guide_id`, `protospacer`, `pam`, `strand`, `seqnames`, `start`, `end`, `cut_site`, `gc_percent`, score columns, `composite_score`, `transcript`, `mature_transcript_position`, `transcript_region`, `exon_number`, `distance_to_tss`, and modality/editing/cloning annotations when applicable.

### Transcript map

The transcript-map CSV contains one row for every guide/transcript combination. It can be much larger than the main guide CSV because each guide is compared with every annotated transcript. Use it to inspect isoform-specific positions.

For an exonic cut, `transcript_position` is the 1-based coordinate on the spliced mature transcript. For an intronic cut, `transcript_position` is `NA` and `transcript_region` is `intron`.

### BED output

The BED file uses zero-based start and one-based-exclusive end conventions expected by genome browsers. Check it against the exact genome assembly and coordinate conventions of the downstream tool.

## Troubleshooting

### The report still has the old gene name

The script may be writing to an old prefix or you may be opening an old HTML file. Check both:

```r
project_name = "NEWGENE_sgRNA_design"
output$prefix = "NEWGENE_sgRNA_design"
target$name = "NEWGENE"
target$gene_id = "NEW_GENE_IDENTIFIER"
```

Then rerun with `Rscript` and open the newly printed `Interactive report:` path. Browser tabs can cache an older local report, so close/reopen the file after regeneration.

### The pairs CSV is blank

A header-only pairs CSV means no pair passed the constraints. Check `paired_guides$enabled`, the minimum and maximum distances, and `opposite_strands_only`. Also check that the normal guide CSV contains at least two guides with finite `cut_site` values. Strict score, GC, TSS, or custom filters can leave too few guides for pairing.

### No guides are returned

Common causes are chromosome names that do not match the genome package (`chr17` versus `17`), target coordinates outside the genome, an incompatible PAM/nuclease, filters that are too strict, the wrong gene identifier system, or a genome/TxDb assembly mismatch. Temporarily relax filters and verify the coordinate system with a known target.

### Scores are missing or all `NA`

Make sure `scoring$add_on_target = TRUE` and that the selected scoring methods support the configured nuclease. If using `mode = "rds"`, confirm the input object already contains usable numeric score columns. Without scores, guides can still be designed, but ranking and score plots are limited.

### There is no transcript position

Transcript positions require `annotation$enabled = TRUE`, a working TxDb, genomic guide coordinates, and a compatible genome/annotation assembly. Custom raw sequences do not automatically map to mature transcripts.

## Reproducibility and biological checks

The pipeline automates guide enumeration and annotation; it does not replace experimental design review. Before ordering or cloning guides:

1. confirm the genome assembly and chromosome naming;
2. verify the selected transcript and exon numbering;
3. inspect the protospacer and PAM sequence in the intended cell context;
4. use a cell-line VCF when sequence variants matter;
5. run off-target analysis with a matching genome index;
6. verify guide orientation, cleavage assumptions, and cloning overhangs; and
7. review guides manually in a genome browser or equivalent tool.

The `cell_context` label is metadata unless a VCF is supplied. A cell-culture name alone cannot alter genomic sequence. A score is a prioritization metric, not a guarantee of activity.

## Current default output

With the default TP53 configuration, output is written under `/home/goro/R/output` using the prefix `TP53_sgRNA_design`. To create a different project, use a separate configuration file or edit the values in `default_config()` before rerunning.

## Package references

This pipeline relies on the following Bioconductor and CRAN packages. Cite the relevant ones when publishing guide designs.

### CRISPR design and nuclease definitions

- **crisprBase** and **crisprDesign** — Hoberecht L, Perampalam P, Lun A, Fortin J. "A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies." _Nature Communications_, 2022, **13**(1), 6568. DOI: [10.1038/s41467-022-34320-7](https://doi.org/10.1038/s41467-022-34320-7).

### Genome ranges, sequence annotation, and infrastructure

The following packages are part of the GenomicRanges/Biostrings ecosystem described in:

- Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan MT, Carey VJ. "Software for computing and annotating genomic ranges." _PLoS Computational Biology_, 2013, **9**(8), e1003118. DOI: [10.1371/journal.pcbi.1003118](https://doi.org/10.1371/journal.pcbi.1003118).

- **BSgenome** — genome sequence representation.
- **Biostrings** — biological string manipulation.
- **GenomicFeatures** — TxDb and genomic interval queries.
- **Rsamtools** — BAM/FASTA file access.
- **AnnotationDbi** — annotation database interface.

### Species-specific annotation data packages

These data packages provide genome sequences and annotations for the target organism. Cite them through the Bioconductor project:

- Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Oleś AK, Pagès H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M. "Orchestrating high-throughput genomic analysis with Bioconductor." _Nature Methods_, 2015, **12**(2), 115–121. DOI: [10.1038/nmeth.3252](https://doi.org/10.1038/nmeth.3252).

- **BSgenome.Hsapiens.UCSC.hg38** — full genome sequences for Homo sapiens (UCSC hg38).
- **TxDb.Hsapiens.UCSC.hg38.knownGene** — transcript annotation for hg38.
- **org.Hs.eg.db** — genome-wide annotation for human.

### Utilities

- **jsonlite** — Ooms J. "The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects." CRAN, 2014. DOI: [10.32614/CRAN.package.jsonlite](https://doi.org/10.32614/CRAN.package.jsonlite).

### Reporting

- **Plotly.js** — Plotly Technologies Inc. Plotly.js: Open-source JavaScript graphing library. https://plotly.com/javascript/

For species other than human, replace the BSgenome, TxDb, and organism packages with the matching annotation packages for that genome assembly.

## Example output files

The repository includes the example output files from a TP53 sgRNA design run (hg38, CDS, SpCas9):

| File | Description |
|---|---|
| `TP53_sgRNA_design.csv` | Full guide table with scores, positions, and annotations |
| `TP53_sgRNA_design.rds` | R serialised data frame of the filtered guides |
| `TP53_sgRNA_design.html` | Interactive HTML report with plots and tables |
| `TP53_sgRNA_design.bed` | BED-formatted guide intervals for genome browsers |
| `TP53_sgRNA_design_pairs.csv` | Paired-guide deletion candidates |
| `TP53_sgRNA_design_transcript_map.csv` | All guide/transcript isoform mappings |