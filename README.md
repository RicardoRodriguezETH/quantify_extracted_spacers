# quantify_extracted_spacers

A pipeline for aligning CRISPR spacer sequences to reference genomes, mapping them to annotated genes, and generating quantitative analyses of spacer integration patterns across samples.

---

## Pipeline overview

```
merged_counts_table_filtered.tsv
          │
          ▼
  extract_spacers          cut -f1 → all_spacers.tsv
          │
          ▼
  align_spacers            Biostrings::matchPattern (± strand, 0–2 mismatches)
                           GenomicRanges::distanceToNearest (gene annotation)
                           → spacers_aligned.tsv
          │
          ▼
  process_alignments       merge counts + alignments
                           filter → unique / multi-mapped / unmapped
                           gene-level detection counts and raw sums
                           log2 CPM normalization (edgeR)
                           diagnostic plots
                           → gene_level/, spacer_level/, plots/
```

---

## Requirements

### Python dependencies

Install with [uv](https://github.com/astral-sh/uv):

```bash
uv pip install piper pyyaml
```

### R packages

Install from Bioconductor and CRAN:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", "BiocParallel", "GenomicRanges",
  "GenomeInfoDb", "rtracklayer", "edgeR", "pwalign"
))
install.packages(c(
  "tidyverse", "data.table", "patchwork",
  "logger", "lubridate", "stringdist", "progressr",
  "pheatmap"   # optional: prettier CPM heatmap
))
```

A full environment spec (conda) is available in `environment.yml`.

---

## Quick start

The repository includes a toy dataset (30k spacers) and reference genomes so you can run the pipeline immediately after cloning:

```bash
python pipeline.py --config config.yaml
```

`config.yaml` is pre-configured to use the toy dataset. To resume after interruption:

```bash
python pipeline.py --config config.yaml --recover
```

To force a clean re-run:

```bash
python pipeline.py --config config.yaml --new-start
```

---

## Configuration

All project-specific settings live in `config.yaml`:

```yaml
project: P3022

# Input data
counts: input/data/P3022/merged_counts_table_filtered_test30k.tsv

# Output
output_dir: output/P3022

# Reference genomes
genome:      input/genomes/chromosomes/MG1655_U00096_DSMZ.fasta
plasmid:     input/genomes/plasmids/pAK0033.fasta
genome_gff:  input/genomes/chromosomes/MG1655_U00096_DSMZ.gff
plasmid_gff: input/genomes/plasmids/pAK0033.gff

# Compute
cores:      32
chunk_size: 15000 # default 100000
```

Relative paths are resolved against the config file's directory, so the pipeline can be run from any working directory.

To switch to the full dataset, change only the `counts` line:

```yaml
counts: input/data/P3022/merged_counts_table_filtered.tsv
```

### Batch size and memory

`chunk_size` and `cores` together control peak memory during alignment. Each worker holds a full copy of the reference genome in memory while processing its slice of the batch, so peak usage scales roughly as:

```
peak memory ≈ chunk_size × cores × (genome size factor)
```

As a reference point, aligning 30k spacers in 2 batches of 15k on 32 cores peaked at **~34 GB**. To reduce memory, lower `cores` or `chunk_size` — both have the same effect. `chunk_size` also controls how often intermediate results are written to disk, which is useful for very large datasets.

---

## Reference genomes

```
input/genomes/
  chromosomes/      ← one .fasta + .gff per bacterial strain
    MG1655_U00096_DSMZ.fasta
    MG1655_U00096_DSMZ.gff
  plasmids/         ← one .fasta + .gff per plasmid
    pAK0033.fasta
    pAK0033.gff
```

`config.yaml` selects which chromosome and plasmid to use for a given run.

---

## Output structure

```
output/<project>/
├── match_type_summary.tsv       # Per-sample unique / multi / unmapped %
│
├── spacer_level/
│   ├── all_spacers.tsv          # target_seq column extracted from counts table
│   ├── spacers_aligned.tsv      # All spacers with genomic coordinates + gene annotations
│   └── spacers_unique_genome.tsv  # Filtered: unique matches on chromosome only
│
├── gene_level/
│   ├── gene_counts.tsv          # Spacers detected (count > 0) per gene per sample
│   ├── gene_count_sums.tsv      # Raw spacer count sums per gene per sample
│   ├── gene_counts_cpm.tsv      # Log2 CPM normalized (from raw sums, edgeR)
│   └── gene_counts_cpm_hits.tsv # High-variance genes (max_cpm > 7, range > 0.5)
│
├── plots/
│   ├── s01_spacer_abundance_distribution.png
│   ├── s02_match_type_composition.png
│   ├── s03_cumulative_match_types.png
│   ├── s04_spacer_length_distribution.png
│   ├── s05_strand_distribution.png
│   ├── s06_distance_to_gene.png
│   ├── s07_genome_coverage.png
│   ├── s08_genic_fraction_vs_abundance.png
│   ├── g01_top_genes_spacer_count.png
│   ├── g02_spacers_per_gene_distribution.png
│   ├── g03_sample_scatter_pairs.png
│   ├── g04_sample_correlation_heatmap.png
│   ├── g05_cpm_variability.png
│   ├── g06_cpm_heatmap.png
│   └── g07_sense_antisense.png
│
└── logs/
    ├── align_spacers_log.txt        # Alignment runtime log (batches, timing)
    └── process_alignments_log.txt   # Processing log (QC stats, file paths)
```

---

## Scripts

| Script | Description |
|--------|-------------|
| `scripts/align_spacers.R` | Streams input in batches (`chunk_size` rows), aligns each batch to chromosome and plasmid on both strands (up to 2 mismatches), maps to nearest gene, appends results to output file |
| `scripts/process_alignments.R` | Merges alignment results with raw counts, computes match-type statistics, builds gene-level tables, normalizes with CPM, generates diagnostic plots |

### `process_alignments.R` optional arguments

```
Rscript scripts/process_alignments.R <aligned_tsv> <counts_tsv> <output_dir> \
    [genome_seqname]    # default: auto-detect (seqname with most unique hits)
    [cpm_max_thresh]    # default: 7
    [cpm_range_thresh]  # default: 0.5 (auto-set to 0 for single-sample runs)
```

---

## Project structure

```
quantify_extracted_spacers/
├── README.md
├── config.yaml              # All project parameters and paths
├── environment.yml          # Minimal environment spec (conda / reference)
├── pipeline.py              # Pipeline runner (reads config.yaml)
├── scripts/
│   ├── align_spacers.R
│   └── process_alignments.R
├── input/
│   ├── data/
│   │   └── P3022/
│   │       └── merged_counts_table_filtered_test30k.tsv   # toy dataset (tracked)
│   └── genomes/                                           # reference genomes (tracked)
│       ├── chromosomes/
│       └── plasmids/
└── output/                  # Pipeline results (not tracked by git)
```
