# Homework 1

This folder contains the deliverables for Homework 1 (GC content analysis of UCSC yeast mRNA FASTA).

## Contents

- `scripts/mrna_gc_content.R`: R script used to parse `mrna.fa.gz`, compute per-record metrics, and generate outputs.
- `results/mrna_metrics.tsv`: TSV table with `accession`, `length`, and `gc_content` (4 decimals), sorted by GC descending.
- `results/gc_content_distribution.png`: GC distribution plot (1600x900 px, 200 dpi) with histogram, density overlay, and mean/median dashed lines.

## Run

From project root:

```bash
conda run -n bch709-r Rscript scripts/mrna_gc_content.R
```
