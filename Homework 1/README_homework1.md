# Homework 1: GC Content Analysis of Yeast mRNA

## Project description

This project analyzes GC content for UCSC yeast (*Saccharomyces cerevisiae*) mRNA FASTA records.
The script reads a gzipped FASTA file directly (no manual unzip), computes per-record sequence metrics, writes a sorted summary table, and generates a GC-content distribution plot.

## Input

- Primary input file:
  - `data/mrna.fa.gz`
- File format:
  - FASTA
  - Header example: `>BC001547 /gb=BC001547 /gi=12654078 /ug=Sc.3456 /len=1254`
  - Sequence is one or more lines after each header

## Output

- Table output:
  - `Homework 1/results/mrna_metrics.tsv`
  - Columns:
    - `accession` (first token after `>` in FASTA header)
    - `length` (total bases across sequence lines)
    - `gc_content` (`(G + C) / length`, rounded to 4 decimals)
  - Sorted by `gc_content` in descending order

- Plot output:
  - `Homework 1/results/gc_content_distribution.png`
  - Specifications:
    - 1600 x 900 pixels
    - 200 dpi
    - Histogram of GC content (x-axis from 0 to 1)
    - Density curve overlay
    - Dashed vertical lines for mean and median
    - Caption with `n`, `mean`, `median`, and `sd` (4 decimals)

## Script

- Script file:
  - `Homework 1/scripts/mrna_gc_content.R`

## Packages and environment details

This work is designed for the `bch709-r` conda environment.

- R package used by this homework script:
  - `ggplot2`

- Environment file (`environment.yml`) also includes:
  - `r-base`
  - `python`
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `seaborn`
  - `biopython`
  - `tqdm`

Only `r-base` and `r-ggplot2` are required to run this Homework 1 GC script.

## How to run

From project root:

```bash
conda run -n bch709-r Rscript "Homework 1/scripts/mrna_gc_content.R"
```

## Method details

1. Open `mrna.fa.gz` with `gzfile()` in text mode.
2. Parse FASTA records by detecting header lines beginning with `>`.
3. Extract accession from header using the first token.
4. Concatenate sequence length logically by summing line lengths after whitespace removal.
5. Count GC bases (`G` and `C`) on uppercase sequence text.
6. Compute per-record GC fraction and sort descending.
7. Write TSV and generate PNG plot with summary statistics.

## Interpretation Q&A

Question: "what is the GC content distribution of yeast mRNA sequences, and are there distinct GC-content subpopulations?''

Answer: "Yeast mRNA GC content is concentrated around ~0.40–0.42 (mean 0.4167, median 0.4011; n=474), producing a largely unimodal distribution. The density curve shows a dominant peak near ~0.39–0.41 with a right-skewed tail toward higher GC values (~0.55–0.70). While there is a slight shoulder at higher GC, the distribution does not show clearly separated peaks, so there is no strong evidence for distinct GC-content subpopulations—rather a main group with a high-GC tail.''
