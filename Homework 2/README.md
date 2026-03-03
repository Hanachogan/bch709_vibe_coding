# Homework 2: Yeast Stress Response (Gasch et al. 2000)

## Project Description

This homework analyzes yeast stress-expression log2 ratios from Gasch et al. (2000) in two stages:

1. Rank genes by coefficient of variation (CV) across stress conditions and select top 200 genes.
2. Build a clustered heatmap for those 200 genes using row-wise Z-scores and k=4 hierarchical clusters.

## Inputs

- Expression matrix:
  - `Homework 2/data/gasch2000.txt`
- Top-200 list (generated in Stage 1):
  - `Homework 2/results/yeast_stress_cv_top200.tsv`

## Scripts

- Stage 1 (CV ranking + heatmap):
  - `Homework 2/scripts/homework2_cv_heatmap.R`
- Stage 2 (row Z-score clustering, k=4):
  - `Homework 2/scripts/cv_top200_cluster_heatmap.R`

## Outputs

- Stage 1 outputs:
  - `Homework 2/results/yeast_stress_cv_top200.tsv`
  - `Homework 2/results/yeast_stress_cv_top200_heatmap.png`
- Stage 2 outputs:
  - `Homework 2/results/cluster_assignment.tsv`
  - `Homework 2/results/cv_top200_cluster_heatmap.pdf`

## Environment

Designed for conda environment with these R packages:

- `r-base`
- `r-data.table`
- `r-ggplot2`
- `r-pheatmap`
- `r-viridislite`
- `r-scales`

Environment file:

- `environment-HW2.yml`

## How to Run

From project root:

```bash
conda run -n bch709_vibe_coding Rscript "Homework 2/scripts/homework2_cv_heatmap.R"
conda run -n bch709_vibe_coding Rscript "Homework 2/scripts/cv_top200_cluster_heatmap.R"
```

If your local environment name is different (for example `bch709_vibe_coding_hw2`), replace the env name in the commands.

## Method Summary

1. Read Gasch matrix and remove metadata columns (`NAME`, `description`, `GWEIGHT`) for numeric-condition analyses.
2. Compute per-gene mean, SD, and CV; filter invalid rows; select top 200 by descending CV.
3. Save top-200 table and render non-clustered top-200 heatmap.
4. For clustering, use the top-200 genes and first 30 condition columns.
5. Row-normalize by Z-score: `(x - row_mean) / row_sd` with safe handling for `row_sd = 0`.
6. Cluster rows using Euclidean distance + `hclust(method = "ward.D2")`.
7. Cut dendrogram into 4 clusters (`cutree(k=4)`), save cluster assignments, and plot annotated heatmap PDF.
