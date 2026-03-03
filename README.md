# bch709_vibe_coding

## Installation

From the project root:

```bash
cd /home/hanachogan/bch709_vibe_coding
conda env create -f environment.yml
```

If the environment already exists:

```bash
conda env update -f environment.yml --prune
```

Activate:

```bash
conda activate bch709-r
```

## Overview

Project workspace for yeast genome/annotation analysis and plotting.

## Verify packages

```bash
python -c "import pandas, numpy, matplotlib, seaborn, Bio, tqdm; print('python ok')"
R -q -e "library(ggplot2); cat('R ok\n')"
```

## Run analyses

Generate the heatmap:

```bash
Rscript data/make_gasch_heatmap.R
```

Count chromosome features:

```bash
python scripts/chr_feature_counts.py
```

## Toy test

Toy inputs are available in `data/toy/`:

- `data/toy/chrom.sizes`
- `data/toy/toy.gff.gz`

Run the feature counter on toy data:

```bash
python scripts/chr_feature_counts.py \
  --gff data/toy/toy.gff.gz \
  --chrom-sizes data/toy/chrom.sizes \
  --out-counts results/toy/chr_feature_counts.tsv \
  --out-dropped results/toy/dropped_seqids.txt
```

## Inputs

- `data/saccharomyces_cerevisiae.gff.gz`
- `data/chrom.sizes`
- `data/gasch2000.txt`

## Outputs

- `data/gasch2000_top10_heatmap.png`
- `results/chr_feature_counts.tsv`
- `results/dropped_seqids.txt`
