# bch709_vibe_coding

This README is dedicated to this repository only.

## Project scope

Yeast genomics workflows in this repo currently include:

- Chromosome-level feature counting from GFF (`scripts/chr_feature_counts.py`)
- mRNA GC content analysis from FASTA (`scripts/mrna_gc_content.R`)

## Environment

Create/update environment from project root:

```bash
conda env create -f environment.yml
# or
conda env update -f environment.yml --prune
```

Run commands with the environment active:

```bash
conda activate bch709-r
```

## Inputs

- `data/saccharomyces_cerevisiae.gff.gz`
- `data/chrom.sizes`
- `data/mrna.fa.gz`
- `data/toy/toy.gff.gz`
- `data/toy/chrom.sizes`

## Workflows

Chromosome feature counts:

```bash
python scripts/chr_feature_counts.py
```

Toy run for chromosome feature counts:

```bash
python scripts/chr_feature_counts.py \
  --gff data/toy/toy.gff.gz \
  --chrom-sizes data/toy/chrom.sizes \
  --out-counts results/toy/chr_feature_counts.tsv \
  --out-dropped results/toy/dropped_seqids.txt
```

mRNA GC content analysis:

```bash
conda run -n bch709-r Rscript scripts/mrna_gc_content.R
```

## Outputs

- `results/chr_feature_counts.tsv`
- `results/dropped_seqids.txt`
- `results/toy/chr_feature_counts.tsv`
- `results/toy/dropped_seqids.txt`
- `results/mrna_metrics.tsv`
- `results/gc_content_distribution.png`

## Notes

The GC workflow computes per-accession GC fraction `(G + C) / length`, sorts records by GC descending, writes a TSV summary, and produces a histogram+density plot with mean and median markers.
