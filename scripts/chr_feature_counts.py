#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GFF_GZ = ROOT / "data" / "saccharomyces_cerevisiae.gff.gz"
DEFAULT_CHROM_SIZES = ROOT / "data" / "chrom.sizes"
DEFAULT_OUT_COUNTS = ROOT / "results" / "chr_feature_counts.tsv"
DEFAULT_OUT_DROPPED = ROOT / "results" / "dropped_seqids.txt"


def load_chrom_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    chrom_order = []
    chrom_lengths = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, length_bp = line.split("\t")
            chrom_order.append(chrom)
            chrom_lengths[chrom] = int(length_bp)
    return chrom_order, chrom_lengths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count chromosome-level gene/exon/tRNA/snoRNA features from a GFF3."
    )
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF_GZ, help="Input GFF3.gz")
    parser.add_argument(
        "--chrom-sizes",
        type=Path,
        default=DEFAULT_CHROM_SIZES,
        help="Input TSV of chrom and length_bp",
    )
    parser.add_argument(
        "--out-counts",
        type=Path,
        default=DEFAULT_OUT_COUNTS,
        help="Output TSV table path",
    )
    parser.add_argument(
        "--out-dropped",
        type=Path,
        default=DEFAULT_OUT_DROPPED,
        help="Output dropped seqids path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    chrom_order, chrom_lengths = load_chrom_sizes(args.chrom_sizes)
    valid_chroms = set(chrom_order)

    n_gene = defaultdict(int)
    n_tRNA = defaultdict(int)
    n_snoRNA = defaultdict(int)
    exon_intervals = defaultdict(set)
    dropped_seqids = set()
    excluded_feature_lines = 0

    with gzip.open(args.gff, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 9:
                continue

            seqid, _, feature_type, start, end, _, strand, _, _ = row

            if seqid not in valid_chroms:
                dropped_seqids.add(seqid)
                excluded_feature_lines += 1
                continue

            if feature_type == "gene":
                n_gene[seqid] += 1
            elif feature_type == "exon":
                exon_intervals[seqid].add((int(start), int(end), strand))
            elif feature_type == "tRNA":
                n_tRNA[seqid] += 1
            elif feature_type == "snoRNA":
                n_snoRNA[seqid] += 1

    args.out_counts.parent.mkdir(parents=True, exist_ok=True)
    args.out_dropped.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for chrom in chrom_order:
        chrom_length_bp = chrom_lengths[chrom]
        chrom_length_mb = chrom_length_bp / 1_000_000
        n_exon_unique = len(exon_intervals[chrom])
        rows.append(
            {
                "chrom": chrom,
                "chrom_length_bp": chrom_length_bp,
                "n_gene": n_gene[chrom],
                "n_exon_unique": n_exon_unique,
                "n_tRNA": n_tRNA[chrom],
                "n_snoRNA": n_snoRNA[chrom],
                "gene_per_Mb": round(n_gene[chrom] / chrom_length_mb, 4),
                "exon_unique_per_Mb": round(n_exon_unique / chrom_length_mb, 4),
                "tRNA_per_Mb": round(n_tRNA[chrom] / chrom_length_mb, 4),
                "snoRNA_per_Mb": round(n_snoRNA[chrom] / chrom_length_mb, 4),
            }
        )

    rows.sort(key=lambda r: r["gene_per_Mb"], reverse=True)

    # Keep 4 decimal places in output file.
    for row in rows:
        row["gene_per_Mb"] = f"{row['gene_per_Mb']:.4f}"
        row["exon_unique_per_Mb"] = f"{row['exon_unique_per_Mb']:.4f}"
        row["tRNA_per_Mb"] = f"{row['tRNA_per_Mb']:.4f}"
        row["snoRNA_per_Mb"] = f"{row['snoRNA_per_Mb']:.4f}"

    with args.out_counts.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "chrom",
                "chrom_length_bp",
                "n_gene",
                "n_exon_unique",
                "n_tRNA",
                "n_snoRNA",
                "gene_per_Mb",
                "exon_unique_per_Mb",
                "tRNA_per_Mb",
                "snoRNA_per_Mb",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    with args.out_dropped.open("w", encoding="utf-8") as handle:
        for seqid in sorted(dropped_seqids):
            handle.write(f"{seqid}\n")

    print(f"number_of_dropped_seqids\t{len(dropped_seqids)}")
    print(f"number_of_excluded_feature_lines\t{excluded_feature_lines}")
    print(
        "chrom\tchrom_length_bp\tn_gene\tn_exon_unique\tn_tRNA\tn_snoRNA\t"
        "gene_per_Mb\texon_unique_per_Mb\ttRNA_per_Mb\tsnoRNA_per_Mb"
    )
    for row in rows[:5]:
        print(
            f"{row['chrom']}\t{row['chrom_length_bp']}\t{row['n_gene']}\t"
            f"{row['n_exon_unique']}\t{row['n_tRNA']}\t{row['n_snoRNA']}\t"
            f"{row['gene_per_Mb']}\t{row['exon_unique_per_Mb']}\t"
            f"{row['tRNA_per_Mb']}\t{row['snoRNA_per_Mb']}"
        )


if __name__ == "__main__":
    main()
