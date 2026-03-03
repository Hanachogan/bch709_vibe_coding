#!/usr/bin/env Rscript
# Prompt:
# TASK: GC content analysis for UCSC yeast mRNA FASTA (Saccharomyces cerevisiae)
# INPUT: mrna.fa.gz FASTA records; header accession is first token after ">".
# REQUIREMENTS: parse gz FASTA directly, compute accession/length/gc_content,
# output sorted TSV (gc_content desc, 4 decimals) to results/mrna_metrics.tsv,
# and save histogram+density plot with mean/median dashed lines and caption stats
# (n, mean, median, sd) to results/gc_content_distribution.png at 1600x900 px, 200 dpi.
#
# Interpretation:
# This script treats sequence length as all characters on sequence lines after
# whitespace removal, and computes GC as (count(G)+count(C))/length using
# uppercase normalized sequence text. Results are written under results/.

suppressPackageStartupMessages({
  library(ggplot2)
})

results_dir <- file.path("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

input_candidates <- c("mrna.fa.gz", file.path("data", "mrna.fa.gz"))
input_path <- input_candidates[file.exists(input_candidates)][1]
if (is.na(input_path)) {
  stop("Could not find mrna.fa.gz in project root or data/ directory.")
}

parse_fasta_metrics <- function(fasta_gz_path) {
  con <- gzfile(fasta_gz_path, open = "rt")
  on.exit(close(con), add = TRUE)

  accessions <- character()
  lengths <- integer()
  gc_contents <- numeric()

  current_accession <- NULL
  current_length <- 0L
  current_gc <- 0L

  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    for (line in lines) {
      if (!nzchar(line)) {
        next
      }

      if (startsWith(line, ">")) {
        if (!is.null(current_accession)) {
          accessions <- c(accessions, current_accession)
          lengths <- c(lengths, current_length)
          gc_contents <- c(
            gc_contents,
            if (current_length > 0L) current_gc / current_length else NA_real_
          )
        }

        current_accession <- sub("^>(\\S+).*$", "\\1", line)
        current_length <- 0L
        current_gc <- 0L
      } else {
        seq_line <- toupper(gsub("\\s+", "", line))
        line_length <- nchar(seq_line)
        current_length <- current_length + line_length

        if (line_length > 0L) {
          gc_count_line <- nchar(gsub("[^GC]", "", seq_line))
          current_gc <- current_gc + gc_count_line
        }
      }
    }
  }

  if (!is.null(current_accession)) {
    accessions <- c(accessions, current_accession)
    lengths <- c(lengths, current_length)
    gc_contents <- c(
      gc_contents,
      if (current_length > 0L) current_gc / current_length else NA_real_
    )
  }

  data.frame(
    accession = accessions,
    length = lengths,
    gc_content = gc_contents,
    stringsAsFactors = FALSE
  )
}

metrics <- parse_fasta_metrics(input_path)

if (nrow(metrics) == 0L) {
  stop("No FASTA records found in input file.")
}

metrics <- metrics[order(metrics$gc_content, decreasing = TRUE, na.last = TRUE), ]

metrics_out <- metrics
metrics_out$gc_content <- sprintf("%.4f", round(metrics_out$gc_content, 4))

write.table(
  metrics_out,
  file = file.path(results_dir, "mrna_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

n <- nrow(metrics)
mean_gc <- mean(metrics$gc_content, na.rm = TRUE)
median_gc <- median(metrics$gc_content, na.rm = TRUE)
sd_gc <- sd(metrics$gc_content, na.rm = TRUE)

caption_txt <- sprintf(
  "n = %d | mean = %.4f | median = %.4f | sd = %.4f",
  n,
  mean_gc,
  median_gc,
  sd_gc
)

p <- ggplot(metrics, aes(x = gc_content)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#5DA5DA", color = "white", alpha = 0.85) +
  geom_density(color = "#D62728", linewidth = 1.1, na.rm = TRUE) +
  geom_vline(xintercept = mean_gc, linetype = "dashed", color = "#1B9E77", linewidth = 0.9) +
  geom_vline(xintercept = median_gc, linetype = "dashed", color = "#E7298A", linewidth = 0.9) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "GC Content Distribution for Saccharomyces cerevisiae mRNA",
    x = "GC content",
    y = "Density",
    caption = caption_txt
  ) +
  theme_minimal(base_size = 13)

png(
  filename = file.path(results_dir, "gc_content_distribution.png"),
  width = 1600,
  height = 900,
  res = 200
)
print(p)
dev.off()

cat(sprintf("Wrote metrics: %s\n", file.path(results_dir, "mrna_metrics.tsv")))
cat(sprintf("Wrote plot: %s\n", file.path(results_dir, "gc_content_distribution.png")))
