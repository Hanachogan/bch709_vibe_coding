#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[[1]] else "data/gasch2000.txt"
output_file <- if (length(args) >= 2) args[[2]] else "data/gasch2000_top10_heatmap.png"
n_top <- if (length(args) >= 3) as.integer(args[[3]]) else 10
n_conditions <- if (length(args) >= 4) as.integer(args[[4]]) else 30

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  comment.char = "#",
  quote = "\"",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(dat) < 2) {
  stop("Input table has too few columns to contain gene IDs and expression values.")
}

gene_col <- colnames(dat)[1]

meta_cols <- c(gene_col, "NAME", "GWEIGHT", "EWEIGHT", "GORDER")
candidate_cols <- setdiff(colnames(dat), meta_cols)

if (length(candidate_cols) == 0) {
  stop("No candidate expression columns found after excluding metadata columns.")
}

numeric_fraction <- vapply(candidate_cols, function(nm) {
  vals <- suppressWarnings(as.numeric(dat[[nm]]))
  mean(!is.na(vals))
}, numeric(1))

condition_cols <- names(numeric_fraction[numeric_fraction >= 0.5])

if (length(condition_cols) == 0) {
  stop("No expression columns look numeric (>=50% parseable values).")
}

if (n_conditions < 1) {
  stop("n_conditions must be >= 1.")
}
condition_cols <- condition_cols[seq_len(min(n_conditions, length(condition_cols)))]

gene_ids <- trimws(as.character(dat[[gene_col]]))
valid_gene_rows <- which(!is.na(gene_ids) & nzchar(gene_ids))

if (length(valid_gene_rows) == 0) {
  stop(sprintf("No valid gene IDs found in first column '%s'.", gene_col))
}

expr <- dat[valid_gene_rows, condition_cols, drop = FALSE]
expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(x)))
expr_mat <- as.matrix(expr)
rownames(expr_mat) <- gene_ids[valid_gene_rows]

has_negative <- any(expr_mat < 0, na.rm = TRUE)
if (!has_negative) {
  expr_mat <- log2(expr_mat + 1)
}

gene_mean <- rowMeans(expr_mat, na.rm = TRUE)
gene_sd <- apply(expr_mat, 1, sd, na.rm = TRUE)
denom <- abs(gene_mean)
denom[denom == 0] <- NA_real_
cv <- gene_sd / denom
cv[is.infinite(cv)] <- NA_real_

if (all(is.na(cv))) {
  stop("All CV values are NA after guarding for near-zero means.")
}

n_top <- max(1, min(n_top, sum(!is.na(cv))))
top_genes <- names(sort(cv, decreasing = TRUE, na.last = NA))[seq_len(n_top)]

plot_mat <- expr_mat[top_genes, , drop = FALSE]
plot_df <- as.data.frame(as.table(plot_mat), stringsAsFactors = FALSE)
colnames(plot_df) <- c("Gene", "Condition", "Expression")

plot_df$Gene <- factor(plot_df$Gene, levels = top_genes)
plot_df$Condition <- factor(plot_df$Condition, levels = rev(colnames(plot_mat)))

p <- ggplot(plot_df, aes(x = Gene, y = Condition, fill = Expression)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_distiller(palette = "PuGn", direction = 1, na.value = "grey90") +
  labs(title = "gene top 10", x = "gene", y = "conditions") +
  theme_minimal(base_family = "Times New Roman", base_size = 11) +
  guides(fill = "none") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = output_file,
  plot = p,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(sprintf("Input file: %s", input_file))
message(sprintf("Gene ID column: %s", gene_col))
message(sprintf("Unique/total genes: %d/%d", length(unique(rownames(expr_mat))), nrow(expr_mat)))
message(sprintf("Condition columns retained: %d", ncol(expr_mat)))
message(sprintf("Missing values in expression matrix: %d", sum(is.na(expr_mat))))
if (has_negative) {
  message("Log-scale check: negative values detected -> treated as already log-scale (no transform).")
} else {
  message("Log-scale check: no negatives detected -> applied log2(x + 1).")
}
message(sprintf("Saved heatmap to %s", output_file))
