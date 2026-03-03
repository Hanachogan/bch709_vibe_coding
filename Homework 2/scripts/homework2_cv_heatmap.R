library(data.table)
library(pheatmap)

input_candidates <- c(
  "Homework 2/data/gasch2000.txt",
  "Homework2/gasch2000.txt",
  "data/gasch2000.txt",
  "gasch2000.txt"
)
input_file <- input_candidates[file.exists(input_candidates)][1]
output_tsv <- "Homework 2/results/yeast_stress_cv_top200.tsv"
output_png <- "Homework 2/results/yeast_stress_cv_top200_heatmap.png"

if (is.na(input_file)) {
  stop("Input file not found. Checked: Homework 2/data/gasch2000.txt, Homework2/gasch2000.txt, data/gasch2000.txt, gasch2000.txt")
}

dir.create("Homework 2/results", showWarnings = FALSE, recursive = TRUE)

# Read full expression matrix.
dt <- fread(input_file, sep = "\t", header = TRUE, data.table = TRUE)

if (!("UID" %in% names(dt))) {
  stop("Expected first column UID was not found in input file.")
}

# Keep UID as gene identifier and remove annotation columns if present.
skip_cols <- c("name", "description", "gweight")
condition_cols <- names(dt)[
  !(tolower(names(dt)) %in% c("uid", skip_cols))
]

if (length(condition_cols) == 0) {
  stop("No condition columns found after removing annotation columns.")
}

expr <- as.matrix(dt[, ..condition_cols])
mode(expr) <- "numeric"

mean_expr <- rowMeans(expr, na.rm = TRUE)
sd_expr <- apply(expr, 1, sd, na.rm = TRUE)
cv <- sd_expr / abs(mean_expr)

summary_dt <- data.table(
  gene_id = dt[["UID"]],
  mean_expr = mean_expr,
  sd_expr = sd_expr,
  cv = cv
)

# Remove rows with all NA means or zero/undefined variance.
summary_dt <- summary_dt[!is.nan(mean_expr) & !is.na(sd_expr) & sd_expr > 0]

setorder(summary_dt, -cv)
top200 <- summary_dt[1:min(200L, .N)]

# Round required columns to 4 decimals for output.
out_dt <- copy(top200)
out_dt[, `:=`(
  mean_expr = round(mean_expr, 4),
  sd_expr = round(sd_expr, 4),
  cv = round(cv, 4)
)]

fwrite(out_dt[, .(gene_id, mean_expr, sd_expr, cv)], output_tsv, sep = "\t")

cat("Top 10 genes by CV:\n")
print(out_dt[1:min(10L, .N), .(gene_id, mean_expr, sd_expr, cv)])

# Heatmap matrix in CV-descending row order and original condition order.
row_idx <- match(top200$gene_id, dt[["UID"]])
heatmap_mat <- expr[row_idx, , drop = FALSE]
rownames(heatmap_mat) <- top200$gene_id

png(filename = output_png, width = 1800, height = 1200, units = "px", res = 200)
pheatmap(
  mat = heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  fontsize_row = 4,
  angle_col = 90,
  main = "Yeast stress response, CV top200 (Gasch et al. 2000)"
)
dev.off()

cat(sprintf("\nSaved table: %s\n", output_tsv))
cat(sprintf("Saved heatmap: %s\n", output_png))
