library(data.table)
library(pheatmap)
library(viridisLite)

# Inputs (run from project root)
top200_candidates <- c(
  "Homework 2/results/yeast_stress_cv_top200.tsv",
  "results/yeast_stress_cv_top200.tsv"
)
expr_candidates <- c(
  "Homework 2/data/gasch2000.txt",
  "data/gasch2000.txt",
  "gasch2000.txt"
)

in_top200 <- top200_candidates[file.exists(top200_candidates)][1]
in_expr <- expr_candidates[file.exists(expr_candidates)][1]

if (is.na(in_top200)) {
  stop("Missing top-200 list. Expected: results/yeast_stress_cv_top200.tsv")
}
if (is.na(in_expr)) {
  stop("Missing expression matrix. Expected: data/gasch2000.txt")
}

dir.create("Homework 2/results", showWarnings = FALSE, recursive = TRUE)
out_pdf <- "Homework 2/results/cv_top200_cluster_heatmap.pdf"
out_tsv <- "Homework 2/results/cluster_assignment.tsv"

# 1) Read top200 gene IDs
top200_dt <- fread(in_top200)
if (!("gene_id" %in% names(top200_dt))) {
  stop("Input top200 file must contain a gene_id column")
}
gene_ids <- unique(top200_dt$gene_id)

# 2) Read gasch data
gasch <- fread(in_expr, sep = "\t", header = TRUE, data.table = TRUE)

# 3) Detect gene identifier column (IDs like YAL001C), then keep only top200 rows
id_pattern <- "^Y[A-Z]{2}[0-9]{3}[CW](?:-[A-Z])?$"
candidate_cols <- intersect(c("NAME", "UID", "gene_id", "ORF"), names(gasch))
if (length(candidate_cols) == 0) {
  stop("No candidate gene ID columns found (expected one of NAME/UID/gene_id/ORF)")
}

match_counts <- sapply(candidate_cols, function(col) {
  vals <- as.character(gasch[[col]])
  sum(grepl(id_pattern, vals), na.rm = TRUE)
})
id_col <- candidate_cols[which.max(match_counts)]

subset_dt <- gasch[get(id_col) %in% gene_ids]
if (nrow(subset_dt) == 0) {
  stop("No matching genes found between top200 list and expression data")
}

# Keep rows in top200 order as much as possible
setkeyv(subset_dt, id_col)
subset_dt <- subset_dt[J(gene_ids), nomatch = 0]

# 4) Remove metadata columns and keep only condition columns
meta_cols <- intersect(c("NAME", "description", "GWEIGHT"), names(subset_dt))
meta_idx <- match(meta_cols, names(subset_dt))
meta_idx <- meta_idx[!is.na(meta_idx)]
if (length(meta_idx) == 0) {
  stop("Could not find metadata columns to define condition block")
}

condition_cols <- names(subset_dt)[(max(meta_idx) + 1):ncol(subset_dt)]
if (length(condition_cols) < 30) {
  stop(sprintf("Need at least 30 condition columns; found %d", length(condition_cols)))
}

# Ensure numeric condition matrix
mat_all <- as.matrix(subset_dt[, ..condition_cols])
mode(mat_all) <- "numeric"

# 5) First 30 conditions only (original order)
mat30 <- mat_all[, 1:30, drop = FALSE]

# Row names are gene IDs (from detected ID column)
row_ids <- as.character(subset_dt[[id_col]])
rownames(mat30) <- row_ids

# 6) Row-wise Z-score: (x - mean) / sd ; handle sd=0 or NA safely
row_means <- rowMeans(mat30, na.rm = TRUE)
row_sds <- apply(mat30, 1, sd, na.rm = TRUE)

z_mat <- mat30
for (i in seq_len(nrow(mat30))) {
  sdi <- row_sds[i]
  if (is.na(sdi) || sdi == 0) {
    z_mat[i, ] <- 0
  } else {
    z_mat[i, ] <- (mat30[i, ] - row_means[i]) / sdi
  }
}

# Replace any residual NA/Inf from missing values with 0 to keep clustering stable
z_mat[!is.finite(z_mat)] <- 0

# 7) Hierarchical clustering on rows (Euclidean + ward.D2)
d_rows <- dist(z_mat, method = "euclidean")
hc_rows <- hclust(d_rows, method = "ward.D2")

# 8) Cut into k=4 clusters
clusters <- cutree(hc_rows, k = 4)
clusters <- as.integer(clusters)
names(clusters) <- rownames(z_mat)

# Assignment table
gene_id_vec <- rownames(z_mat)
assign_dt <- data.table(
  gene_id = gene_id_vec,
  cluster = as.integer(clusters[gene_id_vec])
)
setorder(assign_dt, cluster, gene_id)
fwrite(assign_dt, out_tsv, sep = "\t", quote = FALSE)

# Heatmap row annotation
anno_row <- data.frame(cluster = factor(clusters[gene_id_vec], levels = 1:4))
rownames(anno_row) <- gene_id_vec
anno_colors <- list(cluster = setNames(viridis(4), as.character(1:4)))

# Reorder matrix to clustering leaf order for deterministic display
z_plot <- z_mat[hc_rows$order, , drop = FALSE]

pdf(out_pdf, width = 8, height = 12)
pheatmap(
  mat = z_plot,
  cluster_rows = hc_rows,
  cluster_cols = FALSE,
  annotation_row = anno_row,
  annotation_colors = anno_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "CV top200 yeast stress genes (row Z-score, k=4 clusters)"
)
dev.off()

cat("Wrote:\n")
cat(sprintf("- %s\n", out_pdf))
cat(sprintf("- %s\n", out_tsv))
