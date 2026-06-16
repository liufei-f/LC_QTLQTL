
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments are rquired: output_file_path, gwas_path, qtl_path, gwas_sample_size, qtl_sample_size")
}
output_file_path = args[1]
qtl1_file = args[2]
qtl2_file = args[3]
qtl1_sample_size = as.integer(args[4])
qtl2_sample_size = as.integer(args[5])
qtl1_type = args[6]
qtl2_type = args[7]
qtl1_phenotype_id = args[8]
qtl2_phenotype_id = args[9]
tp1 = as.numeric(args[10]) # prior probability a SNP is associated with trait 1, default 1e-4
tp2 = as.numeric(args[12]) # prior probability a SNP is associated with trait 2, default 1e-4
tp12 = as.numeric(args[13]) # prior probability a SNP is associated with both traits, default 1e-5
LD_path   <- args[14]

# LD_path   <- "/Users/phoebel/Downloads/LD4602.colocboost.tsv"
# eqtl_file <- "/Users/phoebel/Downloads/ENSG00000113504.20.eqtl4602.colocboost.tsv"
# gwas_file <- "/Users/phoebel/Downloads/ENSG00000113504.20.gwas4602.colocboost.tsv"

suppressMessages(library(data.table))
suppressMessages(library(susieR))
suppressMessages(library(coloc))

# ---------- load summary stats ----------
qtl1 <- fread(qtl1_file)
qtl2 <- fread(qtl2_file)

# ---------- load LD matrix ----------
LD_raw <- fread(LD_path)
# first column may be row-name variants or a numeric index
if (is.character(LD_raw[[1]])) {
    variants_LD <- LD_raw[[1]]
    LD_mat      <- as.matrix(LD_raw[, -1, with = FALSE])
} else {
    variants_LD <- colnames(LD_raw)
    LD_mat      <- as.matrix(LD_raw)
}
rownames(LD_mat) <- variants_LD
colnames(LD_mat) <- variants_LD

# remove rows/columns with all NA (e.g. due to missing variants in LD reference)
na_rows <- rowSums(is.na(LD_mat)) == ncol(LD_mat)
LD_mat <- LD_mat[!na_rows, !na_rows]
variants_LD <- variants_LD[!na_rows]

# ---------- align to common variants ----------
common_vars <- Reduce(intersect, list(qtl1$variant, qtl2$variant, rownames(LD_mat)))
if (length(common_vars) == 0) stop("No overlapping variants across QTL1, QTL2, and LD.")

qtl1_sub   <- qtl1[match(common_vars, variant)]
qtl2_sub   <- qtl2[match(common_vars, variant)]
LD_aligned <- LD_mat[common_vars, common_vars]
message(sprintf("Variants: QTL1=%d  QTL2=%d  LD=%d  common=%d",
                nrow(qtl1_sub), nrow(qtl2_sub), nrow(LD_aligned), length(common_vars)))
# ---------- SuSiE-RSS: QTL1 ----------
n_qtl1 <- qtl1_sub$n[1]
message(sprintf("Running susie_rss for QTL1  (n=%d, p=%d)", n_qtl1, length(common_vars)))
res_qtl1 <- susie_rss(qtl1_sub$z, LD_aligned, n = n_qtl1)

# ---------- SuSiE-RSS: QTL2 ----------
n_qtl2 <- qtl2_sub$n[1]
message(sprintf("Running susie_rss for QTL2  (n=%d, p=%d)", n_qtl2, length(common_vars)))
res_qtl2 <- susie_rss(qtl2_sub$z, LD_aligned, n = n_qtl2)

# ---------- colocalization ----------
n_cs_qtl1 <- length(res_qtl1$sets$cs)
n_cs_qtl2 <- length(res_qtl2$sets$cs)
message(sprintf("Credible sets: QTL1=%d  QTL2=%d", n_cs_qtl1, n_cs_qtl2))

if (n_cs_qtl1 == 0 || n_cs_qtl2 == 0) {
    message("No credible sets in one or both traits — skipping coloc.")
} else {
    res_coloc <- coloc.susie(res_qtl1, res_qtl2)
    results_df <- as.data.frame(res_coloc$summary)
    results_df$qtl1_phenotype_id <- qtl1_phenotype_id
    results_df$qtl2_phenotype_id <- qtl2_phenotype_id

    write.table(results_df, output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
    # print(res_coloc$summary)
}


