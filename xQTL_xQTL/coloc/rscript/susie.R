# Assumed GWAS min columns: ["var_id_", "chrom", "position", "beta", "varbeta", "pvalue", "se"]
# Assumed eQTL min columns: ["var_id_", "chrom", "position", "beta", "varbeta", "pvalue", "se", "gene_id", "maf"]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments are rquired: output_file_path, gwas_path, qtl_path, gwas_sample_size, qtl_sample_size")
}
output_file_path = args[1]
gwas_path = '/home/users/nus/e1124850/scratch/lc3casestudy/processed/default/GCST90691928_obesity_europe/Muscle_Skeletal/EUR/eqtl/combined_ld/coloc/input/gwas_ENSG00000090857.14_chr16_69934824_A_G.tsv.gz'
qtl_path = '/home/users/nus/e1124850/scratch/lc3casestudy/processed/default/GCST90691928_obesity_europe/Muscle_Skeletal/EUR/eqtl/combined_ld/coloc/input/qtl_ENSG00000090857.14_chr16_69934824_A_G.tsv.gz'
LD_path = '/home/users/nus/e1124850/scratch/lc3casestudy/processed/default/GCST90691928_obesity_europe/Muscle_Skeletal/EUR/eqtl/combined_ld/susie/input/ENSG00000090857.14_chr16_69934824_A_G.ld'
gwas_path = args[2]
qtl_path = args[3]
LD_path = args[4]
gwas_lead_snp = args[5]
gwas_sample_size = as.integer(args[6])
qtl_sample_size = as.integer(args[7])
gwas_type = args[8]
eqtl_type = args[9]
tp1 = as.numeric(args[10]) # prior probability a SNP is associated with trait 1, default 1e-4
tp2 = as.numeric(args[11]) # prior probability a SNP is associated with trait 2, default 1e-4
tp12 = as.numeric(args[12]) # prior probability a SNP is associated with both traits, default 1e-5
# The threshold to consider overall H4 is true, only when overall H4 is true, the SNP level H4 is considered as relevant, optional
overall_h4_threshold = as.numeric(args[13])

if (is.na(gwas_path) || !file.exists(gwas_path) || !file.info(gwas_path)$size > 0) {
  print("GWAS file does not exist or is empty!")
  q(save = "no")
}
if (is.na(qtl_path) || !file.exists(qtl_path) || !file.info(qtl_path)$size > 0) {
  print("eQTL file does not exist or is empty!")
  q(save = "no")
}
if (is.na(gwas_sample_size)) {
  stop("must supply argument gwas_sample_size")
}
if (is.na(qtl_sample_size)) {
  stop("must supply argument qtl_sample_size")
}
gwas_df = read.table(file = gwas_path, header = T)
eqtl_df = read.table(file = qtl_path, header = T)
if (is.na(gwas_type) || tolower(gwas_type) == 'na' || tolower(gwas_type) == 'none') {
  gwas_type = "cc"
}
if (is.na(eqtl_type) || tolower(eqtl_type) == 'na' || tolower(eqtl_type) == 'none') {
  eqtl_type = "quant"
}
if (is.na(overall_h4_threshold)) {
  overall_h4_threshold = 0
}

if (!require(susieR)) {
  stop("coloc not installed")
}

eqtl_df = eqtl_df[match(gwas_df$var_id_, eqtl_df$var_id_),]

input = merge(gwas_df, eqtl_df, by = "var_id_", all = FALSE, suffixes = c("_gwas", "_eqtl"))
print("input.shape")
print(dim(input))

if (nrow(input) == 0) {
  print("No common snps found in two input dataframes")
  q(save = "no")
}

LD = read.table(LD_path, row.names = 1)
colnames(LD) = row.names(LD)

print("Start susie")
d1 = list(snp = input$var_id_, beta = input$beta_gwas, varbeta = input$varbeta_gwas, position = input$position_gwas, type = gwas_type, N = gwas_sample_size, MAF = as.numeric(input$maf), LD = as.matrix(LD))
d2 = list(snp = input$var_id_, beta = input$beta_eqtl, varbeta = input$varbeta_eqtl, position = input$position_eqtl, type = eqtl_type, N = qtl_sample_size, MAF = as.numeric(input$maf), LD = as.matrix(LD))

# check_prior = FALSE may cause susieR not converge and then keep calculating
s1 = runsusie(d1, check_prior = FALSE, estimate_residual_variance = FALSE, n=gwas_sample_size)
s2 = runsusie(d2, check_prior = FALSE, estimate_residual_variance = FALSE, n=qtl_sample_size)

result = coloc.susie(s1, s2)

# result = susie_rss(bhat = input$beta_gwas,shat = input$se_gwas,R = as.matrix(LD),n=gwas_sample_size, L = 10, estimate_residual_variance = F,coverage = 0.95,niter=10000)

pip_df <- data.frame(result$pip)
cs_index <- result$sets$cs_index

# 初始化：为每个 credible set 创建一列
for (i in seq_along(cs_index)) {
  cs_name <- paste0("cs", i)  # 列名，如 cs1, cs2, ...
  cs_label <- paste0("L", cs_index[i])  # 对应 result$sets$cs$Lx
  
  # 检查这个 credible set 是否存在
  if (!is.null(result$sets$cs[[cs_label]])) {
    snp_idx <- result$sets$cs[[cs_label]]
    pip_df[[cs_name]] <- FALSE
    pip_df[[cs_name]][snp_idx] <- TRUE
  } else {
    # 如果没有对应 Lx，列设为 FALSE
    pip_df[[cs_name]] <- FALSE
  }
}

write.table(pip_df, if (endsWith(output_file_path, ".gz")) gzfile(output_file_path) else output_file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
