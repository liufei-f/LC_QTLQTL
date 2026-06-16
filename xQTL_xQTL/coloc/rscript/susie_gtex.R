
rm(list = ls())
maf_cut <- 0.05
missing_cut <- 0.1
GWAS_hit_cut <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop(paste(
        "Usage: Rscript susie_gtex.R <gene_ensg> <gtex_eqtl_expression_dir>",
        "<gtex_eqtl_covariates_dir> <gtex_genotype_dir> <gene_info_file>",
        "<gwas_file> <gwas_n>",
        "\n  gwas_file: TSV/gz with columns chr, pos, beta, se"
    ))
}
gene_ensg                 <- args[1]
gtex_eqtl_expression_dir <- args[2]
gtex_eqtl_covariates_dir <- args[3]
gtex_genotype_dir         <- args[4]
gene_info_file            <- args[5]
gwas_file                 <- args[6]
gwas_n                    <- as.numeric(args[7])
ld_meta_data <- args[8]

gene_ensg = 'ENSG00000273796.1' # for testing
gtex_eqtl_expression_dir= '/Users/phoebel/Downloads/GTEx_Analysis_v10_eQTL_expression_matrices'
gtex_eqtl_covariates_dir= '/Users/phoebel/Downloads/GTEx_Analysis_v10_eQTL_covariates'
gtex_genotype_dir= '/Users/phoebel/Downloads/Gtex_genotype'
gene_info_file = '/Users/phoebel/github/lc3private/gtf.v39.geneinfo.tsv'
gwas_file = '/Users/phoebel/preprocessed_diagram.mega-meta.chr21.txt.gz'
gwas_n    = 80000
ld_meta_data <- "/Users/phoebel/Downloads/LD_sketch/ld_meta_file_chr21.tsv"

suppressMessages(library(tidyverse))
suppressMessages(library(BEDMatrix))
suppressMessages(library(data.table))
suppressMessages(library(susieR))
suppressMessages(library(coloc))
suppressMessages(library(pecotmr))

# pecotmr bug: file_ext("ADSP.R4.EUR.chr21") returns "chr21" which is not a
# recognised format extension, so the switch hits stop() instead of falling
# through to the stem-file probe.  Patch: only use the switch for known exts.



getQTLByLocus <- function(gene_ensg, gene_info_file,
                          gtex_genotype_dir, gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir,
                          maf_cut, missing_cut,
                          gene_interval = NULL,
                          tissues_considered = NULL,
                          data.table = FALSE,
                          window = 1000000,
                          adjust_cov = TRUE){

    # given a gene locus,
    # get eqtl individual data, correcting for covariates
    # get genotype individual data for all snps within locus, correcting for covariates

    ################ gene information ################
    # gene_info_file: TSV with columns gene_id, gene_name, chr, start, end
    gene_info <- fread(gene_info_file)[gene_id == gene_ensg]
    if (nrow(gene_info) == 0) {
        message("Gene not found in gene info file: ", gene_ensg)
        return(NULL)
    }
    if (nrow(gene_info) > 1) { gene_info <- gene_info[1, ] }
    chr        <- as.numeric(gene_info$chr)
    gene.start <- pmax(0, gene_info$start - window)
    gene.end   <- gene_info$end + window
    if (!is.null(gene_interval)){
        message("Using the provided gene interval: ", gene_interval)
        region_info <- parse_region(gene_interval)
        chr = region_info$chrom
        gene.start = region_info$start
        gene.end = region_info$end
    }
    
    ################ preprocess QTL data ################
    message(paste("------ Processing eQTL genotype data for loci", gene_ensg, "at CHR", chr))
    if (!(chr%in%c(1:22))){ return(NULL) }
    chr_prefix <- paste0(gtex_genotype_dir, "/chr", chr, "/chr", chr)
    genotype <- BEDMatrix(chr_prefix)
    if (data.table){
        bim_data <- fread(paste0(chr_prefix, ".bim"), header = FALSE)
    } else {
        bim_data <- read.delim(paste0(chr_prefix, ".bim"), header = FALSE, stringsAsFactors = FALSE)
    }
    gene.snp <- bim_data %>% 
                    mutate(idx = row_number()) %>% 
                    filter(V4 >= gene.start & V4 <= gene.end)
    genotype <- data.matrix(genotype[, gene.snp$idx])
    # rownames(genotype) <- str_replace(rownames(genotype), "0_", "")
    rownames(genotype) <- str_extract(rownames(genotype), "^[^_]+")
    print(head(rownames(genotype)))
    colnames(genotype) <- paste0("chr", gene.snp$V1, ":", gene.snp$V4)

    
    message(paste("------ Processing eQTL expression data for loci", gene_ensg))
    eQTL.list <- list.files(gtex_eqtl_expression_dir, pattern = ".bed.gz$")
    if (!is.null(tissues_considered)){
        eQTL.list <- eQTL.list[str_detect(eQTL.list, paste(tissues_considered, collapse = "|"))]
    }
    X <- Y <- list()
    phenotypes <- c()
    for (tissue.file in eQTL.list){

        tissue.name <- str_extract(tissue.file, "^[^.]+") 
        message(paste("Start processing eQTL file for tissue", tissue.name))

        # - expression
        if (data.table){
            expression <- fread(paste0(gtex_eqtl_expression_dir, "/", tissue.file))
        } else {
            expression <- read.delim(paste0(gtex_eqtl_expression_dir, "/", tissue.file), stringsAsFactors = FALSE)
        }
        # expression$gene_id <- sub("\\..*", "", expression$gene_id)
        # gene_ensg_base <- sub("\\..*", "", gene_ensg)

        expression <- expression %>%
                        filter(gene_id == gene_ensg) %>%
                        mutate(start = pmax(0, start-window), end = end + window)
        if (nrow(expression) == 0){
            message(paste("No expression for tissue", tissue.name))
            next
        }
        gene.start <- expression$start
        gene.end <- expression$end
        eQTL <- expression[,-c(1:4)]
        sample_expression <- colnames(expression)[-c(1:4)]
        if (!data.table){
            sample_expression <- gsub("\\.", "-", sample_expression)
        }
        # - genotype
        # print(sample_expression)
        # print(head(genotype))
        print(paste("Number of intersecting samples:", length(intersect(rownames(genotype), sample_expression))))
        print(dim(genotype))
        print(dim(sample_expression))
        sample_expression <- sample_expression[sample_expression %in% rownames(genotype)]
        gene.geno <- genotype[sample_expression, , drop = FALSE]
        eQTL <- eQTL[, ..sample_expression]

        gene.geno <- filter_X(gene.geno, missing_rate_thresh = missing_cut, maf_thresh = maf_cut)

        # - covariates
        if (adjust_cov){
            
            if (data.table){
                covariate <- fread(paste0(gtex_eqtl_covariates_dir, "/", tissue.name, ".v10.covariates.txt"))
            } else {
                covariate <- read.delim(paste0(gtex_eqtl_covariates_dir, "/", tissue.name, ".v10.covariates.txt"), stringsAsFactors = FALSE)
                sample_expression <- gsub("-", ".", sample_expression)
            }
            covar <- as.matrix(covariate %>% select(-ID))
            covar <- covar[, sample_expression]
            eQTL.resid <- .lm.fit(x = as.matrix(cbind(1, t(covar))), y = t(as.matrix(eQTL)))$residuals
            geno.resid <- .lm.fit(x = as.matrix(cbind(1, t(covar))), y = as.matrix(gene.geno))$residuals
            X <- c(X, list(geno.resid))
            Y <- c(Y, list(eQTL.resid))
            
        } else {
            X <- c(X, list(as.matrix(gene.geno)))
            Y <- c(Y, list(t(as.matrix(eQTL))))
        }
        phenotypes <- c(phenotypes, tissue.name)

    }
    if (length(phenotypes) == 0){
        message(paste("No eQTL data in GTEx for", gene_ensg))
        return(NULL)
    }
    range_variants <- range(as.numeric(sapply(X, ncol)))
    message(paste("------ Finished preprocess input data for gene", gene_ensg, ".\n",
                  "----- There are", length(phenotypes), "tissues for gene", gene_ensg, 
                  "with number of variants ranging from", range_variants[1], "to", range_variants[2], "."))

    ############ export finalized data ###############
    data <- list(X = X, Y = Y, phenotypes = phenotypes)
    return(data)
    
}


#### this is the function used in line 142 #########
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh = 0, Y=NULL) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, is_zero_variance))
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  if (var_thresh > 0) {
    rm_col <- which(matrixStats::colVars(X) < var_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
  }
  # If Y is provided, remove variants that has zero variance with under the influence of NAs in Y
  if (!is.null(Y) & is.matrix(Y)) {
    drop_snp_indices <- c()
    for (context_idx in 1:ncol(Y)) {
      for (snp_idx in 1:ncol(X)) {
        unique_values <- unique(X[, snp_idx])
        for (value in unique_values) {
          subjects_with_same_genotype <- which(X[, snp_idx] == value)
          subjects_with_na_Y <- which(is.na(Y[, context_idx]))
          if (all(subjects_with_same_genotype %in% subjects_with_na_Y)) {
            # Temporarily remove the specific value and check variance
            temp_X <- X[, snp_idx]
            temp_X[subjects_with_same_genotype] <- NA
            if (sd(temp_X, na.rm = TRUE) == 0) {
              drop_snp_indices <- c(drop_snp_indices, snp_idx)
            } else if (compute_maf(na.omit(temp_X)) <= maf_thresh) {
               drop_snp_indices <- c(drop_snp_indices, snp_idx)
            }
          }
        }
      }
    }
    drop_snp_indices <- unique(drop_snp_indices)
    if (length(drop_snp_indices)) X <- X[, -drop_snp_indices, drop=FALSE]
    message(paste0("Dropped ", length(drop_snp_indices) , " variants with condition of Y subjects, remaining ", ncol(X), "variants. "))
  }
  return(X)
}
compute_maf <- function(geno) {
  f <- mean(geno, na.rm = TRUE) / 2
  return(min(f, 1 - f))
}
compute_missing <- function(geno) {
  miss <- sum(is.na(geno)) / length(geno)
  return(miss)
}
mean_impute <- function(geno) {
  f <- apply(geno, 2, function(x) mean(x, na.rm = TRUE))
  for (i in 1:length(f)) geno[, i][which(is.na(geno[, i]))] <- f[i]
  return(geno)
}
is_zero_variance <- function(x) {
  if (length(unique(x)) == 1) {
    return(T)
  } else {
    return(F)
  }
}



# gene_ensg, gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir,
# gtex_genotype_dir, gene_info_file are all passed in via commandArgs above

# loading all tissues
# QTL_data <- getQTLByLocus(gene_ensg, gene_info_file,
#                           gtex_genotype_dir, gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir,
#                           maf_cut, missing_cut,
#                           data.table = TRUE)

# loading specific tissues




gwas <- fread(gwas_file)
gwas[, z := beta / se]
gwas[, variant_key := paste0("chr", chr, ":", pos)]  # match genotype colname format chr21:pos

geneinfo <- fread(gene_info_file)
gwas_chrs <- unique(gwas$chr)
gene_list <- geneinfo[chr %in% gwas_chrs, gene_id]
message(paste("Scanning", length(gene_list), "genes on chromosomes:", paste(gwas_chrs, collapse = ",")))

coloc_results <- list()

for (g in gene_list) {
    message(paste("===== Gene:", g, "====="))
    tryCatch({

    QTL_data <- getQTLByLocus(g, gene_info_file,
                              gtex_genotype_dir, gtex_eqtl_expression_dir, gtex_eqtl_covariates_dir,
                              maf_cut, missing_cut,
                              tissues_considered = c("Pancreas"),
                              data.table = TRUE)

    QTL_X <- QTL_data$X
    QTL_Y <- QTL_data$Y

    ############ find shared variants between GWAS and eQTL ############

    ld_variants <- colnames(QTL_X[[1]])
    gwas_ld <- gwas[match(ld_variants, variant_key)]   # align to LD variant order
    keep     <- !is.na(gwas_ld$z) & is.finite(gwas_ld$z)

    if (sum(keep) == 0) {
        stop(paste("No overlapping variants between GWAS and eQTL region for gene", g,
                   "\n  eQTL region:", paste(range(ld_variants), collapse = " - "),
                   "\n  GWAS range:  chr21:", min(gwas$pos), "-", max(gwas$pos)))
    }
    message(paste("Overlapping GWAS/eQTL variants:", sum(keep), "/", length(ld_variants)))

    # subset all QTL matrices to the shared variant set
    QTL_X_sub <- lapply(QTL_X, function(x) x[, keep, drop = FALSE])
    LD_sub     <- cor(as.matrix(QTL_X_sub[[1]]))

    ############ run susie on shared variants ############

    # run susie per tissue on subset variants
    res_gtexqtl <- susie(QTL_X_sub[[1]], as.vector(QTL_Y[[1]]))
    names(res_gtexqtl) <- QTL_data$phenotypes

    z_vec       <- gwas_ld$z[keep]
    res_gwassum <- susie_rss(z_vec, LD_sub, n = gwas_n)

    g_region <- paste0("chr", gwas_ld$chr[keep][1], ":", min(gwas_ld$pos[keep]), "-", max(gwas_ld$pos[keep]))
    ld_data <- load_LD_matrix(ld_meta_data, region = g_region, return_genotype = TRUE)

    g_region <- "chr21:10000000-20000000"
    ld_data <- load_LD_matrix(ld_meta_data, region = g_region, return_genotype = TRUE)
    res_coloc <- coloc.susie(res_gwassum, res_gtexqtl)
    print(res_coloc)

    if (!is.null(res_coloc) && !all(is.na(res_coloc$summary$nsnps))) {
        coloc_results[[g]] <- res_coloc
        message(paste("  --> saved result for", g))
    }

    }, error = function(e) message(paste("  ERROR for gene", g, ":", conditionMessage(e))))
}