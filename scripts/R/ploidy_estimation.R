#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
libdir <- args[8]
.libPaths(c(libdir, .libPaths()))

suppressPackageStartupMessages({
  library(vcfR)
  library(future.apply)
  library(ggplot2)
  library(MASS)  # for robust regression
})

if(length(args) < 8){
  stop("Usage: Rscript ploidy_estimation_vcfR_parallel.R <vcf> <chr_or_all> <min_cov> <max_ploidy> <window> <outdir> <threads> <libdir>")
}

vcf_file   <- args[1]
chrom_sel  <- args[2]
min_cov    <- as.numeric(args[3])
max_ploidy <- as.integer(args[4])
window     <- as.integer(args[5])
outdir     <- args[6]
threads    <- as.integer(args[7])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
plan(multisession, workers = threads)

# ------------------------------
# Load VCF
# ------------------------------
vcf <- read.vcfR(vcf_file, verbose = FALSE)
fix <- getFIX(vcf)
chrom <- fix[, "CHROM"]
pos <- as.integer(fix[, "POS"])
samples <- colnames(vcf@gt)[-1]  # drop FORMAT

# Chromosome filter
if(chrom_sel != "all"){
  keep <- chrom == chrom_sel
  vcf <- vcf[keep,]
  fix <- fix[keep,]
  chrom <- chrom[keep]
  pos <- pos[keep]
}

# Biallelic SNPs only
alt <- fix[, "ALT"]
ref <- fix[, "REF"]
is_biallelic <- sapply(strsplit(alt, ","), length) == 1
is_snp <- nchar(ref) == 1 & nchar(alt) == 1
keep_snps <- is_biallelic & is_snp
vcf <- vcf[keep_snps,]
chrom <- chrom[keep_snps]
pos <- pos[keep_snps]

n_snps <- nrow(vcf@gt)
n_samples <- length(samples)

# ------------------------------
# Extract GT and DP
# ------------------------------
gt <- extract.gt(vcf, "GT")
dp <- extract.gt(vcf, "DP", as.numeric = TRUE)

# ------------------------------
# Helper functions
# ------------------------------
count_alt <- function(gt_str){
  if(is.na(gt_str) || gt_str == "." || gt_str == "") return(NA)
  alleles <- unlist(strsplit(gt_str, "[/|]"))
  sum(alleles != "0")
}

smooth_vector <- function(x, window){
  half <- floor(window/2)
  n <- length(x)
  sapply(1:n, function(i){
    idx <- max(1,i-half):min(n,i+half)
    median(x[idx], na.rm=TRUE)
  })
}

# ------------------------------
# Precompute expected AR classes for all ploidies
# ------------------------------
expected_ARs <- lapply(2:max_ploidy, function(P) (1:(P-1))/P)
names(expected_ARs) <- 2:max_ploidy

# ------------------------------
# Compute best-fitting ploidy using weighted RMSE and robust regression
# ------------------------------
compute_best_ploidy <- function(ar_vec, dp_vec, expected_ARs){
  valid <- !is.na(ar_vec) & !is.na(dp_vec) & dp_vec > 0
  if(sum(valid) == 0) return(NA)

  ar_obs <- ar_vec[valid]
  w <- dp_vec[valid]  # weight by coverage

  rmse_per_ploidy <- sapply(expected_ARs, function(ar_exp){
    nearest <- sapply(ar_obs, function(x) ar_exp[which.min(abs(ar_exp - x))])
    fit <- rlm(ar_obs ~ 0 + nearest, weights=w)
    sqrt(sum(w * fit$resid^2)/sum(w))
  })

  as.integer(names(rmse_per_ploidy)[which.min(rmse_per_ploidy)])
}

# ------------------------------
# Per-sample parallel ploidy estimation
# ------------------------------
results_list <- future_lapply(1:n_samples, function(j){

  dp_vec <- as.integer(dp[, j])

  # All heterozygotes
  het_mask <- sapply(gt[, j], function(g){
    alt_c <- count_alt(g)
    total <- length(unlist(strsplit(g, "[/|]")))
    !is.na(alt_c) && total > 0 && alt_c > 0 && alt_c < total
  })

  # Observed AR
  ar_vec <- sapply(gt[, j], function(g){
    alt_c <- count_alt(g)
    total <- length(unlist(strsplit(g, "[/|]")))
    if(!is.na(alt_c) && total > 0) alt_c / total else NA
  })

  # Coverage filter
  ar_vec[dp_vec < min_cov] <- NA
  # Heterozygote filter
  ar_vec[!het_mask] <- NA

  # Compute best-fitting ploidy per SNP/locus
  ploidy_vec <- sapply(1:length(ar_vec), function(i){
    compute_best_ploidy(ar_vec[i], dp_vec[i], expected_ARs)
  })

  # Smooth ploidy along chromosome
  ploidy_smooth <- smooth_vector(ploidy_vec, window)

  data.frame(
    SNP = paste0(chrom, ":", pos),
    CHR = chrom,
    POS = pos,
    Sample = samples[j],
    Ploidy = ploidy_vec,
    SmoothedPloidy = ploidy_smooth
  )

}, future.seed=TRUE)

ploidy_final <- do.call(rbind, results_list)

# ------------------------------
# Write outputs
# ------------------------------
out_base <- paste0(outdir, "/ploidy_", chrom_sel)
write.table(ploidy_final, file=paste0(out_base, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)

chr_summary <- aggregate(SmoothedPloidy ~ Sample + CHR, data=ploidy_final, FUN=function(x) median(x, na.rm=TRUE))
write.table(chr_summary, file=paste0(out_base, "_chr_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

genome_summary <- aggregate(SmoothedPloidy ~ Sample, data=ploidy_final, FUN=function(x) median(x, na.rm=TRUE))
write.table(genome_summary, file=paste0(out_base, "_genome_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ------------------------------
# Visualization
# ------------------------------
plot_df <- ploidy_final[!is.na(ploidy_final$Ploidy), ]
if(nrow(plot_df) > 0){

  # Per-SNP plot with loess smoothing
  p1 <- ggplot(plot_df, aes(x=POS, y=SmoothedPloidy, color=Sample)) +
    geom_point(size=0.8, alpha=0.7) +
    geom_line(aes(group=Sample), alpha=0.5) +
    facet_grid(CHR ~ ., scales="free_x") +
    labs(title="Smoothed ploidy per SNP",
         x="Position (bp)", y="Smoothed Ploidy") +
    theme_bw() +
    theme(legend.position="bottom")
  ggsave(filename=paste0(out_base, "_ploidy_loess.tiff"), plot=p1,
         width=12, height=1 * length(unique(plot_df$CHR)), dpi = 600, compression = "lzw")

  # Boxplot of median ploidy per chromosome
  p2 <- ggplot(chr_summary, aes(x=Sample, y=SmoothedPloidy, fill=Sample)) +
    geom_boxplot(alpha=0.7) +
    labs(x="Sample", y="Ploidy") +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none")
  ggsave(filename=paste0(out_base, "_chr_summary_boxplot.tiff"), plot=p2,
         width=4, height=6, dpi = 600, compression = "lzw")
} else {
  cat("No ploidy points available for plotting.\n")
}
