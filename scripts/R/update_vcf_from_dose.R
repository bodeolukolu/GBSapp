#!/usr/bin/env Rscript

# update_vcf_from_dose.R
# Purpose: Update VCF genotypes using dose file and optionally AR metrics, ploidy-aware
# Works with ploidy 2,4,6,8
# Uses base R (no optparse, no data.table)

suppressPackageStartupMessages({
  library(vcfR)
})

# ----------------------------
# Parse command-line arguments (base R)
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript update_vcf_from_dose.R <vcf> <dose> <ploidy> <outfile> [arfile]")
}

vcf_file   <- args[1]
dose_file  <- args[2]
ploidy     <- as.integer(args[3])
outfile    <- args[4]
ar_file    <- ifelse(length(args) >= 5, args[5], NA)

# ----------------------------
# Helper functions
# ----------------------------
dose_to_gt <- function(dose, ploidy) {
  # Convert numeric dose to VCF GT string
  if (is.na(dose)) {
    return(paste(rep(".", ploidy), collapse="/"))
  }
  alt_count <- as.integer(dose)
  ref_count <- ploidy - alt_count
  gt <- c(rep(0, ref_count), rep(1, alt_count))
  return(paste(gt, collapse="/"))
}

# ----------------------------
# Load dose file
# ----------------------------
dose_dt <- read.table(dose_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = c("NA","NaN"))

# First three columns: SNP, CHROM, POS
snp_col <- names(dose_dt)[1]
chr_col <- names(dose_dt)[2]
pos_col <- names(dose_dt)[3]

samples <- names(dose_dt)[-(1:3)]

# ----------------------------
# Load VCF
# ----------------------------
vcf <- read.vcfR(vcf_file, verbose = FALSE)
vcf_samples <- colnames(vcf@gt)[-1] # skip GT column
vcf_pos <- data.frame(CHROM = vcf@fix[,"CHROM"],
                      POS   = as.integer(vcf@fix[,"POS"]),
                      row_idx = 1:nrow(vcf@fix),
                      stringsAsFactors = FALSE)

# ----------------------------
# Match SNPs
# ----------------------------
dose_pos <- dose_dt[, c(chr_col, pos_col)]
names(dose_pos) <- c("CHROM","POS")
dose_pos$idx <- 1:nrow(dose_pos)

# Merge dose positions with VCF positions
matched <- merge(dose_pos, vcf_pos, by = c("CHROM","POS"))
if (nrow(matched) == 0) stop("No SNPs from dose file match VCF")

# ----------------------------
# Update genotypes
# ----------------------------
for (sample in intersect(samples, vcf_samples)) {
  dose_vals <- dose_dt[matched$idx, sample]
  for (i in seq_len(nrow(matched))) {
    r <- matched$row_idx[i]
    vcf@gt[r, sample] <- dose_to_gt(dose_vals[i], ploidy)
  }
}

# ----------------------------
# Optional AR metrics (placeholder)
# ----------------------------
if (!is.na(ar_file) && file.exists(ar_file)) {
  ar_dt <- read.table(ar_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  # Optional: use AR metrics to correct genotypes
}

# ----------------------------
# Write updated VCF
# ----------------------------
write.vcf(vcf, file = outfile, overwrite = TRUE)
cat("Updated VCF written to:", outfile, "\n")
