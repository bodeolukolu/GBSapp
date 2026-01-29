#!/usr/bin/env Rscript

# update_vcf_from_dose.R
# Purpose: Update VCF genotypes using dose file and optionally AR metrics, ploidy-aware

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(vcfR)
})

option_list <- list(
  make_option(c("--vcf"), type="character", help="Input VCF (.gz supported)"),
  make_option(c("--dose"), type="character", help="Dose file (tab-delimited, first 3 cols: SNP, CHROM, POS, rest samples)"),
  make_option(c("--ploidy"), type="integer", help="Ploidy level (2,4,6,8)"),
  make_option(c("--arfile"), type="character", default=NULL, help="Optional allele ratio metric file"),
  make_option(c("--outfile"), type="character", help="Output VCF path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ----------------------------
# Helper functions
# ----------------------------

dose_to_gt <- function(dose, ploidy) {
  # Convert numeric dose to VCF GT string
  # 0 -> 0/0, 1 -> 0/1 (ploidy=2)
  # 0..ploidy
  if (is.na(dose)) {
    return(paste(rep(".", ploidy), collapse="/"))
  }
  # Generate genotype: reference allele is 0, alternate allele is 1
  alt_count <- as.integer(dose)
  ref_count <- ploidy - alt_count
  gt <- c(rep(0, ref_count), rep(1, alt_count))
  return(paste(gt, collapse="/"))
}

# ----------------------------
# Load dose file
# ----------------------------
dose_dt <- fread(opt$dose, header=TRUE, sep="\t", na.strings=c("NA","NaN"))

# First three columns: SNP, CHROM, POS
snp_col <- names(dose_dt)[1]
chr_col <- names(dose_dt)[2]
pos_col <- names(dose_dt)[3]

samples <- names(dose_dt)[-(1:3)]
dose_dt <- dose_dt[, c(snp_col, chr_col, pos_col, samples), with=FALSE]

# ----------------------------
# Load VCF
# ----------------------------
vcf <- read.vcfR(opt$vcf, verbose=FALSE)

vcf_samples <- colnames(vcf@gt)[-1] # skip GT
vcf_pos <- data.table(CHROM=vcf@fix[,"CHROM"], POS=as.integer(vcf@fix[,"POS"]), row_idx=1:nrow(vcf@fix))

# ----------------------------
# Match SNPs
# ----------------------------
dose_dt[, POS := as.integer(get(pos_col))]
setkeyv(dose_dt, c(chr_col, "POS"))
setkey(vcf_pos, CHROM, POS)
matched <- dose_dt[vcf_pos, nomatch=0] # dose SNPs that exist in VCF

if (nrow(matched)==0) {
  stop("No SNPs from dose file match VCF")
}

# Map row indices of matched SNPs in VCF
row_indices <- matched$row_idx

# ----------------------------
# Update genotypes
# ----------------------------
for (sample in intersect(samples, vcf_samples)) {
  for (i in seq_along(row_indices)) {
    r <- row_indices[i]
    dose_val <- matched[[sample]][i]
    vcf@gt[r, sample] <- dose_to_gt(dose_val, opt$ploidy)
  }
}

# ----------------------------
# Optional: incorporate AR metrics
# ----------------------------
if (!is.null(opt$arfile) && file.exists(opt$arfile)) {
  ar_dt <- fread(opt$arfile, header=TRUE)
  # Optional: could filter or annotate VCF using AR values
  # e.g., force genotype to 0/0 or 1/1 if AR ~0 or 1
  # This can be implemented based on your AR metric format
}

# ----------------------------
# Write output
# ----------------------------
write.vcf(vcf, file=opt$outfile, overwrite=TRUE)
cat("Updated VCF written to:", opt$outfile, "\n")
