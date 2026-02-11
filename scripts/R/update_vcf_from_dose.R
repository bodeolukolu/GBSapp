#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript update_vcf_from_dose.R ",
    "<vcf> <dose> <ploidy> <outfile> [libdir] [threads]"
  )
}

vcf_file  <- args[1]
dose_file <- args[2]
ploidy    <- as.integer(args[3])
outfile   <- args[4]
libdir    <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL

# ----------------------------
# Threads
# ----------------------------
available_cores <- parallel::detectCores()
if (length(args) >= 6 && nzchar(args[6])) {
  threads <- suppressWarnings(as.integer(args[6]))
  threads <- if (is.na(threads)) 1L else max(1L, threads - 2L)
} else {
  threads <- max(1L, available_cores - 1L)
}
cat("Using", threads, "threads\n")

if (!is.null(libdir)) .libPaths(c(.libPaths(), libdir))

suppressPackageStartupMessages({
  library(vcfR)
  library(parallel)
})

# ----------------------------
# Helper function
# ----------------------------
dose_to_gt <- function(dose, ploidy) {
  if (is.na(dose)) return(paste(rep(".", ploidy), collapse = "/"))
  alt <- as.integer(dose)
  ref <- ploidy - alt
  paste(c(rep(0, ref), rep(1, alt)), collapse = "/")
}

# ----------------------------
# Load dose file
# ----------------------------
dose_dt <- read.table(
  dose_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NaN")
)

chr_col <- names(dose_dt)[2]
pos_col <- names(dose_dt)[3]
dose_samples <- names(dose_dt)[-(1:3)]

dose_pos <- dose_dt[, c(chr_col, pos_col)]
names(dose_pos) <- c("CHROM", "POS")
dose_pos$idx <- seq_len(nrow(dose_pos))

# ----------------------------
# Load VCF
# ----------------------------
vcf <- read.vcfR(vcf_file, verbose = FALSE)
vcf_samples <- colnames(vcf@gt)[-1]

vcf_pos <- data.frame(
  CHROM   = vcf@fix[, "CHROM"],
  POS     = as.integer(vcf@fix[, "POS"]),
  row_idx = seq_len(nrow(vcf@fix)),
  stringsAsFactors = FALSE
)

# ----------------------------
# Match SNPs and subset VCF
# ----------------------------
matched <- merge(dose_pos, vcf_pos, by = c("CHROM", "POS"))
if (nrow(matched) == 0)
  stop("No SNPs from dose file match VCF")

vcf <- vcf[matched$row_idx, , drop = FALSE]

# ----------------------------
# Subset samples
# ----------------------------
common_samples <- intersect(dose_samples, vcf_samples)
if (length(common_samples) == 0)
  stop("No overlapping samples between dose file and VCF")

vcf@gt <- vcf@gt[, c("FORMAT", common_samples), drop = FALSE]

# ----------------------------
# Parallel genotype update
# ----------------------------
cl <- makeCluster(threads)
clusterExport(
  cl,
  varlist = c("dose_dt", "matched", "ploidy", "dose_to_gt"),
  envir = environment()
)

gt_updates <- parLapply(cl, common_samples, function(sample) {
  dose_vals <- dose_dt[matched$idx, sample]
  gts <- vapply(dose_vals, dose_to_gt, character(1), ploidy = ploidy)
  list(sample = sample, gt = gts)
})
stopCluster(cl)

# Apply GT updates
format_fields <- strsplit(vcf@gt[, "FORMAT"], ":")[[1]]
gt_index <- which(format_fields == "GT")

if (length(gt_index) == 0) {
  stop("GT field not found in FORMAT column")
}

for (res in gt_updates) {

  sample <- res$sample
  new_gt <- res$gt

  old_vals <- vcf@gt[, sample]

  updated_vals <- mapply(function(old_entry, new_gt_val) {

    if (is.na(old_entry) || old_entry == ".") return(old_entry)

    parts <- strsplit(old_entry, ":")[[1]]

    if (length(parts) >= gt_index) {
      parts[gt_index] <- new_gt_val
    }

    paste(parts, collapse=":")

  }, old_vals, new_gt, USE.NAMES = FALSE)

  vcf@gt[, sample] <- updated_vals
}

# ----------------------------
# Write updated VCF
# ----------------------------
write.vcf(vcf, file = outfile)
cat("Updated VCF written to:", outfile, "\n")
