#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

hap_file   <- args[1]
locus_file <- args[2]
outprefix  <- args[3]
ploidy     <- as.numeric(args[4])
threads    <- as.numeric(args[5])
libdir     <- args[6]

.libPaths(c(libdir, .libPaths()))

suppressPackageStartupMessages({
  library(future.apply)
})

plan(multisession, workers = threads)

cat("collapse_microhaps.R starting\n")
cat("threads:", threads, "\n")
cat("ploidy:", ploidy, "\n")

collapse_file <- function(infile, outfile) {

  cat("Reading:", infile, "\n")

  df <- read.delim(
    infile,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  meta_cols <- df[, 1:3, drop = FALSE]
  dosage_df <- df[, -(1:3), drop = FALSE]

  sample_names <- colnames(dosage_df)

  # Remove ploidy prefix like "2:M19"
  base_names <- sub("^[0-9]+:", "", sample_names)

  unique_samples <- unique(base_names)

  cat("Detected", length(unique_samples), "unique samples\n")

  collapse_sample <- function(s) {
    idx <- which(base_names == s)
    mat <- dosage_df[, idx, drop = FALSE]
    mat[is.na(mat) | mat == "." | mat == ""] <- 0
    rowSums(apply(mat, 2, as.numeric))
  }

  collapsed_list <- future_lapply(unique_samples, collapse_sample)

  collapsed_mat <- do.call(cbind, collapsed_list)
  colnames(collapsed_mat) <- unique_samples

  result <- cbind(meta_cols, collapsed_mat)

  write.table(
    result,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Wrote:", outfile, "\n")
}

hap_out   <- paste0(outprefix, ".haplotype_dosage.collapsed.tsv")
locus_out <- paste0(outprefix, ".locus_dosage.collapsed.tsv")

collapse_file(hap_file, hap_out)
collapse_file(locus_file, locus_out)

cat("collapse_microhaps.R finished\n")
