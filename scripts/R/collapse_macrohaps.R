#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage:",
    "collapse_macrohaps.R <hap_dosage.tsv> <locus_out.tsv>",
    "<ploidy> <prefix> <threads> <libdir>"
  ))
}

hap_file  <- args[1]
locus_out <- args[2]
ploidy    <- as.numeric(args[3])
prefix    <- args[4]
threads   <- as.numeric(args[5])
libdir    <- args[6]

.libPaths(c(libdir, .libPaths()))

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

cat("Threads:", threads, "\n")
plan(multisession, workers = threads)

cat("Reading haplotype dosage file...\n")

hap <- read.table(
  hap_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta_cols <- c("Locus", "Haplotype")

if (!all(meta_cols %in% colnames(hap))) {
  stop("Input must contain columns: Locus and Haplotype")
}

sample_cols <- setdiff(colnames(hap), meta_cols)

cat("Splitting by locus...\n")

hap_split <- split(hap, hap$Locus)

collapse_fun <- function(df) {
  loc <- df$Locus[1]
  vals <- colSums(df[, sample_cols, drop = FALSE])

  vals[vals > ploidy] <- ploidy
  vals[vals < 0] <- 0

  c(Locus = loc, vals)
}

cat("Collapsing loci in parallel...\n")

res_list <- future_lapply(
  hap_split,
  collapse_fun,
  future.seed = TRUE
)

cat("Combining results...\n")

res_mat <- do.call(rbind, res_list)
res_df <- as.data.frame(res_mat, stringsAsFactors = FALSE, check.names = FALSE)

for (j in 2:ncol(res_df)) {
  res_df[[j]] <- as.numeric(res_df[[j]])
}

hap_out_file <- paste0(prefix, ".haplotype_dosage.tsv")
locus_out_file <- paste0(prefix, ".locus_dosage.tsv")

cat("Writing outputs...\n")

write.table(
  hap,
  hap_out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_df,
  locus_out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Macrohap collapse complete.\n")
