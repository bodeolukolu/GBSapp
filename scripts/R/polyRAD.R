#!/usr/bin/env Rscript

# ==========================================================
# polyRAD imputation script
# Stable version for pipeline execution
# ==========================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop("Usage:
  polyRAD.R <input_vcf> <ploidy> <posterior_txt> <hardcall_vcf> <threads> <R_lib_path> <keep_existing>")
}

input_vcf     <- args[1]
ploidy        <- as.numeric(args[2])
posterior_out <- args[3]
hardcall_vcf  <- args[4]
threads       <- as.numeric(args[5])
rlib          <- args[6]
keep_existing <- as.logical(args[7])

cat("polyRAD imputation started\n")
cat(" Input VCF:", input_vcf, "\n")
cat(" Ploidy:", ploidy, "\n")
cat(" Threads:", threads, "\n")
cat(" Posterior output:", posterior_out, "\n")
cat(" Hardcall VCF:", hardcall_vcf, "\n")
cat(" Keep existing:", keep_existing, "\n")

.libPaths(c(rlib, .libPaths()))

# -------------------------
# Load packages
# -------------------------
suppressPackageStartupMessages({
  library(polyRAD)
  library(future)
})

plan(multisession, workers = threads)

# -------------------------
# Read VCF
# -------------------------
cat("Reading VCF into RADdata object...\n")

rad <- polyRAD::VCF2RADdata(
  file = input_vcf,
  min.ind.with.reads = 1
)

cat("Samples:", polyRAD::nSamples(rad), "\n")
cat("Loci:", polyRAD::nLoci(rad), "\n")

# -------------------------
# Imputation
# -------------------------
cat("Running genotype estimation (HWE model)...\n")

rad <- polyRAD::IterateHWE(rad)

# -------------------------
# Posterior mean genotypes
# -------------------------
cat("Extracting posterior mean genotype matrix...\n")

geno <- polyRAD::GetWeightedMeanGenotypes(rad)

cat("Matrix dimensions:",
    dim(geno)[1], "samples ×",
    dim(geno)[2], "loci\n")

write.table(
  geno,
  posterior_out,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# -------------------------
# Export imputed VCF
# -------------------------
cat("Writing imputed hardcall VCF...\n")

polyRAD::ExportVCF(rad, file = hardcall_vcf)

cat("polyRAD imputation completed successfully\n")
