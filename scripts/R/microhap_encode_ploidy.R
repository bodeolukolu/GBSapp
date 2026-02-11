#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript microhap_encode_ploidy.R <microhap.raw.tsv> <ploidy> <out_prefix>")
}

infile    <- args[1]
ploidy    <- as.numeric(args[2])
outprefix <- args[3]

suppressPackageStartupMessages({
  library(data.table)
})

cat("Reading microhap raw table...\n")
raw <- fread(infile, header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# Expect: HapID | CHROM | START | samples...
if (!all(c("HapID", "CHROM", "START") %in% colnames(raw))) {
  stop("Input must contain columns: HapID, CHROM, START")
}

meta <- raw[, .(HapID, CHROM, START)]
geno <- raw[, -(1:3), with = FALSE]

cat("Encoding haplotype strings to numeric dosage...\n")

# Vectorized allele counting: count '1's in each cell
geno_num <- as.data.table(
  lapply(geno, function(col) {
    v <- as.character(col)
    out <- integer(length(v))
    miss <- is.na(v) | v %in% c(".", "./.", "")
    out[miss] <- NA_integer_
    out[!miss] <- lengths(regmatches(v[!miss], gregexpr("1", v[!miss])))
    out
  })
)

# ---- Validation: fail hard if monomorphic ----
vals <- unique(unlist(geno_num, use.names = FALSE))
vals <- vals[!is.na(vals)]

if (length(vals) <= 1) {
  stop("ERROR: Encoded haplotype dosage is monomorphic — check input or encoding logic.")
}

cat("Observed haplotype dosage range:",
    paste(range(vals), collapse = "–"), "\n")

# ---- Write haplotype-level dosage ----
hap_out <- cbind(meta, geno_num)
hap_file <- paste0(outprefix, ".haplotype_dosage.tsv")

cat("Writing haplotype dosage to:", hap_file, "\n")
fwrite(hap_out, hap_file, sep = "\t", quote = FALSE)

# ---- Collapse to locus-level dosage (cap at ploidy) ----
cat("Collapsing to locus-level dosage...\n")

locus_out <- hap_out[
  ,
  lapply(.SD, function(x) pmin(ploidy, sum(x, na.rm = TRUE))),
  by = .(CHROM, START),
  .SDcols = colnames(geno_num)
]

locus_file <- paste0(outprefix, ".locus_dosage.tsv")
fwrite(locus_out, locus_file, sep = "\t", quote = FALSE)

cat("Done.\n")
