#!/usr/bin/env Rscript

# ----------------------------
# Macro-haplotype LD block finder (parallel, no SNPRelate)
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 6){
  stop("Usage: Rscript macrohap_ld_blocks.R <phased.vcf.gz> <ld_window_bp> <ld_r2> <macro_blocks.bed> <R_lib_dir> <threads>")
}

phased_vcf        <- args[1]
ld_window_bp      <- as.numeric(args[2])
ld_r2             <- as.numeric(args[3])
macro_blocks_file <- args[4]
R_LIB_DIR         <- args[5]
threads           <- as.numeric(args[6])

.libPaths(c(R_LIB_DIR, .libPaths()))
suppressPackageStartupMessages({
  library(vcfR)
  library(future.apply)
})

cat("Reading VCF...\n")
vcf <- read.vcfR(phased_vcf, verbose=FALSE)

cat("Extracting genotype matrix...\n")
gt_mat <- extract.gt(vcf, element="GT", as.numeric=FALSE)
# Convert genotypes to 0/1/2 for polyploid counting
gt_mat_num <- apply(gt_mat, 2, function(col){
  sapply(col, function(x){
    if(is.na(x) || x %in% c(".", "./.", ".|.")) return(NA)
    alleles <- unlist(strsplit(x, "[/|]"))
    sum(as.numeric(alleles)!=0)
  })
})

cat("Collecting SNP info...\n")
snp_pos <- as.numeric(vcf@fix[, "POS"])
snp_chr <- vcf@fix[, "CHROM"]
nsnps   <- length(snp_pos)

# Parallel LD block detection per chromosome
cat("Defining LD blocks using future.apply...\n")
plan(multisession, workers = threads)

blocks_list <- future_lapply(unique(snp_chr), function(chr){
  idx <- which(snp_chr == chr)
  chr_pos <- snp_pos[idx]
  chr_gt  <- gt_mat_num[idx, , drop=FALSE]  # rows=SNPs, cols=samples

  blocks <- list()
  i <- 1
  while(i <= length(idx)){
    start_pos <- chr_pos[i]
    j <- i + 1
    while(j <= length(idx) && (chr_pos[j] - chr_pos[j-1]) <= ld_window_bp){
      snp1 <- chr_gt[i, ]
      snp2 <- chr_gt[j, ]
      if(sd(snp1, na.rm=TRUE)==0 || sd(snp2, na.rm=TRUE)==0){
        r2_val <- 1
      } else {
        r2_val <- cor(snp1, snp2, use="complete.obs")^2
      }
      if(is.na(r2_val) || r2_val < ld_r2) break
      j <- j + 1
    }
    end_pos <- chr_pos[j-1]
    blocks[[length(blocks)+1]] <- c(chr, start_pos, end_pos)
    i <- j
  }
  do.call(rbind, blocks)
})

# Combine blocks across chromosomes
blocks_df <- do.call(rbind, blocks_list)
colnames(blocks_df) <- c("CHROM","START","END")
blocks_df <- as.data.frame(blocks_df, stringsAsFactors=FALSE)
blocks_df$START <- as.numeric(blocks_df$START)
blocks_df$END   <- as.numeric(blocks_df$END)

cat("Writing macrohap LD blocks to file...\n")
write.table(blocks_df, macro_blocks_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("Done. Found", nrow(blocks_df), "blocks.\n")
