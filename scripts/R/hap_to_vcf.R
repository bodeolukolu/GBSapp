#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop("Usage: Rscript hap_to_vcf.R <input.tsv> <output.vcf>")
}

infile  <- args[1]
outfile <- args[2]

# Read haplotype dosage table
cat("Reading haplotype file...\n")
hap <- read.table(infile, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Check required columns
if(!all(c("HapID","CHROM","POS") %in% colnames(hap))){
  stop("Input must have columns: HapID, CHROM, POS")
}

sample_cols <- setdiff(colnames(hap), c("HapID","CHROM","POS"))

# VCF header
vcf_header <- c(
  "##fileformat=VCFv4.3",
  paste0("##source=hap_to_vcf.R"),
  paste0("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"),
  paste0("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(sample_cols, collapse="\t"))
)

cat("Writing VCF header...\n")
writeLines(vcf_header, con = outfile)

# Function to convert 0/1/2 to diploid GT
convert_gt <- function(x){
  if(is.na(x) || x == "." || x == "") return("./.")
  x <- as.numeric(x)
  if(x == 0) return("0/0")
  if(x == 1) return("0/1")
  if(x == 2) return("1/1")
  return("./.")  # fallback
}

cat("Processing haplotypes...\n")
for(i in 1:nrow(hap)){
  chrom <- hap$CHROM[i]
  pos   <- hap$POS[i]
  id    <- hap$HapID[i]
  ref   <- "A"   # placeholder REF allele
  alt   <- "T"   # placeholder ALT allele
  qual  <- "."
  filter<- "PASS"
  info  <- "."
  format<- "GT"

  gts <- sapply(hap[i, sample_cols], convert_gt)

  line <- paste(chrom, pos, id, ref, alt, qual, filter, info, format, paste(gts, collapse="\t"), sep="\t")
  write(line, file=outfile, append=TRUE)
}

cat("VCF written to", outfile, "\n")
