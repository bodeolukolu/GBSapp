#!/usr/bin/env Rscript


####################################################################################################################
args <- commandArgs(trailingOnly = TRUE)
args
geno <- read.table(args[1], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
geno$CHROM <- sub(".*\\_","",geno$CHROM)
geno$CHROM <- sub("Chr","",geno$CHROM)
geno$CHROM <- as.numeric(as.character(geno$CHROM))

# convert genotype file to hapmap file
if (colnames(geno[ncol(geno)]) == "pvalue"){geno <- subset(geno, select=-c(pvalue))}

geno[,6:ncol(geno)][is.na(geno[,6:ncol(geno)])] <- "NN"
geno_indel <- subset(geno, nchar(geno$REF) > 1 | nchar(geno$ALT) > 1 )
geno_overlap <- geno  
geno_overlap$ALT <- grepl(pattern = "\\*", geno_overlap$ALT)
geno_overlap <- geno_overlap[geno_overlap$ALT == TRUE, ]
geno_indel <- rbind(geno_indel,geno_overlap)
geno_snp <- subset(geno, nchar(geno$REF) == 1 & nchar(geno$ALT) == 1 )
geno_snp$overlap <- geno_snp$ALT
geno_snp$overlap <- grepl(pattern = "\\*", geno_snp$overlap)
geno_snp <- geno_snp[geno_snp$overlap == FALSE, ]
geno_snp <- subset(geno_snp, select=-c(overlap))

geno <- geno[-c(1:nrow(geno)),]
geno_indel$refcode <- nchar(geno_indel$REF);  geno_indel$refcode[ geno_indel$refcode == 1] <- "-"
geno_indel$refcode[ geno_indel$refcode > 1] <- "+"
geno_indel$altcode <- nchar(geno_indel$ALT);  geno_indel$altcode[ geno_indel$altcode < 2] <- "-"
geno_indel$altcode[ geno_indel$altcode > 1] <- "+"

for (r in c("-","+")){
  for (c in c("-","+")) {
    geno_indelsub <- subset(geno_indel, geno_indel$refcode == r & geno_indel$altcode == c)
    geno_indelsub[,6:ncol(geno_indelsub)][geno_indelsub[,6:ncol(geno_indelsub)] == 0] <- paste(r,r,sep="")
    geno_indelsub[,6:ncol(geno_indelsub)][geno_indelsub[,6:ncol(geno_indelsub)] == 1] <- paste(r,c,sep="")
    geno_indelsub[,6:ncol(geno_indelsub)][geno_indelsub[,6:ncol(geno_indelsub)] == 2] <- paste(c,c,sep="")
    geno_indelsub <- subset(geno_indelsub, select=-c(refcode,altcode))
    geno <- rbind(geno,geno_indelsub)
  }
}
ref_alleles <- unique(geno_snp[,4])
alt_alleles <- unique(geno_snp[,5])
for (r in c(ref_alleles)){
  for (c in c(alt_alleles)) {
    geno_snpsub <- subset(geno_snp, geno_snp$REF == r & geno_snp$ALT == c)
    geno_snpsub[,6:ncol(geno_snpsub)][geno_snpsub[,6:ncol(geno_snpsub)] == 0] <- paste(r,r,sep="")
    geno_snpsub[,6:ncol(geno_snpsub)][geno_snpsub[,6:ncol(geno_snpsub)] == 1] <- paste(r,c,sep="")
    geno_snpsub[,6:ncol(geno_snpsub)][geno_snpsub[,6:ncol(geno_snpsub)] == 2] <- paste(c,c,sep="")
    geno <- rbind(geno,geno_snpsub)
  }
}

geno$ALT <- sub("TRUE","*",geno$ALT)
geno$rs <- geno$SNP
geno$alleles <- paste(geno$REF,"/",geno$ALT,sep="")
geno$chrom <- geno$CHROM
geno$pos <- geno$POS
geno <- geno[,-c(1:5)]
geno$strand <- "+"
geno$assembly <- "NA"
geno$center <- "NA"
geno$protLSID <- "NA"
geno$assayLSID <- "NA"
geno$panelLSID <- "NA"
geno$QCcode <- "NA"
geno <- subset(geno, select=c((ncol(geno)-10):ncol(geno),1:(ncol(geno)-11)))
write.table(geno,"outfile.hmp.txt", row.names=F, quote = FALSE, sep="\t")
