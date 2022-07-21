#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
ploidy <- args[4]
libdir <- args[5]
.libPaths( c( .libPaths(), libdir) )



vcffile <- read.table(args[1], header = TRUE, check.names = FALSE)
GT <- read.table(args[2], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
AR <- read.table(args[3], header = TRUE, check.names = FALSE)

if (colnames(GT[ncol(GT)]) == "pvalue"){GT <- subset(GT, select=-c(pvalue))}
names(AR) <- gsub(x=names(AR), pattern = "_AR", replacement = "")  
samples <- colnames(GT[,1:ncol(GT)])
AR <- AR[,samples]
write.table (AR, file="AR_temp.txt", row.names=F, quote = FALSE, sep = "\t")


if (ploidy == "2x"){
  for (i in 1:nrow(vcffile)){
    for (j in 6:(ncol(GT))) {
      k=j+4
      if (is.na(GT[i,j])) { vcffile[i,k] <- "./.:0,0:.:.:.:.:."}
      if (AR[i,j] > 0 && AR[i,j] < 0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1", replacement = "1/1")}
      if (AR[i,j] > 0 && AR[i,j] < 0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0|1", replacement = "1/1")}
      if (AR[i,j] < 0 && AR[i,j] > -0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1", replacement = "0/0")}
      if (AR[i,j] < 0 && AR[i,j] > -0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0|1", replacement = "0/0")}
    }
    gc()
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "4x"){
  for (i in 1:nrow(vcffile)){
    for (j in 6:(ncol(GT))) {
      k=j+4
      if (is.na(GT[i,j])) { vcffile[i,k] <- "./././.:0,0:.:.:.:.:."}
      if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1", replacement = "1/1/1/1")}
      if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/1", replacement = "0/0/0/0")}
    }
    gc()
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "6x"){
  for (i in 1:nrow(vcffile)){
    for (j in 6:(ncol(GT))) {
      k=j+4
      if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././.:0,0:.:.:.:.:."}
      if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1/1/1", replacement = "1/1/1/1/1/1")}
      if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/0/0/1", replacement = "0/0/0/0/0/0")}
    }
    gc()
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "8x"){
  for (i in 1:nrow(vcffile)){
    for (j in 6:(ncol(GT))) {
      k=j+4
      if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././././.:0,0:.:.:.:.:."}
      if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1/1/1/1/1", replacement = "1/1/1/1/1/1/1/1")}
      if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/0/0/0/0/1", replacement = "0/0/0/0/0/0/0/0")}
    }
    gc()
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}

