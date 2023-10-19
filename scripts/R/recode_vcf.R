#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
ploidy <- args[4]
remove_id_list <- NULL
remove_id_list <- read.table(args[5], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
remove_id_list <- as.vector(remove_id_list[,1])
ExcHet <- args[6]
libdir <- args[7]
.libPaths( c( .libPaths(), libdir) )



vcffile <- read.table(args[1], header = TRUE, check.names = FALSE)
if (length(remove_id_list) > 0) {
  id <- names(vcffile)
  keep_id <- setdiff(id,remove_id_list)
  vcffile <- vcffile[,c(keep_id)]
}
GT <- read.table(args[2], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
AR <- read.table(args[3], header = TRUE, check.names = FALSE)

if (colnames(GT[ncol(GT)]) == "pvalue"){GT <- subset(GT, select=-c(pvalue))}
names(AR) <- gsub(x=names(AR), pattern = "_AR", replacement = "")  
samples <- colnames(GT[,1:ncol(GT)])
cols <- samples[samples %in% names(AR)]
ARnew <- AR[cols]
ARnew[setdiff(samples, cols)] <- NA
AR <- ARnew
write.table (AR, file="AR_temp.txt", row.names=F, quote = FALSE, sep = "\t")

keep_id <- c(colnames(vcffile)[1:9], colnames(AR)[6:ncol(AR)])
vcffile <- vcffile[,c(keep_id)]
AR[6:ncol(AR)][is.na(AR[6:ncol(AR)])] <- 0


if (ploidy == "2x"){
  if (ExcHet == "true") {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./.:0,0:.:."}
        if (AR[i,j] > 0 && AR[i,j] < 0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1", replacement = "1/1")}
        if (AR[i,j] > 0 && AR[i,j] < 0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0|1", replacement = "1/1")}
        if (AR[i,j] < 0 && AR[i,j] > -0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1", replacement = "0/0")}
        if (AR[i,j] < 0 && AR[i,j] > -0.2){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0|1", replacement = "0/0")}
      }
      gc()
    }
  } else {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./.:0,0:.:."}
      }
      gc()
    }
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "3x"){
  if (ExcHet == "true") {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "././.:0,0:.:."}
        if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1", replacement = "1/1/1")}
        if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/1", replacement = "0/0/0")}
      }
      gc()
    }
  } else {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "././.:0,0:.:."}
      }
      gc()
    }
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "4x"){
  if (ExcHet == "true") {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././.:0,0:.:."}
        if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1", replacement = "1/1/1/1")}
        if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/1", replacement = "0/0/0/0")}
      }
      gc()
    }
  } else {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././.:0,0:.:."}
      }
      gc()
    }
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "6x"){
  if (ExcHet == "true") {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././.:0,0:.:."}
        if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1/1/1", replacement = "1/1/1/1/1/1")}
        if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/0/0/1", replacement = "0/0/0/0/0/0")}
      }
      gc()
    }
  } else {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././.:0,0:.:."}
      }
      gc()
    }
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}


if (ploidy == "8x"){
  if (ExcHet == "true") {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././././.:0,0:.:."}
        if (AR[i,j] > 0 && AR[i,j] < 0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/1/1/1/1/1/1/1", replacement = "1/1/1/1/1/1/1/1")}
        if (AR[i,j] < 0 && AR[i,j] > -0.17){vcffile[i,k] <- gsub(x=vcffile[i,k], pattern="0/0/0/0/0/0/0/1", replacement = "0/0/0/0/0/0/0/0")}
      }
      gc()
    }
  } else {
    for (i in 1:nrow(vcffile)){
      for (j in 6:(ncol(GT))) {
        k=j+4
        if (is.na(GT[i,j])) { vcffile[i,k] <- "./././././././.:0,0:.:."}
      }
      gc()
    }
  }
  write.table (vcffile, file="dose_temp.vcf", row.names=F, quote = FALSE, sep = "\t")
}

