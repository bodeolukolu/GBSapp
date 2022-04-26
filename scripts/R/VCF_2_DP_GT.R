#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]
ploidy <- args[2]

if (ploidy == "1x"){
  vcf_to_DP_GT_1x <- function() {
    #Let's load data (information retrieved from vcf files)
    path <- getwd()
    file.names <- dir(path, pattern ="1x_rawSPLIT.*\\.vcf")
    subgenome_1 <- NULL
    for(i in 1:length(file.names)){
      all_content <- readLines(file.names[i])
      vcffile <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
      vcffile <- subset(vcffile, select=-c(ID,QUAL,FILTER,INFO,FORMAT))
      cnames <- colnames(vcffile)
      vcffile_GT <- vcffile; vcffile_DP <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="."] <- NA
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.x", replacement = "_DP")
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.y", replacement = "_GT")
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(vcffile[,5:(((ncol(vcffile)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(5:(((ncol(vcffile)-4)/2)+4), function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] == 0, na.rm = TRUE) <= ((ncol(vcffile)-5)/2)*0.8, ]
      vcffile$maf0 <- rowSums(vcffile == "0", na.rm = TRUE) / (rowSums(!is.na(vcffile))-4)
      vcffile$maf1 <- rowSums(vcffile == "1", na.rm = TRUE) / (rowSums(!is.na(vcffile))-5)
      vcffile <- subset(vcffile, maf0 < 0.99)
      vcffile <- subset(vcffile, maf1 < 0.99)
      vcffile <- subset(vcffile, select=-c(maf0,maf1))
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    write.table (subgenome_1, file=paste(pop,"_1x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*1x_rawSPLIT*",sep=""))
  }
  vcf_to_DP_GT_1x()
}
if (ploidy == "2x"){
  vcf_to_DP_GT_2x <- function() {
    #Let's load data (information retrieved from vcf files)
    path <- getwd()
    file.names <- dir(path, pattern ="2x_rawSPLIT.*\\.vcf")
    subgenome_1 <- NULL
    for(i in 1:length(file.names)){
      all_content <- readLines(file.names[i])
      vcffile <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
      vcffile <- subset(vcffile, select=-c(ID,QUAL,FILTER,INFO,FORMAT))
      cnames <- colnames(vcffile)
      vcffile_GT <- vcffile; vcffile_DP <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./."] <- NA
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.x", replacement = "_DP")
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.y", replacement = "_GT")
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(vcffile[,5:(((ncol(vcffile)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(5:(((ncol(vcffile)-4)/2)+4), function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] == 0, na.rm = TRUE) <= ((ncol(vcffile)-5)/2)*0.8, ]
      vcffile$maf0 <- rowSums(vcffile == "0/0", na.rm = TRUE) / (rowSums(!is.na(vcffile))-4)
      vcffile$maf1 <- rowSums(vcffile == "1/1", na.rm = TRUE) / (rowSums(!is.na(vcffile))-5)
      vcffile <- subset(vcffile, maf0 < 0.99)
      vcffile <- subset(vcffile, maf1 < 0.99)
      vcffile <- subset(vcffile, select=-c(maf0,maf1))
      vcffile[vcffile=="1/0"] <- "0/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    write.table (subgenome_1, file=paste(pop,"_2x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*2x_rawSPLIT*",sep=""))
  }
  vcf_to_DP_GT_2x()
}
if (ploidy == "4x"){
  vcf_to_DP_GT_4x <- function() {
    #Let's load data (information retrieved from vcf files)
    path <- getwd()
    file.names <- dir(path, pattern ="4x_rawSPLIT.*\\.vcf")
    subgenome_1 <- NULL
    for(i in 1:length(file.names)){
      all_content <- readLines(file.names[i])
      vcffile <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
      vcffile <- subset(vcffile, select=-c(ID,QUAL,FILTER,INFO,FORMAT))
      cnames <- colnames(vcffile)
      vcffile_GT <- vcffile; vcffile_DP <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././."] <- NA
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.x", replacement = "_DP")
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.y", replacement = "_GT")
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(vcffile[,5:(((ncol(vcffile)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(5:(((ncol(vcffile)-4)/2)+4), function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] == 0, na.rm = TRUE) <= ((ncol(vcffile)-5)/2)*0.8, ]
      vcffile$maf0 <- rowSums(vcffile == "0/0/0/0", na.rm = TRUE) / (rowSums(!is.na(vcffile))-4)
      vcffile$maf1 <- rowSums(vcffile == "1/1/1/1", na.rm = TRUE) / (rowSums(!is.na(vcffile))-5)
      vcffile <- subset(vcffile, maf0 < 0.99)
      vcffile <- subset(vcffile, maf1 < 0.99)
      vcffile <- subset(vcffile, select=-c(maf0,maf1))
      vcffile[vcffile=="1/0/0/0" | vcffile=="0/1/0/0" | vcffile=="0/0/1/0" ] <- "0/0/0/1"
      vcffile[vcffile=="1/1/0/0" | vcffile=="0/1/1/0"] <- "0/0/1/1"
      vcffile[vcffile=="1/1/1/0"] <- "0/1/1/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    write.table (subgenome_1, file=paste(pop,"_4x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*4x_rawSPLIT*",sep=""))
  }
  vcf_to_DP_GT_4x()
}
if (ploidy == "6x"){
  vcf_to_DP_GT_6x <- function() {
    #Let's load data (information retrieved from vcf files)
    path <- getwd()
    file.names <- dir(path, pattern ="6x_rawSPLIT.*\\.vcf")
    subgenome_1 <- NULL
    for(i in 1:length(file.names)){
      all_content <- readLines(file.names[i])
      vcffile <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
      vcffile <- subset(vcffile, select=-c(ID,QUAL,FILTER,INFO,FORMAT))
      cnames <- colnames(vcffile)
      vcffile_GT <- vcffile; vcffile_DP <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././././."] <- NA
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.x", replacement = "_DP")
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.y", replacement = "_GT")
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(vcffile[,5:(((ncol(vcffile)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(5:(((ncol(vcffile)-4)/2)+4), function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] == 0, na.rm = TRUE) <= ((ncol(vcffile)-5)/2)*0.8, ]
      vcffile$maf0 <- rowSums(vcffile == "0/0/0/0/0/0", na.rm = TRUE) / (rowSums(!is.na(vcffile))-4)
      vcffile$maf1 <- rowSums(vcffile == "1/1/1/1/1/1", na.rm = TRUE) / (rowSums(!is.na(vcffile))-5)
      vcffile <- subset(vcffile, maf0 < 0.99)
      vcffile <- subset(vcffile, maf1 < 0.99)
      vcffile <- subset(vcffile, select=-c(maf0,maf1))
      vcffile[vcffile=="1/0/0/0/0/0" | vcffile=="0/1/0/0/0/0" | 
                    vcffile=="0/0/1/0/0/0" | vcffile=="0/0/0/1/0/0" | vcffile=="0/0/0/0/1/0" ] <- "0/0/0/0/0/1"
      vcffile[vcffile=="1/1/0/0/0/0" | vcffile=="0/1/1/0/0/0" | 
                    vcffile=="0/0/1/1/0/0" | vcffile=="0/0/0/1/1/0"] <- "0/0/0/0/1/1"
      vcffile[vcffile=="1/1/1/0/0/0" | vcffile=="0/1/1/1/0/0" | 
                    vcffile=="0/0/1/1/1/0"] <- "0/0/0/1/1/1"
      vcffile[vcffile=="1/1/1/1/0/0" | vcffile=="0/1/1/1/1/0"] <- "0/0/1/1/1/1"
      vcffile[vcffile=="1/1/1/1/1/0"] <- "0/1/1/1/1/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    write.table (subgenome_1, file=paste(pop,"_6x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*6x_rawSPLIT*",sep=""))
  }
  vcf_to_DP_GT_6x()
}
if (ploidy == "8x"){  
  vcf_to_DP_GT_8x <- function() {
    #Let's load data (information retrieved from vcf files)
    path <- getwd()
    file.names <- dir(path, pattern ="8x_rawSPLIT.*\\.vcf")
    subgenome_1 <- NULL
    for(i in 1:length(file.names)){
      all_content <- readLines(file.names[i])
      vcffile <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
      vcffile <- subset(vcffile, select=-c(ID,QUAL,FILTER,INFO,FORMAT))
      cnames <- colnames(vcffile)
      vcffile_GT <- vcffile; vcffile_DP <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././././././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././././././."] <- NA
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.x", replacement = "_DP")
      names(vcffile) <- gsub(names(vcffile), pattern = "\\.y", replacement = "_GT")
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(vcffile[,5:(((ncol(vcffile)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] <- lapply(5:(((ncol(vcffile)-4)/2)+4), function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:(((ncol(vcffile)-4)/2)+4)] == 0, na.rm = TRUE) <= ((ncol(vcffile)-5)/2)*0.8, ]
      vcffile$maf0 <- rowSums(vcffile == "0/0/0/0/0/0/0/0", na.rm = TRUE) / (rowSums(!is.na(vcffile))-4)
      vcffile$maf1 <- rowSums(vcffile == "1/1/1/1/1/1/1/1", na.rm = TRUE) / (rowSums(!is.na(vcffile))-5)
      vcffile <- subset(vcffile, maf0 < 0.99)
      vcffile <- subset(vcffile, maf1 < 0.99)
      vcffile <- subset(vcffile, select=-c(maf0,maf1))
      vcffile[vcffile=="1/0/0/0/0/0/0/0" | vcffile=="0/1/0/0/0/0/0/0" | 
                    vcffile=="0/0/1/0/0/0/0/0" | vcffile=="0/0/0/1/0/0/0/0" | vcffile=="0/0/0/0/1/0/0/0" |
                    vcffile=="0/0/0/0/0/1/0/0" | vcffile=="0/0/0/0/0/0/1/0"] <- "0/0/0/0/0/0/0/1"
      vcffile[vcffile=="1/1/0/0/0/0/0/0" | vcffile=="0/1/1/0/0/0/0/0" | 
                    vcffile=="0/0/1/1/0/0/0/0" | vcffile=="0/0/0/1/1/0/0/0" | vcffile=="0/0/0/0/1/1/0/0" |
                    vcffile=="0/0/0/0/0/1/1/0"] <- "0/0/0/0/0/0/1/1"
      vcffile[vcffile=="1/1/1/0/0/0/0/0" | vcffile=="0/1/1/1/0/0/0/0" | 
                    vcffile=="0/0/1/1/1/0/0/0" | vcffile=="0/0/0/1/1/1/0/0" | vcffile=="0/0/0/0/1/1/1/0"] <- "0/0/0/0/0/1/1/1"
      vcffile[vcffile=="1/1/1/1/0/0/0/0" | vcffile=="0/1/1/1/1/0/0/0" | 
                    vcffile=="0/0/1/1/1/1/0/0" | vcffile=="0/0/0/1/1/1/1/0"] <- "0/0/0/0/1/1/1/1"
      vcffile[vcffile=="1/1/1/1/1/0/0/0" | vcffile=="0/1/1/1/1/1/0/0" | 
                    vcffile=="0/0/1/1/1/1/1/0"] <- "0/0/0/1/1/1/1/1"
      vcffile[vcffile=="1/1/1/1/1/1/0/0" | vcffile=="0/1/1/1/1/1/1/0"] <- "0/0/1/1/1/1/1/1"
      vcffile[vcffile=="1/1/1/1/1/1/1/0"] <- "0/1/1/1/1/1/1/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    write.table (subgenome_1, file=paste(pop,"_8x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*8x_rawSPLIT*",sep=""))
  }
  vcf_to_DP_GT_8x()
}