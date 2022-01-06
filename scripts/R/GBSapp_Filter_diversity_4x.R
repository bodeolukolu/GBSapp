#!/usr/bin/env Rscript

# pop <-
# gmissingness <-
# smissingness <-
# minRD <-
# hap <-
# MinorAlleleFreq <-
# snpformats <-
# remove_id_list <- NULL
# remove_id_list <- c("NA_1","NA_2")
# remove_id_list <- paste(remove_id_list, "_GT", sep="")
# setwd ("dir")

####################################################################################################################
# set_variables
args <- commandArgs(trailingOnly = TRUE)
args
pop <- args[1]
gmissingness <- args[2]
smissingness <- args[3]
minRD <- args[4]
remove_id_list <- NULL
remove_id_list <- unlist(strsplit(args[5],","))
remove_id_list <- paste(remove_id_list, "_GT", sep="")
libdir <- args[6]
MinorAlleleFreq <- args[7]
snpformats <- args[8]
hap <- args[9]
gmissingness <- as.numeric(gmissingness)
smissingness <- as.numeric(smissingness)
minRD <- as.numeric(minRD)
hap <- as.numeric(hap)
rd <- minRD-1


.libPaths( c( .libPaths(), libdir) )
library(ggplot2)

####################################################################################################################
####################################################################################################################
subgenome_1 <- read.table (file=paste("../../snpcall/",pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
subgenome_1[, 5:(((ncol(subgenome_1)-4)/2)+4)] <- lapply(subgenome_1[,5:(((ncol(subgenome_1)-4)/2)+4)], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
subgenome_1[, 5:(((ncol(subgenome_1)-4)/2)+4)] <- lapply(5:(((ncol(subgenome_1)-4)/2)+4), function(x) as.numeric(subgenome_1[[x]]))
subgenome_1[,(((ncol(subgenome_1)-4)/2)+5):ncol(subgenome_1)] <- lapply(subgenome_1[,(((ncol(subgenome_1)-4)/2)+5):ncol(subgenome_1)], gsub, pattern = "|", replacement = "/", fixed = TRUE)
RD_snpfiltering <- function() {
  #######################################################################################################################################################################################
  # Let's filter the variants based on the following parameters: (1) read depth, (2) gmissingness, and (3) various thresholds for minor allele frequency (maf). Let's plot distribution of maf
  ############################################
  # Filter 4x subgenome
    #remove samples that you want to exclude from the analysis
  if (length(remove_id_list) > 0) {
    remove_id_GT <- remove_id_list
    remove_id_DP <- gsub("_GT", "_DP", remove_id_list)
    remove_id <- c(remove_id_DP, remove_id_GT)
    id <- names(subgenome_1)
    keep_id <- setdiff(id,remove_id)
    subgenome_1 <- subgenome_1[,c(keep_id)]
  }

  subgenome_filtered <- subgenome_1
  subgenome_filtered$no_missing <- apply(subgenome_filtered, 1, function(x) sum(is.na(x)))
  subgenome_filtered <- subset(subgenome_filtered, no_missing <= ((ncol(subgenome_filtered)-5)/2)*gmissingness)
  subgenome_filtered <- subset(subgenome_filtered, select=-c(no_missing))

  subgenome_filtered_AB <- subset(subgenome_filtered, select=c(1:(((ncol(subgenome_filtered)-4)/2)+4)))
  for (i in 5:(((ncol(subgenome_filtered)-4)/2)+4)) {
    j <- i+((ncol(subgenome_filtered)-4)/2)
    subgenome_filtered[,j][subgenome_filtered[,i] < 12] <- NA
    subgenome_filtered[,i][subgenome_filtered[,j] == "0/0/0/0" ] <- minRD
    subgenome_filtered[,i][subgenome_filtered[,j] == "1/1/1/1" ] <- minRD
    subgenome_filtered[,j][subgenome_filtered[,i] <= rd] <- NA
    gc()
  }
  subgenome_filtered_C <- subset(subgenome_filtered, select=c((((ncol(subgenome_filtered)-4)/2)+5):ncol(subgenome_filtered)))
  subgenome_filtered <- cbind(subgenome_filtered_AB, subgenome_filtered_C)
  write.table (subgenome_filtered, file=paste(pop,"_4x_rawRD",rd+1,"_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  subgenome_SD <- subgenome_filtered

  subgenome_SD <- subset(subgenome_SD, select=c(1:4,(((ncol(subgenome_SD)-4)/2)+5):ncol(subgenome_SD)))
  subgenome_SD$no_missing <- apply(subgenome_SD, 1, function(x) sum(is.na(x)))
  subgenome_SD <- subset(subgenome_SD, no_missing <= (ncol(subgenome_SD)-5)*gmissingness)
  subgenome_SD <- subset(subgenome_SD, select=c(-no_missing))
  nonbiallelic <- as.vector(as.matrix(subgenome_SD[,5:ncol(subgenome_SD)]))
  nonbiallelic <- unique(nonbiallelic)
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/0")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/1/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "1/1/1/1")
  if (length(nonbiallelic) > 0) {
    for (i in 1:length(nonbiallelic)){
      subgenome_SD[][subgenome_SD[]==nonbiallelic[i]] <- NA
    }
  }
  subgenome_SD$no_missing <- apply(subgenome_SD, 1, function(x) sum(is.na(x)))
  subgenome_SD <- subset(subgenome_SD, no_missing <= ((ncol(subgenome_SD)-5))*gmissingness)
  subgenome_SD <- subset(subgenome_SD, select=c(-no_missing))
  names(subgenome_SD) <- gsub("X", "", names(subgenome_SD))
  names(subgenome_SD) <- gsub("_GT", "", names(subgenome_SD))

  # remove individuals above missingness threshold
  remove_id_list <- subset(subgenome_SD, select=c(5:ncol(subgenome_SD)))
  remove_id_list <- as.data.frame(t(remove_id_list))
  remove_id_list$percent_missing <- (apply(remove_id_list, 1, function(x) sum(is.na(x))))/ncol(remove_id_list)*100
  remove_id_list <- subset(remove_id_list, select=c(percent_missing))
  remove_id_list$samples <- rownames(remove_id_list)
  remove_id_list <- subset(remove_id_list, select=c("samples","percent_missing"))
  remove_id_list <- remove_id_list[order(remove_id_list$percent_missing),]
  write.table (remove_id_list, file="sample_missing_rate_4x.txt", row.names=F, quote = FALSE, sep = "\t")
  remove_id_list <- as.list(subset(remove_id_list, remove_id_list$percent_missing > smissingness*100))
  write.table (remove_id_list, file="eliminated_samples_4x.txt", row.names=F, quote = FALSE, sep = "\t")
  if (length(remove_id_list) > 0 ) {
    remove_id_list <- remove_id_list[["samples"]]
    remove_id_GT <- remove_id_list
    remove_id_GT <- sub("$", "_GT\\1", remove_id_GT)
    remove_id_DP <- gsub("_GT", "_DP", remove_id_GT)
    remove_id <- c(remove_id_DP, remove_id_GT)
    id <- colnames(subgenome_1[,1:ncol(subgenome_1)])
    keep_id <- setdiff(id,remove_id)
    subgenome_1 <-subset(subgenome_1, select=c(keep_id))
    write.table (subgenome_1, file=paste(pop,"_4x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }

  subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
  subgenome_1[, 5:(((ncol(subgenome_1)-4)/2)+4)] <- lapply(5:(((ncol(subgenome_1)-4)/2)+4), function(x) as.numeric(subgenome_1[[x]]))
  subgenome_filtered <- subgenome_1
  subgenome_filtered$no_missing <- apply(subgenome_filtered, 1, function(x) sum(is.na(x)))
  subgenome_filtered <- subset(subgenome_filtered, no_missing <= ((ncol(subgenome_filtered)-5)/2)*gmissingness)
  subgenome_filtered <- subset(subgenome_filtered, select=-c(no_missing))

  subgenome_filtered_AB <- subset(subgenome_filtered, select=c(1:(((ncol(subgenome_filtered)-4)/2)+4)))
  for (i in 5:(((ncol(subgenome_filtered)-4)/2)+4)) {
    j <- i+((ncol(subgenome_filtered)-4)/2)
    subgenome_filtered[,j][subgenome_filtered[,i] < 12] <- NA
    subgenome_filtered[,i][subgenome_filtered[,j] == "0/0/0/0" ] <- minRD
    subgenome_filtered[,i][subgenome_filtered[,j] == "1/1/1/1" ] <- minRD
    subgenome_filtered[,j][subgenome_filtered[,i] <= rd] <- NA
    gc()
  }
  subgenome_filtered_C <- subset(subgenome_filtered, select=c((((ncol(subgenome_filtered)-4)/2)+5):ncol(subgenome_filtered)))
  subgenome_filtered <- cbind(subgenome_filtered_AB, subgenome_filtered_C)
  subgenome_SD <- subgenome_filtered

  subgenome_SD <- subset(subgenome_SD, select=c(1:4,(((ncol(subgenome_SD)-4)/2)+5):ncol(subgenome_SD)))
  subgenome_SD$no_missing <- apply(subgenome_SD, 1, function(x) sum(is.na(x)))
  subgenome_SD <- subset(subgenome_SD, no_missing <= (ncol(subgenome_SD)-5)*gmissingness)
  subgenome_SD <- subset(subgenome_SD, select=c(-no_missing))
  nonbiallelic <- as.vector(as.matrix(subgenome_SD[,5:ncol(subgenome_SD)]))
  nonbiallelic <- unique(nonbiallelic)
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/0")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/1/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "1/1/1/1")
  if (length(nonbiallelic) > 0) {
    for (i in 1:length(nonbiallelic)){
      subgenome_SD[][subgenome_SD[]==nonbiallelic[i]] <- NA
    }
  }
  subgenome_SD$no_missing <- apply(subgenome_SD, 1, function(x) sum(is.na(x)))
  subgenome_SD <- subset(subgenome_SD, no_missing <= ((ncol(subgenome_SD)-5))*gmissingness)
  subgenome_SD <- subset(subgenome_SD, select=c(-no_missing))
  names(subgenome_SD) <- gsub("X", "", names(subgenome_SD))
  names(subgenome_SD) <- gsub("_GT", "", names(subgenome_SD))
  write.table (subgenome_SD, file=paste(pop,"_4x","_rd",rd+1,".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

  subgenome_SDmaf <- subgenome_SD
  subgenome_SDmaf$freq00 <-rowSums(subgenome_SDmaf == "0/0/0/0", na.rm = TRUE)
  subgenome_SDmaf$freq01 <-rowSums(subgenome_SDmaf == "0/0/0/1", na.rm = TRUE) + rowSums(subgenome_SDmaf == "0/0/1/1", na.rm = TRUE)  + 
                           rowSums(subgenome_SDmaf == "0/1/1/1", na.rm = TRUE)
  subgenome_SDmaf$freq11 <-rowSums(subgenome_SDmaf == "1/1/1/1", na.rm = TRUE)
  subgenome_SDmaf$freq0 <- subgenome_SDmaf$freq00*2 + subgenome_SDmaf$freq01
  subgenome_SDmaf$freq1 <- subgenome_SDmaf$freq11*2 + subgenome_SDmaf$freq01
  subgenome_SDmaf <- subset(subgenome_SDmaf, select=-c(freq00,freq01,freq11))
  maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
  subgenome_SDmaf$min <- apply(subgenome_SDmaf[,(ncol(subgenome_SDmaf)-1):ncol(subgenome_SDmaf)], 1, function(x)x[maxn(2)(x)])
  subgenome_SDmaf$sum <- rowSums(subgenome_SDmaf[,c("freq0","freq1")], na.rm=TRUE)
  subgenome_SDmaf$maf <- subgenome_SDmaf$min/subgenome_SDmaf$sum
  subgenome_SDmaf0.01 <- subset(subgenome_SDmaf, maf > 0.01)
  subgenome_SDmaf0.01 <- subset(subgenome_SDmaf0.01, select=-c(freq0,freq1,min,maf,sum))
  subgenome_SDmaf0.02 <- subset(subgenome_SDmaf, maf > 0.02)
  subgenome_SDmaf0.02 <- subset(subgenome_SDmaf0.02, select=-c(freq0,freq1,min,maf,sum))
  subgenome_SDmaf0.05 <- subset(subgenome_SDmaf, maf > 0.05)
  subgenome_SDmaf0.05 <- subset(subgenome_SDmaf0.05, select=-c(freq0,freq1,min,maf,sum))
  subgenome_SDmaf0.1 <- subset(subgenome_SDmaf, maf > 0.1)
  subgenome_SDmaf0.1 <- subset(subgenome_SDmaf0.1, select=-c(freq0,freq1,min,maf,sum))
  write.table (subgenome_SDmaf, file=paste(pop,"_4x","_rd",rd+1,"_maf0.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  write.table (subgenome_SDmaf0.01, file=paste(pop,"_4x","_rd",rd+1,"_maf0.01.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  write.table (subgenome_SDmaf0.02, file=paste(pop,"_4x","_rd",rd+1,"_maf0.02.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  write.table (subgenome_SDmaf0.05, file=paste(pop,"_4x","_rd",rd+1,"_maf0.05.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  write.table (subgenome_SDmaf0.1, file=paste(pop,"_4x","_rd",rd+1,"_maf0.1.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

  if (MinorAlleleFreq != 0.02) {
    subgenome_SDmafn <- subset(subgenome_SDmaf, maf > MinorAlleleFreq)
    subgenome_SDmafn <- subset(subgenome_SDmafn, select=-c(freq0,freq1,min,maf,sum))
    write.table (subgenome_SDmafn, file=paste(pop,"_2x","_rd",rd+1,"_maf",MinorAlleleFreq,".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }

  if (MinorAlleleFreq != 0.02) {
    maf <- subset(subgenome_SDmaf, maf > MinorAlleleFreq)
    maf <- subset(maf, select="maf")
    mean <- mean(maf$maf, na.rm = TRUE)
    median <- median(maf$maf, na.rm = TRUE)
    plot <- ggplot(data=maf, aes(x=maf)) +
      geom_density(aes(y= ..scaled..), alpha=0.2, fill="cornflowerblue", colour="cornflowerblue") +
      geom_vline(aes(xintercept=mean), color="cornflowerblue", linetype="dashed", size=0.75, alpha=0.5)+
      geom_vline(aes(xintercept=median), color="tomato", linetype="dotted", size=0.75, alpha=0.5)+
      geom_text(aes(x=mean, label=paste("mean = ",round(mean, digits=2),sep=""), y=(0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=3.5) +
      geom_text(aes(x=median, label=paste("median = ",round(median, digits=2),sep=""), y=(0.5)), colour="tomato", angle=90, vjust = 1.2, size=3.5) +
      xlab("Allele Frequency") +
      ylab(paste("Density"))
    ggsave(filename=paste(pop,"_4x","_maf_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=300, compression = "lzw")

    subgenome_SDmafn$SNP <- paste (subgenome_SDmafn$CHROM,"_",subgenome_SDmafn$POS, sep="")
    subgenome_SDmafn <- subgenome_SDmafn[,c(which(colnames(subgenome_SDmafn)=="ALT"),which(colnames(subgenome_SDmafn)!="ALT"))]
    subgenome_SDmafn <- subgenome_SDmafn[,c(which(colnames(subgenome_SDmafn)=="REF"),which(colnames(subgenome_SDmafn)!="REF"))]
    subgenome_SDmafn <- subgenome_SDmafn[,c(which(colnames(subgenome_SDmafn)=="POS"),which(colnames(subgenome_SDmafn)!="POS"))]
    subgenome_SDmafn <- subgenome_SDmafn[,c(which(colnames(subgenome_SDmafn)=="CHROM"),which(colnames(subgenome_SDmafn)!="CHROM"))]
    subgenome_SDmafn <- subgenome_SDmafn[,c(which(colnames(subgenome_SDmafn)=="SNP"),which(colnames(subgenome_SDmafn)!="SNP"))]
    write.table (subgenome_SDmafn, file=paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_binary.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

    if (snpformats == "true") {
      alleles <- unique(subset(subgenome_SDmafn, select=c(4,5)))
      rownames(alleles) <- NULL
      alleles[] <- lapply(alleles, as.character)
      alleles$ref0 <- alleles$REF; alleles$alt0 <- alleles$ALT
      for (i in 1:nrow(alleles)) {
        alleles[i,3] <- gsub(",.*", "", alleles[i,3])
        alleles[i,4] <- gsub(",.*", "", alleles[i,4])
      }
      geno <- NULL
      for (i in 1:nrow(alleles)) {
        ref <- as.vector(alleles[i,3]); REFsub <- as.vector(alleles[i,1])
        alt <- as.vector(alleles[i,4]); ALTsub <- as.vector(alleles[i,2])
        snplen = nchar(ref) + nchar(alt)
        output <- subset(subgenome_SDmafn, REF == REFsub & ALT == ALTsub)
        output[] <- lapply(output, as.character)
        if (snplen == 2) {
          for (j in 1:nrow(output)) {
            output[j,6:ncol(output)] <- gsub("0", ref, output[j,6:ncol(output)])
            output[j,6:ncol(output)] <- gsub("1", alt, output[j,6:ncol(output)])
          }
        }else{
          for (j in 1:nrow(output)) {
            if (nchar(ref) == 1) {
              output[j,6:ncol(output)] <- gsub("0", "-", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("0", "+", output[j,6:ncol(output)])}
            if (nchar(alt) == 1) {
              output[j,6:ncol(output)] <- gsub("1", "-", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("1", "+", output[j,6:ncol(output)])}
            if (nchar(ref) > 1 & nchar(alt) > 1) {
              output[j,6:ncol(output)] <- gsub("0", "+", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("1", "-", output[j,6:ncol(output)])}
          }
        }
        output[] <- lapply(output, as.character)
        for (k in 1:nrow(output)) {
          output[k,6:ncol(output)] <- gsub("/", "", output[k,6:ncol(output)])
        }
        geno <- rbind(geno,output); gc()
      }
      geno <- geno[order(geno$CHROM, geno$POS),]
      write.table (geno, file=paste(pop,"_4x","_rd",rd+1,"_maf0.02_nucleotide.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    }
    
    subgenome_SDmafn <- data.frame(lapply(subgenome_SDmafn, as.character), stringsAsFactors=FALSE, check.names = FALSE)
    subgenome_SDmafn[][subgenome_SDmafn[]=="0/0/0/0"] <- "0"
    subgenome_SDmafn[][subgenome_SDmafn[]=="0/0/0/1"] <- "1"
    subgenome_SDmafn[][subgenome_SDmafn[]=="0/0/1/1"] <- "2"
    subgenome_SDmafn[][subgenome_SDmafn[]=="0/1/1/1"] <- "3"
    subgenome_SDmafn[][subgenome_SDmafn[]=="1/1/1/1"] <- "4"
    write.table (subgenome_SDmafn, file=paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }

  if (MinorAlleleFreq == 0.02) {
    maf <- subset(subgenome_SDmaf, maf > 0.02)
    maf <- subset(maf, select="maf")
    mean <- mean(maf$maf, na.rm = TRUE)
    median <- median(maf$maf, na.rm = TRUE)
    plot <- ggplot(data=maf, aes(x=maf)) +
      geom_density(aes(y= ..scaled..), alpha=0.2, fill="cornflowerblue", colour="cornflowerblue") +
      geom_vline(aes(xintercept=mean), color="cornflowerblue", linetype="dashed", size=0.75, alpha=0.5)+
      geom_vline(aes(xintercept=median), color="tomato", linetype="dotted", size=0.75, alpha=0.5)+
      geom_text(aes(x=mean, label=paste("mean = ",round(mean, digits=2),sep=""), y=(0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=3.5) +
      geom_text(aes(x=median, label=paste("median = ",round(median, digits=2),sep=""), y=(0.5)), colour="tomato", angle=90, vjust = 1.2, size=3.5) +
      xlab("Allele Frequency") +
      ylab(paste("Density"))
    ggsave(filename=paste(pop,"_4x","_maf_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=300, compression = "lzw")

    subgenome_SDmaf0.02$SNP <- paste (subgenome_SDmaf0.02$CHROM,"_",subgenome_SDmaf0.02$POS, sep="")
    subgenome_SDmaf0.02 <- subgenome_SDmaf0.02[,c(which(colnames(subgenome_SDmaf0.02)=="ALT"),which(colnames(subgenome_SDmaf0.02)!="ALT"))]
    subgenome_SDmaf0.02 <- subgenome_SDmaf0.02[,c(which(colnames(subgenome_SDmaf0.02)=="REF"),which(colnames(subgenome_SDmaf0.02)!="REF"))]
    subgenome_SDmaf0.02 <- subgenome_SDmaf0.02[,c(which(colnames(subgenome_SDmaf0.02)=="POS"),which(colnames(subgenome_SDmaf0.02)!="POS"))]
    subgenome_SDmaf0.02 <- subgenome_SDmaf0.02[,c(which(colnames(subgenome_SDmaf0.02)=="CHROM"),which(colnames(subgenome_SDmaf0.02)!="CHROM"))]
    subgenome_SDmaf0.02 <- subgenome_SDmaf0.02[,c(which(colnames(subgenome_SDmaf0.02)=="SNP"),which(colnames(subgenome_SDmaf0.02)!="SNP"))]
    write.table (subgenome_SDmaf0.02, file=paste(pop,"_4x","_rd",rd+1,"_maf0.02_binary.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

    if ( snpformats == "true") {
      alleles <- unique(subset(subgenome_SDmaf0.02, select=c(4,5)))
      rownames(alleles) <- NULL
      alleles[] <- lapply(alleles, as.character)
      alleles$ref0 <- alleles$REF; alleles$alt0 <- alleles$ALT
      for (i in 1:nrow(alleles)) {
        alleles[i,3] <- gsub(",.*", "", alleles[i,3])
        alleles[i,4] <- gsub(",.*", "", alleles[i,4])
      }
      geno <- NULL
      for (i in 1:nrow(alleles)) {
        ref <- as.vector(alleles[i,3]); REFsub <- as.vector(alleles[i,1])
        alt <- as.vector(alleles[i,4]); ALTsub <- as.vector(alleles[i,2])
        snplen = nchar(ref) + nchar(alt)
        output <- subset(subgenome_SDmaf0.02, REF == REFsub & ALT == ALTsub)
        output[] <- lapply(output, as.character)
        if (snplen == 2) {
          for (j in 1:nrow(output)) {
            output[j,6:ncol(output)] <- gsub("0", ref, output[j,6:ncol(output)])
            output[j,6:ncol(output)] <- gsub("1", alt, output[j,6:ncol(output)])
          }
        }
        if (snplen != 2) {
          for (j in 1:nrow(output)) {
            if (nchar(ref) == 1) {
              output[j,6:ncol(output)] <- gsub("0", "-", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("0", "+", output[j,6:ncol(output)])}
            if (nchar(alt) == 1) {
              output[j,6:ncol(output)] <- gsub("1", "-", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("1", "+", output[j,6:ncol(output)])}
            if (nchar(ref) > 1 & nchar(alt) > 1) {
              output[j,6:ncol(output)] <- gsub("0", "+", output[j,6:ncol(output)])
            }else{output[j,6:ncol(output)] <- gsub("1", "-", output[j,6:ncol(output)])}
          }
        }
        output[] <- lapply(output, as.character)
        for (k in 1:nrow(output)) {
          output[k,6:ncol(output)] <- gsub("/", "", output[k,6:ncol(output)])
        }
        geno <- rbind(geno,output); gc()
      }
      geno <- geno[order(geno$CHROM, geno$POS),]
      write.table (geno, file=paste(pop,"_4x","_rd",rd+1,"_maf0.02_nucleotide.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    }
    
    subgenome_SDmaf0.02 <- data.frame(lapply(subgenome_SDmaf0.02, as.character), stringsAsFactors=FALSE, check.names = FALSE)
    subgenome_SDmaf0.02[][subgenome_SDmaf0.02[]=="0/0/0/0"] <- "0"
    subgenome_SDmaf0.02[][subgenome_SDmaf0.02[]=="0/0/0/1"] <- "1"
    subgenome_SDmaf0.02[][subgenome_SDmaf0.02[]=="0/0/1/1"] <- "2"
    subgenome_SDmaf0.02[][subgenome_SDmaf0.02[]=="0/1/1/1"] <- "3"
    subgenome_SDmaf0.02[][subgenome_SDmaf0.02[]=="1/1/1/1"] <- "4"
    write.table (subgenome_SDmaf0.02, file=paste(pop,"_4x","_rd",rd+1,"_maf0.02_dose.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }

  sumfreq <- read.table(paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = FALSE)
  sumfreq <- subset(sumfreq, select=-c(1:5))
  SNP <- sumfreq
  SNP$percent <- (apply(SNP, 1, function(x) sum(is.na(x))))/ncol(SNP)*100
  SNP <- subset(SNP, select=c(percent))
  plot <- ggplot(data=SNP, aes(x=percent)) +
    geom_bar(color="gray60", fill="white", alpha=0.1) +
    geom_vline(aes(xintercept=mean(percent)),color="cornflowerblue", linetype="dashed", size=0.5, alpha=0.5) +
    geom_vline(aes(xintercept=median(percent)),color="tomato", linetype="dotted", size=0.5, alpha=0.5) +
    geom_text(aes(x=mean(percent), label=paste("mean = ",round(mean(percent)),sep=""), y=(max(table(SNP))*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=3) +
    geom_text(aes(x=median(percent), label=paste("median = ",round(median(percent)),sep=""), y=(max(table(SNP))*0.5)), colour="tomato", angle=90, vjust = 1.2, size=3) +
    xlim(-5,105) +
    labs(title="missing rate per variant",x="Percent", y = "Count")
  ggsave(filename=paste(pop,"_4x","_Variant_missing_rate_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=300, compression = "lzw")
  sample <- as.data.frame(t(sumfreq))
  sample$percent <- (apply(sample, 1, function(x) sum(is.na(x))))/ncol(sample)*100
  sample <- subset(sample, select=c(percent))
  plot <- ggplot(data=sample, aes(x=percent)) +
    geom_bar(color="gray60", fill="white", alpha=0.1) +
    geom_vline(aes(xintercept=mean(percent)),color="cornflowerblue", linetype="dashed", size=0.5, alpha=0.5) +
    geom_vline(aes(xintercept=median(percent)),color="tomato", linetype="dotted", size=0.5, alpha=0.5) +
    geom_text(aes(x=mean(percent), label=paste("mean = ",round(mean(percent)),sep=""), y=(max(table(sample))*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=3) +
    geom_text(aes(x=median(percent), label=paste("median = ",round(median(percent)),sep=""), y=(max(table(sample))*0.5)), colour="tomato", angle=90, vjust = 1.2, size=3) +
    xlim(-5,105) +
    labs(title="missing rate per sample",x="Percent", y = "Count")
  ggsave(filename=paste(pop,"_4x","_sample_missing_rate_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=300, compression = "lzw")
}
RD_snpfiltering()

####################################################################################################################
subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
subgenome_1[, 5:(((ncol(subgenome_1)-4)/2)+4)] <- lapply(5:(((ncol(subgenome_1)-4)/2)+4), function(x) as.numeric(subgenome_1[[x]]))
rd_boxplot <- function() {
  #######################################################################################################################################################################################
  # Now, let's make boxplots use raw data set. It reflects read depth distribution across sample IDs
  for (t in c(0,rd)) {
    subgenome_1_plots <- subgenome_1[c(5:(((ncol(subgenome_1)-4)/2)+4))]
    subgenome_1_plots[][subgenome_1_plots[]=="0"] <- NA
    subgenome_1_plots$no_missing <- apply(subgenome_1_plots, 1, function(x) sum(is.na(x)))
    subgenome_1_plots <- subset(subgenome_1_plots, no_missing <= (ncol(subgenome_1_plots)-5)*gmissingness)
    subgenome_1_plots <- subset(subgenome_1_plots, select=-c(no_missing))
    subgenome_1_plots[][subgenome_1_plots[] <= t] <- NA
    names(subgenome_1_plots) <- gsub(paste("_DP",sep=""), "", names(subgenome_1_plots))
    subgenome_1_plots$no_missing <- apply(subgenome_1_plots, 1, function(x) sum(is.na(x)))
    subgenome_1_plots <- subset(subgenome_1_plots, no_missing <= (ncol(subgenome_1_plots)-5)*gmissingness)
    subgenome_1_plots <- subset(subgenome_1_plots, select=c(-no_missing))
    subgenome_1_plots <- as.matrix(subgenome_1_plots)
    subgenome_1_plots <- as.data.frame(as.table(subgenome_1_plots))
    subgenome_1_plots <- na.omit(subgenome_1_plots)
    subgenome_1_boxplot <- subset(subgenome_1_plots, select=c(2,3))
    names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Var2"] <- "samples"
    names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Freq"] <- "DP"

    subgenome_1_boxplot$DP <- as.numeric(as.character(subgenome_1_boxplot$DP))
    subgenome_1_boxplot <- na.omit(subgenome_1_boxplot)
    quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.95), na.rm= TRUE)
    nsamples <- length(unique(subgenome_1_boxplot$samples))
    if (nsamples <= 20) { nsamples <- 36 } else { nsamples <- sqrt(22500/nsamples) }
    boxplot <- ggplot(subgenome_1_boxplot, aes(x = reorder(samples,DP, na.rm = TRUE), y=DP), stat='identity')+
      geom_boxplot(fill="white", colour="cornflowerblue",
                   outlier.alpha = 0.01, outlier.colour="tomato", outlier.size=1.0)+
      coord_flip()+
      scale_y_continuous(expand = c(0,0), limits = c(0,quantile999)) +
      theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
            axis.text.y=element_text(colour="cornflowerblue", size=nsamples),
            axis.title=element_text(size=36)) +
      xlab(paste(pop," Diversity Population",sep="")) +
      ylab("Read Depth (4x Genotypes)")
    ggsave(filename= paste(pop,"_4x","_boxplot_rd",t+1,".tiff",sep=""), plot=boxplot, width=15, height= 25, dpi=300, compression = "lzw")

    meanDP <- mean(subgenome_1_boxplot$DP, na.rm=T)
    medianDP <- median(subgenome_1_boxplot$DP, na.rm=T)
    maxDP <- max(table(subgenome_1_boxplot$DP))
    quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.95), na.rm= TRUE)
    subgenome_1_dist <- as.data.frame(table(subgenome_1_boxplot$DP))
    subgenome_1_dist$Var1 <- as.numeric(as.character(subgenome_1_dist$Var1))
    subgenome_1_distm <- subset(subgenome_1_dist, Freq > quantile999)
    if (nrow(subgenome_1_distm) == 0) {maxX <- max(subgenome_1_dist[,1])} else {maxX <- max(subgenome_1_distm[,1])}
    if (maxX >= 100) { intervalm <- round(maxX/10,-1) } else { intervalm <- round(maxX/10,0) } 
    histogram <- ggplot(subgenome_1_dist, aes(x=Var1, y=Freq)) +
      geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="cornflowerblue", fill="white")+
      geom_vline(aes(xintercept=meanDP), color="cornflowerblue", linetype="dashed", size=3, alpha=0.5)+
      geom_vline(aes(xintercept=medianDP), color="tomato", linetype="dotted", size=3, alpha=0.5)+
      geom_text(aes(x=meanDP, label=paste("mean = ",round(meanDP),sep=""), y=(maxDP*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=7.5) +
      geom_text(aes(x=medianDP, label=paste("median = ",round(medianDP),sep=""), y=(maxDP*0.5)), colour="tomato", angle=90, vjust = 1.2, size=7.5) +
      scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",scientific = T)) +
      scale_x_continuous(breaks=seq(0,maxX,intervalm)) +
      theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
            axis.text.y=element_text(colour="cornflowerblue", size=24),
            axis.title=element_text(size=30)) +
      ylab("Frequency") +
      xlab(paste("Read Depth Distribution (", "Diversity Population)",sep=""))
    ggsave(filename= paste(pop,"_4x","_histogram_rd",t+1,".tiff",sep=""), plot=histogram, width=25, height= 15, dpi=300, compression = "lzw")

    boxplot <- NULL
    subgenome_1_boxplot <- NULL
    subgenome_1_plots <- NULL
    gc()
  }
  # Extract read depth values specifically for filtered SNPs, then plot boxplot
  # Also, plot histogram of read depth across data set
  subgenome_1_plots <- subgenome_1
  subgenome_final <- read.table (file=paste(pop,"_4x_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), header=T, sep="\t", check.names = FALSE)
  subgenome_final <- as.data.frame(subgenome_final[,2:3])
  subgenome_1_plots <- merge(subgenome_1_plots, subgenome_final, by=c("CHROM","POS"), all.y=TRUE)
  subgenome_1_plots[][subgenome_1_plots[] < 12] <- NA
  for (i in 5:(((ncol(subgenome_1_plots)-4)/2)+4)) {
    j <- i+((ncol(subgenome_1_plots)-4)/2)
    subgenome_1_plots[,i][subgenome_1_plots[,j] != "0/0/0/0" && subgenome_1_plots[,j] != "1/1/1/1" &&
                            subgenome_1_plots[,i] <= rd] <- NA
    gc()
  }
  subgenome_1_plots <- subgenome_1_plots[c(5:(((ncol(subgenome_1_plots)-4)/2)+4))]
  names(subgenome_1_plots) <- gsub(paste("_DP",sep=""), "", names(subgenome_1_plots))
  subgenome_1_plots <- as.matrix(subgenome_1_plots)
  subgenome_1_plots <- as.data.frame(as.table(subgenome_1_plots))
  subgenome_1_plots <- na.omit(subgenome_1_plots)
  subgenome_1_boxplot <- subset(subgenome_1_plots, select=c(2,3))
  names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Var2"] <- "samples"
  names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Freq"] <- "DP"

  subgenome_1_boxplot$DP <- as.numeric(as.character(subgenome_1_boxplot$DP))
  subgenome_1_boxplot <- na.omit(subgenome_1_boxplot)
  quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.95), na.rm= TRUE)
  nsamples <- length(unique(subgenome_1_boxplot$samples))
  if (nsamples <= 20) { nsamples <- 36 } else { nsamples <- sqrt(22500/nsamples) }
  boxplot <- ggplot(subgenome_1_boxplot, aes(x = reorder(samples,DP, na.rm = TRUE), y=DP), stat='identity')+
    geom_boxplot(fill="white", colour="cornflowerblue",
                 outlier.alpha = 0.01, outlier.colour="tomato", outlier.size=1.0)+
    coord_flip()+
    scale_y_continuous(expand = c(0,0), limits = c(0,quantile999)) +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
          axis.text.y=element_text(colour="cornflowerblue", size=nsamples),
          axis.title=element_text(size=36)) +
    xlab(paste(pop," Diversity Population",sep="")) +
    ylab("Read Depth (4x Genotypes)")
  ggsave(filename= paste(pop,"_4x","_boxplot_filtered.tiff",sep=""), plot=boxplot, width=15, height= 25, dpi=300, compression = "lzw")

  meanDP <- mean(subgenome_1_boxplot$DP, na.rm=T)
  medianDP <- median(subgenome_1_boxplot$DP, na.rm=T)
  maxDP <- max(table(subgenome_1_boxplot$DP))
  quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.95), na.rm= TRUE)
  subgenome_1_dist <- as.data.frame(table(subgenome_1_boxplot$DP))
  subgenome_1_dist$Var1 <- as.numeric(as.character(subgenome_1_dist$Var1))
  subgenome_1_distm <- subset(subgenome_1_dist, Freq > quantile999)
  if (nrow(subgenome_1_distm) == 0) {maxX <- max(subgenome_1_dist[,1])} else {maxX <- max(subgenome_1_distm[,1])}
  if (maxX >= 100) { intervalm <- round(maxX/10,-1) } else { intervalm <- round(maxX/10,0) } 
  histogram <- ggplot(subgenome_1_dist, aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="cornflowerblue", fill="white")+
    geom_vline(aes(xintercept=meanDP), color="cornflowerblue", linetype="dashed", size=3, alpha=0.5)+
    geom_vline(aes(xintercept=medianDP), color="tomato", linetype="dotted", size=3, alpha=0.5)+
    geom_text(aes(x=meanDP, label=paste("mean = ",round(meanDP),sep=""), y=(maxDP*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=7.5) +
    geom_text(aes(x=medianDP, label=paste("median = ",round(medianDP),sep=""), y=(maxDP*0.5)), colour="tomato", angle=90, vjust = 1.2, size=7.5) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",scientific = T)) +
    scale_x_continuous(breaks=seq(0,maxX,intervalm)) +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
          axis.text.y=element_text(colour="cornflowerblue", size=24),
          axis.title=element_text(size=30)) +
    ylab("Frequency") +
    xlab(paste("Read Depth Distribution (", "Diversity Population)",sep=""))
  ggsave(filename= paste(pop,"_4x","_histogram_filtered.tiff",sep=""), plot=histogram, width=25, height= 15, dpi=300, compression = "lzw")

  boxplot <- NULL
  subgenome_1_boxplot <- NULL
  subgenome_1_plots <- NULL
  gc()
}
rd_boxplot()
raw_alleles <- function() {
  #######################################################################################################################################################################################
  # Let's plot the distribution of multi-allelic variants

  subgenome_1_plots <- read.table (file=paste(pop,"_4x_","rawRD",rd+1,"_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
  subgenome_1_plots <- subset(subgenome_1_plots, select=c((((ncol(subgenome_1_plots)-4)/2)+5):ncol(subgenome_1_plots)))
  Multiallelic <- as.data.frame(table(as.matrix(subgenome_1_plots)))
  names(Multiallelic)[names(Multiallelic) == "Var1"] <- "Genotype"
  Multiallelic <- subset(Multiallelic, Genotype != "./././.")
  Multiallelic[][Multiallelic[]=="NA"] <- "0"
  Multiallelic$Freq <- round((Multiallelic$Freq)/((ncol(subgenome_1)-4)/2),0)
  sum <- sum(Multiallelic$Freq)
  max <- max(Multiallelic$Freq)
  max <- max*1.3
  Multiallelic$percentage <- ((Multiallelic$Freq)/sum)*100
  Multiallelic <- subset(Multiallelic, percentage >= 0.01)
  Multiallelic[,3] <- round(Multiallelic[,3], 2)
  Multiallelic$length0 <- lengths(regmatches(Multiallelic$Genotype, gregexpr("0", Multiallelic$Genotype)))
  Multiallelic$length1 <- lengths(regmatches(Multiallelic$Genotype, gregexpr("1", Multiallelic$Genotype)))
  Multiallelic <- Multiallelic[order(Multiallelic$length0, Multiallelic$length1, Multiallelic$Genotype),]
  Multiallelic$Genotype <- factor(Multiallelic$Genotype, levels = Multiallelic$Genotype)
  plot <- ggplot(Multiallelic, aes(x = Genotype, y=Freq, fill=Genotype, group=Genotype)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9,
             colour="black")+
    geom_text(aes(x=Genotype, y=Freq, label = paste(percentage, "%"), group=Genotype),
              position=position_dodge(0.95), hjust = -0.1, size=4, color="black", fontface="italic")+
    theme(axis.text.x=element_text(colour="cornflowerblue", size=12),
          axis.text.y=element_text(colour="cornflowerblue", size=12),
          axis.title=element_text(size=14)) +
    coord_flip()+
    theme(legend.text=element_text(size=12)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.4, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab("Genotypes") +
    ylab(paste("Proportion of Genotypes (",sum,")", sep=""))
  ggsave(filename=paste(pop,"_4x_rawRD",rd+1,"_Variants.tiff",sep=""), plot=plot, width=7.5, height= 5, dpi=300, compression = "lzw")
}
raw_alleles()

####################################################################################################################
copy_filter <- function(){
  dir.create(unique_mapped)
  snpid <- read.table(paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = F)
  snpid <- subset(snpid, select=c("CHROM","POS"))
  snpid$POS <- as.numeric(as.character(snpid$POS))
  snpid$position <- lapply(snpid$POS / 100, as.integer)
  snpidN <- read.table("../../alignment_summaries/refgenome_paralogs.txt", header=T, sep="\t", quote="", check.names=FALSE, fill=F)
  snpidN <- snpidN[!(snpidN$CHROM=="" | snpidN$POS=="" | snpidN$nloci==""), ]
  snpidN$POS <- as.numeric(as.character(snpidN$POS))
  snpidN$position <- lapply(snpidN$POS / 100, as.integer); snpidN1 <- snpidN; snpidN2 <- snpidN
  snpidN1$position <- as.numeric(as.character(snpidN1$position)); snpidN2$position <- as.numeric(as.character(snpidN2$position))
  snpidN1$position <- snpidN1$position + 1; snpidN2$position <- snpidN2$position - 1
  snpidN <- rbind(snpidN,snpidN1); snpidN <- rbind(snpidN,snpidN2)
  snpid$match <- paste(snpid$CHROM,snpid$position,sep="_")
  snpidN$match <- paste(snpidN$CHROM,snpidN$position,sep="_")
  snpid$nloci <- snpidN$nloci[ match(snpid$match, snpidN$match)]
  snpidM <- subset(snpid, select=-c(position,match))
  write.table (snpidM, file=paste("./unique_mapped/",pop,"_4x","refgenome_nloci_matched.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  
  snpidM <- subset(snpidM, snpidM$nloci <= hap)
  snpidM <- subset(snpidM, select=c("CHROM","POS"))
  snpid <- read.table(paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = F)
  subgenome_uniqmap <- merge(snpid, snpidM, by = c("CHROM","POS"), all.y = TRUE)
  write.table (subgenome_uniqmap, file=paste("./unique_mapped/",pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose_unique_mapped.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
}
copy_filter()

####################################################################################################################
pop_struc <- function() {
  for (i in c(0.01,0.05,0.1)) {
    pop_data <- read.table(paste(pop,"_4x","_rd",rd+1,"_maf",MinorAlleleFreq,"_dose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = FALSE)
    pop_data <- subset(pop_data, select=-c(1:5))
    pop_data$no_missing <- apply(pop_data, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
    pop_data <- subset(pop_data, no_missing <= ncol(pop_data)*i)
    pop_data <- subset(pop_data, select=-c(no_missing))
    if (nrow(pop_data) >= 1000) {
      print(paste("missing rate = ",i))
      break
    }
  }
  normalize_kinmat <- function(kinmat){
    #normalize kinship so that Kij \in [0,1]
    tmp=kinmat - min(kinmat)
    tmp=tmp/max(tmp)
    tmp[1:9,1:9]
    #fix eigenvalues to positive
    diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
    tmp[1:9,1:9]  
    return(tmp)
  }
  if (nrow(pop_data) >= 100) {
    pop_data <- as.matrix(t(pop_data))
    
    #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
    Gmatrix <- function(SNPmatrix = pop_data, method = "VanRaden", missingValue = NA, 
                        maf = 0, thresh.missing = 0.1, verify.posdef = FALSE, ploidy = 4, 
                        pseudo.diploid = FALSE, integer = FALSE, ratio = FALSE, impute.method = TRUE, 
                        ratio.check = FALSE) {
      Time = proc.time()
      
      if(ratio){ #This allows to enter in the scaled crossprod condition
        method="VanRaden"
        ploidy=8
      }
      
      if (!is.na(missingValue)) {
        m <- match(SNPmatrix, missingValue, 0)
        SNPmatrix[m > 0] <- NA
      }
      
      # Internal function to check input Gmatrix arguments
      check_Gmatrix_data <- function(SNPmatrix,ploidy,method, ratio=FALSE, integer=TRUE){
        if (is.null(SNPmatrix)) {
          stop(deparse("Please define the variable SNPdata"))
        }
        if (all(method != c("Yang", "VanRaden", "Slater", "Su", "Vitezica", "MarkersMatrix","Endelman"))) {
          stop("Method to build Gmatrix has to be either `Yang` or `VanRaden` for marker-based additive relationship matrix, or `Su` or `Vitezica` or `Endelman` for marker-based dominance relationship matrx, or `MarkersMatrix` for matrix with amount of shared-marks by individuals pairs")
        }
        
        #  if( method=="Yang" && ploidy>2)
        #    stop("Change method to 'VanRaden' for ploidies higher than 2 for marker-based additive relationship matrix")
        
        if( method=="Su" && ploidy>2)
          stop("Change method to 'Slater' for ploidies higher than 2 for marker-based non-additive relationship matrix")
        
        if( method=="Vitezica" && ploidy>2)
          stop("Change method to 'Slater' for ploidies higher than 2 for marker-based non-additive relationship matrix")
        
        if(class(SNPmatrix)!="matrix"){
          cat("SNPmatrix class is:",class(SNPmatrix),"\n")
          stop("SNPmatrix class must be matrix. Please verify it.")
        }
        
        if(!ratio){
          if( ploidy > 20 | (ploidy %% 2) != 0)
            stop(deparse("Only even ploidy from 2 to 20"))
          
          t <- max(SNPmatrix,na.rm = TRUE)
          if( t > ploidy )
            stop(deparse("Check your data, it has values above ploidy number"))
          
          t <- min(SNPmatrix,na.rm=TRUE)
          if( t < 0 )
            stop(deparse("Check your data, it has values under 0"))
          
          if(integer)
            if(prod(SNPmatrix == round(SNPmatrix),na.rm = TRUE)==0)
              stop(deparse("Check your data, it has not integer values"))
        }
        
        if(ratio){
          t <- max(SNPmatrix,na.rm = TRUE)
          if( t > 1)
            stop(deparse("Check your data, it has values above 1. It is expected a ratio values [0;1]."))
          
          t <- min(SNPmatrix,na.rm=TRUE)
          if( t < 0 )
            stop(deparse("Check your data, it has values under 0. It is expected a ratio values [0;1]."))
        }
      }
      check_Gmatrix_data(SNPmatrix=SNPmatrix,method=method,ploidy=ploidy,ratio=ratio,integer=integer)
      
      NumberMarkers <- ncol(SNPmatrix)
      nindTotal <- colSums(!is.na(SNPmatrix))
      nindAbs <- max(nindTotal)
      cat("Initial data: \n")
      cat("\tNumber of Individuals:", max(nindTotal), "\n")
      cat("\tNumber of Markers:", NumberMarkers, "\n")
      
      # Function by Luis F. V. Ferrao
      # Internal function to maf cutoff and impute data
      snp.check = function(M = NULL,
                           ploidy=4,
                           thresh.maf = 0.05,
                           thresh.missing = 0.9,
                           impute.method = "mean"){
        # SNP missing data
        ncol.init <- ncol(M)
        
        missing <- apply(M, 2, function(x) sum(is.na(x))/nrow(M))
        missing.low = missing <= thresh.missing
        cat("\nMissing data check: \n")
        if(any(missing.low)){
          cat("\tTotal SNPs:", ncol(M),"\n")
          cat("\t",ncol(M) - sum(missing.low), "SNPs dropped due to missing data threshold of", thresh.missing,"\n")
          cat("\tTotal of:",sum(missing.low), " SNPs \n")
          idx.rm <- which(missing.low)
          M <- M[, idx.rm, drop=FALSE]
        } else{
          cat("\tNo SNPs with missing data, missing threshold of = ", thresh.missing,"\n")
        }
        
        # Minor alele frequency
        MAF <- apply(M, 2, function(x) {
          AF <- mean(x, na.rm = T)/ploidy
          MAF <- ifelse(AF > 0.5, 1 - AF, AF) # Minor allele freq can be ref allele or not
        })
        snps.low <- MAF < thresh.maf
        cat("MAF check: \n")
        if(any(snps.low)){
          cat("\t",sum(snps.low), "SNPs dropped with MAF below", thresh.maf,"\n")
          cat("\tTotal:",ncol(M) - sum(snps.low), " SNPs \n")
          idx.rm <- which(snps.low)
          M <- M[, -idx.rm, drop=FALSE]
        } else{
          cat("\tNo SNPs with MAF below", thresh.maf,"\n")
        }
        
        # SNPs monomorficos
        mono <- apply(M, 2, function(x) {
          equal <- isTRUE(all.equal(x, rep(x[1], length(x))))
        })
        cat("Monomorphic check: \n")
        if(any(mono)){
          cat("\t",sum(mono), "monomorphic SNPs \n")
          cat("\tTotal:",ncol(M) - sum(mono), "SNPs \n")
          idx.rm <- which(mono)
          M <- M[, -idx.rm, drop=FALSE]
        } else{
          cat("\tNo monomorphic SNPs \n")
        }
        
        # Imputing by mode
        if(impute.method=="mean"){
          ix <- which(is.na(M))
          if (length(ix) > 0) {
            M[ix] <- mean(M,na.rm = TRUE)
          }
        }
        
        if(impute.method=="mode"){
          ix <- which(is.na(M))
          if (length(ix) > 0) {
            M[ix] <- as.integer(names(which.max(table(M))))
          }
        }
        
        # Total of SNPs
        cat("Summary check: \n")
        cat("\tInitial: ", ncol.init, "SNPs \n")
        cat("\tFinal: ", ncol(M), " SNPs (", ncol.init - ncol(M), " SNPs removed) \n \n")
        return(M)
      }
      
      if(ratio==FALSE){
        SNPmatrix <- snp.check(SNPmatrix,
                               ploidy = ploidy, 
                               thresh.maf = maf, 
                               thresh.missing = thresh.missing,
                               impute.method = impute.method)
      }
      
      ## Testing ratio check function: not final!
      if(ratio && ratio.check){
        SNPmatrix <- snp.check(SNPmatrix,
                               ploidy = ploidy, 
                               thresh.maf = maf, 
                               thresh.missing = thresh.missing,
                               impute.method = impute.method)
      }
      
      ## Internal Functions ##
      # Coding SNPmatrix as Slater (2016) Full autotetraploid model including non-additive effects (Presence/Absence per Genotype per Marker)
      slater_par <- function(X,ploidy){
        prime.index <- c(3,5,7,11,13,17,19,23,29,31,37,
                         41,43,47,53,59,61,67,71,73,79)
        
        NumberMarkers <- ncol(X)
        nindTotal <- nrow(X)
        X <- X+1
        
        ## Breaking intervals to use less RAM
        temp <- seq(1,NumberMarkers,10000)
        temp <- cbind(temp,temp+9999)
        temp[length(temp)] <- NumberMarkers
        prime.index <- prime.index[1:(ploidy+1)]
        
        ## Uses Diagonal (which is Sparse mode, uses less memmory)
        for(i in 1:nrow(temp)){
          X.temp <- X[,c(temp[i,1]:temp[i,2])]
          NumberMarkers <- ncol(X.temp)
          X.temp <- X.temp %*% t(kronecker(diag(NumberMarkers),prime.index))
          X.temp[which(as.vector(X.temp) %in%
                         c(prime.index*c(1:(ploidy+1))))] <- 1
          X.temp[X.temp!=1] <- 0
          if(i==1){
            X_out <- X.temp
          }else{
            X_out <- cbind(X_out,X.temp)
          }   
        }
        gc()
        return(X_out)
      }
      
      if(method=="Slater"){
        P <- colSums(SNPmatrix,na.rm = TRUE)/nrow(SNPmatrix)
        SNPmatrix[,which(P>ploidy/2)] <- ploidy-SNPmatrix[,which(P>(ploidy/2))]
        SNPmatrix <- slater_par(SNPmatrix,ploidy=ploidy)
        NumberMarkers <- ncol(SNPmatrix)
        Frequency <- colSums(SNPmatrix,na.rm=TRUE)/nrow(SNPmatrix)
        FreqP <- matrix(rep(Frequency, each = nrow(SNPmatrix)), 
                        ncol = ncol(SNPmatrix))
      }
      
      if(ploidy==2){
        alelleFreq <- function(x, y) {
          (2 * length(which(x == y)) + length(which(x == 1)))/(2 * 
                                                                 length(which(!is.na(x))))
        }
        Frequency <- cbind(apply(SNPmatrix, 2, function(x) alelleFreq(x,0))
                           , apply(SNPmatrix, 2, function(x) alelleFreq(x, 2)))
        
        #   if (any(Frequency[, 1] <= maf) & maf != 0) {
        #      cat("\t", length(which(Frequency[, 1] <= maf)), "markers dropped due to maf cutoff of", maf, "\n")
        #      SNPmatrix <- SNPmatrix[,-which(Frequency[, 1] <= maf)]
        #      cat("\t", ncol(SNPmatrix), "markers kept \n")
        #      Frequency <- as.matrix(Frequency[-which(Frequency[,1] <= 
        #                                                maf), ])
        #      NumberMarkers <- ncol(SNPmatrix)
        #    }
        FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
                        ncol = ncol(SNPmatrix))
      }
      
      if(ploidy>2 && pseudo.diploid){## Uses Pseudodiploid model
        P <- colSums(SNPmatrix,na.rm = TRUE)/nrow(SNPmatrix)
        SNPmatrix[,which(P>ploidy/2)] <- ploidy-SNPmatrix[,which(P>(ploidy/2))]
        Frequency <- colSums(SNPmatrix,na.rm=TRUE)/(ploidy*nrow(SNPmatrix))
        Frequency <- cbind(1-Frequency,Frequency)
        FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
                        ncol = ncol(SNPmatrix))
        SNPmatrix[SNPmatrix %in% c(1:(ploidy-1))] <- 1
        SNPmatrix[SNPmatrix==ploidy] <- 2
      }
      
      if (method == "MarkersMatrix") {
        Gmatrix <- !is.na(SNPmatrix)
        Gmatrix <- tcrossprod(Gmatrix, Gmatrix)
        return(Gmatrix)
      }
      
      ## VanRaden ##
      if (method == "VanRaden") {
        if(ploidy==2){
          TwoPQ <- 2 * t(Frequency[, 1]) %*% Frequency[, 2]
          SNPmatrix <- SNPmatrix- 2 * FreqP
          SNPmatrix[is.na(SNPmatrix)] <- 0
          Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(TwoPQ)
        }else{
          if(ploidy>2){
            SNPmatrix<-scale(SNPmatrix,center=TRUE,scale=FALSE) 
            K<-sum(apply(X=SNPmatrix,FUN=var,MARGIN=2,na.rm=TRUE))
            SNPmatrix[which(is.na(SNPmatrix))] <- 0
            Gmatrix<-tcrossprod(SNPmatrix)/K
          }
        }
      }
      
      if (method == "Yang") {
        FreqPQ <- matrix(rep(2 * Frequency[, 1] * Frequency[, 
                                                            2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
        G.all <- (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 
                    2 * (FreqP^2))/FreqPQ
        G.ii <- as.matrix(colSums(t(G.all), na.rm = T))
        SNPmatrix <- (SNPmatrix - (2 * FreqP))/sqrt(FreqPQ)
        G.ii.hat <- 1 + (G.ii)/NumberMarkers
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
        diag(Gmatrix) <- G.ii.hat
      }
      
      if (method == "Su"){
        TwoPQ <- 2*(FreqP)*(1-FreqP)
        SNPmatrix[SNPmatrix==2 | SNPmatrix==0] <- 0
        SNPmatrix <- SNPmatrix - TwoPQ
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/
          sum(TwoPQ[1,]*(1-TwoPQ[1,]))        
      }
      
      if (method == "Vitezica"){
        TwoPQ <- 2*(FreqP[1,])*(1-FreqP[1,])
        SNPmatrix[is.na(SNPmatrix)] <- 0
        SNPmatrix <- (SNPmatrix==0)*-2*(FreqP^2) +
          (SNPmatrix==1)*2*(FreqP)*(1-FreqP) +
          (SNPmatrix==2)*-2*((1-FreqP)^2)
        Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/sum(TwoPQ^2)
      }
      
      if (method == "Slater"){
        drop.alleles <- which(Frequency==0)
        if(length(drop.alleles)>0){
          Frequency <- Frequency[-drop.alleles]
          SNPmatrix <- SNPmatrix[,-drop.alleles]
          FreqP <- FreqP[,-drop.alleles]
        }
        FreqPQ <- matrix(rep(Frequency * (1-Frequency),
                             each = nrow(SNPmatrix)),
                         ncol = ncol(SNPmatrix))
        SNPmatrix[which(is.na(SNPmatrix))] <- 0
        G.ii <- (SNPmatrix^2 - (2 * FreqP) * SNPmatrix + FreqP^2)/FreqPQ
        G.ii <- as.matrix(colSums(t(G.ii), na.rm = T))
        G.ii <- 1 + (G.ii)/NumberMarkers
        SNPmatrix <- (SNPmatrix - (FreqP))/sqrt(FreqPQ)
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
        diag(Gmatrix) <- G.ii
      }
      
      if( method == "Endelman" ){
        if( ploidy != 4 ){
          cat( stop( "'Endelman' method is just implemented for ploidy=4" ))
        }
        Frequency <- colSums(SNPmatrix)/(nrow(SNPmatrix)*ploidy)
        Frequency <- cbind(Frequency,1-Frequency)
        SixPQ <- 6 * t((Frequency[, 1]^2)) %*% (Frequency[, 2]^2)
        SNPmatrix <- 6 * t((Frequency[, 1]^2)%*%t(rep(1,nrow(SNPmatrix)))) - 
          3*t((Frequency[, 1])%*%t(rep(1,nrow(SNPmatrix))))*SNPmatrix + 0.5 * SNPmatrix*(SNPmatrix-1)
        Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(SixPQ)
      }
      
      if (verify.posdef) {
        e.values <- eigen(Gmatrix, symmetric = TRUE)$values
        indicator <- length(which(e.values <= 0))
        if (indicator > 0) 
          cat("\t Matrix is NOT positive definite. It has ", indicator, 
              " eigenvalues <= 0 \n \n")
      }
      
      Time = as.matrix(proc.time() - Time)
      cat("Completed! Time =", Time[3], " seconds \n")
      gc()
      return(Gmatrix)
    }
    
    G_matrix <- Gmatrix()
    Gmat<- normalize_kinmat(as.matrix(G_matrix)); Gmat[is.na(Gmat)] <- 0
    write.table (Gmat, file=paste(pop,"_4x","_rd",rd+1,"_Kinship_Matrix.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    tiff(paste(pop,"_",ncol(pop_data),"markers_relatedness_heatmap_dendogram_4x.tiff",sep=""), width=30, height=30, units = 'in', res = 300, compression = 'lzw')
    heatmap(as.matrix(Gmat))
    dev.off()
    
    pca <- prcomp(Gmat)
    pca_out <- as.data.frame(pca$x)
    percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
    percentage <- paste( colnames(pca_out), "(", paste( as.character(percentage), "%", ")", sep="") )
    plot <- ggplot(pca_out,aes(x=PC1,y=PC2, label=row.names(pca_out))) +
      geom_point(color="gray") +
      geom_text(size=5, color="cornflowerblue") +
      xlab(percentage[1]) +
      ylab(percentage[2])
    ggsave(filename=paste(pop,"_",ncol(pop_data),"markers_2D_pca_4x.tiff",sep=""), plot=plot, width=15, height= 15, dpi=300, compression = "lzw")
    
    # Run this outside of the pipeline if required:
    # library(rgl)
    # x1=floor(min(pca_out$PC1, na.rm = TRUE))
    # x2=ceiling(max(pca_out$PC1, na.rm = TRUE))
    # y1=floor(min(pca_out$PC2, na.rm = TRUE))
    # y2=ceiling(max(pca_out$PC2, na.rm = TRUE))
    # z1=floor(min(pca_out$PC3, na.rm = TRUE))
    # z2=ceiling(max(pca_out$PC3, na.rm = TRUE))
    # plot3d(pca_out[,1:3], col=c(1:4), size=7.5, type='p',
    #        xlim = c(x1,x2), ylim=c(y1,y2), zlim=c(z1,z2))
    # text3d(pca_out[,1]+0, pca_out[,2]+0, pca_out[,3]+0,
    #        texts=c(rownames(pca_out)), cex= 1, pos=3)
    # # Wait! Adjust (moving image around with mouse) the 3D image before saving below.
    # rgl.snapshot(paste(pop,"_3D_pca_4x.png",sep=""), "png")
    
  } else {
    print ("Not enough markers to compute Gmatrix (i.e. threshold of 100 markers)")
  }
}
pop_struc()


