#!/usr/bin/env Rscript

# pop <- 
# p1 <- 
# p2 <- 
# gmissingness <- 
# smissingness <- 
# minRD <- 
# remove_id_list <- NULL
# remove_id_list <- c("NA_1","NA_2")
# remove_id_list <- paste(remove_id_list, "_GT", sep="")
# setwd ("dir")

####################################################################################################################
# set_variables
args <- commandArgs(trailingOnly = TRUE)
args
pop <- args[1]
p1 <- args[2]
p2 <- args[3]
gmissingness <- args[4]
smissingness <- args[5]
minRD <- args[6]
remove_id_list <- NULL
remove_id_list <- unlist(strsplit(args[7],","))
remove_id_list <- paste(remove_id_list, "_GT", sep="")
libdir <- args[8]
snpformats <- args[9]
gmissingness <- as.numeric(gmissingness)
smissingness <- as.numeric(smissingness)
minRD <- as.numeric(minRD)
rd <- minRD-1


.libPaths( c( .libPaths(), libdir) )
library(ggplot2)

####################################################################################################################
####################################################################################################################
subgenome_1 <- read.table (file=paste("../../snpcall/",pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
subgenome_1[,(((ncol(subgenome_1)-4)/2)+5):ncol(subgenome_1)] <- lapply(subgenome_1[,(((ncol(subgenome_1)-4)/2)+5):ncol(subgenome_1)], gsub, pattern = "|", replacement = "/", fixed = TRUE)
RD_snpfiltering <- function(){
  #remove samples theat you want to exclude from the analysis
  if (length(remove_id_list) > 0) {
    remove_id_GT <- remove_id_list
    remove_id_DP <- gsub("_GT", "_DP", remove_id_list)
    remove_id <- c(remove_id_DP, remove_id_GT)
    id <- names(subgenome_1)
    keep_id <- setdiff(id,remove_id)
    subgenome_1 <- subgenome_1[,c(keep_id)]
  }
  write.table (subgenome_1, file=paste(pop,"_4x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  
  #Filter based on DP >= rd in Parents
  subgenome_1_filtered <- subgenome_1
  subgenome_1_filtered <- subset(subgenome_1_filtered, subgenome_1_filtered[,paste(p1,"_DP",sep="")] >= 12) 
  subgenome_1_filtered <- subset(subgenome_1_filtered, subgenome_1_filtered[,paste(p2,"_DP",sep="")] >= 12)
  subgenome_1_filtered_1 <- subgenome_1_filtered
  subgenome_1_filtered_1$no_missing <- rowSums(is.na(subgenome_1_filtered_1))
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, no_missing <= ((ncol(subgenome_1_filtered_1)-5)/2)*gmissingness)
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, select=-c(no_missing))

  #remove monomorphic and non-segregating markers
  subgenome_1_filtered_1$noseg <- paste(subgenome_1_filtered_1[,paste(p1,"_GT",sep="")],"_",subgenome_1_filtered_1[,paste(p2,"_GT",sep="")], sep="")
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, noseg != "0/0/0/0_0/0/0/0")
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, noseg != "0/0/0/0_1/1/1/1")
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, noseg != "1/1/1/1_0/0/0/0")
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, noseg != "1/1/1/1_1/1/1/1")
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, select=-c(noseg))
  subgenome_1_filtered_0 <- subgenome_1_filtered_1
    
  #change genotype to missing if read depth is < rd in complete dataset (including progenies)
  subgenome_1_filtered_1_AB <- subset(subgenome_1_filtered_1, select=c(1:(((ncol(subgenome_1_filtered_1)-4)/2)+4)))
  for (i in 5:(((ncol(subgenome_1_filtered_1)-4)/2)+4)) {
    j <- i+((ncol(subgenome_1_filtered)-4)/2)
    subgenome_filtered[,j][subgenome_filtered[,i] < 12] <- NA
    subgenome_1_filtered_1[,i][subgenome_1_filtered_1[,j] == "0/0/0/0/0/0" ] <- minRD
    subgenome_1_filtered_1[,i][subgenome_1_filtered_1[,j] == "1/1/1/1/1/1" ] <- minRD
    subgenome_1_filtered_1[,j][subgenome_1_filtered_1[,i] <= rd] <- NA
    gc()
  }
  subgenome_1_filtered_1_C <- subset(subgenome_1_filtered_1, select=c((((ncol(subgenome_1_filtered_1)-4)/2)+5):ncol(subgenome_1_filtered_1)))
  subgenome_1_filtered_1 <- cbind(subgenome_1_filtered_1_AB, subgenome_1_filtered_1_C)
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, !is.na(subgenome_1_filtered_1[,paste(p1,"_DP",sep="")])) 
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, !is.na(subgenome_1_filtered_1[,paste(p2,"_DP",sep="")]))
  subgenome_1_filtered_1$no_missing <- rowSums(is.na(subgenome_1_filtered_1))
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, no_missing <= ((ncol(subgenome_1_filtered_1)-5)/2)*gmissingness)
  subgenome_1_filtered_1 <- subset(subgenome_1_filtered_1, select=-c(no_missing))

    
  subgenome_1_SD_1 <- subgenome_1_filtered_1
  subgenome_1_SD_1 <- subset(subgenome_1_SD_1, select=c(1:4,(((ncol(subgenome_1_SD_1)-4)/2)+5):ncol(subgenome_1_SD_1)))
  # check for non-biallelic genotypes and change to missing (NA)
  nonbiallelic <- as.vector(as.matrix(subgenome_1_SD_1[,5:ncol(subgenome_1_SD_1)]))
  nonbiallelic <- unique(nonbiallelic)
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/0")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/1/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "1/1/1/1")
  if (length(nonbiallelic > 0)){
    for (i in 1:length(nonbiallelic)){
      subgenome_1_SD_1[][subgenome_1_SD_1[]==nonbiallelic[i]] <- NA
    }
  }
  subgenome_1_SD_1$no_missing <- rowSums(is.na(subgenome_1_SD_1))
  subgenome_1_SD_1 <- subset(subgenome_1_SD_1, no_missing <= (ncol(subgenome_1_SD_1)-5)*gmissingness)
  subgenome_1_SD_1 <- subset(subgenome_1_SD_1, select=-c(no_missing))
  col_idx <- grep(paste(p1,"_GT",sep=""), names(subgenome_1_SD_1))
  subgenome_1_SD_1 <- subgenome_1_SD_1[, c((1:ncol(subgenome_1_SD_1))[-col_idx],col_idx)]
  col_idx <- grep(paste(p2,"_GT",sep=""), names(subgenome_1_SD_1))
  subgenome_1_SD_1 <- subgenome_1_SD_1[, c((1:ncol(subgenome_1_SD_1))[-col_idx],col_idx)]
  write.table (subgenome_1_SD_1, file=paste(pop,"_4x","_SD_rd",rd+1,".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

  #change genotype to missing if read depth is < 12 in complete dataset (including progenies)
  for (i in 5:(((ncol(subgenome_1_filtered_0)-4)/2)+4)) {
    j <- i+((ncol(subgenome_1_filtered_0)-4)/2)
    subgenome_1_filtered_0[,j][subgenome_1_filtered_0[,i] < 12] <- NA
    gc
  }
  subgenome_1_filtered_0$no_missing <- rowSums(is.na(subgenome_1_filtered_0))
  subgenome_1_filtered_0 <- subset(subgenome_1_filtered_0, no_missing <= ((ncol(subgenome_1_filtered_0)-5)/2)*gmissingness)
  subgenome_1_filtered_0 <- subset(subgenome_1_filtered_0, select=-c(no_missing))
  
  subgenome_1_SD_0 <- subgenome_1_filtered_0
  subgenome_1_SD_0 <- subset(subgenome_1_SD_0, select=c(1:4,(((ncol(subgenome_1_SD_0)-4)/2)+5):ncol(subgenome_1_SD_0)))
  # check for non-biallelic genotypes and change to missing (NA)
  nonbiallelic <- as.vector(as.matrix(subgenome_1_SD_0[,5:ncol(subgenome_1_SD_0)]))
  nonbiallelic <- unique(nonbiallelic)
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/0")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/0/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/0/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "0/1/1/1")
  nonbiallelic <- subset(nonbiallelic, nonbiallelic != "1/1/1/1")
  if (length(nonbiallelic > 0)){
    for (i in 1:length(nonbiallelic)){
      subgenome_1_SD_0[][subgenome_1_SD_0[]==nonbiallelic[i]] <- NA
    }
  }
  subgenome_1_SD_0$no_missing <- rowSums(is.na(subgenome_1_SD_0))
  subgenome_1_SD_0 <- subset(subgenome_1_SD_0, no_missing <= (ncol(subgenome_1_SD_0)-5)*gmissingness)
  subgenome_1_SD_0 <- subset(subgenome_1_SD_0, select=-c(no_missing))
  col_idx <- grep(paste(p1,"_GT",sep=""), names(subgenome_1_SD_0))
  subgenome_1_SD_0 <- subgenome_1_SD_0[, c((1:ncol(subgenome_1_SD_0))[-col_idx],col_idx)]
  col_idx <- grep(paste(p2,"_GT",sep=""), names(subgenome_1_SD_0))
  subgenome_1_SD_0 <- subgenome_1_SD_0[, c((1:ncol(subgenome_1_SD_0))[-col_idx],col_idx)]
  write.table (subgenome_1_SD_0, file=paste(pop,"_4x","_SD_rd12",".txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }
RD_snpfiltering()

subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
SD_snpfiltering <- function() {
  #######################################################################################################################################################################################
  #Filter for segregation distortion at unmethlyated loci
  #############################################################
  # Every thing at Read Depth of 20 except for configurations with duplex
  #	simplex x other configuration
  # 0/0/0/1 x 0/0/0/0   ->   0/0/0/0   0/0/0/1                                   -> 1:1
  # 0/0/0/1 x 0/0/0/1   ->   0/0/0/0   0/0/0/1   0/0/1/1                         -> 1:2:1
  # 0/0/0/1 x 0/0/1/1   ->   0/0/0/0   0/0/0/1   0/0/1/1   0/1/1/1               -> 1:5:5:1
  # 0/0/0/1 x 0/1/1/1   ->   0/0/0/1   0/0/1/1   0/1/1/1                         -> 1:2:1
  # 0/0/0/1 x 1/1/1/1   ->   0/0/1/1   0/1/1/1                                   -> 1:1
  #
  # Duplex and other configurations
  # 0/0/1/1 x 0/0/0/0   ->    0/0/0/0   0/0/0/1   0/0/1/1                         -> 1:4:1
  # 0/0/1/1 x 0/0/1/1   ->    0/0/0/0   0/0/0/1   0/0/1/1   0/1/1/1   1/1/1/1     -> 1:8:18:8:1
  # 0/0/1/1 x 0/1/1/1   ->    0/0/0/1   0/0/1/1   0/1/1/1   1/1/1/1               -> 1:5:5:1
  # 0/0/1/1 x 1/1/1/1   ->    0/0/1/1   0/1/1/1   1/1/1/1                         -> 1:4:1
  #
  # Triplex x other configurations
  # 0/1/1/1 x 0/0/0/0   ->    0/0/0/1   0/0/1/1                                   -> 1:1
  # 0/1/1/1 x 0/1/1/1   ->    0/0/1/1   0/1/1/1   1/1/1/1                         -> 1:2:1
  # 0/1/1/1 x 1/1/1/1   ->    0/1/1/1   1/1/1/1                                   -> 1:1
  #
  #############################################################
  #############################################################
  # Simplex x Nulliplex
  #	0/0/0/0 x 0/0/0/1   ->    0/0/0/0   0/0/0/1       ->	1:1
  #	0/0/0/1 x 0/0/0/0   ->    0/0/0/0   0/0/0/1       ->  1:1
  ## File Code: "subgenome_1_SD_1_G0001G0000"
  
  subgenome_1_SD_0 <- read.table (file=paste(pop,"_4x","_SD_rd12",".txt",sep=""), header=T, sep="\t", check.names = FALSE)
  subgenome_1_SD_1 <- read.table (file=paste(pop,"_4x","_SD_rd",rd+1,".txt",sep=""), header=T, sep="\t", check.names = FALSE)

  subgenome_1_SD_1a <- subgenome_1_SD_0
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/1/1"] <- "0/0/0/1"
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/1/1/1"] <- "0/0/0/1"
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/0/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/0")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/0/0" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/1"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/0"))
  }
  if ( nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0000"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0" #
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0001=as.numeric(GT0001), GT0000=as.numeric(GT0000))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:4],p=c(0.5,0.5))
      subgenome_1_SD_1a_count[i,5] <- X$p.value
    }
    subgenome_1_SD_1_G0001G0000 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0001G0000 <- subset(subgenome_1_SD_1_G0001G0000, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6)))
    write.table (subgenome_1_SD_1_G0001G0000, file=paste(pop,"_4x","_SD_1_G0001G0000_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0001G0000 <- subset(subgenome_1_SD_1_G0001G0000, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0001G0000, file=paste(pop,"_4x","_SD_1_G0001G0000.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0001G0000 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Simplex x Simplex
  #	0/0/0/1 x 0/0/0/1   ->    0/0/0/0   0/0/0/1   0/0/1/1   ->  1:2:1
  #File Code: "subgenome_1_SD_1_G0001G0001"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/1/1/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/0/0/1")
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/0/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="0/0/0/1")
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT=="0/0/0/0")
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0000"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0011=as.numeric(GT0011), GT0001=as.numeric(GT0001), GT0000=as.numeric(GT0000))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:5],p=c(0.25,0.5,0.25))
      subgenome_1_SD_1a_count[i,6] <- X$p.value
    }
    subgenome_1_SD_1_G0001G0001 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0001G0001 <- subset(subgenome_1_SD_1_G0001G0001, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7)))
    write.table (subgenome_1_SD_1_G0001G0001, file=paste(pop,"_4x","_SD_1_G0001G0001_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0001G0001 <- subset(subgenome_1_SD_1_G0001G0001, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0001G0001, file=paste(pop,"_4x","_SD_1_G0001G0001.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0001G0001 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Simplex x Duplex
  #	0/0/0/1 x 0/0/1/1   ->    0/0/0/0   0/0/0/1   0/0/1/1   0/1/1/1   ->  1:5:5:1
  #	0/0/1/1 x 0/0/0/1   ->    0/0/0/0   0/0/0/1   0/0/1/1   0/1/1/1   ->  1:5:5:1
  #File Code: "subgenome_1_SD_1_G0001G0011"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/0/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/0/0/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1aw <- subset(subgenome_1_SD_1aa, GT=="0/1/1/1")
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/0/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="0/0/0/1")
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT=="0/0/0/0")
  }
  if (nrow(subgenome_1_SD_1aw)>0 && nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1aw_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1aw, FUN=length)
    names(subgenome_1_SD_1aw_count)[names(subgenome_1_SD_1aw_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0000"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1aw_count, subgenome_1_SD_1ax_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0111=as.numeric(GT0111), GT0011=as.numeric(GT0011), GT0001=as.numeric(GT0001), GT0000=as.numeric(GT0000))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:6],p=c(0.0833333333333333,0.4166666666666670,
                                                         0.4166666666666670,0.0833333333333333), simulate.p.value = TRUE)
      subgenome_1_SD_1a_count[i,7] <- X$p.value
    }
    subgenome_1_SD_1_G0001G0011 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0001G0011 <- subset(subgenome_1_SD_1_G0001G0011, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7),-(((ncol(subgenome_1)-4)/2)+8)))
    write.table (subgenome_1_SD_1_G0001G0011, file=paste(pop,"_4x","_SD_1_G0001G0011_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0001G0011 <- subset(subgenome_1_SD_1_G0001G0011, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0001G0011, file=paste(pop,"_4x","_SD_1_G0001G0011.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0001G0011 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Simplex x Triplex
  #	0/0/0/1 x 0/1/1/1   ->    0/0/0/1   0/0/1/1   0/1/1/1   ->  1:2:1
  #	0/1/1/1 x 0/0/0/1   ->    0/0/0/1   0/0/1/1   0/1/1/1   ->  1:2:1
  #File Code: "subgenome_1_SD_1_G0001G0111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/1/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/0/0/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/0/0/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="0/0/1/1")
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT=="0/1/1/1")
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0001=as.numeric(GT0001), GT0011=as.numeric(GT0011), GT0111=as.numeric(GT0111))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:5],p=c(0.25,0.5,0.25))
      subgenome_1_SD_1a_count[i,6] <- X$p.value
    }
    subgenome_1_SD_1_G0001G0111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0001G0111 <- subset(subgenome_1_SD_1_G0001G0111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7)))
    write.table (subgenome_1_SD_1_G0001G0111, file=paste(pop,"_4x","_SD_1_G0001G0111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0001G0111 <- subset(subgenome_1_SD_1_G0001G0111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0001G0111, file=paste(pop,"_4x","_SD_1_G0001G0111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0001G0111 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Simplex x Quadruplex
  #	0/0/0/1 x 1/1/1/1   ->    0/0/1/1   0/1/1/1       ->	1:1
  #	1/1/1/1 x 0/0/0/1   ->    0/0/1/1   0/1/1/1       ->  1:1
  # File Code: "subgenome_1_SD_1_G0001G1111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/0/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "1/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="1/1/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("0/0/1/1"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/1/1/1"))
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0011=as.numeric(GT0011), GT0111=as.numeric(GT0111))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:4],p=c(0.5,0.5))
      subgenome_1_SD_1a_count[i,5] <- X$p.value
    }
    subgenome_1_SD_1_G0001G1111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0001G1111 <- subset(subgenome_1_SD_1_G0001G1111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6)))
    write.table (subgenome_1_SD_1_G0001G1111, file=paste(pop,"_4x","_SD_1_G0001G1111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0001G1111 <- subset(subgenome_1_SD_1_G0001G1111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0001G1111, file=paste(pop,"_4x","_SD_1_G0001G1111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0001G1111 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Duplex x Nulliplex
  #	0/0/1/1 x 0/0/0/0   ->    0/0/0/0   0/0/0/1   0/0/1/1       ->	1:4:1
  #	0/0/0/0 x 0/0/1/1   ->    0/0/0/0   0/0/0/1   0/0/1/1       ->	1:4:1
  # File Code: "subgenome_1_SD_1_G0011G0000"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/1/1/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/0")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/0/0" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/1/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/0"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/1"))
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT==c("0/0/1/1"))
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0000"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0000=as.numeric(GT0000), GT0001=as.numeric(GT0001), GT0011=as.numeric(GT0011))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:5],p=c(0.166666666666667,0.666666666666667,0.166666666666667))
      subgenome_1_SD_1a_count[i,6] <- X$p.value
    }
    subgenome_1_SD_1_G0011G0000 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0011G0000 <- subset(subgenome_1_SD_1_G0011G0000, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7)))
    write.table (subgenome_1_SD_1_G0011G0000, file=paste(pop,"_4x","_SD_1_G0011G0000_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0011G0000 <- subset(subgenome_1_SD_1_G0011G0000, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0011G0000, file=paste(pop,"_4x","_SD_1_G0011G0000.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0011G0000 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Duplex x Duplex
  #	0/0/1/1 x 0/0/1/1   ->    0/0/0/0   0/0/0/1   0/0/1/1   0/1/1/1   1/1/1/1     -> 1:8:18:8:1
  # File Code: "subgenome_1_SD_1_G0011G0011"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/1/1")
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1av <- subset(subgenome_1_SD_1aa, GT==c("1/1/1/1"))
    subgenome_1_SD_1aw <- subset(subgenome_1_SD_1aa, GT==c("0/1/1/1"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("0/0/1/1"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/1"))
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/0"))
  }
  if (nrow(subgenome_1_SD_1av)>0 && nrow(subgenome_1_SD_1aw)>0 && nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1av_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1av, FUN=length)
    names(subgenome_1_SD_1av_count)[names(subgenome_1_SD_1av_count) == "GT"] <- "GT1111"
    subgenome_1_SD_1aw_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1aw, FUN=length)
    names(subgenome_1_SD_1aw_count)[names(subgenome_1_SD_1aw_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0000"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1av_count, subgenome_1_SD_1aw_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1ax_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT1111=as.numeric(GT1111), GT0111=as.numeric(GT0111), GT0011=as.numeric(GT0011), GT0001=as.numeric(GT0001), GT0000=as.numeric(GT0000))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:7],p=c(0.0277777777777778,0.2222222222222220,0.5000000000000000,
                                                         0.2222222222222220,0.0277777777777778), simulate.p.value = TRUE)
      subgenome_1_SD_1a_count[i,8] <- X$p.value
    }
    subgenome_1_SD_1_G0011G0011 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0011G0011 <- subset(subgenome_1_SD_1_G0011G0011, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7),-(((ncol(subgenome_1)-4)/2)+8),-(((ncol(subgenome_1)-4)/2)+9)))
    write.table (subgenome_1_SD_1_G0011G0011, file=paste(pop,"_4x","_SD_1_G0011G0011_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0011G0011 <- subset(subgenome_1_SD_1_G0011G0011, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0011G0011, file=paste(pop,"_4x","_SD_1_G0011G0011.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0011G0011 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Duplex x Triplex
  #	0/0/1/1 x 0/1/1/1   ->    0/0/0/1   0/0/1/1   0/1/1/1   1/1/1/1       ->	1:5:5:1
  #	0/1/1/1 x 0/0/1/1   ->    0/0/0/1   0/0/1/1   0/1/1/1   1/1/1/1       ->	1:5:5:1
  # File Code: "subgenome_1_SD_1_G0011G0111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/1/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/1/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1aw <- subset(subgenome_1_SD_1aa, GT==c("0/0/0/1"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("0/0/1/1"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/1/1/1"))
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT==c("1/1/1/1"))
  }
  if (nrow(subgenome_1_SD_1aw)>0 && nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1aw_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1aw, FUN=length)
    names(subgenome_1_SD_1aw_count)[names(subgenome_1_SD_1aw_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT1111"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1aw_count, subgenome_1_SD_1ax_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0001=as.numeric(GT0001), GT0011=as.numeric(GT0011), GT0111=as.numeric(GT0111), GT1111=as.numeric(GT1111))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:6],p=c(0.0833333333333333,0.4166666666666670,
                                                         0.4166666666666670,0.0833333333333333), simulate.p.value = TRUE)
      subgenome_1_SD_1a_count[i,7] <- X$p.value
    }
    subgenome_1_SD_1_G0011G0111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0011G0111 <- subset(subgenome_1_SD_1_G0011G0111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7),-(((ncol(subgenome_1)-4)/2)+8)))
    write.table (subgenome_1_SD_1_G0011G0111, file=paste(pop,"_4x","_SD_1_G0011G0111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0011G0111 <- subset(subgenome_1_SD_1_G0011G0111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0011G0111, file=paste(pop,"_4x","_SD_1_G0011G0111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0011G0111 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Duplex x  Quadruplex
  #	0/0/1/1 x 1/1/1/1   ->    0/0/1/1   0/1/1/1   1/1/1/1       ->	1:4:1
  #	1/1/1/1 x 0/0/1/1   ->    0/0/1/1   0/1/1/1   1/1/1/1       ->	1:4:1
  # File Code: "subgenome_1_SD_1_G0011G1111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="0/0/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "1/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  =="1/1/1/1" & subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/1/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT==c("1/1/1/1"))
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT==c("0/1/1/1"))
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT==c("0/0/1/1"))
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT1111"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT1111=as.numeric(GT1111), GT0111=as.numeric(GT0111), GT0011=as.numeric(GT0011))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:5],p=c(0.166666666666667,0.666666666666667,0.166666666666667))
      subgenome_1_SD_1a_count[i,6] <- X$p.value
    }
    subgenome_1_SD_1_G0011G1111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0011G1111 <- subset(subgenome_1_SD_1_G0011G1111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7)))
    write.table (subgenome_1_SD_1_G0011G1111, file=paste(pop,"_4x","_SD_1_G0011G1111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0011G1111 <- subset(subgenome_1_SD_1_G0011G1111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0011G1111, file=paste(pop,"_4x","_SD_1_G0011G1111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0011G1111 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Triplex x Nulliplex
  #	0/1/1/1 x 0/0/0/0   ->    0/0/0/1   0/0/1/1   ->  1:1
  #	0/0/0/0 x 0/1/1/1   ->    0/0/0/1   0/0/1/1   ->  1:1
  #File Code: "subgenome_1_SD_1_G0111G0000"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="1/1/1/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/1/1/1"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/1/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/0/0/0")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/0/0/0" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/1/1/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/0/0/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="0/0/1/1")
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0001"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0001=as.numeric(GT0001), GT0011=as.numeric(GT0011))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:4],p=c(0.5,0.5))
      subgenome_1_SD_1a_count[i,5] <- X$p.value
    }
    subgenome_1_SD_1_G0111G0000 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0111G0000 <- subset(subgenome_1_SD_1_G0111G0000, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6)))
    write.table (subgenome_1_SD_1_G0111G0000, file=paste(pop,"_4x","_SD_1_G0111G0000_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0111G0000 <- subset(subgenome_1_SD_1_G0111G0000, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0111G0000, file=paste(pop,"_4x","_SD_1_G0111G0000.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0111G0000 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Triplex x Triplex
  #	0/1/1/1 x 0/1/1/1   ->    0/0/1/1   0/1/1/1   1/1/1/1   ->  1:2:1
  #File Code: "subgenome_1_SD_1_G0111G0111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_1
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/1"] <- NA
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/1/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/1/1/1")
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/0/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="0/1/1/1")
    subgenome_1_SD_1az <- subset(subgenome_1_SD_1aa, GT=="1/1/1/1")
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0 && nrow(subgenome_1_SD_1az)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0011"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1az_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1az, FUN=length)
    names(subgenome_1_SD_1az_count)[names(subgenome_1_SD_1az_count) == "GT"] <- "GT1111"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1a_count, subgenome_1_SD_1az_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0011=as.numeric(GT0011), GT0111=as.numeric(GT0111), GT1111=as.numeric(GT1111))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:5],p=c(0.25,0.5,0.25))
      subgenome_1_SD_1a_count[i,6] <- X$p.value
    }
    subgenome_1_SD_1_G0111G0111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0111G0111 <- subset(subgenome_1_SD_1_G0111G0111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6),-(((ncol(subgenome_1)-4)/2)+7)))
    write.table (subgenome_1_SD_1_G0111G0111, file=paste(pop,"_4x","_SD_1_G0111G0111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0111G0111 <- subset(subgenome_1_SD_1_G0111G0111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0111G0111, file=paste(pop,"_4x","_SD_1_G0111G0111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0111G0111 configuration failed segregation distortion test")
  }
  #############################################################
  #############################################################
  # Triplex x Quadruplex
  #	0/1/1/1 x 1/1/1/1   ->    0/1/1/1   1/1/1/1   ->  1:1
  #	1/1/1/1 x 0/1/1/1   ->    0/1/1/1   1/1/1/1   ->  1:1
  #File Code: "subgenome_1_SD_1_G0111G1111"
  
  subgenome_1_SD_1a <- subgenome_1_SD_0
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/0"] <- NA
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/0/1"] <- "0/1/1/1"
  subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2][subgenome_1_SD_1a[, 5:ncol(subgenome_1_SD_1a)-2]=="0/0/1/1"] <- "0/1/1/1" 
  subgenome_1_SD_1a$no_missing <- rowSums(is.na(subgenome_1_SD_1a))
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, no_missing <= ((ncol(subgenome_1)-5)/2)*gmissingness)
  subgenome_1_SD_1a <- subset(subgenome_1_SD_1a, select=c(-(((ncol(subgenome_1)-4)/2)+5)))
  if (nrow(subgenome_1_SD_1a) > 0) {
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "0/1/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "1/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1a, subgenome_1_SD_1a[paste(p2,"_GT",sep="")]  == "1/1/1/1" & subgenome_1_SD_1a[paste(p1,"_GT",sep="")]  == "0/1/1/1")
    subgenome_1_SD_1a <- rbind(subgenome_1_SD_1ax, subgenome_1_SD_1ay)
  }
  subgenome_1_SD_1aa <- data.frame()
  subgenome_1_SD_1av <- data.frame()
  subgenome_1_SD_1aw <- data.frame()
  subgenome_1_SD_1ax <- data.frame()
  subgenome_1_SD_1ay <- data.frame()
  subgenome_1_SD_1az <- data.frame()
  if (nrow(subgenome_1_SD_1a) >= 1) {
    subgenome_1_SD_1aa <- reshape(subgenome_1_SD_1a, direction ="long",
                                  varying=list(c(5:(((ncol(subgenome_1)-4)/2)+2))),
                                  v.names=c("GT"))
    subgenome_1_SD_1ax <- subset(subgenome_1_SD_1aa, GT=="0/1/1/1")
    subgenome_1_SD_1ay <- subset(subgenome_1_SD_1aa, GT=="1/1/1/1")
  }
  if (nrow(subgenome_1_SD_1ax)>0 && nrow(subgenome_1_SD_1ay)>0) {
    subgenome_1_SD_1ax_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ax, FUN=length)
    names(subgenome_1_SD_1ax_count)[names(subgenome_1_SD_1ax_count) == "GT"] <- "GT0111"
    subgenome_1_SD_1ay_count <- aggregate(GT ~ CHROM + POS, subgenome_1_SD_1ay, FUN=length)
    names(subgenome_1_SD_1ay_count)[names(subgenome_1_SD_1ay_count) == "GT"] <- "GT1111"
    subgenome_1_SD_1a_count <- merge(subgenome_1_SD_1ax_count, subgenome_1_SD_1ay_count, all=TRUE)
    subgenome_1_SD_1a_count[][is.na(subgenome_1_SD_1a_count[])] <- "0"
    subgenome_1_SD_1a_count <- transform(subgenome_1_SD_1a_count, GT0111=as.numeric(GT0111), GT1111=as.numeric(GT1111))
    subgenome_1_SD_1a_count$pvalue <- NA
    for (i in 1:nrow(subgenome_1_SD_1a_count)) {
      X <- chisq.test(subgenome_1_SD_1a_count[i,3:4],p=c(0.5,0.5))
      subgenome_1_SD_1a_count[i,5] <- X$p.value
    }
    subgenome_1_SD_1_G0111G1111 <- merge(subgenome_1_SD_1a, subgenome_1_SD_1a_count, by=c("CHROM", "POS"), all=TRUE)
    subgenome_1_SD_1_G0111G1111 <- subset(subgenome_1_SD_1_G0111G1111, select=c(-(((ncol(subgenome_1)-4)/2)+5),-(((ncol(subgenome_1)-4)/2)+6)))
    write.table (subgenome_1_SD_1_G0111G1111, file=paste(pop,"_4x","_SD_1_G0111G1111_plusSD.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    subgenome_1_SD_1_G0111G1111 <- subset(subgenome_1_SD_1_G0111G1111, pvalue >=0.001)
    write.table (subgenome_1_SD_1_G0111G1111, file=paste(pop,"_4x","_SD_1_G0111G1111.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  } else {
    print ("all SNPs in G0111G1111 configuration failed segregation distortion test")
  }
  
}
SD_snpfiltering()

# set functions for RD and SD filtering..... remove samples with more than threshold (e.g. 20%) missing data 
sample_subsample <- function(){
  subgenome_1_noSD <- NULL
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0000 <- read.table(paste(pop,"_4x","_SD_1_G0001G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0001G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0001.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0001 <- read.table(paste(pop,"_4x","_SD_1_G0001G0001.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0001G0001) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0011.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0011 <- read.table(paste(pop,"_4x","_SD_1_G0001G0011.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0001G0011) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0111 <- read.table(paste(pop,"_4x","_SD_1_G0001G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0001G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0001G1111 <- read.table(paste(pop,"_4x","_SD_1_G0001G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0001G1111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0000 <- read.table(paste(pop,"_4x","_SD_1_G0011G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0011G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0011.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0011 <- read.table(paste(pop,"_4x","_SD_1_G0011G0011.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0011G0011) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0111 <- read.table(paste(pop,"_4x","_SD_1_G0011G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0011G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0011G1111 <- read.table(paste(pop,"_4x","_SD_1_G0011G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0011G1111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0111G0000 <- read.table(paste(pop,"_4x","_SD_1_G0111G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0111G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0111G0111 <- read.table(paste(pop,"_4x","_SD_1_G0111G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0111G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0111G1111 <- read.table(paste(pop,"_4x","_SD_1_G0111G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_noSD <- rbind(subgenome_1_noSD,subgenome_1_SD_1_G0111G1111) }
  remove_id_list <- subset(subgenome_1_noSD, select=-c(1:4,pvalue))
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
    remove_id_DP <- gsub("_GT", "_DP", remove_id_list)
    remove_id <- c(remove_id_DP, remove_id_GT)
    id <- colnames(subgenome_1[,1:ncol(subgenome_1)])
    keep_id <- setdiff(id,remove_id)
    subgenome_1 <-subset(subgenome_1, select=c(keep_id))
    write.table (subgenome_1, file=paste(pop,"_4x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }
}
sample_subsample()

subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
RD_snpfiltering()

subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
SD_snpfiltering()

####################################################################################################################
final_summary <- function(){
  ############################################################################################################################################################################################################
  # Plot all markers (regular) with segregation distortion
  parents1 <- read.table(paste(pop,"_4x","_SD_rd",rd+1,".txt",sep=""), header=T, sep="\t", check.names = FALSE)
  parents1 <- subset(parents1, select=c(paste(p1,"_GT",sep=""),paste(p2,"_GT",sep="")))
  parents1 <- na.omit(parents1)
  parents1$cross <- paste (parents1[,paste(p1,"_GT",sep="")], parents1[,paste(p2,"_GT",sep="")], sep=" x ")
  SNP4x <- parents1
  SNP4x <- subset(SNP4x, select="cross")
  SNP4x <- as.data.frame(table(SNP4x$cross))
  sum <- sum(SNP4x$Freq, na.rm=T)
  max <- max(SNP4x$Freq, na.rm=T)
  
  max <- max*1.2
  SNP4x$percentage <- ((SNP4x$Freq)/sum)*100
  SNP4x[,3] <- round(SNP4x[,3], 2)
  plot_SNP <- ggplot(SNP4x,aes(x=Var1, y=Freq, color=Var1)) +
    geom_col(fill="white") +
    scale_fill_hue() +
    coord_flip()+
    geom_text(aes(label=paste(Freq," = ",percentage,"%")),hjust=-0.2, size=3, color="black", fontface="italic") +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=10),
          axis.text.y=element_text(colour="cornflowerblue", size=10),
          axis.title=element_text(size=12)) +
    theme(legend.position="none")+
    theme(legend.text=element_text(colour="cornflowerblue", size=10)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    theme(legend.title=element_blank()) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
    ylab(paste("Number of Variants (",sum," loci)"))
  ggsave(filename="plot_SNP4x_wSD.tiff", plot=plot_SNP, width=10, height= 5, dpi=600, compression = "lzw")
  
  
  ########################################################################################################################################
  # Plot all markers (regular, restriction polymorphism-based):
  # noSD (no segregation distorted markers)
  ########################################################################################################################################
  subgenome_1_SD_1 <- NULL
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0000 <- read.table(paste(pop,"_4x","_SD_1_G0001G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0001G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0001.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0001 <- read.table(paste(pop,"_4x","_SD_1_G0001G0001.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0001G0001) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0011.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0011 <- read.table(paste(pop,"_4x","_SD_1_G0001G0011.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0001G0011) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0001G0111 <- read.table(paste(pop,"_4x","_SD_1_G0001G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0001G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0001G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0001G1111 <- read.table(paste(pop,"_4x","_SD_1_G0001G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0001G1111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0000 <- read.table(paste(pop,"_4x","_SD_1_G0011G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0011G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0011.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0011 <- read.table(paste(pop,"_4x","_SD_1_G0011G0011.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0011G0011) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0011G0111 <- read.table(paste(pop,"_4x","_SD_1_G0011G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0011G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0011G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0011G1111 <- read.table(paste(pop,"_4x","_SD_1_G0011G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0011G1111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G0000.txt",sep=""))) {
    subgenome_1_SD_1_G0111G0000 <- read.table(paste(pop,"_4x","_SD_1_G0111G0000.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0111G0000) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G0111.txt",sep=""))) {
    subgenome_1_SD_1_G0111G0111 <- read.table(paste(pop,"_4x","_SD_1_G0111G0111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0111G0111) }
  if (file.exists(paste(pop,"_4x","_SD_1_G0111G1111.txt",sep=""))) {
    subgenome_1_SD_1_G0111G1111 <- read.table(paste(pop,"_4x","_SD_1_G0111G1111.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE,colClasses=c("character"), check.names = FALSE)
    subgenome_1_SD_1 <- rbind(subgenome_1_SD_1,subgenome_1_SD_1_G0111G1111) }
  subgenome_1_noSD <- subgenome_1_SD_1
  
  subgenome_1_SD_1 <- subset(subgenome_1_SD_1, select=c(paste(p1,"_GT",sep=""),paste(p2,"_GT",sep="")))
  subgenome_1_SD_1$cross <- paste (subgenome_1_SD_1[,paste(p1,"_GT",sep="")], subgenome_1_SD_1[,paste(p2,"_GT",sep="")], sep=" x ")
  subgenome_1_SD_1 <- na.omit(subgenome_1_SD_1)
  SNP4x <- subgenome_1_SD_1
  SNP4x <- subset(SNP4x, select="cross")
  SNP4x <- as.data.frame(table(SNP4x$cross))
  sum <- sum(SNP4x$Freq, na.rm=T)
  max <- max(SNP4x$Freq, na.rm=T)
  max <- max*1.2
  SNP4x$percentage <- ((SNP4x$Freq)/sum)*100
  SNP4x[,3] <- round(SNP4x[,3], 2)
  
  plot_SNP <- ggplot(SNP4x,aes(x=Var1, y=Freq, color=Var1)) +
    geom_col(fill="white") +
    coord_flip()+
    geom_text(aes(label=paste(Freq," = ",percentage,"%")),hjust=-0.2, size=3, color="black", fontface="italic") +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=10),
          axis.text.y=element_text(colour="cornflowerblue", size=10),
          axis.title=element_text(size=12)) +
    theme(legend.position="none")+
    theme(legend.text=element_text(colour="cornflowerblue", size=10)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    theme(legend.title=element_blank()) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
    ylab(paste("Number of Variants (",sum," loci)"))
  ggsave(filename="plot_SNP4x_noSD.tiff", plot=plot_SNP, width=10, height= 6, dpi=600, compression = "lzw")
  
  
  ########################################################################################################################################
  #Generate final filtered genotype file for "Regular Markers"
  names(subgenome_1_noSD) <- gsub("_GT", "", names(subgenome_1_noSD))
  subgenome_1_noSD$SNP <- paste (subgenome_1_noSD$CHROM,"_",subgenome_1_noSD$POS, sep="")
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)==paste(p2,sep="")),which(colnames(subgenome_1_noSD)!=paste(p2,sep="")))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)==paste(p1,sep="")),which(colnames(subgenome_1_noSD)!=paste(p1,sep="")))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)=="ALT"),which(colnames(subgenome_1_noSD)!="ALT"))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)=="REF"),which(colnames(subgenome_1_noSD)!="REF"))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)=="POS"),which(colnames(subgenome_1_noSD)!="POS"))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)=="CHROM"),which(colnames(subgenome_1_noSD)!="CHROM"))]
  subgenome_1_noSD <- subgenome_1_noSD[,c(which(colnames(subgenome_1_noSD)=="SNP"),which(colnames(subgenome_1_noSD)!="SNP"))]
  write.table (subgenome_1_noSD, file=paste(pop,"_4x","_rd",rd+1,"_noSDbinary.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

  if (snpformats == "true") {  
    alleles <- unique(subset(subgenome_1_noSD, select=c(4,5)))
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
      output <- subset(subgenome_1_noSD, REF == REFsub & ALT == ALTsub)
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
    write.table (geno, file=paste(pop,"_4x","_rd",rd+1,"_noSDnucleotide.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
  }
  
  subgenome_1_noSD[][subgenome_1_noSD[]=="0/0/0/0"] <- "0"
  subgenome_1_noSD[][subgenome_1_noSD[]=="0/0/0/1"] <- "1"
  subgenome_1_noSD[][subgenome_1_noSD[]=="0/0/1/1"] <- "2"
  subgenome_1_noSD[][subgenome_1_noSD[]=="0/1/1/1"] <- "3"
  subgenome_1_noSD[][subgenome_1_noSD[]=="1/1/1/1"] <- "4"
  write.table (subgenome_1_noSD, file=paste(pop,"_4x","_rd",rd+1,"_noSDdose.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")

  sumfreq <- read.table(paste(pop,"_4x","_rd",rd+1,"_noSDdose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = FALSE)
  sumfreq <- subset(sumfreq, select=-c(1:5))
  sumfreq <- subset(sumfreq, select=-c(pvalue))
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
    labs(title="missing rate per Variant",x="Percent", y = "Count")
  ggsave(filename=paste(pop,"_4x","_Variant_missing_rate_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=600, compression = "lzw")
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
  ggsave(filename=paste(pop,"_4x","_sample_missing_rate_rd",rd+1,".tiff",sep=""), plot=plot, width=5, height= 5, dpi=600, compression = "lzw")
}
final_summary()

pop_struc <- function() {
  for (i in c(0,0.01,0.05,0.1)) {
    pop_data <- read.table(paste(pop,"_4x","_rd",rd+1,"_noSDdose.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE, check.names = FALSE)
    pop_data <- subset(pop_data, select=-c(1:5))
    pop_data <- subset(pop_data, select=-c(pvalue))
    pop_data$no_missing <- apply(pop_data, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )
    pop_data <- subset(pop_data, no_missing <= ncol(pop_data)*i)
    pop_data <- subset(pop_data, select=-c(no_missing))
    if (nrow(pop_data) >= 1000) {
      print(paste("missing rate = ",i))
      break
    }
  }

  if (nrow(pop_data) >= 1000) {
    pop_data <- as.matrix(t(pop_data))
    
    #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
    Gmatrix <- function (SNPmatrix = pop_data, method = "VanRaden", missingValue = NA, 
                         maf = 0, thresh.missing = 0.1, verify.posdef = FALSE, ploidy = 4, 
                         pseudo.diploid = FALSE, integer = FALSE, ratio = FALSE, impute.method = FALSE, 
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
    
    Gmat <- Gmatrix()
    write.table (Gmat, file=paste(pop,"_4x","_rd",rd+1,"_Kinship_Matrix.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    tiff(paste(pop,"_",ncol(pop_data),"markers_relatedness_heatmap_dendogram_4x.tiff",sep=""), width=30, height=30, units = 'in', res = 300, compression = 'lzw')
    heatmap(as.matrix(Gmat))
    dev.off()

    pedigree <- as.data.frame(Gmat[,c(p1,p2)])
    colnames(pedigree) <- c("Distance_p1","Distance_p2")
    ped <- ggplot(pedigree, aes(x=Distance_p1, y=Distance_p2, label=row.names(Gmat))) +
      geom_text(size=5, color="cornflowerblue", nudge_y = -0.02)+
      geom_point(color="gray20", size=3) +
      geom_point(data=pedigree[c(p1,p2),], aes(x=Distance_p1,y=Distance_p2, label=row.names(pedigree[c(p1,p2),]), colour=c(p1,p2)), color="tomato") +
      xlab("Kinship Coefficient Relative to Maternal Parent (P1)") + 
      ylab("Kinship Coefficient Relative to Paternal Parent (P2)")
    ggsave(filename=paste(pop,"_",ncol(pop_data),"markers_Parentage_Inference_4x.tiff",sep=""), plot=ped, width=15, height= 15, dpi=600, compression = "lzw")
    
  } else {
    print ("Not enough markers to compute Gmatrix (i.e. threshold of 100 markers)")
  }
}
pop_struc()

####################################################################################################################
subgenome_1 <- read.table (file=paste(pop,"_4x","_DP_GT.txt",sep=""), header=T, sep="\t", check.names = FALSE)
rd_boxplot <- function(){
  #######################################################################################################################################################################################
  # Now, let's make boxplots use raw data set. It reflects read depth distribution across sample IDs
  for (t in c(0,rd)) {
    subgenome_1_plots <- subgenome_1[c(5:(((ncol(subgenome_1)-4)/2)+4))]
    subgenome_1_plots[][subgenome_1_plots[]=="0"] <- NA
    subgenome_1_plots$no_missing <- rowSums(is.na(subgenome_1_plots))
    subgenome_1_plots <- subset(subgenome_1_plots, no_missing <= (ncol(subgenome_1_plots)-5)*gmissingness)
    subgenome_1_plots <- subset(subgenome_1_plots, select=-c(no_missing))
    subgenome_1_plots[][subgenome_1_plots[] <= t] <- NA
    names(subgenome_1_plots) <- gsub(paste("_DP",sep=""), "", names(subgenome_1_plots))
    subgenome_1_plots$no_missing <- rowSums(is.na(subgenome_1_plots))
    subgenome_1_plots <- subset(subgenome_1_plots, no_missing <= (ncol(subgenome_1_plots)-5)*gmissingness)
    subgenome_1_plots <- subset(subgenome_1_plots, select=-c(no_missing))
    names(subgenome_1_plots)[names(subgenome_1_plots) == p1] <- paste("P1_",p1,sep="")
    names(subgenome_1_plots)[names(subgenome_1_plots) == p2] <- paste("P2_",p2,sep="")
    subgenome_1_plots <- as.matrix(subgenome_1_plots)
    subgenome_1_plots <- as.data.frame(as.table(subgenome_1_plots))
    subgenome_1_plots <- na.omit(subgenome_1_plots)
    subgenome_1_boxplot <- subset(subgenome_1_plots, select=c(2,3))
    names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Var2"] <- "samples"
    names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Freq"] <- "DP"
    
    subgenome_1_boxplot$DP <- as.numeric(as.character(subgenome_1_boxplot$DP))
    subgenome_1_boxplot <- na.omit(subgenome_1_boxplot)
    quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.999), na.rm= TRUE)
    nsamples <- length(unique(subgenome_1_boxplot$samples))
    if (nsamples <= 20) { nsamples <- 36 } else { nsamples <- sqrt(22500/nsamples) }
    boxplot <- ggplot(subgenome_1_boxplot, aes(x = reorder(samples,DP,na.rm = TRUE), y=DP), stat='identity')+
      geom_boxplot(fill="white", colour="cornflowerblue",
                   outlier.alpha = 0.01, outlier.colour="tomato", outlier.size=1.0)+
      coord_flip()+
      scale_y_continuous(expand = c(0,0), limits = c(0,quantile999)) +
      theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
            axis.text.y=element_text(colour="cornflowerblue", size=nsamples),
            axis.title=element_text(size=36)) +
      xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
      ylab("Read Depth (4x Genotypes)")
    ggsave(filename= paste(pop,"_4x","_boxplot_rd",t+1,".tiff",sep=""), plot=boxplot, width=15, height= 25, dpi=600, compression = "lzw")
    
    meanDP <- mean(subgenome_1_boxplot$DP, na.rm=T)
    medianDP <- median(subgenome_1_boxplot$DP, na.rm=T)
    maxDP <- max(table(subgenome_1_boxplot$DP))
    quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.999), na.rm= TRUE)
    subgenome_1_dist <- as.data.frame(table(subgenome_1_boxplot$DP))
    subgenome_1_dist$Var1 <- as.numeric(as.character(subgenome_1_dist$Var1))
    subgenome_1_dist <- subset(subgenome_1_dist, Freq > quantile999)
    maxX <- max(subgenome_1_dist[,1])
    histogram <- ggplot(subgenome_1_dist, aes(x=Var1, y=Freq)) +
      geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="cornflowerblue", fill="white")+
      geom_vline(aes(xintercept=meanDP), color="cornflowerblue", linetype="dashed", size=3, alpha=0.5)+
      geom_vline(aes(xintercept=medianDP), color="tomato", linetype="dotted", size=3, alpha=0.5)+
      geom_text(aes(x=meanDP, label=paste("mean = ",round(meanDP),sep=""), y=(maxDP*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=7.5) +
      geom_text(aes(x=medianDP, label=paste("median = ",round(medianDP),sep=""), y=(maxDP*0.5)), colour="tomato", angle=90, vjust = 1.2, size=7.5) +
      scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",scientific = T)) +
      scale_x_continuous(breaks=seq(0,maxX,round(maxX/10,-1))) +
      theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
            axis.text.y=element_text(colour="cornflowerblue", size=24),
            axis.title=element_text(size=30)) +
      ylab("Frequency") + 
      xlab(paste("Read Depth Distribution (", p1," x ",p2," Biparental Mapping Population)",sep=""))
    ggsave(filename= paste(pop,"_4x","_histogram_rd",t+1,".tiff",sep=""), plot=histogram, width=25, height= 15, dpi=600, compression = "lzw")

    boxplot <- NULL
    subgenome_1_boxplot <- NULL
    subgenome_1_plots <- NULL
    gc()
  }
  # Extract read depth values specifically for filtered SNPs, then plot boxplot 
  # Also, plot histogram of read depth across data set
  subgenome_1_plots <- subgenome_1
  subgenome_final <- read.table (file=paste(pop,"_4x","_rd",rd+1,"_noSDdose.txt",sep=""), header=T, sep="\t", check.names = FALSE)
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
  subgenome_1_plots <- subgenome_1[c(5:(((ncol(subgenome_1)-5)/2)+4))]
  for (i in 0:rd) {
    subgenome_1_plots[][subgenome_1_plots[]==i] <- NA
  }
  names(subgenome_1_plots) <- gsub(paste("_DP",sep=""), "", names(subgenome_1_plots))
  names(subgenome_1_plots)[names(subgenome_1_plots) == p1] <- paste("P1_",p1,sep="")
  names(subgenome_1_plots)[names(subgenome_1_plots) == p2] <- paste("P2_",p2,sep="")
  subgenome_1_plots <- as.matrix(subgenome_1_plots)
  subgenome_1_plots <- as.data.frame(as.table(subgenome_1_plots))
  subgenome_1_plots <- na.omit(subgenome_1_plots)
  subgenome_1_boxplot <- subset(subgenome_1_plots, select=c(2,3))
  names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Var2"] <- "samples"
  names(subgenome_1_boxplot)[names(subgenome_1_boxplot) == "Freq"] <- "DP"
  
  subgenome_1_boxplot$DP <- as.numeric(as.character(subgenome_1_boxplot$DP))
  subgenome_1_boxplot <- na.omit(subgenome_1_boxplot)
  quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.999), na.rm= TRUE)
  nsamples <- length(unique(subgenome_1_boxplot$samples))
  if (nsamples <= 20) { nsamples <- 36 } else { nsamples <- sqrt(22500/nsamples) }
  boxplot <- ggplot(subgenome_1_boxplot, aes(x = reorder(samples,DP,na.rm = TRUE), y=DP), stat='identity')+
    geom_boxplot(fill="white", colour="cornflowerblue",
                 outlier.alpha = 0.01, outlier.colour="tomato", outlier.size=1.0)+
    coord_flip()+
    scale_y_continuous(expand = c(0,0), limits = c(0,quantile999)) +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
          axis.text.y=element_text(colour="cornflowerblue", size=nsamples),
          axis.title=element_text(size=36)) +
    xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
    ylab("Read Depth (4x Genotypes)")
  ggsave(filename= paste(pop,"_4x","_boxplot_filtered.tiff",sep=""), plot=boxplot, width=15, height= 25, dpi=600, compression = "lzw")
  
  meanDP <- mean(subgenome_1_boxplot$DP, na.rm=T)
  medianDP <- median(subgenome_1_boxplot$DP, na.rm=T)
  maxDP <- max(table(subgenome_1_boxplot$DP))
  quantile999 <- quantile(subgenome_1_boxplot$DP, probs = c(0.999), na.rm= TRUE)
  subgenome_1_dist <- as.data.frame(table(subgenome_1_boxplot$DP))
  subgenome_1_dist$Var1 <- as.numeric(as.character(subgenome_1_dist$Var1))
  subgenome_1_dist <- subset(subgenome_1_dist, Freq > quantile999)
  maxX <- max(subgenome_1_dist[,1])
  histogram <- ggplot(subgenome_1_dist, aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="cornflowerblue", fill="white")+
    geom_vline(aes(xintercept=meanDP), color="cornflowerblue", linetype="dashed", size=3, alpha=0.5)+
    geom_vline(aes(xintercept=medianDP), color="tomato", linetype="dotted", size=3, alpha=0.5)+
    geom_text(aes(x=meanDP, label=paste("mean = ",round(meanDP),sep=""), y=(maxDP*0.25)), colour="cornflowerblue", angle=90, vjust = 1.2, size=7.5) +
    geom_text(aes(x=medianDP, label=paste("median = ",round(medianDP),sep=""), y=(maxDP*0.5)), colour="tomato", angle=90, vjust = 1.2, size=7.5) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",scientific = T)) +
    scale_x_continuous(breaks=seq(0,maxX,round(maxX/10,-1))) +
    theme(axis.text.x=element_text(colour="cornflowerblue", size=24),
          axis.text.y=element_text(colour="cornflowerblue", size=24),
          axis.title=element_text(size=30)) +
    ylab("Frequency") + 
    xlab(paste("Read Depth Distribution (", p1," x ",p2," Biparental Mapping Population)",sep=""))
  ggsave(filename= paste(pop,"_4x","_histogram_filtered.tiff",sep=""), plot=histogram, width=25, height= 15, dpi=600, compression = "lzw")
  boxplot <- NULL
  subgenome_1_boxplot <- NULL
  subgenome_1_plots <- NULL
  gc()
}
rd_boxplot()

raw_alleles <- function(){
  #######################################################################################################################################################################################
  #plot all markers classes based on number of alleles present
  #change DP >0 and < rd to DP=1, consider this as missing and use this to remove genotypes greater than 20% missing (progenies)
  subgenome_1_plots <- subgenome_1
  subgenome_1_plots <- subset(subgenome_1_plots, select=c("CHROM", "POS", paste(p1,"_DP",sep=""), paste(p2,"_DP",sep=""),
                                                          paste(p1,"_GT",sep=""), paste(p2,"_GT",sep="")))
  subgenome_1_plots[,paste(p1,"_DP",sep="")] <- as.numeric(as.character(subgenome_1_plots[,paste(p1,"_DP",sep="")]))
  subgenome_1_plots[,paste(p2,"_DP",sep="")] <- as.numeric(as.character(subgenome_1_plots[,paste(p2,"_DP",sep="")]))
  subgenome_1_plots$DiffMethyl <- subgenome_1_plots[,paste(p1,"_DP",sep="")] * subgenome_1_plots[,paste(p2,"_DP",sep="")]
  subgenome_1_plots <- subset (subgenome_1_plots, DiffMethyl > 0)
  subgenome_1_plots <- subset(subgenome_1_plots, select=-c(DiffMethyl))
  for (i in 3:(((ncol(subgenome_1_plots)-2)/2)+2)) {
    j <- i+((ncol(subgenome_1_plots)-2)/2)
    subgenome_1_plots[,j][subgenome_1_plots[,i] <= rd] <- NA
  }
  subgenome_1_plots <-na.omit(subgenome_1_plots)
  #If progeny read depth="0" and both parents have DP>0
  for (i in 3:4){
    subgenome_1_plots[,i] <- as.numeric(as.character(subgenome_1_plots[,i]))
  }
  subgenome_1_plots$DiffMethyl <- subgenome_1_plots[,paste(p1,"_DP",sep="")] * subgenome_1_plots[,paste(p2,"_DP",sep="")]
  subgenome_1_plots <- subset (subgenome_1_plots, DiffMethyl > 0)
  subgenome_1_plots$DiffMethyl <- subgenome_1_plots[,paste(p1,"_DP",sep="")] + subgenome_1_plots[,paste(p2,"_DP",sep="")]
  subgenome_1_plots <- subset (subgenome_1_plots, DiffMethyl != 1)
  subgenome_1_plots <- subset (subgenome_1_plots, DiffMethyl != 2)
  subgenome_1_plots <- subset(subgenome_1_plots, select=c(-7))
  subgenome_1_plots <- subset(subgenome_1_plots, select=c(5,6))
  
  Cross <- subgenome_1_plots
  Cross <- na.omit(Cross)
  Cross$Cross <- paste (Cross[,paste(p1,"_GT",sep="")], Cross[,paste(p2,"_GT",sep="")], sep=" x ")
  Cross <- subset(Cross, select="Cross")
  Cross <- as.data.frame(table(Cross))
  sum <- sum(Cross$Freq)
  Cross$percentage <- ((Cross$Freq)/sum)*100
  Cross <- subset(Cross, percentage >= 0.05)
  Cross[,3] <- round(Cross[,3], 2)
  sum <- sum(Cross$Freq)
  max <- max(Cross$Freq)
  max <- max *1.2
  plot <- ggplot(Cross, aes(x=Cross, y=Freq, fill=Cross, group=Cross)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="black")+
    geom_text(aes(x=Cross, y=Freq, label = paste(Freq,"=", percentage, "%"), group=Cross),
              position=position_dodge(0.95), hjust = -0.1, size=5, color="black", fontface="italic")+
    coord_flip()+
    theme(axis.text.x=element_text(colour="cornflowerblue", size=14),
          axis.text.y=element_text(colour="cornflowerblue", size=14),
          axis.title=element_text(size=20)) +
    theme(legend.position="none")+
    theme(legend.position="none")+
    theme(legend.text=element_text(colour="cornflowerblue",size=14)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
    ylab(paste("Proportion of Genotypic Configuration (",sum,")", sep=""))
  ggsave(filename="Multiallelic_cross_4x.tiff", plot=plot, width=15, height= 10, dpi=600, compression = "lzw")
  
  Cross <- subset(Cross, Cross!="0/0/0/0 x 0/0/0/0")
  Cross <- subset(Cross, Cross!="1/1/1/1 x 1/1/1/1")
  Cross <- subset(Cross, Cross!="2/2/2/2 x 2/2/2/2")
  Cross <- subset(Cross, Cross!="3/3/3/3 x 3/3/3/3")
  Cross <- subset(Cross, Cross!="4/4/4/4 x 4/4/4/4")
  Cross <- subset(Cross, Cross!="5/5/5/5 x 5/5/5/5")
  Cross <- subset(Cross, select=c(-3))
  sum <- sum(Cross$Freq)
  Cross$percentage <- ((Cross$Freq)/sum)*100
  Cross <- subset(Cross, percentage >= 0.05)
  Cross[,3] <- round(Cross[,3], 2)
  sum <- sum(Cross$Freq)
  max <- max(Cross$Freq)
  max <- max *1.2
  plot <- ggplot(Cross, aes(x=Cross, y=Freq, fill=Cross, group=Cross)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="black")+
    geom_text(aes(x=Cross, y=Freq, label = paste(Freq,"=", percentage, "%"), group=Cross),
              position=position_dodge(0.95), hjust = -0.1, size=5, color="black", fontface="italic")+
    coord_flip()+
    theme(axis.text.x=element_text(colour="cornflowerblue", size=14),
          axis.text.y=element_text(colour="cornflowerblue", size=14),
          axis.title=element_text(size=20)) +
    theme(legend.position="none")+
    theme(legend.text=element_text(colour="cornflowerblue",size=14)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab(paste(p1," x ",p2," Biparental Mapping Population",sep="")) +
    ylab(paste("Proportion of Genotypic Configuration (",sum,")", sep=""))
  ggsave(filename="Multiallelic_cross_subset_4x.tiff", plot=plot, width=15, height= 10, dpi=600, compression = "lzw")
  plot <- NULL
  
  Multiallelic_p1_4x <- as.data.frame(table(subgenome_1[,paste(p1,"_GT",sep="")]))
  names(Multiallelic_p1_4x)[names(Multiallelic_p1_4x) == "Var1"] <- "Genotype"
  Multiallelic_p1_4x <- subset(Multiallelic_p1_4x, Genotype != "./././.")
  sum <- sum(Multiallelic_p1_4x$Freq)
  Multiallelic_p1_4x$percentage <- ((Multiallelic_p1_4x$Freq)/sum)*100
  Multiallelic_p1_4x <- subset(Multiallelic_p1_4x, percentage >= 0.05)
  Multiallelic_p1_4x[,3] <- round(Multiallelic_p1_4x[,3], 2)
  Multiallelic_p1_4x$Parents <- rep(paste(p1),nrow(Multiallelic_p1_4x))
  Multiallelic_p2_4x <- as.data.frame(table(subgenome_1[,paste(p2,"_GT",sep="")]))
  names(Multiallelic_p2_4x)[names(Multiallelic_p2_4x) == "Var1"] <- "Genotype"
  Multiallelic_p2_4x <- subset(Multiallelic_p2_4x, Genotype !="./././.")
  sum <- sum(Multiallelic_p2_4x$Freq)
  Multiallelic_p2_4x$percentage <- ((Multiallelic_p2_4x$Freq)/sum)*100
  Multiallelic_p2_4x <- subset(Multiallelic_p2_4x, percentage >= 0.05)
  Multiallelic_p2_4x[,3] <- round(Multiallelic_p2_4x[,3], 2)
  Multiallelic_p2_4x$Parents <- rep(paste(p2),nrow(Multiallelic_p2_4x))
  Multiallelic_4x <- rbind(Multiallelic_p1_4x, Multiallelic_p2_4x)
  sum <- sum(Multiallelic_4x$Freq)
  max <- max(Multiallelic_4x$Freq)
  max <- max *1.2
  plot <- ggplot(Multiallelic_4x, aes(x=Parents, y=Freq, fill=Genotype, group=Genotype)) +
    geom_bar(stat="identity", position=position_dodge(0.95), width=0.9, colour="black")+
    geom_text(aes(x=Parents, y=Freq, label = paste(Genotype," (", Freq, "=", percentage, "%)", sep=""), group=Genotype),
              position=position_dodge(0.95), hjust = -0.1, size=5, color="black", fontface="italic")+
    coord_flip()+
    theme(axis.text.x=element_text(colour="cornflowerblue", size=14),
          axis.text.y=element_text(colour="cornflowerblue", size=14),
          axis.title=element_text(size=20)) +
    theme(legend.position="none")+
    theme(legend.text=element_text(colour="cornflowerblue",size=14)) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    guides(fill=guide_legend(ncol=1))+
    scale_y_continuous(expand = c(0, 0), labels = function(x) format(x, big.mark = ",",
                                                                     scientific = FALSE)) +
    expand_limits(y = c(0, max))+
    xlab("Genotypes") +
    ylab(paste("Proportion of Genotypes (",sum,")", sep=""))
  ggsave(filename="Multiallelic_4x.tiff", plot=plot, width=15, height= 10, dpi=600, compression = "lzw")
  plot <- NULL
}
raw_alleles()

####################################################################################################################


