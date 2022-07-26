#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]
ploidy <- args[2]
libdir <- args[3]
.libPaths( c( .libPaths(), libdir) )
library(ggplot2)


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
      names(vcffile_GT)[5:length(cnames)] <- paste(colnames(vcffile_GT[,5:length(cnames)]), "_GT", sep="")
      names(vcffile_DP)[5:length(cnames)] <- paste(colnames(vcffile_DP[,5:length(cnames)]), "_DP", sep="")
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="."] <- NA
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/2)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      DPstart <- 5; DPend <- (((ncol(vcffile)-4)/2)+4); increments <- DPend - DPstart
      GTstart <- DPend + 1; GTend <-  GTstart +increments
      sample_size <- increments + 1
      vcffile[, 5:DPend] <- lapply(vcffile[,5:DPend], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:DPend] <- lapply(5:DPend, function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:DPend] == 0, na.rm = TRUE) <= ((ncol(vcffile)-4)/3)*0.8, ]
      vcffile$GTsum <- (rowSums(!is.na(vcffile))-(4 + sample_size + sample_size))
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
      vcffile_GT <- vcffile; vcffile_DP <- vcffile; vcffile_AR <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_AR[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile_AR_label <- vcffile_AR[,1:4]
      nosplit <- vcffile_AR_label
      for (k in c(5:(ncol(vcffile_AR)))) {
        splitAR <- strsplit(as.character(vcffile_AR[,k]),',')
        suppressWarnings(splitAR <- as.data.frame(do.call(rbind, splitAR)))
        suppressWarnings(splitAR <- as.data.frame(as.numeric(splitAR$V1) / as.numeric(splitAR$V2)))
        splitAR[][splitAR[]=="NaN"] <- NA
        splitAR[][splitAR[]=="Inf"] <- "0"
        splitAR <- cbind(vcffile_AR_label,splitAR)
        names(splitAR)[5] <- colnames(vcffile_AR[k])
        splitAR[,5] <- as.numeric(splitAR[,5])
        splitAR_na <- splitAR[is.na(splitAR[,5]),]
        splitAR_na[splitAR_na=='NA'] <- NA
        splitAR_pos <- subset(splitAR, splitAR[,5] <= 1)
        splitAR_neg <- subset(splitAR, splitAR[,5] > 1)
        splitAR_neg[,5] <- -(1/as.numeric(splitAR_neg[,5]))
        splitAR <- rbind(splitAR_pos, splitAR_neg)
        splitAR <- rbind(splitAR, splitAR_na)
        nosplit <- merge(nosplit,splitAR, by=c("CHROM","POS","REF","ALT"))
      }
      vcffile_AR <- nosplit
      
      
      names(vcffile_GT)[5:length(cnames)] <- paste(colnames(vcffile_GT[,5:length(cnames)]), "_GT", sep="")
      names(vcffile_DP)[5:length(cnames)] <- paste(colnames(vcffile_DP[,5:length(cnames)]), "_DP", sep="")
      names(vcffile_AR)[5:length(cnames)] <- paste(colnames(vcffile_AR[,5:length(cnames)]), "_AR", sep="")
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile <- merge(vcffile,vcffile_AR, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./."] <- NA
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/3)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      DPstart <- 5; DPend <- (((ncol(vcffile)-4)/3)+4); increments <- DPend - DPstart
      GTstart <- DPend + 1; GTend <-  GTstart +increments
      ARstart <- GTend + 1; ARend <- ARstart + increments
      sample_size <- increments + 1
      vcffile[, 5:DPend] <- lapply(vcffile[,5:DPend], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:DPend] <- lapply(5:DPend, function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:DPend] == 0, na.rm = TRUE) <= ((ncol(vcffile)-4)/3)*0.8, ]
      vcffile$GTsum <- (rowSums(!is.na(vcffile))-(4 + sample_size + sample_size))
      vcffile[vcffile=="1/0"] <- "0/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    
    GTincre <- (ARend-ARstart) + 1
    for (m in c(ARstart:ARend)) {
      n <- m-GTincre
      subgenome_1[,n][subgenome_1[,m] > 0 && subgenome_1[,m] < 0.2] <- "1/1"
      subgenome_1[,n][subgenome_1[,m] < 0 && subgenome_1[,m] > -0.2] <- "0/0"
      gc()
    }
    AR <- subgenome_1
    AR1 <- subgenome_1[,c(1:4,ARstart:ARend)]
    subgenome_1 <- subgenome_1[,-c(ARstart:ARend)]
    write.table (subgenome_1, file=paste(pop,"_2x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    write.table (AR1, file=paste(pop,"_2x","_AR.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    AR1 <- NULL
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*2x_rawSPLIT*",sep=""))
    
    
    AR$propHet <- (rowSums(AR == "0/1" | AR == "0|1", na.rm = TRUE)) / (GTincre - (rowSums(AR == "./.", na.rm = TRUE)))
    AR <- AR[,c(ARstart:ARend,ncol(AR))]
    AR <- subset(AR, AR$propHet != 0)
    AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
    AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
    suppressWarnings(AR <- na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
    AR <- AR[,1:2]; names(AR) <- c("propHet","Allele_Ratio")
    
    ARplot <- ggplot(AR,aes(x=propHet,y=Allele_Ratio))+
      geom_point(aes(x = propHet, y = Allele_Ratio, color = Allele_Ratio), size = 1, pch=19, alpha=0.1)+
      geom_density_2d(bins=50)+
      annotate("rect", xmin=0, xmax=1, ymin=-0.2, ymax=0.2, alpha=0.2, fill="tomato")+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      scale_colour_gradient2(low="darkorange3", mid="darkgoldenrod1", high ="cornflowerblue",
                             breaks=c(0.75,-0.75), limits=c(-1,1), 
                             labels=c("Minor Allele: Ref","Minor Allele: Alt"))+
      theme(legend.title = element_blank())+
      geom_hline(yintercept = 0.2, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.2, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 1, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -1, color="grey20", size=0.5, linetype="dashed")+ 
      annotate("text", x=0.92, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=1, label="Heterozygote (balanced allele ratio)", vjust=1.0, fontface="italic")+
      annotate("text", x=0.85, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1.0, fontface="italic")+
      annotate("text", x=0.86, y=0.2, label="Allele Ratio filter threshold (> 0.2)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.86, y=-0.2, label="Allele Ratio filter threshold (< -0.2)", vjust=1.2, fontface="italic")+
      xlab("Proportion of Heterozygote per Locus (diploid)") +
      ylab("Allele Read Depth Ratio per Genotype")
    ggsave(file=paste("raw2x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
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
      vcffile_GT <- vcffile; vcffile_DP <- vcffile; vcffile_AR <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_AR[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile_AR_label <- vcffile_AR[,1:4]
      nosplit <- vcffile_AR_label
      for (k in c(5:(ncol(vcffile_AR)))) {
        splitAR <- strsplit(as.character(vcffile_AR[,k]),',')
        suppressWarnings(splitAR <- as.data.frame(do.call(rbind, splitAR)))
        suppressWarnings(splitAR <- as.data.frame(as.numeric(splitAR$V1) / as.numeric(splitAR$V2)))
        splitAR[][splitAR[]=="NaN"] <- NA
        splitAR[][splitAR[]=="Inf"] <- "0"
        splitAR <- cbind(vcffile_AR_label,splitAR)
        names(splitAR)[5] <- colnames(vcffile_AR[k])
        splitAR[,5] <- as.numeric(splitAR[,5])
        splitAR_na <- splitAR[is.na(splitAR[,5]),]
        splitAR_na[splitAR_na=='NA'] <- NA
        splitAR_pos <- subset(splitAR, splitAR[,5] <= 1)
        splitAR_neg <- subset(splitAR, splitAR[,5] > 1)
        splitAR_neg[,5] <- -(1/as.numeric(splitAR_neg[,5]))
        splitAR <- rbind(splitAR_pos, splitAR_neg)
        splitAR <- rbind(splitAR, splitAR_na)
        nosplit <- merge(nosplit,splitAR, by=c("CHROM","POS","REF","ALT"))
      }
      vcffile_AR <- nosplit
      
      
      names(vcffile_GT)[5:length(cnames)] <- paste(colnames(vcffile_GT[,5:length(cnames)]), "_GT", sep="")
      names(vcffile_DP)[5:length(cnames)] <- paste(colnames(vcffile_DP[,5:length(cnames)]), "_DP", sep="")
      names(vcffile_AR)[5:length(cnames)] <- paste(colnames(vcffile_AR[,5:length(cnames)]), "_AR", sep="")
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile <- merge(vcffile,vcffile_AR, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././."] <- NA
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/3)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      DPstart <- 5; DPend <- (((ncol(vcffile)-4)/3)+4); increments <- DPend - DPstart
      GTstart <- DPend + 1; GTend <-  GTstart +increments
      ARstart <- GTend + 1; ARend <- ARstart + increments
      sample_size <- increments + 1
      vcffile[, 5:DPend] <- lapply(vcffile[,5:DPend], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:DPend] <- lapply(5:DPend, function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:DPend] == 0, na.rm = TRUE) <= ((ncol(vcffile)-4)/3)*0.8, ]
      vcffile$GTsum <- (rowSums(!is.na(vcffile))-(4 + sample_size + sample_size))
      vcffile[vcffile=="1/0/0/0" | vcffile=="0/1/0/0" | vcffile=="0/0/1/0" ] <- "0/0/0/1"
      vcffile[vcffile=="1/1/0/0" | vcffile=="0/1/1/0"] <- "0/0/1/1"
      vcffile[vcffile=="1/1/1/0"] <- "0/1/1/1"
      subgenome_1 <- rbind(subgenome_1,vcffile)
      gc()
    }
    GTincre <- (ARend-ARstart) + 1
    for (m in c(ARstart:ARend)) {
      n <- m-GTincre
      subgenome_1[,n][subgenome_1[,m] > 0 && subgenome_1[,m] < 0.17] <- "1/1/1/1"
      subgenome_1[,n][subgenome_1[,m] < 0 && subgenome_1[,m] > -0.17] <- "0/0/0/0"
      gc()
    }
    AR <- subgenome_1
    AR1 <- subgenome_1[,c(1:4,ARstart:ARend)]
    subgenome_1 <- subgenome_1[,-c(ARstart:ARend)]
    write.table (subgenome_1, file=paste(pop,"_4x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    write.table (AR1, file=paste(pop,"_4x","_AR.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    AR1 <- NULL
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*4x_rawSPLIT*",sep=""))
    
    
    AR$propHet <- (rowSums(AR == "0/0/0/1" | AR == "0/0/1/1" | AR == "0/1/1/1" | 
                           AR == "0|0|0|1" | AR == "0|0|1|1" | AR == "0|1|1|1", na.rm = TRUE)) / (GTincre - (rowSums(AR == "./././.", na.rm = TRUE)))
    AR <- AR[,c(ARstart:ARend,ncol(AR))]
    AR <- subset(AR, AR$propHet != 0)
    AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
    AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
    suppressWarnings(AR <- na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
    AR <- AR[,1:2]; names(AR) <- c("propHet","Allele_Ratio")
    
    ARplot <- ggplot(AR,aes(x=propHet,y=Allele_Ratio))+
      geom_point(aes(x = propHet, y = Allele_Ratio, color = Allele_Ratio), size = 1, pch=19, alpha=0.1)+
      geom_density_2d(bins=50)+
      annotate("rect", xmin=0, xmax=1, ymin=-0.2, ymax=0.2, alpha=0.2, fill="tomato")+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      scale_colour_gradient2(low="darkorange3", mid="darkgoldenrod1", high ="cornflowerblue",
                             breaks=c(0.75,-0.75), limits=c(-1,1), 
                             labels=c("Minor Allele: Ref","Minor Allele: Alt"))+
      theme(legend.title = element_blank())+
      geom_hline(yintercept = 0.17, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.17, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 1, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0.33, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.33, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -1, color="grey20", size=0.5, linetype="dashed")+ 
      annotate("text", x=0.92, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.89, y=0.33, label="Heterozygote (0/0/0/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
      annotate("text", x=0.85, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
      annotate("text", x=0.89, y=-0.33, label="Heterozygote (0/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.85, y=0.17, label="Allele Ratio filter threshold (> 0.17)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=-0.17, label="Allele Ratio filter threshold (< -0.17)", vjust=1.2, fontface="italic")+
      xlab("Proportion of Heterozygote per Locus (tetraploid)") +
      ylab("Allele Read Depth Ratio per Genotype")
    ggsave(file=paste("raw4x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
    
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
      vcffile_GT <- vcffile; vcffile_DP <- vcffile; vcffile_AR <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_AR[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile_AR_label <- vcffile_AR[,1:4]
      nosplit <- vcffile_AR_label
      for (k in c(5:(ncol(vcffile_AR)))) {
        splitAR <- strsplit(as.character(vcffile_AR[,k]),',')
        suppressWarnings(splitAR <- as.data.frame(do.call(rbind, splitAR)))
        suppressWarnings(splitAR <- as.data.frame(as.numeric(splitAR$V1) / as.numeric(splitAR$V2)))
        splitAR[][splitAR[]=="NaN"] <- NA
        splitAR[][splitAR[]=="Inf"] <- "0"
        splitAR <- cbind(vcffile_AR_label,splitAR)
        names(splitAR)[5] <- colnames(vcffile_AR[k])
        splitAR[,5] <- as.numeric(splitAR[,5])
        splitAR_na <- splitAR[is.na(splitAR[,5]),]
        splitAR_na[splitAR_na=='NA'] <- NA
        splitAR_pos <- subset(splitAR, splitAR[,5] <= 1)
        splitAR_neg <- subset(splitAR, splitAR[,5] > 1)
        splitAR_neg[,5] <- -(1/as.numeric(splitAR_neg[,5]))
        splitAR <- rbind(splitAR_pos, splitAR_neg)
        splitAR <- rbind(splitAR, splitAR_na)
        nosplit <- merge(nosplit,splitAR, by=c("CHROM","POS","REF","ALT"))
      }
      vcffile_AR <- nosplit
      
      
      names(vcffile_GT)[5:length(cnames)] <- paste(colnames(vcffile_GT[,5:length(cnames)]), "_GT", sep="")
      names(vcffile_DP)[5:length(cnames)] <- paste(colnames(vcffile_DP[,5:length(cnames)]), "_DP", sep="")
      names(vcffile_AR)[5:length(cnames)] <- paste(colnames(vcffile_AR[,5:length(cnames)]), "_AR", sep="")
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile <- merge(vcffile,vcffile_AR, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././././."] <- NA
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/3)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      DPstart <- 5; DPend <- (((ncol(vcffile)-4)/3)+4); increments <- DPend - DPstart
      GTstart <- DPend + 1; GTend <-  GTstart +increments
      ARstart <- GTend + 1; ARend <- ARstart + increments
      sample_size <- increments + 1
      vcffile[, 5:DPend] <- lapply(vcffile[,5:DPend], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:DPend] <- lapply(5:DPend, function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:DPend] == 0, na.rm = TRUE) <= ((ncol(vcffile)-4)/3)*0.8, ]
      vcffile$GTsum <- (rowSums(!is.na(vcffile))-(4 + sample_size + sample_size))
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
    GTincre <- (ARend-ARstart) + 1
    for (m in c(ARstart:ARend)) {
      n <- m-GTincre
      subgenome_1[,n][subgenome_1[,m] > 0 && subgenome_1[,m] < 0.14] <- "1/1/1/1/1/1"
      subgenome_1[,n][subgenome_1[,m] < 0 && subgenome_1[,m] > -0.14] <- "0/0/0/0/0/0"
      gc()
    }
    AR <- subgenome_1
    AR1 <- subgenome_1[,c(1:4,ARstart:ARend)]
    subgenome_1 <- subgenome_1[,-c(ARstart:ARend)]
    write.table (subgenome_1, file=paste(pop,"_6x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    write.table (AR1, file=paste(pop,"_6x","_AR.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    AR1 <- NULL
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*6x_rawSPLIT*",sep=""))
    
    AR$propHet <- (rowSums(AR == "0/0/0/0/0/1" | AR == "0/0/0/0/1/1" | AR == "0/0/0/1/1/1" | AR == "0/0/0/1/1/1" | AR == "0/0/1/1/1/1" | AR == "0/1/1/1/1/1" | 
                             AR == "0|0|0|0|0|1" | AR == "0|0|0|0|1|1" | AR == "0|0|0|1|1|1" | AR == "0|0|0|0|0|1" | AR == "0|0|0|0|1|1" | AR == "0|0|0|1|1|1", na.rm = TRUE)) / (GTincre - (rowSums(AR == "./././././.", na.rm = TRUE)))
    AR <- AR[,c(ARstart:ARend,ncol(AR))]
    AR <- subset(AR, AR$propHet != 0)
    AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
    AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
    suppressWarnings(AR <- na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
    AR <- AR[,1:2]; names(AR) <- c("propHet","Allele_Ratio")
    
    ARplot <- ggplot(AR,aes(x=propHet,y=Allele_Ratio))+
      geom_point(aes(x = propHet, y = Allele_Ratio, color = Allele_Ratio), size = 1, pch=19, alpha=0.1)+
      geom_density_2d(bins=50)+
      annotate("rect", xmin=0, xmax=1, ymin=-0.2, ymax=0.2, alpha=0.2, fill="tomato")+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      scale_colour_gradient2(low="darkorange3", mid="darkgoldenrod1", high ="cornflowerblue",
                             breaks=c(0.75,-0.75), limits=c(-1,1), 
                             labels=c("Minor Allele: Ref","Minor Allele: Alt"))+
      theme(legend.title = element_blank())+
      geom_hline(yintercept = 0.14, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.14, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 1, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0.2, color="grey20", size=0.5, linetype="dashed")+
      geom_hline(yintercept = 0.5, color="grey20", size=0.5, linetype="dashed")+
      geom_hline(yintercept = -0.2, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.5, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -1, color="grey20", size=0.5, linetype="dashed")+ 
      annotate("text", x=0.92, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.88, y=0.2, label="Heterozygote (0/0/0/0/0/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.88, y=0.5, label="Heterozygote (0/0/0/0/1/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
      annotate("text", x=0.85, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
      annotate("text", x=0.88, y=-0.5, label="Heterozygote (0/0/1/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.88, y=-0.2, label="Heterozygote (0/1/1/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.85, y=0.14, label="Allele Ratio filter threshold (< 0.14)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=-0.14, label="Allele Ratio filter threshold (> -0.14)", vjust=1.2, fontface="italic")+
      xlab("Proportion of Heterozygote per Locus (hexaploid)") +
      ylab("Allele Read Depth Ratio per Genotype")
    ggsave(file=paste("raw6x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")

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
      vcffile_GT <- vcffile; vcffile_DP <- vcffile; vcffile_AR <- vcffile
      all_content <- NULL
      for (j in 5:length(cnames)) {
        vcffile_GT[,j] <- gsub(":.*", "", vcffile_GT[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_AR[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("^[^:]*:", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub(":.*", "", vcffile_DP[,j])
        vcffile_DP[,j] <- gsub("./././././././.", "0", vcffile_DP[,j], fixed = TRUE)
        vcffile_DP[,j] <- gsub(".", "0", vcffile_DP[,j], fixed = TRUE)
      }
      vcffile_AR_label <- vcffile_AR[,1:4]
      nosplit <- vcffile_AR_label
      for (k in c(5:(ncol(vcffile_AR)))) {
        splitAR <- strsplit(as.character(vcffile_AR[,k]),',')
        suppressWarnings(splitAR <- as.data.frame(do.call(rbind, splitAR)))
        suppressWarnings(splitAR <- as.data.frame(as.numeric(splitAR$V1) / as.numeric(splitAR$V2)))
        splitAR[][splitAR[]=="NaN"] <- NA
        splitAR[][splitAR[]=="Inf"] <- "0"
        splitAR <- cbind(vcffile_AR_label,splitAR)
        names(splitAR)[5] <- colnames(vcffile_AR[k])
        splitAR[,5] <- as.numeric(splitAR[,5])
        splitAR_na <- splitAR[is.na(splitAR[,5]),]
        splitAR_na[splitAR_na=='NA'] <- NA
        splitAR_pos <- subset(splitAR, splitAR[,5] <= 1)
        splitAR_neg <- subset(splitAR, splitAR[,5] > 1)
        splitAR_neg[,5] <- -(1/as.numeric(splitAR_neg[,5]))
        splitAR <- rbind(splitAR_pos, splitAR_neg)
        splitAR <- rbind(splitAR, splitAR_na)
        nosplit <- merge(nosplit,splitAR, by=c("CHROM","POS","REF","ALT"))
      }
      vcffile_AR <- nosplit
      
      
      names(vcffile_GT)[5:length(cnames)] <- paste(colnames(vcffile_GT[,5:length(cnames)]), "_GT", sep="")
      names(vcffile_DP)[5:length(cnames)] <- paste(colnames(vcffile_DP[,5:length(cnames)]), "_DP", sep="")
      names(vcffile_AR)[5:length(cnames)] <- paste(colnames(vcffile_AR[,5:length(cnames)]), "_AR", sep="")
      vcffile <- merge(vcffile_DP,vcffile_GT, by=c("CHROM","POS","REF","ALT"))
      vcffile <- merge(vcffile,vcffile_AR, by=c("CHROM","POS","REF","ALT"))
      vcffile[][vcffile[]=="./././././././."] <- NA
      vcffile$no_missing <- rowSums(is.na(vcffile))
      vcffile <- subset(vcffile, no_missing <= ((ncol(vcffile)-5)/3)*0.8)
      vcffile <- subset(vcffile, select=-c(no_missing))
      DPstart <- 5; DPend <- (((ncol(vcffile)-4)/3)+4); increments <- DPend - DPstart
      GTstart <- DPend + 1; GTend <-  GTstart +increments
      ARstart <- GTend + 1; ARend <- ARstart + increments
      sample_size <- increments + 1
      vcffile[, 5:DPend] <- lapply(vcffile[,5:DPend], gsub, pattern = "0,0", replacement = "0", fixed = TRUE)
      vcffile[, 5:DPend] <- lapply(5:DPend, function(x) as.numeric(vcffile[[x]]))
      vcffile <- vcffile[rowSums(vcffile[, 5:DPend] == 0, na.rm = TRUE) <= ((ncol(vcffile)-4)/3)*0.8, ]
      vcffile$GTsum <- (rowSums(!is.na(vcffile))-(4 + sample_size + sample_size))
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
    AR <- subgenome_1
    AR1 <- subgenome_1[,c(1:4,ARstart:ARend)]
    subgenome_1 <- subgenome_1[,-c(ARstart:ARend)]
    write.table (subgenome_1, file=paste(pop,"_8x","_DP_GT.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    write.table (AR1, file=paste(pop,"_8x","_AR.txt",sep=""), row.names=F, quote = FALSE, sep = "\t")
    AR1 <- NULL
    vcffile <- NULL
    vcffile_DP <- NULL
    vcffile_GT <- NULL
    unlink(paste("*8x_rawSPLIT*",sep=""))
    
    
    AR$propHet <- (rowSums(AR == "0/0/0/0/0/0/0/1" | AR == "0/0/0/0/0/0/1/1" | AR == "0/0/0/0/0/1/1/1" | AR == "0/0/0/0/1/1/1/1" | AR == "0/0/0/1/1/1/1/1" | AR == "0/0/1/1/1/1/1/1" | AR == "0/1/1/1/1/1/1/1" | 
                             AR == "0|0|0|0|0|0|0|1" | AR == "0|0|0|0|0|0|1|1" | AR == "0|0|0|0|0|1|1|1" | AR == "0|0|0|0|1|1|1|1" | AR == "0|0|0|1|1|1|1|1" | AR == "0|0|1|1|1|1|1|1" | AR == "0|1|1|1|1|1|1|1", na.rm = TRUE)) / (GTincre - (rowSums(AR == "./././././.", na.rm = TRUE)))
    AR <- AR[,c(ARstart:ARend,ncol(AR))]
    AR <- subset(AR, AR$propHet != 0)
    AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
    AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
    suppressWarnings(AR <- na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
    AR <- AR[,1:2]; names(AR) <- c("propHet","Allele_Ratio")
    
    ARplot <- ggplot(AR,aes(x=propHet,y=Allele_Ratio))+
      geom_point(aes(x = propHet, y = Allele_Ratio, color = Allele_Ratio), size = 1, pch=19, alpha=0.1)+
      geom_density_2d(bins=50)+
      annotate("rect", xmin=0, xmax=1, ymin=-0.2, ymax=0.2, alpha=0.2, fill="tomato")+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      scale_colour_gradient2(low="darkorange3", mid="darkgoldenrod1", high ="cornflowerblue",
                             breaks=c(0.75,-0.75), limits=c(-1,1), 
                             labels=c("Minor Allele: Ref","Minor Allele: Alt"))+
      theme(legend.title = element_blank())+
      geom_hline(yintercept = 0.09, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.09, color="tomato", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 1, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = 0.6, color="grey20", size=0.5, linetype="dashed")+
      geom_hline(yintercept = 0.33, color="grey20", size=0.5, linetype="dashed")+
      geom_hline(yintercept = 0.14, color="grey20", size=0.5, linetype="dashed")+
      geom_hline(yintercept = -0.6, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.33, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -0.14, color="grey20", size=0.5, linetype="dashed")+ 
      geom_hline(yintercept = -1, color="grey20", size=0.5, linetype="dashed")+ 
      annotate("text", x=0.92, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.87, y=0.14, label="Heterozygote (0/0/0/0/0/0/0/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.87, y=0.33, label="Heterozygote (0/0/0/0/0/0/1/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.87, y=0.6, label="Heterozygote (0/0/0/0/0/1/1/1)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
      annotate("text", x=0.85, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
      annotate("text", x=0.87, y=-0.6, label="Heterozygote (0/0/0/0/1/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.87, y=-0.33, label="Heterozygote (0/0/1/1/1/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.87, y=-0.14, label="Heterozygote (0/1/1/1/1/1/1/1)", vjust=0.5, fontface="italic")+
      annotate("text", x=0.85, y=0.09, label="Allele Ratio filter threshold (< 0.09)", vjust=-0.5, fontface="italic")+
      annotate("text", x=0.85, y=-0.09, label="Allele Ratio filter threshold (> -0.09)", vjust=1.2, fontface="italic")+
      xlab("Proportion of Heterozygote per Locus (octaploid)") +
      ylab("Allele Read Depth Ratio per Genotype")
    ggsave(file=paste("raw8x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
    
  }
  vcf_to_DP_GT_8x()
}
