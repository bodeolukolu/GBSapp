#!/usr/bin/env Rscript

####################################################################################################################
# set_variables
args <- commandArgs(trailingOnly = TRUE)
cat(args,"\n")
p1 <- NULL
p2 <- NULL
rd <- args[1]
rd <- as.numeric(as.character(rd))
libdir <- args[2]
dir.create("ploidy_aneuploidy_test")
dir.create("genotype_accuracy")

.libPaths( c( .libPaths(), libdir) )
library(ggplot2)


####################################################################################################################
####################################################################################################################
AD <- read.table ("snp_allele_depth.txt", header=T, sep="\t", comment.char="", check.names = FALSE)
snpno <- ncol(AD) * 0.01
AD$SNP <- paste(AD$CHROM,AD$POS,sep="_")
AD <- subset(AD, select=c(ncol(AD),1:(ncol(AD)-1)))
AD <- subset(AD, select=-c(2:10)); gc()
hold <- NULL
for (i in 2:ncol(AD)){
  transit <- AD[,c(1,i)]; transit$samples <- colnames(AD[i])
  colnames(transit) <- c("SNP","value","samples"); transit <- subset(transit, select=c("SNP","samples","value"))
  transit <- subset(transit, transit$value != "./." | transit$value != "./././." | transit$value != "./././././." | transit$value != "./././././././.")
  hold <- rbind(hold,transit); gc()
}
transit <- NULL; AD <- hold; hold <- NULL
AD$Genotype <- gsub(":.*", "", AD$value)
AD$AD <- gsub("^[^:]*:", "", AD$value)
AD$AD <- gsub(":.*", "", AD$AD)
AD$ref <- gsub(",.*", "", AD$AD)
AD$alt <- gsub(".*,", "", AD$AD)
AD[][AD[]=="."] <- NA
AD <- na.omit(AD)
AD$ref <- as.numeric(as.character(AD$ref))
AD$alt <- as.numeric(as.character(AD$alt))

AD <- subset(AD, select=-c(value,AD))
AD <- subset(AD, Genotype=="0/0" | Genotype=="0/1" | Genotype=="1/1" |
               Genotype=="0/0/0/0" | Genotype=="0/0/0/1" | Genotype=="0/0/1/1" | Genotype=="0/1/1/1" | Genotype=="1/1/1/1" | 
               Genotype=="0/0/0/0/0/0" | Genotype=="0/0/0/0/0/1" | Genotype=="0/0/0/0/1/1" | Genotype=="0/0/0/1/1/1" | 
               Genotype=="0/0/1/1/1/1" | Genotype=="0/1/1/1/1/1" | Genotype=="1/1/1/1/1/1" | 
               Genotype=="0/0/0/0/0/0/0/0" | Genotype=="0/0/0/0/0/0/0/1" | Genotype=="0/0/0/0/0/0/1/1" | Genotype=="0/0/0/0/0/1/1/1" | 
               Genotype=="0/0/0/0/1/1/1/1" | Genotype=="0/0/0/1/1/1/1/1" | Genotype=="0/0/1/1/1/1/1/1" | Genotype=="0/1/1/1/1/1/1/1" |
               Genotype=="1/1/1/1/1/1/1/1")

AD$chrom <- AD$SNP; AD$chrom <- gsub("^[^_]*_", "", AD$chrom)
AD$chrom <- gsub("_.*", "", AD$chrom)
AD$chrom <- gsub("[^0-9.-]", "", AD$chrom)
AD$chrom <- as.numeric(as.character(AD$chrom))
AD$pos <- AD$SNP; AD$pos <- gsub(".*_", "", AD$pos)

ADerror <- AD
ADerror$error <- ADerror$ref * ADerror$alt
ADerror0 <- subset(ADerror, error == 0)
ADerror0  <- subset(ADerror0, Genotype=="0/1" | Genotype=="0/0/0/1" | Genotype=="0/0/1/1" | Genotype=="0/1/1/1" | 
                      Genotype=="0/0/0/0/0/1" | Genotype=="0/0/0/0/1/1" | Genotype=="0/0/0/1/1/1" | Genotype=="0/0/1/1/1/1" | 
                      Genotype=="0/1/1/1/1/1" | Genotype=="0/0/0/0/0/0/0/1" | Genotype=="0/0/0/0/0/0/1/1" | 
                      Genotype=="0/0/0/0/0/1/1/1" | Genotype=="0/0/0/0/1/1/1/1" | Genotype=="0/0/0/1/1/1/1/1" | Genotype=="0/0/1/1/1/1/1/1" | 
                      Genotype=="0/1/1/1/1/1/1/1")
ADerror0 <- as.data.frame(table(ADerror0$SNP))


ADerror <- ADerror0
remove_id <- unique(as.vector(ADerror$Var1))
id <- as.vector(unique(AD$SNP))
keep_id <- setdiff(id,remove_id)
keep_id_anchored <- gsub("^[^_]*_","",keep_id)
keep_id <- as.data.frame(keep_id); names(keep_id) <- "SNP"
keep_id_anchored <- as.data.frame(keep_id_anchored); names(keep_id_anchored) <- "SNP"
keep <- rbind(keep_id, keep_id_anchored)
path <- getwd()
path <- paste(path,"/",sep="")
file.names <- dir(path,pattern ="maf.*_.*\\.txt$",full.names=TRUE)
for(i in 1:length(file.names)){
  all_content <- readLines(paste(file.names[i],sep=""))
  mafdiv <- read.table(textConnection(all_content), header = TRUE, check.names = FALSE)
  mafdiv <- merge(mafdiv, keep, by="SNP")
  write.table (mafdiv, file=paste(file.names[i]), row.names=F, quote = FALSE, sep = "\t")
}

AD <- merge(AD, keep_id, by="SNP")

ADp <- AD
AD[][AD[]=="0/0"] <- 0
AD[][AD[]=="0/1"] <- 1
AD[][AD[]=="1/1"] <- 2
AD[][AD[]=="./."] <- NA
AD[][AD[]=="0/0/0/0"] <- 0
AD[][AD[]=="0/0/0/1"] <- 1
AD[][AD[]=="0/0/1/1"] <- 2
AD[][AD[]=="0/1/1/1"] <- 3
AD[][AD[]=="1/1/1/1"] <- 4
AD[][AD[]=="./././."] <- NA
AD[][AD[]=="0/0/0/0/0/0"] <- 0
AD[][AD[]=="0/0/0/0/0/1"] <- 1
AD[][AD[]=="0/0/0/0/1/1"] <- 2
AD[][AD[]=="0/0/0/1/1/1"] <- 3
AD[][AD[]=="0/0/1/1/1/1"] <- 4
AD[][AD[]=="0/1/1/1/1/1"] <- 5
AD[][AD[]=="1/1/1/1/1/1"] <- 6
AD[][AD[]=="./././././."] <- NA
AD[][AD[]=="0/0/0/0/0/0/0/0"] <- 0
AD[][AD[]=="0/0/0/0/0/0/0/1"] <- 1
AD[][AD[]=="0/0/0/0/0/0/1/1"] <- 2
AD[][AD[]=="0/0/0/0/0/1/1/1"] <- 3
AD[][AD[]=="0/0/0/0/1/1/1/1"] <- 4
AD[][AD[]=="0/0/0/1/1/1/1/1"] <- 5
AD[][AD[]=="0/0/1/1/1/1/1/1"] <- 6
AD[][AD[]=="0/1/1/1/1/1/1/1"] <- 7
AD[][AD[]=="1/1/1/1/1/1/1/1"] <- 8
AD[][AD[]=="./././././././."] <- NA
AD <- na.omit(AD)
maxdose <- max(AD$Genotype)
ploidy <- ADp[1,3]; ploidy <- unlist(strsplit(ploidy,"/")); ploidy <- length(ploidy) + 4


####################################################################################################################
####################################################################################################################

ADp[][ADp[]=="0/0"] <- NA
ADp[][ADp[]=="0/1"] <- 1
ADp[][ADp[]=="1/1"] <- NA
ADp[][ADp[]=="./."] <- NA
ADp[][ADp[]=="0/0/0/0"] <- NA
ADp[][ADp[]=="0/0/0/1"] <- 1
ADp[][ADp[]=="0/0/1/1"] <- 2
ADp[][ADp[]=="0/1/1/1"] <- 3
ADp[][ADp[]=="1/1/1/1"] <- NA
ADp[][ADp[]=="./././."] <- NA
ADp[][ADp[]=="0/0/0/0/0/0"] <- NA
ADp[][ADp[]=="0/0/0/0/0/1"] <- 1
ADp[][ADp[]=="0/0/0/0/1/1"] <- 2
ADp[][ADp[]=="0/0/0/1/1/1"] <- 3
ADp[][ADp[]=="0/0/1/1/1/1"] <- 4
ADp[][ADp[]=="0/1/1/1/1/1"] <- 5
ADp[][ADp[]=="1/1/1/1/1/1"] <- NA
ADp[][ADp[]=="./././././."] <- NA
ADp[][ADp[]=="0/0/0/0/0/0/0/0"] <- NA
ADp[][ADp[]=="0/0/0/0/0/0/0/1"] <- 1
ADp[][ADp[]=="0/0/0/0/0/0/1/1"] <- 2
ADp[][ADp[]=="0/0/0/0/0/1/1/1"] <- 3
ADp[][ADp[]=="0/0/0/0/1/1/1/1"] <- 4
ADp[][ADp[]=="0/0/0/1/1/1/1/1"] <- 5
ADp[][ADp[]=="0/0/1/1/1/1/1/1"] <- 6
ADp[][ADp[]=="0/1/1/1/1/1/1/1"] <- 7
ADp[][ADp[]=="1/1/1/1/1/1/1/1"] <- NA
ADp[][ADp[]=="./././././././."] <- NA
ADp <- na.omit(ADp)

chromtrim <- unique(subset(ADp, select=c("chrom","pos")))
chromtrim <- aggregate(pos ~ chrom, data = chromtrim, max)
chromtrim$pos <- as.numeric(as.character(chromtrim$pos))
chromtrim_sum <- sum(chromtrim[,2])
if (nrow(chromtrim) > 50){
  chromtrim[,2][chromtrim[,2]<1000000] <- NA
  chromtrim <- na.omit(chromtrim)
  chromtrim <- chromtrim[order(-(chromtrim$pos)),]
  if (nrow(chromtrim) > 50) {
    chromtrim <- chromtrim[1:50,]
  }
  chromtrim_sum <- sum(chromtrim[,2])
  chromtrim <- chromtrim[,1]
  catADp <- NULL
  for (i in 1:length(chromtrim)){
    j <- paste("^",chromtrim[i],"$",sep="")
    subADp <- ADp[grep(j, ADp$chrom),]
    catADp <- rbind(catADp, subADp)
  }
  ADp <- catADp
}
chromtrim_sum <- sum(chromtrim[,2])

ADp$ref <- as.numeric(as.character(ADp$ref))
ADp$alt <- as.numeric(as.character(ADp$alt))
ADp$chrom <- as.numeric(as.character(ADp$chrom))
ADp$pos <- as.numeric(as.character(ADp$pos))
ADp$DP <- ADp$ref+ADp$alt
ADp <- subset(ADp, ADp$DP >= rd)
quantile5 <- quantile(ADp$DP, probs = c(0.05), na.rm= TRUE)
quantile95 <- quantile(ADp$DP, probs = c(0.95), na.rm= TRUE)
ADp <- subset(ADp, ADp$DP > quantile5 & ADp$DP < quantile95)
quantile5 <- quantile(ADp$alt, probs = c(0.05), na.rm= TRUE)
quantile95 <- quantile(ADp$alt, probs = c(0.95), na.rm= TRUE)
ADp <- subset(ADp, ADp$alt > quantile5 & ADp$alt < quantile95)
quantile5 <- quantile(ADp$ref, probs = c(0.05), na.rm= TRUE)
quantile95 <- quantile(ADp$ref, probs = c(0.95), na.rm= TRUE)
ADp <- subset(ADp, ADp$ref > quantile5 & ADp$ref < quantile95)

samples <- unique(subset(ADp, select=c(samples))); names(samples)[1] <- "sample"
for (k in 1:nrow(samples)) {
  sample <- samples[k,]
  rdpos <- subset(ADp, samples==sample)
  rdpos <- subset(rdpos, select=c(chrom,pos,ref,alt))
  rdpos$ref <- as.numeric(as.character(rdpos$ref))
  rdpos$alt <- as.numeric(as.character(rdpos$alt))
  rdpos$DP <- as.numeric(rdpos$ref + rdpos$alt)
  chromlist <- unique(rdpos$chrom)
  sumaggr <- NULL
  for (m in 1:length(chromlist)) {
    aggr <- subset(rdpos, chrom == chromlist[m]); aggr <- mean(aggr$DP); aggr <- t(as.data.frame(c(chromlist[m],aggr)))
    sumaggr <- as.data.frame(rbind(sumaggr, aggr[1:2]))
  }
  colnames(sumaggr) <- c("chrom","mean")
  rdpos <- merge(rdpos, sumaggr, by="chrom")
  rdpos$norm <- rdpos$mean / mean(rdpos$DP)
  sumaggr <- NULL
  rdpos$min <- do.call(pmin, rdpos[,3:4]); rdpos$max <- do.call(pmax, rdpos[,3:4])
  rdpos$ratio <- (rdpos$DP / rdpos$min) * rdpos$norm
  rdpos0 <- subset(rdpos, (rdpos$DP / rdpos$min) >= 2.5)
  rdpos1 <- subset(rdpos, (rdpos$DP / rdpos$min) < 2.5)
  rdpos1a <- rdpos1; rdpos1a$ratio <- ((rdpos1$DP / rdpos1$min) * rdpos1$norm)
  rdpos1b <- rdpos1; rdpos1b$ratio <- ((rdpos1$DP / rdpos1$max) * rdpos1$norm)
  rdpos1 <- rbind(rdpos1a, rdpos1b)
  rdpos <- rbind(rdpos0, rdpos1)
  rdpos <- na.omit(rdpos); rdpos <- subset(rdpos, ratio != "Inf" | ratio != 0)
  rdpos$pos <- as.numeric(as.character(rdpos$pos))
  rdpos <- rdpos[order(rdpos$chrom, rdpos$pos),]
  rdpos <- subset(rdpos, chrom != 0)
  quantile5 <- quantile(rdpos$ratio, probs = c(0.05), na.rm= TRUE)
  quantile95 <- quantile(rdpos$ratio, probs = c(0.95), na.rm= TRUE)
  rdpos <- subset(rdpos, rdpos$ratio > quantile5 & rdpos$ratio < quantile95)
  
  if (nrow(rdpos) >= 100 && chromtrim_sum >= 50000000) {
    genomeploidy <- ggplot(rdpos, aes(x=interaction(chrom,pos), y=ratio, group=chrom)) +
      geom_point(aes(colour=factor(chrom)), pch=19, size=1.5, alpha=0.05) +
      geom_density_2d(aes(colour=factor(chrom)), alpha=1.0) +
      theme(text = element_text(size=10)) + 
      facet_grid(. ~ chrom, scales="free_x", space="free_x") +
      scale_color_manual(values=c('cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange',
                                  'cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange','cornflowerblue','orange')) +
      xlab("\nChromosomes and Physical Map Position (bp)")+ ylab("Estimated Ploidy") +
      theme(panel.spacing = unit(0.1, "lines"), panel.grid.major = element_line(colour ="grey50"), panel.background = element_blank(),
            panel.border = element_blank(), axis.text.x=element_blank(),strip.text.x = element_text(size=10,color="grey20"), 
            strip.background = element_blank(), strip.text.y = element_text(size=10,color="black"), axis.line=element_line(colour="grey50")) +
      scale_x_discrete(breaks = NULL) +
      coord_cartesian(ylim = c(0,ploidy), expand=FALSE) +
      theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
    ggsave(file=paste("./ploidy_aneuploidy_test/",sample,"_rd",rd,"_ploidytest.tiff",sep=""), plot=genomeploidy, width=12, height=3, units=("in"), dpi=240, compression = "lzw")
    gc()
    Sys.sleep(1)
  }
}

if (nrow(rdpos) < 100 ) {
  print("Not enough markers to estimate ploidy level and aneuploidy")
}

####################################################################################################################
####################################################################################################################
AD$ref <- as.numeric(as.character(AD$ref))
AD$alt <- as.numeric(as.character(AD$alt))
AD$chrom <- as.numeric(as.character(AD$chrom))
AD$pos <- as.numeric(as.character(AD$pos))
snplist <- unique(subset(AD, select=c(1)))

if (maxdose==2) {
  AD <- subset(AD, AD$ref+AD$alt >= rd)
  for (i in 1:nrow(snplist)) {
    genoaccu <- subset(AD, SNP==snplist[i,])
    h1 <- subset(genoaccu, samples == p1); h1 <- c(h1$ref,h1$alt)
    h2 <- subset(genoaccu, samples == p2); h2 <- c(h2$ref,h2$alt)
    maxY <- max(genoaccu$alt, na.rm=T); maxX <- max(genoaccu$ref, na.rm=T)
    max <- max(maxY,maxX)
    marker <- gsub(".*::", "", snplist[i,])
    plot <- ggplot(genoaccu, aes(x=ref, y=alt), group="Genotype") +
      annotate("point", x=h1[1], y=h1[2], size=3, color="grey50", pch=1) + coord_equal() +
      geom_text(label=paste("\np1: ",p1,"\n",sep=""), x=h1[1], y=h1[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      annotate("point", x=h2[1], y=h2[2], size=3, color="grey50", pch="1") + 
      geom_text(label=paste("\np2: ",p2,"\n",sep=""), x=h2[1], y=h2[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      geom_point(aes(colour=Genotype), size=2.5, pch=19, alpha=0.5) +
      scale_colour_manual(values=c('black','tomato','cornflowerblue','orange')) +
      geom_hline(aes(yintercept = 0), color="black", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=0, label="0", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 1/1, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(max*1.1), label="1", color="black", size=3.5) +
      geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=0, y=(max*1.1), label="2", color="black", size=3.5) +
      xlim(0,(max*1.1)) + ylim(0,(max*1.1)) + theme(legend.title = element_text(colour="black", size=9, face="italic")) +
      theme(legend.text = element_text(colour="black", size=8, face="italic")) +
      guides(color = guide_legend(override.aes = list(size=4))) +
      theme(strip.text.x = element_text(size=40,color="black"), strip.text.y = element_text(size=40,color="black"), 
            axis.line=element_line(colour="white")) + theme(plot.title = element_text(hjust = 0.0)) +
      labs(title= paste("Genotype call accuracy (",marker,")\nExpected (dashed lines) vs. Estimated (colored dots) dose",sep="")) +
      theme(plot.title = element_text(color = "grey20", size =10, face="italic")) +
      xlab("Reference Allele Read Depth")+ ylab("Alternate Allele Read Depth")
    ggsave(file=paste("./genotype_accuracy/",marker,"_rd",rd,".tiff",sep=""), plot=plot, width=6, height=6, units=("in"), dpi=120, compression = "lzw")
    gc()
    Sys.sleep(1)
  }
}
if (maxdose==4) {
  AD <- subset(AD, AD$ref+AD$alt >= rd)
  for (i in 1:nrow(snplist)) {
    genoaccu <- subset(AD, SNP==snplist[i,])
    h1 <- subset(genoaccu, samples == p1); h1 <- c(h1$ref,h1$alt)
    h2 <- subset(genoaccu, samples == p2); h2 <- c(h2$ref,h2$alt)
    maxY <- max(genoaccu$alt, na.rm=T); maxX <- max(genoaccu$ref, na.rm=T)
    max <- max(maxY,maxX)
    marker <- gsub(".*::", "", snplist[i,])
    plot <- ggplot(genoaccu, aes(x=ref, y=alt), group="Genotype") +
      annotate("point", x=h1[1], y=h1[2], size=10, color="grey50", pch="*") + coord_equal() +
      geom_text(label=paste("\np1: ",p1,"\n",sep=""), x=h1[1], y=h1[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      annotate("point", x=h2[1], y=h2[2], size=10, color="grey50", pch="*") + 
      geom_text(label=paste("\np2: ",p2,"\n",sep=""), x=h2[1], y=h2[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      geom_point(aes(colour=Genotype), size=2.5, pch=19, alpha=0.5) +
      scale_colour_manual(values=c('black','tomato','cornflowerblue','orange','darkolivegreen3','mediumorchid')) +
      geom_hline(aes(yintercept = 0), color="black", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=0, label="0", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 1/3, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(1/3)*(max*1.1), label="1", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 2/2, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(max*1.1), label="2", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 3/1, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(1/3)*(max*1.1), y=(max*1.1), label="3", color="black", size=3.5) +
      geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=0, y=(max*1.1), label="4", color="black", size=3.5) +
      xlim(0,(max*1.1)) + ylim(0,(max*1.1)) + theme(legend.title = element_text(colour="black", size=9, face="italic")) +
      theme(legend.text = element_text(colour="black", size=8, face="italic")) +
      guides(color = guide_legend(override.aes = list(size=4))) +
      theme(strip.text.x = element_text(size=40,color="black"), strip.text.y = element_text(size=40,color="black"), 
            axis.line=element_line(colour="white")) + theme(plot.title = element_text(hjust = 0.0)) +
      labs(title= paste("Genotype call accuracy (",marker,")\nExpected (dashed lines) vs. Estimated (colored dots) dose",sep="")) +
      theme(plot.title = element_text(color = "grey20", size =10, face="italic")) +
      xlab("Reference Allele Read Depth")+ ylab("Alternate Allele Read Depth")
    ggsave(file=paste("./genotype_accuracy/",marker,"_rd",rd,".tiff",sep=""), plot=plot, width=6, height=6, units=("in"), dpi=120, compression = "lzw")
    gc()
    Sys.sleep(1)
  }
}
if (maxdose==6) {
  AD <- subset(AD, AD$ref+AD$alt >= rd)
  for (i in 1:nrow(snplist)) {
    genoaccu <- subset(AD, SNP==snplist[i,])
    h1 <- subset(genoaccu, samples == p1); h1 <- c(h1$ref,h1$alt)
    h2 <- subset(genoaccu, samples == p2); h2 <- c(h2$ref,h2$alt)
    maxY <- max(genoaccu$alt, na.rm=T); maxX <- max(genoaccu$ref, na.rm=T)
    max <- max(maxY,maxX)
    marker <- gsub(".*::", "", snplist[i,])
    plot <- ggplot(genoaccu, aes(x=ref, y=alt), group="Genotype") +
      annotate("point", x=h1[1], y=h1[2], size=10, color="grey50", pch="*") + coord_equal() +
      geom_text(label=paste("\np1: ",p1,"\n",sep=""), x=h1[1], y=h1[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      annotate("point", x=h2[1], y=h2[2], size=10, color="grey50", pch="*") + 
      geom_text(label=paste("\np2: ",p2,"\n",sep=""), x=h2[1], y=h2[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      geom_point(aes(colour=Genotype), size=2.5, pch=19, alpha=0.5) +
      scale_colour_manual(values=c('black','tomato','cornflowerblue','orange','darkolivegreen3','mediumorchid','dodgerblue4','gold4')) +
      geom_hline(aes(yintercept = 0), color="black", linetype="dashed", size=0.5) +
      annotate("text", x=max*1.1, y=0, label="0", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 1/5, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(1/5)*(max*1.1), label="1", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 2/4, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(2/4)*(max*1.1), label="2", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 3/3, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(max*1.1), label="3", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 4/2, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(2/4)*(max*1.1), y=(max*1.1), label="4", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 5/1, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(1/5)*(max*1.1), y=(max*1.1), label="5", color="black", size=3.5) +
      geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=0, y=(max*1.1), label="6", color="black", size=3.5) +
      xlim(0,max*1.1) + ylim(0,max*1.1) + theme(legend.title = element_text(colour="black", size=9, face="italic")) +
      theme(legend.text = element_text(colour="black", size=8, face="italic")) +
      guides(color = guide_legend(override.aes = list(size=4))) +
      theme(strip.text.x = element_text(size=40,color="black"), strip.text.y = element_text(size=40,color="black"), 
            axis.line=element_line(colour="white")) + theme(plot.title = element_text(hjust = 0.0)) +
      labs(title= paste("Genotype call accuracy (",marker,")\nExpected (dashed lines) vs. Estimated (colored dots) dose",sep="")) +
      theme(plot.title = element_text(color = "grey20", size =10, face="italic")) +
      xlab("Reference Allele Read Depth")+ ylab("Alternate Allele Read Depth")
    ggsave(file=paste("./genotype_accuracy/",marker,"_rd",rd,".tiff",sep=""), plot=plot, width=6, height=6, units=("in"), dpi=120, compression = "lzw")
    gc()
    Sys.sleep(1)
  }
}
if (maxdose==8) {
  AD <- subset(AD, AD$ref+AD$alt >= rd)
  for (i in 1:nrow(snplist)) {
    genoaccu <- subset(AD, SNP==snplist[i,])
    h1 <- subset(genoaccu, samples == p1); h1 <- c(h1$ref,h1$alt)
    h2 <- subset(genoaccu, samples == p2); h2 <- c(h2$ref,h2$alt)
    maxY <- max(genoaccu$alt, na.rm=T); maxX <- max(genoaccu$ref, na.rm=T)
    max <- max(maxY,maxX)
    marker <- gsub(".*::", "", snplist[i,])
    plot <- ggplot(genoaccu, aes(x=ref, y=alt), group="Genotype") +
      annotate("point", x=h1[1], y=h1[2], size=10, color="grey50", pch="*") + coord_equal() +
      geom_text(label=paste("\np1: ",p1,"\n",sep=""), x=h1[1], y=h1[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      annotate("point", x=h2[1], y=h2[2], size=10, color="grey50", pch="*") + 
      geom_text(label=paste("\np2: ",p2,"\n",sep=""), x=h2[1], y=h2[2], vjust ="inward", hjust="inward", color="grey50", size=4, fontface='italic') +
      geom_point(aes(colour=Genotype), size=2.5, pch=19, alpha=0.5) +
      scale_colour_manual(values=c('black','tomato','cornflowerblue','orange','darkolivegreen3','mediumorchid','dodgerblue4','gold4','green3','bisque4')) +
      geom_hline(aes(yintercept = 0), color="black", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=0, label="0", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 1/7, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(1/7)*(max*1.1), label="1", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 2/6, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(2/6)*(max*1.1), label="2", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 3/5, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(3/5)*(max*1.1), label="3", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 4/4, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(max*1.1), y=(max*1.1), label="4", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 5/3, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(3/5)*(max*1.1), y=(max*1.1), label="5", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 6/2, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(2/6)*(max*1.1), y=(max*1.1), label="6", color="black", size=3.5) +
      geom_abline(intercept = 0, slope = 7/1, color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=(1/7)*(max*1.1), y=(max*1.1), label="7", color="black", size=3.5) +
      geom_vline(aes(xintercept = 0), color="grey50", linetype="dashed", size=0.5) +
      annotate("text", x=0, y=(max*1.1), label="8", color="black", size=3.5) +
      xlim(0,(max*1.1)) + ylim(0,(max*1.1)) + theme(legend.title = element_text(colour="black", size=9, face="italic")) +
      theme(legend.text = element_text(colour="black", size=8, face="italic")) +
      guides(color = guide_legend(override.aes = list(size=4))) +
      theme(strip.text.x = element_text(size=40,color="black"), strip.text.y = element_text(size=40,color="black"), 
            axis.line=element_line(colour="white")) + theme(plot.title = element_text(hjust = 0.0)) +
      labs(title= paste("Genotype call accuracy (",marker,")\nExpected (dashed lines) vs. Estimated (colored dots) dose",sep="")) +
      theme(plot.title = element_text(color = "grey20", size =10, face="italic")) +
      xlab("Reference Allele Read Depth")+ ylab("Alternate Allele Read Depth")
    ggsave(file=paste("./genotype_accuracy/",marker,"_rd",rd,".tiff",sep=""), plot=plot, width=6, height=6, units=("in"), dpi=120, compression = "lzw")
    gc()
    Sys.sleep(1)
  }
}


