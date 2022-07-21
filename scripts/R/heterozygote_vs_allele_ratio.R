#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
ploidy <- args[3]
subgenome <- args[4]
libdir <- args[5]
.libPaths( c( .libPaths(), libdir) )
library(ggplot2)


GT <- read.table(args[1], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
GT <- subset(GT, select=-c(SNP))
if (colnames(GT[ncol(GT)]) == "pvalue"){GT <- subset(GT, select=-c(pvalue))}
AR <- read.table(args[2], header=T, sep="\t", check.names=FALSE,stringsAsFactors=FALSE)
if (subgenome == 1) { AR$CHROM <-  gsub(".*_", "", AR$CHROM) }
GT_AR <- merge(GT, AR, by=c("CHROM","POS","REF","ALT"), all.x=TRUE)
keepAR <- subset(GT_AR, select=-c(CHROM,POS,REF,ALT))
keepAR <- keepAR[, -grep(pattern = "_AR$", colnames(keepAR))]
keepAR <- paste(colnames(keepAR[]),"_AR",sep="")
keepAR <- c(keepAR,"propHet")
AR <- GT_AR; GTlen <- length(keepAR)-1


if (ploidy == "2x"){
  AR$propHet <- (rowSums(AR == "1", na.rm = TRUE)) / (GTlen - (apply(AR[,1:GTlen], 1, function(x) sum(is.na(x)))))
  AR <- AR[,c(keepAR)]
    AR <- subset(AR, AR$propHet != 0)
  AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
  AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
  AR <- suppressWarnings(na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
  AR <- AR[,1:2]; names(AR) <- c("propHet","Allele_Ratio")
  AR2 <- AR
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
    annotate("text", x=0.94, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.87, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
    annotate("text", x=0.87, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
    annotate("text", x=0.88, y=0.2, label="Allele Ratio threshold (> 0.2): 1/1", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.88, y=-0.2, label="Allele Ratio threshold (< -0.2): 0/0", vjust=1.2, fontface="italic")+
    xlab("Proportion of Heterozygote per Locus (diploid)") +
    ylab("Allele Read Depth Ratio per Genotype")
  ggsave(file=paste("./visualizations/","filtered2x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
}

if (ploidy == "4x"){
  AR$propHet <- (rowSums(AR == "1" | AR == "2" | AR == "3", na.rm = TRUE)) / (GTlen - (apply(AR[,1:GTlen], 1, function(x) sum(is.na(x)))))
  AR <- AR[,c(keepAR)]
  AR <- subset(AR, AR$propHet != 0)
  AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
  AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
  AR <- suppressWarnings(na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
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
    annotate("text", x=0.94, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.91, y=0.33, label="Heterozygote (0/0/0/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.87, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
    annotate("text", x=0.87, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
    annotate("text", x=0.91, y=-0.33, label="Heterozygote (0/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.87, y=0.17, label="Allele Ratio threshold (> 0.17): 1/1/1/1", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.87, y=-0.17, label="Allele Ratio threshold (< -0.17): 0/0/0/0", vjust=1.2, fontface="italic")+
    xlab("Proportion of Heterozygote per Locus (tetraploid)") +
    ylab("Allele Read Depth Ratio per Genotype")
  ggsave(file=paste("./visualizations/","filtered4x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
}

if (ploidy == "6x"){
  AR$propHet <- (rowSums(AR == "1" | AR == "2" | AR == "3" | AR == "4" | AR == "5", na.rm = TRUE)) / (GTlen - (apply(AR[,1:GTlen], 1, function(x) sum(is.na(x)))))
  AR <- AR[,c(keepAR)]
  AR <- subset(AR, AR$propHet != 0)
  AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
  AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
  AR <- suppressWarnings(na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
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
    annotate("text", x=0.94, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.9, y=0.2, label="Heterozygote (0/0/0/0/0/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.9, y=0.5, label="Heterozygote (0/0/0/0/1/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.87, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
    annotate("text", x=0.87, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
    annotate("text", x=0.9, y=-0.5, label="Heterozygote (0/0/1/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.9, y=-0.2, label="Heterozygote (0/1/1/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.86, y=0.14, label="Allele Ratio threshold (> 0.14): 1/1/1/1/1/1", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.86, y=-0.14, label="Allele Ratio threshold (< -0.14): 0/0/0/0/0/0", vjust=1.2, fontface="italic")+
    xlab("Proportion of Heterozygote per Locus (hexaploid)") +
    ylab("Allele Read Depth Ratio per Genotype")
  ggsave(file=paste("./visualizations/","filtered6x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
}

if (ploidy == "8x"){  
  AR$propHet <- (rowSums(AR == "1" | AR == "2" | AR == "3" | AR == "4" | AR == "5" | AR == "6" | AR == "7", na.rm = TRUE)) / (GTlen - (apply(AR[,1:GTlen], 1, function(x) sum(is.na(x)))))
  AR <- AR[,c(keepAR)]
  AR[,1:(ncol(AR)-1)][AR[,1:(ncol(AR)-1)] == 0] <- NA
  AR <- AR[,c(which(colnames(AR)=="propHet"),which(colnames(AR)!="propHet"))]
  AR <- suppressWarnings(na.omit(cbind(AR[1], stack(AR[-1]), row.names = NULL)))
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
    annotate("text", x=0.94, y=0, label="Homozygote", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.89, y=0.14, label="Heterozygote (0/0/0/0/0/0/0/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.89, y=0.33, label="Heterozygote (0/0/0/0/0/0/1/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.89, y=0.6, label="Heterozygote (0/0/0/0/0/1/1/1)", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.87, y=1, label="Heterozygote (balanced allele ratio)", vjust=1, fontface="italic")+
    annotate("text", x=0.87, y=-1, label="Heterozygote (balanced allele ratio)", vjust=-1, fontface="italic")+
    annotate("text", x=0.89, y=-0.6, label="Heterozygote (0/0/0/0/1/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.89, y=-0.33, label="Heterozygote (0/0/1/1/1/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.89, y=-0.14, label="Heterozygote (0/1/1/1/1/1/1/1)", vjust=0.5, fontface="italic")+
    annotate("text", x=0.85, y=0.09, label="Allele Ratio threshold (0.09): 1/1/1/1/1/1/1/1", vjust=-0.5, fontface="italic")+
    annotate("text", x=0.85, y=-0.09, label="Allele Ratio threshold (0.09): 0/0/0/0/0/0/0/0", vjust=1.2, fontface="italic")+
    xlab("Proportion of Heterozygote per Locus (octaploid)") +
    ylab("Allele Read Depth Ratio per Genotype")
  ggsave(file=paste("./visualizations/","filtered8x_Allele_Ratio_Heterozygosity_plot",".tiff",sep=""), plot=ARplot, width=12, height=4, units=("in"), dpi=300, compression = "lzw")
}
