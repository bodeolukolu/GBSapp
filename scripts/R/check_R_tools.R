#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
args
libdir <- args[1]
.libPaths( c( .libPaths(), libdir) )

if(!require("reshape2")) install.packages("reshape2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")
if(!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")
if(!require("CMplot")) install.packages("CMplot", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")
