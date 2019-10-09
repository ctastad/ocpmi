#!/usr/bin/env Rscript

library(DNAcopy)

args <- commandArgs(TRUE)
sampleName <- gsub(pattern = "\\.sortedByCoord.out.called$", "", args)

cn <- read.table(args,header=T)
cna.object <- CNA(cbind(cn$normal_depth, cn$tumor_depth), cn$chrom, cn$chr_start, data.type="logratio", sampleid=sampleName)

# Smoothing object
cna.object.smooth <- smooth.CNA(cna.object)

#Segmentation at default parameters
cna.object.smooth.segd <- segment(cna.object.smooth, verbose=1)

#Plot whole studies

# pdf('rplot2.pdf')
# plot(cna.object.smooth.segd, plot.type="w")
# dev.off()

pdf(paste0(sampleName, ".pdf"))
plot(cna.object.smooth.segd, plot.type="w")
dev.off()
