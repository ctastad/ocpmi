#!/usr/bin/env Rscript

################################################################################
#
# This scriopt creates a clustered graphical heatmap from the same dataset as
# geneAssociationGraph.R using the pheatmap package.
#
################################################################################

# Pull arguments passed from original command
args <- commandArgs(TRUE)
dir <- paste0("~/ocpmi/results/excavator2/", args[1], "/Results")
setwd(dir)

library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(scales)

geneAssociationTable <- read.table("geneAssociationTable",
                                   header = F,
                                   sep = "\t")
message("created OG table")
colnames(geneAssociationTable) <- c("seqId", "geneName", "call", "sampleId")
tableWide_lab <- dcast(geneAssociationTable,
                   seqId + geneName ~ sampleId,
                   value.var = "call",
                   fun = median,
                   fill = 0)
message("reshaped table to wide format")

tableWide <- tableWide_lab[,-2]
rownames(tableWide) <- tableWide[,1]
tableWide[,1] <- NULL
mat_unmelt <- as.matrix(tableWide)
message("created matrix")

# Call for column colorbar by patient subtyping
subtypes <- read.csv("~/ocpmi/reference/wes/patient_subtyping.csv")
rownames(subtypes) <- subtypes$sample
oubtypes[1] <- NULL

chromosomes <- read.table("~/ocpmi/data/gene_ref/tim/master_gene_list.bed")
chromosomes <- data.frame(chromosomes[4], chromosomes[1])
rownames(chromosomes) <- chromosomes$V4
chromosomes[1] <- NULL


message("working on heatmap")
pheatmap(mat_unmelt,
         show_rownames = F,
         #color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu"))),
         color = brewer_pal(palette = "RdBu", direction = -1)(5),
         annotation_col = subtypes,
         annotation_row = chromosomes,
         filename = "clustered_plot.pdf",
         main = args[2])
