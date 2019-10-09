#!/usr/bin/env Rscript

################################################################################
#
# This script is used to create the graphical output of the geneAssociation
# function. The first argument should be the experiment name, and the second
# should be the figure title.
#
################################################################################

library(ggplot2)

args <- commandArgs(TRUE)
dir <- paste0("~/ocpmi/results/excavator2/", args[1], "/Results")
setwd(dir)

geneAssociationTable <- read.table("geneAssociationTable", header = F, sep = "\t")
colnames(geneAssociationTable) <- c("seqID", "call", "sampleId")

g <- ggplot(geneAssociationTable) +
    geom_tile(aes(sampleId, seqID, fill = call)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(title = args[2]) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_blank())

ggsave("plot.png", plot = g)
