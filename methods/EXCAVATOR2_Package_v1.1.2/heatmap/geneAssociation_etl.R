#!/usr/bin/env Rscript

################################################################################
#
# This is the helper script for the gene_association process. It takes the
# output bedfile subsets from bedtools intersect and assembles them into a
# workable dataframe.
#
################################################################################

# Draw file names from local dir
callFile <- Sys.glob("*call_subset")
bedFile <- Sys.glob("*bed_subset.txt")

bedTable <- read.table(bedFile, sep = "\t")
callTable <- read.table(callFile, sep = ":")
df <- data.frame(geneId = bedTable$V4,
                 geneName = bedTable$V5,
                 call = callTable$V8,
                 sampleId = c(rep(callFile, nrow(bedTable))))
write.table(df,
            file = paste0("etl_output_", callFile, ".txt"),
            row.names = F, sep = "\t",
            col.names = F,
            quote = F)
