#!/usr/bin/env Rscript

callFile <- Sys.glob("*call_subset")
bedFile <- Sys.glob("*bed_subset.txt")

bedTable <- read.table(bedFile, sep = "\t")
callTable <- read.table(callFile, sep = ":")
df <- data.frame(geneId = bedTable$V4,
                 call = callTable$V8,
                 sampleId = c(rep(callFile, nrow(bedTable))))
write.table(df,
            file = paste0("etl_output_", callFile, ".txt"),
            row.names = F, sep = "\t",
            col.names = F,
            quote = F)
