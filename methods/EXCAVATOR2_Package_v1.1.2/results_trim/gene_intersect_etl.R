#!/usr/bin/env Rscript

################################################################################
#
# This is the helper script for the gene_association process. It takes the
# output bedfile subsets from bedtools intersect and assembles them into a
# long form dataframe.
#
################################################################################

library(tidyverse)


# draw file names from local dir
callFile <- Sys.glob("*_intersect")
callTable <- read.table(callFile, sep = "\t")

# transform
callTable <- callTable %>%
    separate(into = c("GT", "CN", "CNF", "FCL", "FCP"), V15, sep = ":") %>%
    separate(into = c("IMP", "TYPE", "END", "SVLEN"), V13, sep = ";") %>%
    select(V1:V5, V7, END, SVLEN, FCL) %>%
    mutate(PT = callFile)

# write table to disk
write.table(callTable,
  file = paste0("etl_output_", callFile, ".csv"),
  row.names = F,
  sep = ",",
  col.names = F,
  quote = F
)
