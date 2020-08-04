library(tidyverse)


df <- read.table("region_intersect.bed")
colnames(df) <- c("chrom", "gene_s", "gene_e", "ensembl", "gene", "cnv_s", "cnv_e", "cnv_len", "cnv_call", "pt_id", "chrom1", "region_s", "region_e", "region")

largeTable <- df %>%
    select(chrom, ensembl, gene, cnv_call, pt_id, region)

focalTable <- df %>%
    select(chrom, ensembl, gene, cnv_s:pt_id)

geneTable <- df %>%
    select(chrom, ensembl, gene, cnv_call, pt_id)

write.table(largeTable,
    file = "large_table.csv",
    row.names = F,
    sep = ",",
    quote = F
)

write.table(focalTable,
    file = "focal_table.csv",
    row.names = F,
    sep = ",",
    quote = F
)

write.table(geneTable,
    file = "gene_table.csv",
    row.names = F,
    sep = ",",
    quote = F
)
