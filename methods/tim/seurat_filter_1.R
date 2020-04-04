# First Seurat analysis to determine min & max for
#   genes per cell, UMI per cell, % mito genes per cell, and 
#   min cells per gene
# Author: Tim Starr
# Last updated: 11/15/18 using Seurat v 3.0.0.9000

# Set the working directory
setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Clear the workspace global environment
rm(list=ls())
# Load required libraries
library(Seurat)
library(dplyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(rgl)
library(reshape2)

# Save the version of Seurat
seurat_version <- as.character(packageVersion("Seurat"))

#######################################################################
# Load data
patient.data <- Read10X("cell_ranger_input_files/")
write.csv(patient.data@Dim, "output/sparse_matrix_dimensions.csv")
########################################################################
# CreateSeuratObject

# CreateSeuratObject parameters
assay_parameter <- "RNA"
min_cells_per_gene <- 1 # Include features/genes detected in this many cells
min_genes_per_cell <- 1 # Include cells with at least this many features/genes
############################
# Create Seurat object 
#   and print out cell violin plots for gene_count, UMI_count and % Mito
patient <- CreateSeuratObject(patient.data, project = "Patient", 
                              assay = assay_parameter, min.cells = min_cells_per_gene, 
                              min.features = min_genes_per_cell, names.field = 1,
                              names.delim = "-", meta.data = NULL)

vln_genes <- VlnPlot(patient, features = 'nFeature_RNA', 
      cols = NULL, pt.size = 1, idents = NULL,
      sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL,
      adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE,
      ncol = NULL, combine = TRUE, slot = "data")

png('output/violin_genes_all.png', width = 360, height = 360)
vln_genes + 
  theme(panel.background = element_rect(fill = "lightgray"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  labs(title="Genes per cell", 
       y="# of genes", x="single cells") +
  guides(fill=FALSE)
dev.off() 

vln_UMI <- VlnPlot(patient, features = 'nCount_RNA', 
                   cols = NULL, pt.size = 1, idents = NULL,
                   sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL,
                   adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE,
                   ncol = NULL, combine = TRUE, slot = "data")

png('output/violin_UMI_all.png', width = 360, height = 360)
vln_UMI + 
  theme(panel.background = element_rect(fill = "lightgray"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  labs(title="UMI total per cell", 
       y="# of UMI", x="single cells") +
  guides(fill=FALSE)
dev.off()

# Create a metadata entry containing % mitochondrial genes
# Select the appropriate pattern for matching mito gene names
# mito.genes <- grep(pattern = "^GRch38_MT-", x = rownames(x = patient@data), value = TRUE)
mito.genes <- grep(pattern = "^MT-", x = rownames(patient@assays$RNA), value = TRUE)
percent.mito <- Matrix::colSums(patient@assays$RNA@counts[mito.genes, ])/Matrix::colSums(patient@assays$RNA@counts)
patient <- AddMetaData(object = patient, metadata = percent.mito, col.name = "percent.mito")

vln_mito <- VlnPlot(patient, features = "percent.mito")
png('output/violin_mito_all.png', width = 360, height = 360)
vln_mito + 
  theme(panel.background = element_rect(fill = "lightgray"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  labs(title="% Mitochondrial genes per cell", 
       y="% Mitochondrial genes", x="single cells") +
  guides(fill=FALSE)
dev.off()

# Calculate cell count per gene: For each gene calculate the number of cells with UMI > 0
cells.per.gene <- rowSums(patient@assays$RNA@data !=0)
sorted.cells.per.gene <- data.frame(sort(cells.per.gene))

# Place the cells.per.gene in the meta.features slot
# metadata.cells.per.gene <- data.frame(cells.per.gene)
# patient@assays$RNA@meta.features <- cbind(patient@assays$RNA@meta.features, metadata.cells.per.gene)

# Write the cells.per.gene data to an output file
write.csv(cells.per.gene, "output/cells_per_gene_all.csv")

# Graph the entire vector of cells.per.gene
png('output/cells_per_gene_all.png', width = 720, height = 360)
plot(sort(cells.per.gene), xlab = 'Genes in rank order',
     ylab = 'Cells per gene', main = 'Cells per gene')
dev.off()

# Graph the bottom 50% of genes
png('output/cells_per_gene_bottom_all.png', width = 720, height = 360)
plot(sorted.cells.per.gene[1:(patient@assays$RNA@data@Dim[1] * 0.5),], xlab = 'Genes in rank order',
     ylab = 'Cells per gene', main = 'Bottom 50% of genes')
dev.off()

# Graph the top 10% of genes
# png('output/cells_per_gene_top_all.png', width = 1200, height = 400)
# plot(sorted.cells.per.gene[(patient@assays$RNA@data@Dim[1] - (patient@assays$RNA@data@Dim[1] * 0.1)):patient@assays$RNA@data@Dim[1],], 
#     xlab = 'Genes in rank order',
#     ylab = 'Cells per gene', main = 'Top 500 genes')
# dev.off()

# Write a table of cells with gene counts, UMI counts and percent.mito
write.csv(patient@meta.data, "output/barcode_nGene_nUMI_mito.csv")

# Save these files to output_filter_1 directory
# Use the graphs produced above to determine filtering parameters.
# Enter the filtering parameters in the next script, seurat_filter_2.R

