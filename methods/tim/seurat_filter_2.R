# Second Seurat analysis to determine 

# Author: Tim Starr
# Last updated: 11/15/18 using Seurat v 3.0.0.9000

# Set the working directory
#   Use this for working on my own laptop
# setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Clear the workspace global environment
rm(list=ls())
#   Use this for running on MSI and change the patient name
setwd('~/ovarian/ssrnaseq/cr_starr_msi')

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

# Filtered parameters
min_cells_per_gene <- 20

min_genes_per_cell <- 100
max_genes_per_cell <- 10000

min_umi_cell <- 1
max_umi_cell <- 200000

min_percent_mito <- 0.0
max_percent_mito <- 0.75
############################
# Create Seurat object 
#   and print out cell violin plots for gene_count, UMI_count and % Mito
patient <- CreateSeuratObject(patient.data, project = "Patient", 
                              assay = assay_parameter, min.cells = min_cells_per_gene, 
                              min.features = min_genes_per_cell, names.field = 1,
                              names.delim = "-", meta.data = NULL)

# Create a metadata entry containing % mitochondrial genes
# Select the appropriate pattern for matching mito gene names
# mito.genes <- grep(pattern = "^GRch38_MT-", x = rownames(x = patient@data), value = TRUE)
mito.genes <- grep(pattern = "^MT-", x = rownames(patient@assays$RNA), value = TRUE)

percent.mito <- Matrix::colSums(patient@assays$RNA@counts[mito.genes, ])/Matrix::colSums(patient@assays$RNA@counts)
patient <- AddMetaData(object = patient, metadata = percent.mito, col.name = "percent.mito")

########################################################################
# Filter the patient object based on parameters above
filtered_patient <- subset(patient, 
    subset = nFeature_RNA < max_genes_per_cell &
    nCount_RNA > min_umi_cell & nCount_RNA < max_umi_cell &
    percent.mito > min_percent_mito & percent.mito < max_percent_mito)

# Replace original patient object with filtered object
#   and remove the filtered_patient object
patient <- filtered_patient
rm(filtered_patient)

#####################################################################

# Calculate cell count per gene: For each gene calculate the number of cells with UMI > 0
cells.per.gene <- rowSums(patient@assays$RNA@data !=0)
sorted.cells.per.gene <- data.frame(sort(cells.per.gene))

# Place the cells.per.gene in the meta.features slot
# metadata.cells.per.gene <- data.frame(cells.per.gene)
# patient@assays$RNA@meta.features <- cbind(patient@assays$RNA@meta.features, metadata.cells.per.gene)

########################################################################
# Write a table of cells with gene counts, UMI counts and percent.mito
write.csv(patient@meta.data, "output/barcode_nGene_nUMI.csv")
#########################################################################
# NormalizeData parameters
normalization_method <- "LogNormalize" # other option is CLR
scale_factor <- 10000 # According to Rachel, the larger the scale factor, smaller UMI counts with small changes will be emphasized more

# Normalize the data (Maybe try doing this once without normalizing?)
patient <- NormalizeData(object = patient, 
          normalization.method = normalization_method, 
         scale.factor = scale_factor)
#########################################################################
# FindVariableFeatures parameters
# Play around with the parameters and view the graph
#   go with your gut instinct

selection_method <- "vst" # options are dispersion, mean.var.plot or vst
# Dispersion favors high expression genes compared to vst, while max.var.plot can be manipulated by setting minimums and maximums on the mean and dispersion
loess_span <- 0.3 #only applicable to vst selection_method
clip_max <- "auto" #only applicable to vst selection_method
num_bin <- 20
binning_method <- "equal_width"
nfeatures_parameter <- 0.067 #only applicable to dispersion and vst selection_method.
# I set the nFeatures at 6.7% of genes, which will result in ~1000 genes if there are 15,000 genes
min_mean_cutoff <- 0.1 #only applicable to mean.var.plot
max_mean_cutoff <- 10 #only applicable to mean.var.plot
min_dispersion_cutoff <- 1 #only applicable to mean.var.plot
max_dispersion_cutoff <- 20 #only applicable to mean.var.plot

# Find variable genes
patient <- FindVariableFeatures(patient, assay = assay_parameter,
          selection.method = selection_method, loess.span = loess_span, 
          clip.max = clip_max,
          num.bin = num_bin, binning.method = binning_method, 
          nfeatures = round(nfeatures_parameter * patient@assays$RNA@data@Dim[1]), 
          mean.cutoff = c(min_mean_cutoff, max_mean_cutoff), 
          dispersion.cutoff = c(min_dispersion_cutoff, max_dispersion_cutoff))

VariableFeaturePlot(patient, cols = c("black", "red"), pt.size = 1,
                    log = NULL, assay = assay_parameter)


png('output/variable_gene_plot.png', width = 1200, height = 400)
VariableFeaturePlot(patient, cols = c("black", "red"), pt.size = 1,
                    log = NULL, assay = assay_parameter)
dev.off()
write.csv(patient@assays$RNA@var.features, "output/variable_genes.csv")

#########################################################################
# ScaleData parameters
# I don't know how this affects things, so I'm going with default
vars_to_regress <- c('nCount_RNA', 'percent.mito')
model_use <- "linear" # Other options are poisson and negbinom. Josh says linear is probably best
use_umi <- FALSE # Default is FALSE for linear modeling, but automatically set to TRUE of model.use is 'negbinom' or 'poisson'
do_scale <- TRUE
do_center <- TRUE
scale_max <- 50
block_size <- 1000
min_cells_to_block <- 1000

# Scale the data
patient <- ScaleData(patient, features = NULL, assay = assay_parameter,
      vars.to.regress = vars_to_regress, model.use = model_use, 
      use.umi = use_umi, do.scale = do_scale, do.center = do_center, 
      scale.max = scale_max, block.size = block_size, 
      min.cells.to.block = min_cells_to_block)

#########################################################################
# RunPCA and Jackstraw parameters
#   Start with 50 PCs and view jackstraw graph, graph smaller values
#   until the jackstraw plots are well above the dotted line.
number_pcs <- 50
score_thresh <- 1e-05
num_replicate <- 100
prop_freq <- 0.01
maxit_parameter <- 1000

# Add 10 to the number of pcs for calculating above
npcs_plus_10 <- number_pcs + 10

# Run PCA
patient <- RunPCA(patient, assay = assay_parameter, features = NULL,
        npcs = npcs_plus_10, rev.pca = FALSE, weight.by.var = TRUE, 
        verbose = TRUE, ndims.print = 1:5, nfeatures.print = 30, 
        reduction.name = "pca", reduction.key = "PC_", seed.use = 42)

patient <- ProjectDim(patient, reduction = "pca", assay = NULL, 
          dims.print = 1:5, nfeatures.print = 20, overwrite = FALSE, 
          do.center = FALSE, verbose = TRUE)

write.csv(patient@reductions$pca@cell.embeddings, file = "output/pc_values.csv")

png('output/pca_plot.png', width = 360, height = 360)
DimPlot(patient, dims = 1:2, cells = NULL, cols = NULL,
        pt.size = NULL, reduction = NULL, group.by = NULL,
        split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,
        label.size = 4, repel = FALSE, cells.highlight = NULL,
        cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",
        combine = TRUE)
dev.off()

png('output/pca_heatmap.png')
DimHeatmap(patient, dims = 1:6, nfeatures = 30, cells = NULL,
           reduction = "pca", disp.min = -2.5, disp.max = NULL,
           balanced = TRUE, projected = FALSE, ncol = NULL, combine = TRUE,
           fast = TRUE, slot = "scale.data", assays = assay_parameter)
dev.off()

png('output/PCElbowPlot.png', width = 360, height = 360)
ElbowPlot(patient, ndims = npcs_plus_10, reduction = "pca")
dev.off()

patient <- JackStraw(patient, reduction = "pca",
        assay = assay_parameter, dims = number_pcs,
        num.replicate = num_replicate, prop.freq = prop_freq,
        verbose = TRUE, maxit = maxit_parameter)

patient <- ScoreJackStraw(patient, reduction = "pca",
        dims = 1:number_pcs,
        score.thresh = score_thresh, do.plot = FALSE)

png('output/jackstraw.png', width = 1200, height = 400)
  JackStrawPlot(patient, dims = 1:number_pcs, reduction = "pca",
              xmax = 0.1, ymax = 0.3)
dev.off()

write.csv(patient@reductions$pca@jackstraw@overall.p.values,  
          file = "output/jackstraw_pvalues.csv")

save(patient, file = "output/Patient.Robj")


