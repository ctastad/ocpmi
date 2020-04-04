# Final seurat pipeline to be used after all parameters are set 

# Author: Tim Starr
# Last updated: 1/5/19 using Seurat v 3.0.0.9000

# Set the working directory
#   Use this one for running on my mac
# setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Clear the workspace global environment
rm(list=ls())
#   Use this one for running at MSI (change the patient folder)
setwd('~/ovarian/ssrnaseq/pat_#')
# FindNeighbors parameters
dims_max <- #
k_param <- #
prune_SNN <- 0.08
# FindClusters parameters
resolution_parameter <- 1.4


# Load required libraries
library(Seurat)
library(dplyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(rgl)
library(reshape2)
library(gplots)
library(RColorBrewer)

# Save the version of Seurat
seurat_version <- as.character(packageVersion("Seurat"))

# CreateSeuratObject
assay_parameter <- "RNA"
min_cells_per_gene <- 20
min_genes_per_cell <- 100
max_genes_per_cell <- 10000
min_umi_cell <- 1
max_umi_cell <- 200000
min_percent_mito <- 0.0
max_percent_mito <- 0.75

# NormalizeData parameters
normalization_method <- "LogNormalize"
scale_factor <- 10000

# FindVariableFeatures parameters
selection_method <- "vst"
loess_span <- 0.3
clip_max <- "auto"
num_bin <- 20
binning_method <- "equal_width"
nfeatures_parameter <- 0.067
min_mean_cutoff <- 0.1
max_mean_cutoff <- 10
min_dispersion_cutoff <- 1
max_dispersion_cutoff <- 20

# ScaleData parameters
vars_to_regress <- c('nCount_RNA', 'percent.mito')
model_use <- "linear"
use_umi <- FALSE
do_scale <- TRUE
do_center <- TRUE
scale_max <- 50
block_size <- 1000
min_cells_to_block <- 1000

# RunPCA and Jackstraw parameters
number_pcs <- 50
score_thresh <- 1e-05
num_replicate <- 100
prop_freq <- 0.01
maxit_parameter <- 1000

# FindNeighbors fixed parameters
nn_eps <- 0

# FindClusters parameters fixed parameters
modularity_fxn <- 1
algorithm_parameter <- 1
n_start <- 10
n_iter <- 10

# RunTSNE parameters
add_iter <- 0
dim_embed <- 2

# FindAllMarkers and FindMarkers parameters
logfc_threshold <- 0.25
test_use <- "wilcox"
min_pct <- 0.1
min_diff_pct <- -Inf
only_pos <- FALSE
max_cells_per_ident <- Inf
min_cells_feature <- 3
min_cells_group <- 3
return_thresh <- 0.01


#######################################################################
# Load data
patient.data <- Read10X("cell_ranger_input_files/")
write.csv(patient.data@Dim, "output/sparse_matrix_dimensions.csv")
########################################################################
# Create Seurat object 
patient <- CreateSeuratObject(patient.data, project = "Patient", 
                              assay = assay_parameter, min.cells = min_cells_per_gene, 
                              min.features = min_genes_per_cell, names.field = 1,
                              names.delim = "-", meta.data = NULL)

mito.genes <- grep(pattern = "^MT-", x = rownames(patient@assays$RNA), value = TRUE)
percent.mito <- Matrix::colSums(patient@assays$RNA@counts[mito.genes, ])/Matrix::colSums(patient@assays$RNA@counts)
patient <- AddMetaData(object = patient, metadata = percent.mito, col.name = "percent.mito")

########################################################################
# Filter the patient object based on parameters above
filtered_patient <- subset(patient, 
    subset = nFeature_RNA < max_genes_per_cell &
    nCount_RNA > min_umi_cell & nCount_RNA < max_umi_cell &
    percent.mito > min_percent_mito & percent.mito < max_percent_mito)

patient <- filtered_patient
rm(filtered_patient)

#####################################################################
# Print out post-filter graphics and statistics
vln_genes <- VlnPlot(patient, features = 'nFeature_RNA', 
                     cols = NULL, pt.size = 1, idents = NULL,
                     sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL,
                     adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE,
                     ncol = NULL, combine = TRUE, slot = "data")

png('output/violin_genes.png', width = 360, height = 360)
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
png('output/violin_UMI.png', width = 360, height = 360)
vln_UMI + 
  theme(panel.background = element_rect(fill = "lightgray"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  labs(title="UMI total per cell", 
       y="# of UMI", x="single cells") +
  guides(fill=FALSE)
dev.off()

vln_mito <- VlnPlot(patient, features = "percent.mito")
png('output/violin_mito.png', width = 360, height = 360)
vln_mito + 
  theme(panel.background = element_rect(fill = "lightgray"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  labs(title="% Mitochondrial genes per cell", 
       y="% Mitochondrial genes", x="single cells") +
  guides(fill=FALSE)
dev.off()

# Create a single graphic with three violin plots for genes, UMI and % mito
# png('output/violin_gene_UMI_mito.png')
# VlnPlot(object = patient, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
# dev.off()

# Calculate cell count per gene: For each gene calculate the number of cells with UMI > 0
cells.per.gene <- rowSums(patient@assays$RNA@data !=0)
sorted.cells.per.gene <- data.frame(sort(cells.per.gene))

# Place the cells.per.gene in the meta.features slot
metadata.cells.per.gene <- data.frame(cells.per.gene)
patient@assays$RNA@meta.features <- cbind(patient@assays$RNA@meta.features, metadata.cells.per.gene)

# Write the cells.per.gene data to an output file
write.csv(cells.per.gene, "output/cells_per_gene.csv")

# Graph the entire vector of cells.per.gene
png('output/cells_per_gene.png', width = 720, height = 360)
plot(sort(cells.per.gene), xlab = 'Genes in rank order',
     ylab = 'Cells per gene', main = 'Cells per gene')
dev.off()

# Graph the bottom 50% of genes
png('output/cells_per_gene_bottom.png', width = 720, height = 360)
plot(sorted.cells.per.gene[1:(patient@assays$RNA@data@Dim[1] * 0.5),], xlab = 'Genes in rank order',
     ylab = 'Cells per gene', main = 'Bottom 50% of genes')
dev.off()

########################################################################
# Write the list of genes, cells and matrix dimensions
write.csv(patient@assays$RNA@data@Dimnames[1], file = "output/genes.csv")
write.csv(patient@assays$RNA@data@Dimnames[2], file = "output/cells.csv")
write.csv(patient@assays$RNA@data@Dim, "output/matrix_dimensions.csv")

# Write a table of cells with gene counts, UMI counts and percent.mito
write.csv(patient@meta.data, "output/barcode_nGene_nUMI.csv")
#########################################################################
# Normalize the data (Maybe try doing this once without normalizing?)
patient <- NormalizeData(object = patient, 
          normalization.method = normalization_method, 
         scale.factor = scale_factor)

#########################################################################
# Find variable genes
patient <- FindVariableFeatures(patient, assay = assay_parameter,
          selection.method = selection_method, loess.span = loess_span, 
          clip.max = clip_max,
          num.bin = num_bin, binning.method = binning_method, 
          nfeatures = round(nfeatures_parameter * patient@assays$RNA@data@Dim[1]), 
          mean.cutoff = c(min_mean_cutoff, max_mean_cutoff), 
          dispersion.cutoff = c(min_dispersion_cutoff, max_dispersion_cutoff))

png('output/variable_gene_plot.png', width = 1200, height = 400)
VariableFeaturePlot(patient, cols = c("black", "red"), pt.size = 1,
                    log = NULL, assay = assay_parameter)
dev.off()
write.csv(patient@assays$RNA@var.features, "output/variable_genes.csv")

#########################################################################
# Scale the data
patient <- ScaleData(patient, features = NULL, assay = assay_parameter,
      vars.to.regress = vars_to_regress, model.use = model_use, 
      use.umi = use_umi, do.scale = do_scale, do.center = do_center, 
      scale.max = scale_max, block.size = block_size, 
      min.cells.to.block = min_cells_to_block)

#########################################################################
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

png('output/jackstraw_all.png', width = 1200, height = 400)
  JackStrawPlot(patient, dims = 1:number_pcs, reduction = "pca",
              xmax = 0.1, ymax = 0.3)
dev.off()

png('output/jackstraw_significant.png', width = 1200, height = 400)
JackStrawPlot(patient, dims = 1:dims_max, reduction = "pca",
              xmax = 0.1, ymax = 0.3)
dev.off()

write.csv(patient@reductions$pca@jackstraw@overall.p.values,  
          file = "output/jackstraw_pvalues.csv")
#########################################################################
# Find Neighbors
patient <- FindNeighbors(patient, reduction = "pca", 
          dims = 1:dims_max,
          assay = assay_parameter, features = NULL, 
          k.param = k_param, prune.SNN = prune_SNN,
          nn.eps = nn_eps, verbose = TRUE, force.recalc = FALSE,
          do.plot = FALSE, graph.name = NULL)

#########################################################################
# Find Clusters
patient <- FindClusters(patient, graph.name = NULL,
                        modularity.fxn = modularity_fxn, 
                        resolution = resolution_parameter, 
                        algorithm = algorithm_parameter,
                        n.start = n_start, n.iter = n_iter, random.seed = 0,
                        temp.file.location = NULL, edge.file.name = NULL, 
                        verbose = TRUE)

write.csv(patient@active.ident, "output/cluster_assignments.csv")
write.csv(table(patient@active.ident), "output/cluster_size.csv")
png('output/cluster_size.png', width = 600, height = 400)
barplot(table(patient@active.ident), xlab = "Cluster", 
        ylab = "Number of cells", main = "Cluster size",
        cex.lab = 1.5, cex.main = 1.5)
dev.off()

#########################################################################
# RunTSNE
patient <- RunTSNE(patient, reduction = "pca", cells = NULL,
                   dims = 1:dims_max, features = NULL, seed.use = 1, 
                   tsne.method = "Rtsne",
                   add.iter = add_iter, dim.embed = dim_embed, 
                   distance.matrix = NULL,
                   reduction.name = "tsne", reduction.key = "tSNE_")
write.csv(patient@reductions$tsne@cell.embeddings, "output/tsne_cell_embeddings.csv")

png('output/TSNE_2D.png', width = 360, height = 360)
DimPlot(patient, dims = c(1, 2), cells = NULL, cols = NULL,
        pt.size = 0.1, reduction = "tsne", group.by = "ident",
        split.by = NULL, shape.by = NULL, order = NULL, label = TRUE,
        label.size = 4, repel = FALSE, cells.highlight = NULL,
        cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",
        combine = TRUE)
dev.off()

#########################################################################
# FindAllMarkers and FindMarkers
all_cluster_markers <- FindAllMarkers(patient, assay = assay_parameter, 
    features = NULL, logfc.threshold = logfc_threshold, 
    test.use = test_use, min.pct = min_pct,
    min.diff.pct = min_diff_pct, verbose = TRUE, 
    only.pos = only_pos, max.cells.per.ident = max_cells_per_ident, 
    random.seed = 1, latent.vars = NULL,
    min.cells.feature = min_cells_feature, 
    min.cells.group = min_cells_group, pseudocount.use = 1,
    return.thresh = return_thresh)
write.csv(all_cluster_markers, 'output/cluster_markers_all.csv')

# FindMarkers
cluster_number <- nlevels(patient@active.ident)
for (i in 1:cluster_number)
{
  cluster_id = i-1
  custom_cluster_markers <- FindMarkers(patient, ident.1 = cluster_id, 
      ident.2 = NULL, assay = assay_parameter, features = NULL, 
      logfc.threshold = logfc_threshold, test.use = test_use, 
      min.pct = min_pct, min.diff.pct = min_diff_pct,
      verbose = TRUE, only.pos = only_pos, 
      max.cells.per.ident = max_cells_per_ident, random.seed = 1, 
      latent.vars = NULL, min.cells.feature = min_cells_feature,
      min.cells.group = min_cells_group, pseudocount.use = 1)
  file_name = paste("output/custom_markers_",i-1, ".csv", sep="")
  write.csv(custom_cluster_markers, file_name)
}

############################
# Print out the parameters
parameter_headers <- c("seurat_version",
                       "assay_parameter",
                       "min_cells_per_gene",
                       "min_genes_per_cell",
                       "max_genes_per_cell",
                       "min_umi_cell",
                       "max_umi_cell",
                       "min_percent_mito",
                       "max_percent_mito",
                       "normalization_method",
                       "scale_factor",
                       "selection_method",
                       "loess_span",
                       "clip_max",
                       "num_bin",
                       "binning_method",
                       "nfeatures_parameter",
                       "min_mean_cutoff",
                       "max_mean_cutoff",
                       "min_dispersion_cutoff",
                       "max_dispersion_cutoff",
                       "vars_to_regress_nCount","vars_to_regress_percent_mito", 
                       "model_use",
                       "use_umi",
                       "do_scale",
                       "do_center",
                       "scale_max",
                       "block_size",
                       "min_cells_to_block",
                       "number_pcs",
                       "score_thresh",
                       "num_replicate",
                       "prop_freq",
                       "maxit_parameter",
                       "dims_max",
                       "k_param",
                       "prune_SNN",
                       "nn_eps",
                       "modularity_fxn",
                       "resolution_parameter",
                       "algorithm_parameter",
                       "n_start",
                       "n_iter",
                       "add_iter",
                       "dim_embed",
                       "logfc_threshold",
                       "test_use",
                       "min_pct",
                       "min_diff_pct",
                       "only_pos",
                       "max_cells_per_ident",
                       "min_cells_feature",
                       "min_cells_group",
                       "return_thresh")
parameters <- c(seurat_version,
                assay_parameter,
                min_cells_per_gene,
                min_genes_per_cell,
                max_genes_per_cell,
                min_umi_cell,
                max_umi_cell,
                min_percent_mito,
                max_percent_mito,
                normalization_method,
                scale_factor,
                selection_method,
                loess_span,
                clip_max,
                num_bin,
                binning_method,
                nfeatures_parameter,
                min_mean_cutoff,
                max_mean_cutoff,
                min_dispersion_cutoff,
                max_dispersion_cutoff,
                vars_to_regress[1], vars_to_regress[2],
                model_use,
                use_umi,
                do_scale,
                do_center,
                scale_max,
                block_size,
                min_cells_to_block,
                number_pcs,
                score_thresh,
                num_replicate,
                prop_freq,
                maxit_parameter,
                dims_max,
                k_param,
                prune_SNN,
                nn_eps,
                modularity_fxn,
                resolution_parameter,
                algorithm_parameter,
                n_start,
                n_iter,
                add_iter,
                dim_embed,
                logfc_threshold,
                test_use,
                min_pct,
                min_diff_pct,
                only_pos,
                max_cells_per_ident,
                min_cells_feature,
                min_cells_group,
                return_thresh)
parameter_output <- cbind(parameter_headers,parameters)
write.csv(parameter_output, "output/parameters.csv")
##############################

save(patient, file = "output/Patient.Robj")

#######################################################
# Calculate overlap of custom_markers with annotated gene lists
#   Uses the "cell_type_marker_analysis_multiple_files.R" script
logfc <- 0.25
adj_pval <- 0.05
# Read in the file names of the cluster marker files
output_file_names <- list.files('output/')
custom_marker_indices <- startsWith(output_file_names, "custom_markers")
custom_marker_file_names <- subset(output_file_names, custom_marker_indices)

annotated_files <- list.files("~/ovarian/ssrnaseq/annotated_gene_lists/")
# Create an empty matrix with annotated list columns and cluster name rows
matrix <- matrix(0, dimnames = list(custom_marker_file_names, annotated_files), 
                   nrow = length(custom_marker_file_names), ncol = length(annotated_files))
  
  ###############################################################
  # Loop through each cluster marker file
  #   Identify genes with positive fold change and adj p-val > 0.05
  #   Calculate overlap of these genes with each annotated gene list
  #   Store the percentage overlap in the matrix
  
for (i in 1:length(custom_marker_file_names))
{
  # i=1
  # Input the custom_markers.csv file
  temp_filename_cluster_markers <- paste("output/",
                                         custom_marker_file_names[i],sep = "")
  custom_cluster_markers <- read.csv(temp_filename_cluster_markers, 
                                       header = TRUE, sep = ",", quote = "\"",
                                       dec = ".", fill = TRUE)
  ########################################################################
  # Select markers that are positive and have an adj p val < p value set above
  cluster_subset <- custom_cluster_markers[(custom_cluster_markers$avg_logFC > logfc) & (custom_cluster_markers$p_val_adj < adj_pval), ]
    
 ########################################################################
  # Loop through the annotated gene lists and calculate the
  #   percentage of genes in cluster_subset that overlap
  for (j in 1:length(annotated_files))
  {
    # j = 1
    temp_filename_annotated_files <- paste ("~/ovarian/ssrnaseq/annotated_gene_lists/",
                                            annotated_files[j], sep = "")
    temp_annotated_list <- scan(temp_filename_annotated_files, 
                                what = "character")
    genes_in_list <- cluster_subset[cluster_subset$X %in% temp_annotated_list,]
    matrix[custom_marker_file_names[i],annotated_files[j]] <- nrow(genes_in_list)/nrow(cluster_subset)
  }
  write.csv(matrix, file = "output/cell_type_percentages.csv")
  
  # Create a heatmap from the output file
  my_palette <- colorRampPalette(c("white","gray","red"))(n = 100)
  png(filename = "output/cell_type_heatmap.png", width = 800, height = 600)
  heatmap.2(matrix, dendrogram='none', Rowv = FALSE, 
            Colv = FALSE, trace = 'none', density.info = "none",
            margins = c(28,16), col = my_palette,
            sepcolor = "black", keysize = 0.9,
            colsep = 1:ncol(matrix),
            rowsep = 1:nrow(matrix),
            offsetCol = 0.1, offsetRow = 0.1,
            cexRow = 1.5, cexCol = 2, lhei = c(0.5,2.5),
            main = "Cell Type Heatmap")
  dev.off()
  
  # Create a scaled heatmap
  my_palette <- colorRampPalette(c("white","white","red"))(n = 100)
  png(filename = "output/cell_type_heatmap_scaled.png", width = 800, height = 600)
  heatmap.2(matrix, dendrogram='none', Rowv = FALSE, scale = "row",
            Colv = FALSE, trace = 'none', density.info = "none",
            margins = c(28,16), col = my_palette,
            sepcolor = "black", keysize = 0.9,
            colsep = 1:ncol(matrix),
            rowsep = 1:nrow(matrix),
            offsetCol = 0.1, offsetRow = 0.1,
            cexRow = 1.5, cexCol = 2, lhei = c(0.5,2.5),
            main = "Scaled Cell Type Heatmap")
  dev.off()
  
  
}