# Run a Seurat analysis on my macbook
# Author: Tim Starr
# Last updated: 11/9/18 using Seurat v 3.0.0.9000

# This script will perform various analyses of Seurat output

# 1) 2D TSNE using tsne_cell_embeddings
# 2) 2D TNSE plot using annotation_table_seurat.csv
# 3) 2D TSNE plot using Seurat Patient.Robj
# 4) 3D TSNE plot using Seurat Patient.Robj
# 5) 2D TSNE plot colored with heatmap of single gene
# 6) FindMarkers comparing selected groups
# 7) Violin plots using barcode_nGene_nUMI_mito.csv files
# 8) Return a table of JackStraw PC p-values
# 9) Find Clusters using patient.Robj
# 10) 2D UMAP using patient.Robj
# 11) Custom Marker Analysis one patient




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
####################################################
# Load previously generate seurat object from a patient folder
patient_folder <- 'Patient_59'
main_directory <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/'
sub_directory <- '/TenX_Analysis/seurat/current/output_final/Patient.Robj'
load(file = paste(main_directory,patient_folder,sub_directory, sep = ""))

# Load the annotation_table_seurat.csv from the input folder
annotation_table_seurat <- read.csv(file = 'input/annotation_table_seurat.csv')

# Load a previously generated seurat object from the input folder
load(file = "input/Patient.Robj")
###################################################
# Set all the default parameters
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

# FindNeighbors parameters
dims_max <- 38
k_param <- 33
prune_SNN <- 0.08
nn_eps <- 0

# FindClusters parameters
modularity_fxn <- 1
resolution_parameter <- 1.4
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

####################################################
# Generate a list of default colors used by ggplot to color clusters
# in the tSNE plot. Enter the number of clusters
cluster_number <- 25
cluster_colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
write.csv(cluster_colors(cluster_number), 'output/cluster_colors.csv')

####################################################
# 1) 2D TSNE using tsne_cell_embeddings
# Create a 2D TSNE plot using tsne_cell_embeddings.csv placed in input folder
#   Use the colors.txt file also in input folder
setwd('/Users/star0044old/Documents/Programming/R/Seurat')
data_matrix <- read.csv("input/tsne_cell_embeddings.csv", row.names = 1, header = TRUE)
color_list <- scan("input/colors.txt", what ="character")
png('output/tsne_colored.png', width = 1200, height = 1200)
plot(x=data_matrix$tSNE_1 , y=data_matrix$tSNE_2, type = "p", 
     main = "TSNE Plot", 
     xlab = "tSNE-1", ylab = "tSNE-2", 
     cex = 2.5, pch = 20,
     col = color_list)
dev.off()

# Use this for testing 2D plot
color_list <- scan("input/colors.txt", what ="character")
plot(x=data_matrix$tSNE_1 , y=data_matrix$tSNE_2, type = "p", 
     main = "TSNE Plot", 
     xlab = "tSNE-1", ylab = "tSNE-2", 
     cex = 0.6, pch = 20,
     col = color_list)

####################################################
# 2) 2D TNSE plot using annotation_table_seurat.csv
png('output/tsne_annot_table.png', width = 360, height = 360)
plot(x=annotation_table_seurat$tsne1 , y=annotation_table_seurat$tsne2, 
     type = "p", main = "TSNE Plot", 
     xlab = "tSNE-1", ylab = "tSNE-2", cex = 0.1, pch = 20,
     col = as.vector(annotation_table_seurat$original_colors))
dev.off()

##################################################
# 3) 2D TSNE plot using Seurat Patient.Robj

png('output/TSNE_2D_seurat_colors.png', width = 700, height = 500)
DimPlot(patient, dims = c(1, 2), cells = NULL, cols = NULL,
        pt.size = 3, reduction = "tsne", group.by = "ident",
        split.by = NULL, shape.by = NULL, order = NULL, label = TRUE,
        label.size = 12, repel = FALSE, cells.highlight = NULL,
        cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",
        combine = TRUE)
dev.off()

# Create TSNE 2D plot above with custom colors
color_list <- scan("input/colors.txt", what ="character")
png('output/TSNE_2D_colors.png', width = 500, height = 500)
DimPlot(patient, dims = c(1, 2), cells = NULL, cols = color_list,
        pt.size = 0.3, reduction = "tsne", group.by = "ident",
        split.by = NULL, shape.by = NULL, order = NULL, label = TRUE,
        label.size = 4, repel = FALSE, cells.highlight = NULL,
        cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",
        combine = TRUE)
dev.off()

##################################################
# 4) 3D TSNE plot using Seurat Patient.Robj

patient <- RunTSNE(object = patient, dim.embed = 3, 
                   dims.use = 1:significant_PCs, do.fast = TRUE)

# Create a data matrix for plotting the three dimensions.
tSNE_3D_matrix <- as.data.frame(patient@reductions$tsne@cell.embeddings)
# Save the TSNE coordinates
write.csv(patient@reductions$tsne@cell.embeddings, file = "output/tsne_3D_embeddings.csv")

# Plot the 3D matrix in black (this requires loading the library rgl)
# plot3d(tSNE_3D_matrix$tSNE_1, tSNE_3D_matrix$tSNE_2, tSNE_3D_matrix$tSNE_3, type="s", size=0.5, xlab = "tSNE-1", ylab = "tSNE-2", zlab = "tSNE-3")

# Plot the 3D matrix in colors (this requires loading the library rgl)
color_list <- scan("input/colors.txt", what ="character")

plot3d(tSNE_3D_matrix$tSNE_1, tSNE_3D_matrix$tSNE_2, 
       tSNE_3D_matrix$tSNE_3, 
       type="s", col = color_list, size=0.5, 
       xlab = "tSNE-1", ylab = "tSNE-2", zlab = "tSNE-3")

plot3d(tSNE_3D_matrix$tSNE_1, tSNE_3D_matrix$tSNE_2, tSNE_3D_matrix$tSNE_3, 
       type="s", col = color_list, size=0.5, xlab = "", ylab = "", zlab = "",
       axes(edges = bbox, ticks=FALSE))


# Create 90 sequential frames of the 3d image and place them into 
# the animation folder
for (i in 1:90) {
  view3d(userMatrix=rotationMatrix(2*pi * i/90, 1, -1, -1))
  rgl.snapshot(filename=paste("animation/frame-",
                              sprintf("%03d", i), ".png", sep=""))
}
# The frames are there, in the animation directory. You can quit R and get ready to assemble the frames in a .gif file. At the shell command line, you then need to use convert from imagemagick to create the animated gif.
# convert -delay 10 -loop 0 frame*.png animated.gif
# The option -delay 10 indicates the time between the frames and -loop 0 indicates to create a looping animation.

##################################################
# 5) 2D TSNE plot colored with heatmap of single gene 
#   no cut-offs
gene <- "UNC45A"
min_cutoff <- NA
max_cutoff <- NA

png(filename = paste("output/",gene,"_min_cutoff_", min_cutoff,
                     "_tsne.png", sep = ""), width = 500, height = 500)
FeaturePlot(patient, features = gene, dims = c(1, 2), cells = NULL,
            cols = c("yellow", "red"), pt.size = 0.4, 
            min.cutoff = min_cutoff,
            max.cutoff = max_cutoff, 
            reduction = "tsne", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = FALSE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()

# Create TSNE plot with heatmap of single gene with a min-cutoff
gene <- "GAPDH"
min_cutoff <- 0.9
max_cutoff <- min_cutoff + 0.00001

png(filename = paste("output/",gene,"_min_cutoff_", min_cutoff,
                      "_tsne.png", sep = ""), width = 500, height = 500)
FeaturePlot(patient, features = gene, dims = c(1, 2), cells = NULL,
            cols = c("yellow", "red"), pt.size = 0.5, 
            min.cutoff = min_cutoff,
            max.cutoff = max_cutoff, 
            reduction = "tsne", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = FALSE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()

# Create a series of TSNE plot heatmaps from a list of genes
gene_list <- scan("input/genes.txt", what ="character")
for (i in 1:length(gene_list))
{
  # i=2
  png(filename = paste("output/",gene_list[i],"_tsne_gradient.png", sep = ""), 
      width = 500, height = 500)
  temp_ggplot <- FeaturePlot(patient, features = gene_list[i], dims = c(1, 2), cells = NULL,
              cols = c("yellow", "red"), pt.size = 0.5, min.cutoff = NA,
              max.cutoff = NA, reduction = "tsne", split.by = NULL,
              shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
              order = NULL, label = FALSE, label.size = 4, ncol = NULL,
              combine = TRUE, coord.fixed = FALSE)
  print(temp_ggplot)
  dev.off()
}

##################################################################
# 6) FindMarkers comparing selected groups
cluster_group_1 <- c(0,1,3,6,8)
cluster_group_2 <- c(7) # Note: enter "NULL" if comparing against all other clusters
min_diff_pct <- 0.1

custom_cluster_markers <- FindMarkers(patient, ident.1 = cluster_group_1, 
            ident.2 = cluster_group_2, assay = assay_parameter, 
            features = NULL, 
            logfc.threshold = logfc_threshold, test.use = test_use, 
            min.pct = min_pct, min.diff.pct = min_diff_pct,
            verbose = TRUE, only.pos = only_pos, 
            max.cells.per.ident = max_cells_per_ident, random.seed = 1, 
            latent.vars = NULL, min.cells.feature = min_cells_feature,
            min.cells.group = min_cells_group, pseudocount.use = 1)

group_1_name <- paste(cluster_group_1, collapse = "-")
group_2_name <- paste(cluster_group_2, collapse = "-")
file_name = paste("output/group_1=", group_1_name,
                  "_group_2=", group_2_name, ".csv", sep="")
write.csv(custom_cluster_markers, file_name)


#################################################################
# 7) Violin plots using barcode_nGene_nUMI_mito.csv files
# Use the barcode_nGene_nUMI.csv file to create violin plots with custom
#   colors given to the different clusters.

setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Read in the data file
barcode_nGene_nUMI_mito <- read.csv('input/barcode_nGene_nUMI.csv', 
                                    header = TRUE)
# Read in the cluster assignments
cluster_assignments <- read.csv('input/cluster_assignments.csv', 
                                header = TRUE)
colnames(cluster_assignments) <- c("barcode", "cluster_number")
# Read in the cluster colors (a list of colors that is the same length
#   as the number of clusters present)
color_list <- scan("input/colors.txt", what ="character")

# Add the cluster_assignments to the table
barcode_nGene_nUMI_mito$cluster <- as.factor(cluster_assignments$cluster_number)

# print out the genes/cell violin plot
png(filename = "output/genes_per_cell_violin_custom_color.png")
ggplot(data = barcode_nGene_nUMI_mito, aes(x="cells",
          y=barcode_nGene_nUMI_mito$nFeature_RNA))  +
    geom_violin(fill = "white") +
    guides(fill=FALSE) +
  geom_jitter(cex = 0.3, aes(colour = factor(cluster))) +
  scale_colour_manual(name="colour", values = color_list) +
    ylab("Number of genes") +
    lims(y = c(0,10000)) +
    xlab(NULL) +
    ggtitle("Genes/cell")
  dev.off()  

# Print the umi violin plot
png(filename = "output/umi_per_cell_violin_custom_color.png")
  ggplot(data = barcode_nGene_nUMI_mito, aes(x="cells",
                                       y=nCount_RNA)) +
    geom_violin(fill = "white") +
    guides(fill=FALSE) +
    geom_jitter(cex = 0.3, aes(colour = factor(cluster))) +
    scale_colour_manual(name="colour", values = color_list) +
    ylab("Number of UMI") +
    lims(y = c(0,300000)) +
    xlab(NULL) +
    ggtitle("UMI/cell")
dev.off()  

# Print the mito violin plot
png(filename = "output/mito_per_cell_violin_custom_color.png")
ggplot(data = barcode_nGene_nUMI_mito, aes(x="cells",
                                           y=barcode_nGene_nUMI_mito$percent.mito)) +
  geom_violin(fill = "white") +
  guides(fill=FALSE) +
  geom_jitter(cex = 0.3, aes(colour = factor(cluster))) +
  scale_colour_manual(name="colour", values = color_list) +
  ylab("Percent Mito genes") +
  lims(y = c(0,1)) +
  xlab(NULL) +
  ggtitle("% Mito/cell")
dev.off()  
########################################################################
# 8) Return a table of JackStraw PC p-values
write.csv(patient@reductions$pca@jackstraw@overall.p.values,  
          file = "output/jackstraw_pvalues.csv")
#########################################################################
# Find Neighbors
dims_max <- 39
k_param <- 31
prune_SNN <- 0.16
resolution_parameter <- 0.8

patient <- FindNeighbors(patient, reduction = "pca", 
                         dims = 1:dims_max,
                         assay = assay_parameter, features = NULL, 
                         k.param = k_param, prune.SNN = prune_SNN,
                         nn.eps = nn_eps, verbose = TRUE, force.recalc = FALSE,
                         do.plot = FALSE, graph.name = NULL)

#########################################################################
# 9) Find Clusters using patient.Robj
patient <- FindClusters(patient, graph.name = NULL,
                        modularity.fxn = modularity_fxn, 
                        resolution = resolution_parameter, 
                        algorithm = algorithm_parameter,
                        n.start = n_start, n.iter = n_iter, random.seed = 0,
                        temp.file.location = NULL, edge.file.name = NULL, 
                        verbose = TRUE)
file_name_prefix <- paste(dims_max,"-",k_param,"-",prune_SNN,"-",resolution_parameter,"-", sep = "")
write.csv(patient@active.ident, paste("output/",file_name_prefix,"cluster_assignments.csv", sep = ""))
write.csv(table(patient@active.ident), paste("output/",file_name_prefix,"cluster_size.csv", sep = ""))
png(filename = paste("output/",file_name_prefix,"cluster_size.png", sep = ""), width = 600, height = 400)
barplot(table(patient@active.ident), xlab = "Cluster", 
        ylab = "Number of cells", main = "Cluster size",
        cex.lab = 1.5, cex.main = 1.5)
dev.off()

#########################################################################
# Create a pie chart of major cell types using the 
#   annotations_table.csv file

# Read in the file
patient_id <- "Pat_88_per_02"
annotations_table <- read.csv("input/annotation_table.csv",
                              header = TRUE)
pie_chart <- ggplot(annotations_table, aes(x="",
                              y=major_color,
                              fill = major_color)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggtitle(patient_id) +
  scale_fill_manual(values = c("blue","red","yellow")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()
  )

png(filename = paste("output/major_type_pie_chart.png", sep = ""), 
    width = 600, height = 400)
pie_chart
dev.off()

#########################################################################
# Fetch data

# Fetch the normalized UMI count for a single gene
fetched_data <- FetchData(patient, vars = 'UNC45A')
write.csv(file = 'output/fetched_data.csv', fetched_data)

#########################################################################
# 10) 2D UMAP using patient.Robj
patient <- RunUMAP(patient, dims = 1:5, reduction = "pca",
        features = NULL, assay = assay_parameter, 
        nneighbors = 30L, max.dim = 2L,
        min.dist = 0.3, reduction.name = "umap", 
        reduction.key = "UMAP_",
        metric = "correlation", seed.use = 42)

# Plot the UMAP
png('output/umap.png', width = 700, height = 500)
DimPlot(patient, dims = c(1, 2), cells = NULL, cols = NULL,
        pt.size = 0.5, reduction = "umap", group.by = "ident",
        split.by = NULL, shape.by = NULL, order = NULL, 
        label = TRUE, label.size = 8, repel = TRUE, 
        cells.highlight = NULL, cols.highlight = "red", 
        sizes.highlight = 1, na.value = "grey50",
        combine = TRUE)
dev.off()

# Save the UMAP cell embeddings
write.csv(patient@reductions$umap@cell.embeddings, file = "output/umap_cell_embeddings.csv")

# Create a 2D plot using umap_cell_embeddings.csv file
data_matrix <- read.csv("input/umap_cell_embeddings.csv", row.names = 1, header = TRUE)
color_list <- scan("input/colors.txt", what ="character")
png('output/umap_colored.png', width = 1200, height = 1200)
plot(x=data_matrix$UMAP_1 , y=data_matrix$UMAP_2, type = "p", 
     main = "TSNE Plot", 
     xlab = "tSNE-1", ylab = "tSNE-2", 
     cex = 2.5, pch = 20,
     col = color_list)
dev.off()
#########################################################################
# 11) Custom Marker Analysis one patient
# Analyze the custom marker analysis files to do the following
#   Create a table of # of significant up and down regulated genes 
#   per cluster

# Set the working directory
setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Clear the workspace global environment
rm(list=ls())

# Set the patient folder name (folder should be in the analyze folder)
patient_folder_name <- "ava_5903"

# Set the desired cutoff parameters
max_p_val_adj <- .05
min_avg_logFC <- 0.25
max_avg_logFC <- 0
min_high_avg_logFC <- 1
max_low_avg_logFC <- -1

# Set the data directories
main_directory <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/'
sub_directory <- '/TenX_Analysis/seurat/current/output_final/'

# Retrieve the names of the custom marker files
output_file_names <- list.files(paste(main_directory,patient_folder_name, 
                                      sub_directory, sep = ""))
custom_marker_indices <- startsWith(output_file_names, "custom_markers")
custom_marker_file_names <- subset(output_file_names, custom_marker_indices)

# Create an empty output table
custom_marker_output_table <- as.data.frame(custom_marker_file_names)
# Use custom marker file names as rownames
row.names(custom_marker_output_table) <- custom_marker_output_table$custom_marker_file_names
# Initialize output columns
custom_marker_output_table$num_sig_markers <-0
custom_marker_output_table$num_pos_markers <- 0
custom_marker_output_table$num_hi_pos_markers <- 0
custom_marker_output_table$num_neg_markers <- 0
custom_marker_output_table$num_low_neg_markers <-0

# Loop through each custom_marker file and extract information
for (i in 1:length(custom_marker_file_names))
{
  # i=1
  # Open the file
  temp_file_name <- paste(main_directory,patient_folder_name,
                          sub_directory,custom_marker_file_names[i], sep = "")
  temp_cluster_markers <- read.csv(temp_file_name, header = TRUE)

  # Identify the row index in custom_marker_output_table for 
  #   the cluster being analyzed
  temp_row_index <- which(rownames(custom_marker_output_table) %in% custom_marker_file_names[i])

  # Select markers that are significant (adj p val < max_p_val_adj)
  temp_sig_subset <- temp_cluster_markers[temp_cluster_markers$p_val_adj < max_p_val_adj, ]
  # Add the number of significant markers to the output table
  custom_marker_output_table$num_sig_markers[temp_row_index] <- nrow(temp_sig_subset)
  
    
  # Select markers that are positive (avg_logFC > min_avg_logFC) 
  #   and have an adj p val < max_p_val_adj
  temp_pos_subset <- temp_cluster_markers[(temp_cluster_markers$avg_logFC > min_avg_logFC) & (temp_cluster_markers$p_val_adj < max_p_val_adj), ]
  # Add the number of positive markers to the output table
  custom_marker_output_table$num_pos_markers[temp_row_index] <- nrow(temp_pos_subset)
 
  # Select markers that are highly positive min_high_avg_logFC)
  temp_hi_pos_subset <- temp_cluster_markers[(temp_cluster_markers$avg_logFC > min_high_avg_logFC) & (temp_cluster_markers$p_val_adj < max_p_val_adj), ]
  # Add the number of positive markers to the output table
  custom_marker_output_table$num_hi_pos_markers[temp_row_index] <- nrow(temp_hi_pos_subset)
  
  # Select markers that are negative and have an adj p val < p value set above
  temp_neg_subset <- temp_cluster_markers[(temp_cluster_markers$avg_logFC < max_avg_logFC) & (temp_cluster_markers$p_val_adj < max_p_val_adj), ]
  # Add the number of negative markers to the output table
  custom_marker_output_table$num_neg_markers[temp_row_index] <- nrow(temp_neg_subset)
  
  # Select markers that are highly negative
  temp_low_neg_subset <- temp_cluster_markers[(temp_cluster_markers$avg_logFC < max_low_avg_logFC) & (temp_cluster_markers$p_val_adj < max_p_val_adj), ]
  # Add the number of negative markers to the output table
  custom_marker_output_table$num_low_neg_markers[temp_row_index] <- nrow(temp_low_neg_subset)
  
  
  
}
# Write the output file
write.csv(custom_marker_output_table, file = "output/custom_marker_analysis.csv")




