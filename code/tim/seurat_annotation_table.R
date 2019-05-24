# Create an annotation_table of barcodes with their cluster assignment,
#    cell type assignment, TSNE values and cluster colors

# Author: Tim Starr
# Last updated: 2/26/19

# Purpose: Generate an annotation table to be used for uploading into 
#   mysql and for generating colored TSNE plots

# Input file 1: cluster_assignments.csv: 
#   Table of barcodes and their cluster assignment created by seurat_final.R
#   and stored in output_final folder for each sample:
#   /Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/Patient_59/TenX_Analysis/seurat/current/output_final/cluster_assignments.csv

# Input file 2: tsne_cell_embeddings.csv: 
#   Table of TSNE coordinates for each barcode. Stored in same location as cluster_assignments.csv

# Note: Both cluster_assignments.csv and tsne_cell_embeddings are sorted
#   in the same order on barcodes.

# Input file 3: cell_type_by_cluster.csv
#   Table listing all clusters and their cell type annotation. This file
#   is generated using the R script cell_type_annotation_by_matrix.R.
#   This file is located in '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/ss_rna_seq_compilation/cell_types_seurat/'

# Input file 4: cell_type_dictionary.csv: 
#   A table of all possible cell annotations, their description and associated colors

# Output files:
#   annotation_table_seurat.csv: A table of barcodes and their cluster assignment, 
#     tsne embeddings, and cell type annotations (and colors)

#############################################################
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
# Define the file locations
main_directory <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/'
sub_directory <- '/TenX_Analysis/seurat/current/output_final/'
dictionary_filename <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/cell_type_dictionary.csv'
cell_type_filename <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/ss_rna_seq_compilation/cell_types_seurat/cell_type_by_cluster.csv'

# Read in the patient folder names
patient_folder_names <- list.files(main_directory)

# Load the cell_type_dictionary.csv
cell_type_dictionary <- read.csv(file = dictionary_filename)
cell_type_dictionary$cell_type_code <- as.character(cell_type_dictionary$cell_type_code)

# Load the cell_type_by_cluster.csv
cell_type_by_cluster <- read.csv(file = cell_type_filename)

# Manipulate cell_type_by_cluster to 
# Delete first column
cell_type_by_cluster$X <- NULL
# Rename cluster_id column
colnames(cell_type_by_cluster)[1] <- "sample_id"
# Copy the first column
cell_type_by_cluster$cluster_number <- cell_type_by_cluster$sample_id
# Clean up sample_id  
cell_type_by_cluster$sample_id <- gsub("_custom_markers_*.*", "", cell_type_by_cluster$sample_id)
# Clean up cluster_number
cell_type_by_cluster$cluster_number <- gsub(".*_custom_markers_", "", cell_type_by_cluster$cluster_number)
cell_type_by_cluster$cluster_number <- gsub(".csv", "", cell_type_by_cluster$cluster_number)

####################################################
# Loop through each patient folder, retrieve files, extract data, generate
#   output annotation table

for (i in 1:length(patient_folder_names))
{
  # i=1
  # Read in input files: cluster_assignments.csv and tsne_cell_embeddings.csv:
  cluster_assignments <- read.csv(file = paste(main_directory,patient_folder_names[i],
                                               sub_directory, "cluster_assignments.csv", sep = ""))
  tsne_cell_embeddings <- read.csv(file = paste(main_directory,patient_folder_names[i],
                                                sub_directory, "tsne_cell_embeddings.csv", sep = ""))

  
  # select the cluster assignments from cell_type_by_cluster based on sample_id
  temp_cluster_annotations <- subset(cell_type_by_cluster, sample_id == patient_folder_names[i])
  temp_cluster_annotations$cluster_number <- as.integer(as.character(temp_cluster_annotations$cluster_number))

  ####################################################
  # Generate a list of original colors used by Seurat (this is based on a
  #   color pallet that is constructed using the following function)
  cluster_colors <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  original_colors <- cluster_colors(nrow(temp_cluster_annotations))
  original_colors <- as.data.frame(original_colors)
  original_colors$cluster_number <- seq(from = 0, to = (nrow(temp_cluster_annotations)-1))
  ####################################################
  # Create the output annotation_table
  # Combine the cluster_assignments and tsne_cell_embeddings
  annotation_table <- cbind(cluster_assignments,tsne_cell_embeddings)
  colnames(annotation_table) <- c("barcode","cluster_number","barcode2","tsne1","tsne2")
  annotation_table$barcode2 = NULL

  # Add the annotations
  annotation_table <- left_join(x = annotation_table, y = temp_cluster_annotations,
                            by = c("cluster_number" = "cluster_number"))
  annotation_table$cell_type_code <- as.character(annotation_table$cell_type_code)
  # Add the colors and descriptions
  annotation_table <- left_join(x = annotation_table, y = cell_type_dictionary,
                                by = c("cell_type_code" = "cell_type_code"))
  annotation_table <- left_join(x = annotation_table, y = original_colors,
                                by = c("cluster_number" = "cluster_number"))
  
  write.csv(annotation_table, paste(main_directory,patient_folder_names[i],
                                    "/annotation_table_seurat.csv", sep = ""))

}