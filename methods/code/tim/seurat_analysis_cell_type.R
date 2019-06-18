# Create annotation_table and TSNE plots based on cell type assignments
# Author: Tim Starr
# Last updated: 1/5/19 using Seurat v 3.0.0.9000

# Purpose: Generate TSNE plots colored by cell type assignments and
#   generate annotation file listing barcodes with their cluster and cell 
#   type assignments

# Input files
#   The first two files are created by seurat_final.R and placed in the output folder
#     cluster_assignments.csv: Table of barcodes and their cluster assignment
#       Note: this file contains only filtered barcodes and is sorted
#     tsne_cell_embeddings.csv: Table of TSNE coordinates for each barcode
#       Note: this file is sorted to match cluster_assignments.csv
#   The third and fourth files are created manually and stored in the analysis folder
#     cell_type_dictionary.csv: A table of all possible cell annotations, their description
#       and associated colors
#     cluster_annotations.csv: A table of cluster annotations
#       This is generated manually and lists the cluster number and its assigned
#       annotation based on judgement calls

# Output files:
#   annotation_table.csv: A table of barcodes and their cluster assignment, 
#     tsne embeddings, and cell type annotations (and colors)
#   tsne_original_colors.png
#   tsne_major_colors.png
#   tsne_minor_colors.png

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
patient_folder_name <- "Patient_111"
main_directory <- '/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/'
sub_directory_output <- '/TenX_Analysis/seurat/current/output_final/'
sub_directory_analysis <- '/TenX_Analysis/seurat/current/analysis/cell_types/'
sub_directory_tsne <- '/TenX_Analysis/seurat/current/analysis/tsne/'


# Load the input files
cluster_annotations <- read.csv(file = paste(main_directory,patient_folder_name,
                      sub_directory_analysis, "cluster_annotations.csv", sep = ""), 
                      header = TRUE, sep = ",")
cell_type_dictionary <- read.csv(file = paste(main_directory,patient_folder_name,
                        sub_directory_analysis, "cell_type_dictionary.csv", sep = ""), 
                        header = TRUE, sep = ",")
cluster_assignments <- read.csv(file = paste(main_directory,patient_folder_name,
                         sub_directory_output, "cluster_assignments.csv", sep = ""), 
                         header = TRUE, sep = ",")
tsne_cell_embeddings <- read.csv(file = paste(main_directory,patient_folder_name,
                        sub_directory_output, "tsne_cell_embeddings.csv", sep = ""), 
                        header = TRUE, sep = ",")

####################################################
# Generate a list of original colors used by Seurat (this is based on a
#   color pallet that is constructed using the following function)
cluster_colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
original_colors <- cluster_colors(nrow(cluster_annotations))
original_colors <- as.data.frame(original_colors)
original_colors$cluster_number <- seq(from = 0, to = (nrow(cluster_annotations)-1))
####################################################
# Create the output annotation_table
# Combine the cluster_assignments and tsne_cell_embeddings
annotation_table <- cbind(cluster_assignments,tsne_cell_embeddings)
colnames(annotation_table) <- c("barcode","cluster_number","barcode2","tsne1","tsne2")
annotation_table$barcode2 = NULL

# Add the annotations
annotation_table <- left_join(x = annotation_table, y = cluster_annotations,
                          by = c("cluster_number" = "cluster_number"))
# Add the colors and descriptions
annotation_table <- left_join(x = annotation_table, y = cell_type_dictionary,
                              by = c("cell_type_code" = "cell_type_code"))
annotation_table <- left_join(x = annotation_table, y = original_colors,
                              by = c("cluster_number" = "cluster_number"))

write.csv(annotation_table, paste(main_directory,patient_folder_name,
                                  sub_directory_analysis, "annotation_table.csv", sep = ""))
####################################################
# Create 2D TSNE plots using tsne_cell_embeddings.csv and the different colors

# Set the point size (cex)
cex_value <- 0.3

# Create 2D TSNE plot using original seurat colors
png(paste(main_directory,patient_folder_name,
          sub_directory_tsne, "tsne_original_colors.png", sep = ""), 
    width = 360, height = 360)
plot(x=annotation_table$tsne1 , y=annotation_table$tsne2, type = "p", 
     main = paste(patient_folder_name,"Original Seurat Colors"), 
     xlab = "TSNE 1", ylab = "TSNE 2", cex = cex_value, pch = 20,
     col = as.character(annotation_table$original_colors))
dev.off()

# Create 2D TSNE plot using major cell type colors
png(paste(main_directory,patient_folder_name,
        sub_directory_tsne, "tsne_major_colors.png", sep = ""), 
        width = 360, height = 360)
plot(x=annotation_table$tsne1 , y=annotation_table$tsne2, type = "p", 
     main = paste(patient_folder_name,"Major Cell Type Colors"), 
     xlab = "TSNE 1", ylab = "TSNE 2", cex = cex_value, pch = 20,
     col = as.character(annotation_table$major_color))
dev.off()

# Create 2D TSNE plot using minor cell type colors
png(paste(main_directory,patient_folder_name,
          sub_directory_tsne, "tsne_minor_colors.png", sep = ""), 
    width = 360, height = 360)
plot(x=annotation_table$tsne1 , y=annotation_table$tsne2, type = "p", 
     main = paste(patient_folder_name,"Minor Cell Type Colors"), 
     xlab = "TSNE 1", ylab = "TSNE 2", cex = cex_value, pch = 20,
     col = as.character(annotation_table$minor_color))
dev.off()