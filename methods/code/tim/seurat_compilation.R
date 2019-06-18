# Analyze the seurat output files
# by Tim Starr
# Last updated 5/23/18

# Use this script to do various things with the Seurat output
#   files

# 1) Compile cluster sizes
# 2) Create pie charts
# 3) Color tsne charts by major cell type
# 4) Retrive number of clusters per sample and create a spreadsheet and graph

# Input files:
#   Patient-specific input files are all stored in a similar directory format:
#     /Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyzed/Patient_54/TenX_Analysis/seurat/current/output/
#     With the only difference being the patient # (in the example above, this would be Patient_54)
#     This script will look in the directory: /Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyzed/
#       and extract the number of patient folders
#  cluster_size.csv = a file listing each cluster and the number of cells in the cluster
#     There will be one of these for each patient in the output folder
#   cell_type.csv = a manually generated file listing all 

# Output files:
#   cluster_size_compilation.csv = a mat


library(dplyr)
library(gplots)
library(RColorBrewer)
setwd('/Users/star0044old/Documents/Programming/R/Seurat')
# Clear the workspace global environment
rm(list=ls())

# Define the input directories
analyzed_dir <- "/Users/star0044old/Documents/Projects/Ovarian/Ovarian_Avatars/patients/analyze/"
output_final_dir <- "/TenX_Analysis/seurat/current/output_final/"

# Create a vector of patient folder names
patient_folder_names <- list.files(analyzed_dir)
number_of_patients <- length(patient_folder_names)

# Alternatively, define the number of samples  based on manual input
# number_of_patients <- 11



#################################################################
# 1) Compile cluster sizes

# Create output data.frame
cluster_size_compilation <- data.frame(patient = character(),
                     cluster=integer(),
                     frequency = integer()
                      )
for (i in 1:number_of_patients)
{
  # Read in the cluster_size.csv file for each patient
  filename <- paste(analyzed_dir, patient_folder_names[i], output_final_dir, 
                    "cluster_size.csv", sep = "")
  cluster_size <- read.csv(filename, header = TRUE)
  # Remove the X column from the imported data frame
  cluster_size[,1] <- NULL
  # Add in a column with the patient name
  cluster_size <- cbind(patient_folder_names[i], cluster_size)
  # Add a column with the total cells for the patient
  cluster_size <- cbind(cluster_size,sum(cluster_size$Freq))
  # Add a column with the percentage
  cluster_size <- cbind(cluster_size, cluster_size$Freq/sum(cluster_size$Freq))
  # Add the data to the output dataframe
  cluster_size_compilation <- rbind(cluster_size_compilation,cluster_size)
}
# Rename the column headers and save the table
names(cluster_size_compilation) <- c("Patient", "Cluster_number", 
      "Number_of_cells", "Total_cells", "Percent_of_total")
write.csv(cluster_size_compilation, file = "output/cluster_size_compilation.csv")

# Save a table of cluster numbers per sample
sample_cluster_count <- data.frame(table(cluster_size_compilation$Patient))
colnames(sample_cluster_count) <- c('sample_id','cluster_count')
write.csv(sample_cluster_count, file = "output/sample_cluster_count.csv")


#################################################################
# 2) Calculate % of major cell types and create pie charts
#   These are created using a manually curated cell_type.csv file
#   To generate this file see the worksheet "Cluster Assignment" in
#   the "Single_cell_Compilation" spreadsheet

# The cell_type.csv file should look like this:
#   patient	cluster_number	number_of_cells	total_cells	percent_of_total	cell_type_major	cell_type_minor	high_cycle
#   Patient_88_dia_01	0	485	1083	0.447830102	immune		
#   Patient_88_dia_01	1	218	1083	0.201292705	stroma		
#   Patient_88_dia_01	2	131	1083	0.120960295	epithelial		
#   Patient_88_dia_01	3	127	1083	0.117266851	immune	b-cells
# Note: this script only uses the "cell_type_major" annotation 
#   column and ignores the other two.

# Input the cell_type.csv file
cell_type <- read.csv("input/cell_type.csv", header = TRUE)

# Create a table of patients
patient_table <- table(cell_type$patient)
patient_table_length <- length(patient_table)

# Create output dataframe 
#   (Note: need to use exact words used in creating cell_type.csv file)
major_cell_type <- data.frame(c("epithelial","stroma","immune",
              "endothelial","unknown"))
colnames(major_cell_type) <- "cell_type_major"

# Create an output table that lists each patient and the percentage of
#   each cell type'

for (i in 1:number_of_patients)
{
# Extract each patient's data from the table and create a new table
temp_patient <- subset(cell_type, patient==names(patient_table[i]),
                    select = c(percent_of_total,cell_type_major))

# Create a new table by summing all percentages for each cell type
aggregated_table <- aggregate(percent_of_total~cell_type_major, 
                  data = temp_patient, FUN = "sum")

# Add the aggregated data to the major_cell_type output dataframe
major_cell_type <- merge(major_cell_type, aggregated_table, 
                         by = intersect(names(major_cell_type), names(aggregated_table)),
                         all = TRUE)
colnames(major_cell_type)[which(colnames(major_cell_type) == 'percent_of_total')] <- patient_folder_names[i]

# Create a color lookup table
color_table <- cbind(c("epithelial","immune","stroma","endothelial","unknown"),
                     c("blue","yellow","red","green","gray"))
color_table <- as.data.frame(color_table)
colnames(color_table) <- c("cell_type_major","color")

# Add the corresponding color to the aggregated_table
aggregated_table <- merge(aggregated_table,color_table, 
             by = intersect(names(aggregated_table), names(color_table)))

# Create the pie chart
filename <- paste("output/",patient_folder_names[i],
                  "_major_pie.png", sep = "")

png(filename, width = 500, height = 500)
pie(aggregated_table$percent_of_total, 
    labels = aggregated_table$cell_type_major, 
    main=names(patient_table[i]),
    col = as.vector(aggregated_table$color),
    cex.main = 2, cex = 1.5)
dev.off()
}

# Write the major_cell_type to a file
# First change the NA values to 0
major_cell_type[is.na(major_cell_type)] <- 0
write.csv(major_cell_type, file = "output/major_cell_type.csv")

#################################################################
# 3) Color tsne charts by major cell type

# Input the cell_type.csv file
cell_type <- read.csv("input/cell_type.csv", header = TRUE)

# Create a table of patients
patient_table <- table(cell_type$patient)
patient_table_length <- length(patient_table)

for (i in 1:number_of_patients)
{
# Extract each patient's data and create a new table
temp_patient <- subset(cell_type, patient==names(patient_table[i]),
                       select = c(cell_type_major,cluster_number))

# Create a color lookup table
color_table <- cbind(c("epithelial","immune","stroma","endothelial","unknown"),
                     c("blue","yellow","red","green","gray"))
color_table <- as.data.frame(color_table)
colnames(color_table) <- c("cell_type_major","color")

# Add the corresponding color to the temp_patient table
temp_patient <- merge(temp_patient,color_table, 
                          by = intersect(names(temp_patient), names(color_table)))

# Read in the cluster_assignments.csv file for each patient
filename <- paste(dir_1, patient_folder_names[i], dir_2, 
                  "cluster_assignments.csv", sep = "")
cluster_assignments <- read.csv(filename, header = TRUE)
colnames(cluster_assignments) <- c("barcode","cluster_number")

# Add the corresponding color to the cluster_assignments table
cluster_assignments <- merge(cluster_assignments,temp_patient, 
                      by = intersect(names(cluster_assignments), 
                                     names(temp_patient)))

# Read in the tsne_cell_embeddings.csv file for each patient
filename <- paste(dir_1, patient_folder_names[i], dir_2, 
                  "tsne_cell_embeddings.csv", sep = "")
tsne_embeddings <- read.csv(filename, header = TRUE)
colnames(tsne_embeddings) <- c("barcode","tSNE_1","tSNE_2")

# Merge the tsne_embeddings table and the cluster_assignments table
tsne_embeddings <- merge(tsne_embeddings,cluster_assignments, 
    by = intersect(names(tsne_embeddings),names(cluster_assignments)))

# Print out the tSNE plot
filename <- paste("output/",patient_folder_names[i],
                  "_major_tsne.png", sep = "")

chartname <- paste(patient_folder_names[i],
                  " major cell type tSNE", sep = "")

png(filename, width = 500, height = 500)
plot(tsne_embeddings$tSNE_1,tsne_embeddings$tSNE_2,
     xlab = "tSNE-1", ylab = "tSNE-2", pch = 20,
     col = as.vector(tsne_embeddings$color),
     main = chartname, cex.lab = 1.5, cex.main = 2)
dev.off()
}