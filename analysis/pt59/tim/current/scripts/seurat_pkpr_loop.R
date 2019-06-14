# Analyze prune, resolution and number_of_pcs parameters for 
#   their effect on cluster numbers
# Author: Tim Starr
# Last updated: 1/5/19

# Purpose: Analyze the effect of changing PCs, K-param, Prune and Resolution
#   (pcpr) on cluster number using the Seurat analysis package

# Input files: 
#   Patient.Robj. A Seurat S4 object created by the Seurat analysis package
#   I am using the Patient.Robj created by my seurat_filter_2.R script

# Input variables:
#   Resolution_variable: Using a fixed set of resolution values, 0.4 - 2.2 by 0.2 (n=10)
#   PruneSNN: Using a fixed set of values, 0.02 to 0.16 by 0.02 (n = 8)
#   dims_max: Significant PCs. Using a range of 8 consecutive significnat PCs 
#     based on the Jackstraw plot p-values returned seurat_filter_2.R
#   k-param: Using a range of 8 consecutvie values, with the minimum value 
#     based on estimated number of cells (see compilation spreadsheet).
#   Note: 10 x 8 x 8 x 8 variables = 5,120 tests


# Output files
#   At the end of each triple loop, three files will be writtin:
#     #-dims_#-k_res_prune.png: A cluster number x resolution graph with lines
#       for each prune value. One graph per dims_max x k-param value pair.
#     res_prune_matrix_dims_#_k_#.csv: The matrix of cluster number values for
#       a given dims_max x k-param value pair.
#     res_prune_count_dims_#_k_#.csv: A table of cluster numbers and their 
#       frequency for each dims_max x k-param value pair
#   After the looping has concluded, four summary files will be written:
#     cluster_number_graphs.gif: an animated .gif file of the .png files above
#     cluster_number_matrix: A summary table of the res_prune_matrix files
#     cluster_number_count: A summary table of the res_prune_count files
#     cluster_freq_graph.png: A line graph of the cluster_number_count table
#     

# Pseudo Code:
# Enter the four input vectors listed above
# Load the Patient.Robj
# Run FindNeighbors and FindClusters in a triple loop, each time varying
#   the dims_max, k-param, and prune values.  The Resolution values are
#   all computed in a single Seurat run
#   The final step within the loop will create three output files (see above)
# After the triple loop is complete, create three summary files (see above)

##################################################################
# Set the working directory
#   Use this working directory for running on my mac
# setwd('/Users/star0044old/Documents/Programming/R/Seurat')
#   Use this working directory for running on MSI (change patient)
setwd('~/ovarian/ssrnaseq/pat_59_b')
# Clear the workspace global environment
rm(list=ls())
# Load the required libraries
library(Seurat)
library(dplyr)
library(tibble)
library(Matrix)
library(rgl)
library(reshape2)
library(ggplot2)

##################################################################
# Set the dims_max
max_dims_max <- 50
dims_max <- seq(from = (max_dims_max-7), to = max_dims_max, by = 1)

# Set the k-param
min_k_param <- 38
k_param <- seq(from = min_k_param, to = (min_k_param + 5), by =1)

# Set the prune factors to use
prune_SNN <- seq(from = 0.02, to = 0.16, by = 0.02)

# Set the number of resolution factors to use (Note, cannot change
#   this unless you change script below also)
resolution_parameter <- seq(from = 0.4, to = 2.2, by = 0.2)


##################################################################
# Save the version of Seurat
seurat_version <- as.character(packageVersion("Seurat"))

####################################################
# Load data
load(file = "input/Patient.Robj")

###################################################
# Set the parameters (See Seurat Help for explanations)

# FindNeighbors parameters
nn_eps <- 0

# FindClusters parameters
modularity_fxn <- 1
algorithm_parameter <- 1
n_start <- 10
n_iter <- 10
############################
# Print out the parameters
parameter_headers <- c("seurat_version",
                       "nn_eps",
                       "modularity_fxn",
                       "algorithm_parameter",
                       "n_start",
                       "n_iter"
                       )
parameters <- c(seurat_version,
                nn_eps,
                modularity_fxn,
                algorithm_parameter,
                n_start,
                n_iter
                )
parameter_output <- cbind(parameter_headers,parameters)
write.csv(parameter_output, "output/loop_parameters.csv")

############################
# Loop through the different k-params

for(k in 1:length(k_param))
{

for(j in 1:length(dims_max))
  {
  # Create the output matrix
  output_matrix <- matrix(0, dimnames = list(prune_SNN,resolution_parameter), 
                          nrow = length(prune_SNN), ncol = length(resolution_parameter))
 
  for(i in 1:length(prune_SNN))
    {
    patient <- FindNeighbors(patient, reduction = "pca", 
                  dims = 1:dims_max[j], 
                  k.param = k_param[k], prune.SNN = prune_SNN[i],
                  nn.eps = nn_eps, verbose = FALSE, force.recalc = TRUE,
                  do.plot = FALSE)

    patient <- FindClusters(patient,
                          modularity.fxn = modularity_fxn, 
                          resolution = resolution_parameter, 
                          algorithm = algorithm_parameter,
                          n.start = n_start, n.iter = n_iter, random.seed = 0,
                          verbose = FALSE)

    resolution_comparison <- cbind(
                               patient@meta.data$RNA_snn_res.0.4,
                               patient@meta.data$RNA_snn_res.0.6,
                               patient@meta.data$RNA_snn_res.0.8,
                               patient@meta.data$RNA_snn_res.1,
                               patient@meta.data$RNA_snn_res.1.2,
                               patient@meta.data$RNA_snn_res.1.4,
                               patient@meta.data$RNA_snn_res.1.6,
                               patient@meta.data$RNA_snn_res.1.8,
                               patient@meta.data$RNA_snn_res.2,
                               patient@meta.data$RNA_snn_res.2.2
                              )
    # write.csv(resolution_comparison, "output_res_loop_2/resolution_comparison.csv")

    # Use the nlevels from each resolution to get cluster size
    cluster_numbers <- cbind(
                               nlevels(patient@meta.data$RNA_snn_res.0.4),
                               nlevels(patient@meta.data$RNA_snn_res.0.6),
                               nlevels(patient@meta.data$RNA_snn_res.0.8),
                               nlevels(patient@meta.data$RNA_snn_res.1),
                               nlevels(patient@meta.data$RNA_snn_res.1.2),
                               nlevels(patient@meta.data$RNA_snn_res.1.4),
                               nlevels(patient@meta.data$RNA_snn_res.1.6),
                               nlevels(patient@meta.data$RNA_snn_res.1.8),
                               nlevels(patient@meta.data$RNA_snn_res.2),
                               nlevels(patient@meta.data$RNA_snn_res.2.2)
                              )
    # Change to a data frame and add column and row names
    cluster_numbers <- as.data.frame(cluster_numbers)
    colnames(cluster_numbers) <- resolution_parameter
    rownames(cluster_numbers) <- prune_SNN[i]
    cluster_numbers <- as.numeric(cluster_numbers)
    output_matrix[i,] <- cluster_numbers
    }

  output_filename_2 <- paste("output_res_loop_2/res_prune_matrix_dims_", dims_max[j], "_k_", k_param[k], ".csv", sep = "")
  write.csv(output_matrix, output_filename_2)

  table_output <- as.data.frame(table(output_matrix))
  colnames(table_output) <- c("cluster_number","freq")
  output_filename_3 <- paste("output_res_loop_2/res_prune_count_dims_",dims_max[j],"_k_",k_param[k],".csv", sep = "")
  write.csv(table_output, output_filename_3)

  output_filename_4 <- paste("output_res_loop_2/", dims_max[j], "-dims_", k_param[k], "-k_res_prune.png", sep = "")
  png(filename = output_filename_4, width = 360, height = 288)
    # Change the output_matrix to a data frame
    plot_matrix <- as.data.frame(output_matrix)
    #  Create a first column with rownames
    plot_matrix <- tibble::rownames_to_column(plot_matrix,var = "prune")
    plot_matrix <- melt(plot_matrix, id.vars = "prune",
                             variable.name = "resolution")
    # Change prune to factor
    plot_matrix$prune <- as.factor(plot_matrix$prune)
    # Change resolution to numbers
    plot_matrix$resolution <- as.numeric(as.character(plot_matrix$resolution))
    # plot the data
    print(ggplot(plot_matrix, aes(x = resolution, y = value, colour = prune)) +
      geom_line(position = position_dodge(width = 0.2)) +
      ylab("Number of Clusters") +
      ylim(0,40) +
      ggtitle(output_filename_4) +
      theme_bw())
  dev.off()
  }
}

###########################################################################
# Create summary files of the three sets of files created above
# Create an empty output table to populate with the res_prune_count... files
output_table <- data.frame(cluster_number = integer(), freq = integer(),
                           dims = integer(), k_params = integer())

# Loop through the res_prune_count files and append to the output table
for (i in 1:length(dims_max))
{
  for (k in 1:length(k_param))
  {
    # Establish filename
    temp_filename <- paste('output_res_loop_2/res_prune_count_dims_', 
                           dims_max[i], '_k_', k_param[k],'.csv', sep = "")
    # Import file
    temp_table <- read.csv(temp_filename, header = TRUE)
    # Create the dim and k_param column and remove the X column
    temp_table$dim <- dims_max[i]
    temp_table$k_param <- k_param[k]
    temp_table$X <- NULL
    
    # append temp_table to output_table
    output_table <- rbind(output_table,temp_table)
  }
}

cluster_num_freq_totals <- aggregate(freq~cluster_number, output_table, sum)

# Save the output_table and the cluster_num_freq_totals
write.csv(output_table, "output_res_loop_2/cluster_number_matrix.csv")
write.csv(cluster_num_freq_totals, "output_res_loop_2/cluster_number_count.csv")

# Create a line graph
png(filename = "output_res_loop_2/cluster_freq_graph.png", 
    width = 360, height = 288)
print(ggplot(cluster_num_freq_totals, aes(x=cluster_num_freq_totals$cluster_number,
                                          y=cluster_num_freq_totals$freq)) +
        geom_line() +
        lims(x = c(0,40), y= c(0,1600)) +
        xlab("Cluster Number") +
        ylab("Cluster Frequency") +
        ggtitle("Cluster Number Frequency"))
dev.off()


