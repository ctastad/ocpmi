################################################################################
# This is the config file for full_send.sh
#
# Use of this config file should be paired with full_send.sh.
# Both should be placed in the Excavator2 working directory and executed there.
#
# After setting all variables, the function can be executed with the following
# command in the Excavator2 working directory
#
# ./full_send.sh
#
# In addition to the output generated in the results directory, this function
# will move all config files to a directory it creates in the Excavator2 working
# directory labeled config_files.
################################################################################


################################################################################
# Analysis variables
################################################################################
# Experiment name
expName="ocpmi_v_1kgp_full_50k1"
# Window size
window="50000"
# Number of processors for parallelization
numProc="8"
# Step 3 EXCAVATORDataAnalysis.pl sample analysis mode [pooling, paired]
anMode="pooling"

################################################################################
# Gene reference target files
################################################################################
# Reference assembly to use [hg19,hg38]
refAssm="hg38"
# Location of genome reference fasta
genRef="/ocpmi/data/gene_ref/ucsc/hg38/hg38.fa"
# Location of bed file
bedFile="/ocpmi/data/gene_ref/ucsc/hg38/hg38_sortmerge.bed"
# Location of BigWig file
bigWig="/ocpmi/methods/EXCAVATOR2_Package_v1.1.2/data/GCA_000001405.15_GRCh38.bw"

################################################################################
# Location of working and data directories
################################################################################
# Location of the Excavator working directory
excv2Dir="/ocpmi/methods/EXCAVATOR2_Package_v1.1.2"
# Location of the bam files
bamDir="/ocpmi/data/bam/wes"
# Location of results directory [format will be /path/to/dir/$expName]
resDir="/ocpmi/results/excavator2"

