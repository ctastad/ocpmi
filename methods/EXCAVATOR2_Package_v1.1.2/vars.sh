# This is the config file for full_send.sh
#

# Experiment name
expName="ocpmi_1kgp_compare"

# Choose to run step 1 TargetPerla.pl [1,0] (not required if desired targets already built)
runTp="1"

# The target name as output from the TargetPerla.sh script
target="hg38_target"

# Reference assembly to use [hg19,hg38]
refAssm="hg38"

# Window size
window="50000"

# Location of bed file
bedFile="/ocpmi/data/gene_ref/ucsc/hg38/hg38_sortmerge.bed"

# Location of BigWig file
bigWig="/ocpmi/methods/EXCAVATOR2_Package_v1.1.2/data/GCA_000001405.15_GRCh38.bw"

# Location of genome reference fasta
genRef="/ocpmi/data/gene_ref/ucsc/hg38/hg38.fa"

# Number of processors
numProc="8"

# Step 3 EXCAVATORDataAnalysis.pl sample analysis mode [pooling, paired]
anMode="pooling"

# Location of the Excavator working directory
excv2Dir="/ocpmi/methods/EXCAVATOR2_Package_v1.1.2"

# Location of the bam files
bamDir="/ocpmi/data/bam/wes/pilot/ocmpi_1kgp_compare"

# Location of results directory [format will be /path/to/dir/$expName]
resDir="/ocpmi/results/excavator2"
