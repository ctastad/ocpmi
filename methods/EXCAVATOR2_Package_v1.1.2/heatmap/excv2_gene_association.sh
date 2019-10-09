#!/bin/bash

################################################################################
#
# The Results dir should be used as the primary argument when executing this
# script. The name of the figure should be the second argument.
#
################################################################################

# module load bedtools R/3.6.0

# Establish path for function working dir
dir=$HOME/ocpmi/results/excavator2/$1/Results
cd $dir
# geneAssociation scripts location
scriptDir=$HOME/ocpmi/methods/EXCAVATOR2_Package_v1.1.2/heatmap

for i in */; do
    cd $i
    sampleName=$(pwd | awk -F'/' '{print $NF}')
    bedtools intersect \
        -a EXCAVATORRegionCall_* \
        -b $HOME/ocpmi/data/gene_ref/biomart/hg38/ensembl_genes.bed \
        > ${sampleName}_call_subset
        # -b $HOME/ocpmi/data/gene_ref/ucsc/hg38/hg38.bed \
    echo "Created call subset file for" "$sampleName"
    bedtools intersect \
        -a $HOME/ocpmi/data/gene_ref/biomart/hg38/ensembl_genes.bed \
        -b EXCAVATORRegionCall_* \
        > bed_subset.txt
        # -a $HOME/ocpmi/data/gene_ref/ucsc/hg38/hg38.bed \
    echo "Created bed subset file for" "$sampleName"
    echo "Running ETL"
    $scriptDir/geneAssociation_etl.R
    echo "Doing cleanup"
    mv etl_output* $dir
    rm *_call_subset
    rm bed_subset.txt
    cd $dir
done

echo "Finished creating sample files"

echo "Creating final table"
cat *.txt > geneAssociationTable
rm *.txt
sed -i -e 's/_Aligned.sortedByCoord.out_call_subset//g' geneAssociationTable

echo "Graphing output"
$scriptDir/geneAssociationGraph.R $1 $2

echo "All done"
