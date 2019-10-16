#!/bin/bash

################################################################################
#
# The Results dir name should be used as the primary argument when executing
# this script. The name of the figure should be the second argument.
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
        -b $HOME/ocpmi/data/gene_ref/tim/master_gene_list.bed \
        > ${sampleName}_call_subset
    echo "Created call subset file for" "$sampleName"
    bedtools intersect \
        -a $HOME/ocpmi/data/gene_ref/tim/master_gene_list.bed \
        -b EXCAVATORRegionCall_* \
        > bed_subset.txt
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

echo "Graphing unclustered output"
$scriptDir/geneAssociationGraph.R $1 $2
echo "Graphing clustered output"
$scriptDir/geneAssociationGraph_clust.R $1 $2

echo "All done"
