#!/bin/bash

################################################################################
#
# The Results dir name should be used as the primary argument when executing
# this script.
#
################################################################################

module load bedtools R/3.6.0

# establish path for function working dir
resDir=$HOME/ocpmi/results/excavator2/$1/Results
cd $resDir

# gene intersect scripts location
scriptDir=$HOME/ocpmi/methods/EXCAVATOR2_Package_v1.1.2/results_trim

for i in */; do
    cd $i
    sampleName=$(pwd | awk -F'/' '{print $NF}')
    bedtools intersect -wb \
        -a $HOME/ocpmi/data/gene_ref/tim/master_gene_list.bed \
        -b EXCAVATORRegionCall_* \
        > ${sampleName}_intersect
    echo "Created call intersect file for" "$sampleName"
    echo "Running ETL"
    Rscript $scriptDir/gene_intersect_etl.R
    echo "Doing cleanup"
    mv etl_output* $resDir
    rm *_intersect
    cd $resDir
done

echo "Finished creating sample files"

echo "Creating final table"
echo "chrom,gene_s,gene_e,ensembl_id,gene,cnv_s,cnv_e,cnv_length,cnv_call,pt_id" \
    > header_output_.csv
cat header_output_.csv etl_output* > gene_intersect.csv
rm *_output_*.csv
# doing text clean up
sed -i -e 's/_Aligned.sortedByCoord.out_intersect//g' \
    -i -e 's/SVLEN=//g'\
    -i -e 's/END=//g' \
    gene_intersect.csv
cat gene_intersect.csv | tr "," "\\t" > gene_intersect.bed

echo "Creating intersection with regions"
bedtools intersect -wb \
    -a gene_intersect.bed \
    -b $HOME/ocpmi/data/gene_ref/tim/cytoBand_hg38.bed \
    > region_intersect.bed

echo "All done"
