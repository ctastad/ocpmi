#!/bin/bash
# This script is meant to prep the configuration and execution of the Excavator2 pipeline. It first auto generates two text files that will serve as configs for the 2nd and 3rd pipeline scripts. It then creates directories for each sample at the prepDir site.

source vars.sh
cd $HOME$bamDir

# Overwrite any existing config files
> ExperimentalFilePrepare.txt
> ExperimentalFileAnalysis.txt

# Loop to iterate through each bam to create string for FilePrepare and FileAnalysis configs
count=1
for i in *.{bam,cram}; do
    echo "$PWD/$i" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFilePrepare.txt
    echo "T$count" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFileAnalysis.txt
    mkdir -p $HOME$excv2Dir/data/prepare/$expName/${i%.*}
    count=$[count+1]
done

# Function to remove leading path from sshfs local mount path and create control
sed -i -e 's/\/home\/chris/~/g' ./ExperimentalFilePrepare.txt
sed -i '$s/T[[:digit:]]*/C1/g' ./ExperimentalFileAnalysis.txt

# Relocate configs to Excavator working dir
mv ExperimentalFile* $HOME$excv2Dir
mkdir $HOME$resDir/$expName

# Generate SourceTarget.txt and execute step 1 TargetPerla.pl
cd $HOME$excv2Dir
if [ "$runTp"="1" ]; then
    echo "$HOME$bigWig" "$HOME$genRef" > SourceTarget.txt
    perl TargetPerla.pl SourceTarget.txt $HOME$bedFile $target $window $refAssm
fi

# Executing steps 2 and 3 of the Excavator pipeline
perl EXCAVATORDataPrepare.pl ExperimentalFilePrepare.txt --processors $numProc --target $target --assembly $refAssm
perl EXCAVATORDataAnalysis.pl ExperimentalFileAnalysis.txt --processors $numProc --target $target --assembly $refAssm --output $HOME$resDir/$expName --mode $anMode



