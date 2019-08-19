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
if [[ $runTp1 = 1 ]]; then
    echo "$HOME$bigWig" "$HOME$genRef" > SourceTarget.txt
    perl TargetPerla.pl SourceTarget.txt $HOME$bedFile $target $window $refAssm;
else
    echo "Skipping step 1 TargetPerla.pl"
fi

# Executing step 2 EXCAVATORDataPrepare.pl
if [[ $runDp2 = 1 ]]; then
    perl EXCAVATORDataPrepare.pl -v ExperimentalFilePrepare.txt --processors $numProc --target $target --assembly $refAssm
else
    echo "Skipping step 2 EXCAVATORDataPrepare.pl"
fi

# Executing step 3 EXCAVATORDataAnalysis.pl
if [[ $runDa3 = 1 ]]; then
    perl EXCAVATORDataAnalysis.pl -v ExperimentalFileAnalysis.txt --processors $numProc --target $target --assembly $refAssm --output $HOME$resDir/$expName --mode $anMode
else
    echo "Skipping step 3 EXCAVATORDataAnalysis.pl"
fi

