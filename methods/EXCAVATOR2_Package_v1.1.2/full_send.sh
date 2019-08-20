#!/bin/bash
# This script is meant to prep the configuration and execution of the Excavator2 pipeline. It first auto generates two text files that will serve as configs for the 2nd and 3rd pipeline scripts. It then creates directories for each sample at the prepDir site.

source full_send_vars.sh
cd $HOME$bamDir

# Overwrite any existing config files
rm -f ExperimentalFile*
target="${expName}_target"
runTp1=1 # step 1 TargetPerla.pl
runDp2=1 # step 2 EXCAVATORDataPrepare.pl
runDa3=1 # step 3 EXCAVATORDataAnalysis.pl

# Loop to iterate through each bam to create string for FilePrepare and FileAnalysis configs
count=1
for i in *.{bam,cram}; do
    echo "$PWD/$i" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFilePrepare.$expName.txt
    echo "T$count" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFileAnalysis.$expName.txt
    mkdir -p $HOME$excv2Dir/data/prepare/$expName/${i%.*}
    count=$[count+1]
done

# Function to remove leading path from sshfs local mount path and create control
sed -i -e 's/\/home\/chris/~/g' ./ExperimentalFilePrepare.$expName.txt
sed -i '$s/T[[:digit:]]*/C1/g' ./ExperimentalFileAnalysis.$expName.txt

# Relocate configs to Excavator working dir
mv ExperimentalFile* $HOME$excv2Dir
mkdir -p $HOME$resDir/$expName
mkdir -p $HOME$excv2Dir/config_files/$expName

# Generate SourceTarget.txt and execute step 1 TargetPerla.pl
cd $HOME$excv2Dir
if [[ $runTp1 = 1 ]]; then
    echo "$HOME$bigWig" "$HOME$genRef" > SourceTarget.$expName.txt
    perl TargetPerla.pl SourceTarget.$expName.txt $HOME$bedFile $target $window $refAssm;
else
    echo "Skipping step 1 TargetPerla.pl"
fi

# Executing step 2 EXCAVATORDataPrepare.pl
if [[ $runDp2 = 1 ]]; then
    perl EXCAVATORDataPrepare.pl -v ExperimentalFilePrepare.$expName.txt --processors $numProc --target $target --assembly $refAssm
else
    echo "Skipping step 2 EXCAVATORDataPrepare.pl"
fi

# Executing step 3 EXCAVATORDataAnalysis.pl
if [[ $runDa3 = 1 ]]; then
    perl EXCAVATORDataAnalysis.pl -v ExperimentalFileAnalysis.$expName.txt --processors $numProc --target $target --assembly $refAssm --output $HOME$resDir/$expName --mode $anMode
else
    echo "Skipping step 3 EXCAVATORDataAnalysis.pl"
fi

# Copy config files for record
cp ParameterFile.txt ParameterFile.$expName.txt
cp full_send_vars.sh full_send_vars.$expName.sh
mv *.{$expname}.* $HOME$excv2Dir/config_files/$expName

# Cleanup bad file naming convention
cd $HOME$resDir/$expName/Plots
for i in */ ; do
    cd $i
    for j in {1..9}; do
        mv PlotResults_chr$j.pdf PlotResults_chr0$j.pdf
    done
    cd ../
done
