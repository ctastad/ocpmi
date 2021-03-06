#!/bin/bash
# This script is meant to prep the configuration and execution of the Excavator2 pipeline. It auto generates two text files that serve as configs for the pipeline scripts. Proper use of this script is paired with the full_send_vars.txt variable file to indicate pointers and variable names. Both should be placed in the Excavator2 working directory with this script being executed from there.

source full_send_vars.txt
cd $HOME$bamDir

## Overwrite existing config files and updating some vars
rm -f ExperimentalFile*
expName="${expName}_${window}_`date +"%Y-%m-%d_%H%M"`"
target="${expName}_target"
# runTp1=1 # step 1 TargetPerla.pl
# runDp2=1 # step 2 EXCAVATORDataPrepare.pl
# runDa3=1 # step 3 EXCAVATORDataAnalysis.pl

## Loop to iterate through each bam to create string for FilePrepare and FileAnalysis configs
count=1
for i in *.bam; do
    echo "$PWD/$i" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFilePrepare.$expName.txt
    echo "T$count" "data/prepare/$expName/${i%.*}" "${i%.*}" >> ExperimentalFileAnalysis.$expName.txt
    mkdir -p $HOME$excv2Dir/data/prepare/$expName/${i%.*}
    count=$[count+1]
done

## Function to remove leading path from sshfs local mount path and create control
sed -i -e 's/\/home\/chris/~/g' ./ExperimentalFilePrepare.$expName.txt
sed -i -e '$s/T[0-9]* /C1 /g' ./ExperimentalFileAnalysis.$expName.txt

## Relocate configs to Excavator working dir
mv ExperimentalFile* $HOME$excv2Dir
mkdir -p $HOME$resDir/$expName
mkdir -p $HOME$excv2Dir/config_files/$expName

## Generate SourceTarget.txt and execute step 1 TargetPerla.pl

# if [[ $runTp1 = 1 ]]; then
#     echo "$HOME$bigWig" "$HOME$genRef" > SourceTarget.$expName.txt
#     perl TargetPerla.pl SourceTarget.$expName.txt $HOME$bedFile $target $window $refAssm;
# else
#     echo "Skipping step 1 TargetPerla.pl"
# fi

cd $HOME$excv2Dir
cp ParameterFile.txt ParameterFile.$expName.txt
cp full_send_vars.txt full_send_vars.$expName.txt

## Generate SourceTarget file and run step 1
echo "$HOME$bigWig" "$HOME$genRef" > SourceTarget.$expName.txt
perl TargetPerla.pl -v SourceTarget.$expName.txt $HOME$bedFile $target $window $refAssm;

## Executing step 2 EXCAVATORDataPrepare.pl
perl EXCAVATORDataPrepare.pl -v ExperimentalFilePrepare.$expName.txt --processors $numProc --target $target --assembly $refAssm

## Executing step 3 EXCAVATORDataAnalysis.pl
perl EXCAVATORDataAnalysis.pl -v ExperimentalFileAnalysis.$expName.txt --processors $numProc --target $target --assembly $refAssm --output $HOME$resDir/$expName --mode $anMode

## Copy config files for record
mv *.$expName.* $HOME$excv2Dir/config_files/$expName

## Cleanup bad file naming convention
cd $HOME$resDir/$expName/Plots
for i in */ ; do
    cd $i
    for j in {1..9}; do
        mv PlotResults_chr$j.pdf PlotResults_chr0$j.pdf
    done
    cd ../
done
