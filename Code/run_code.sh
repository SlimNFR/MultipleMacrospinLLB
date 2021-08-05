#This script is used to run the multispin code, save its output & sim files, plot and clean the current dir

#Global variables

##Paths
	plotScripts_Dir="/home/slimnfr/Desktop/Multiple-Macrospins_LLB-code/GnuplotScripts"
	outDir="ForceDW_by_fixingEdgespinsOn_+me_z_-me_z_40_spins_300K_Ni_3200000_steps"
	outTxt="$outDir/Output"
	inTxt="$outDir/Input"
	outFigs="$outDir/Graphs"

##Material/Simulation
	N=$1 #Number of macrospins
	timeStep=$2 #plot points every timeStep 

#!/bin/bash

#Create output directory

if [ -d $outDir ]; then
echo "$outDir exists!"
exit
else
mkdir $outDir
mkdir "$outFigs"
mkdir "$outTxt"
mkdir "$inTxt"
echo "$outDir and its subfolders were created!"
fi

#Run executable
./a.out


#Plot

gnuplot -c "$plotScripts_Dir/plot_M_and_T_time.gnu" "$N" "$timeStep"

#Save output and clean

cp matInput.txt $inTxt
cp simInput.txt $inTxt
mv output_* $outTxt
mv *.png $outFigs


