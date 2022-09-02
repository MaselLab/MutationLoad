#!/bin/bash

# set library path for shared libraries
LD_LIBRARY_PATH=/home/wmawass/gsl/lib
export LD_LIBRARY_PATH

# run Makefile
make

#Mutationload variables
timeSteps=120000
initialPopsize=6000
mud=2.0
chromosomesize=50
numberofchromosomes=23
bentodelratio=0.001
sb=0.010
#0 for point; 1 for exponential; and 2 for uniform
bendist=1
#0 for root; 1 for single
typeofrun=1
slope=0
seed=24
K=7500
#0 for relative; 1 for absolute
fitnesstype=1
d0=0.2
r=0.98
sdmin=0.001
#0 for runs without modular epistasis; 1 for runs with modular epistasis
modularepis=0
elementsperl=1
#0 for no tskit; 1 for tskit on; 2 for tskit on after burnin
tskitstatus=2
SdtoSbratio=0.5
#0 for point; 1 for exponential
deldist=1
#rawdata file size in datapoints
rawdatafilesize=100

if [ $fitnesstype -eq 0 ]
then
	fitnessstring="relative_"
elif [ $fitnesstype -eq 1 ]
then
	fitnessstring="absolute_"
fi

if [ $bendist -eq 0 ]
then
	bendiststring="point_"
elif [ $bendist -eq 1 ]
then
	bendiststring="exponential_"
elif [ $bendist -eq 2 ]
then
	bendiststring="uniform_"
fi

if [ $tskitstatus -eq 0 ]
then
	tskitstatusstring="OFF_"
elif [ $tskitstatus -eq 1 ]
then
	tskitstatusstring="ON_"
elif [ $tskitstatus -eq 2 ]
then
	tskitstatusstring="ON_AFTER_BURNIN_"
fi

if [ $deldist -eq 0 ]
then
	deldiststring="point_"
elif [ $deldist -eq 1 ]
then
	deldiststring="exponential_"
fi

#mub is written as a formated double in mutation load program
mub=$(echo "$mud * $bentodelratio" | bc -l)
mub=$(printf "%.4f" $mub)

#creates 2 strings; directory refers to the folder where data for the specified parameters will be stored; file is the snapshot of the simulation at its end.
if [ $modularepis -eq 0 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstatusstring"d0_"$d0"r_"$r"sdmin_"$sdmin"K_"$K$deldiststring"mud_"$mud$bendiststring"mub_"$mub"popsize_"$initialPopsize"chromosomes_"$numberofchromosomes"L_"$chromosomesize"seed_"$seed"/"
elif [ $modularepis -eq 1 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstatusstring"elementsperlb_"$elementsperl"d0_"$d0"r_"$r"sdmin_"$sdmin"K_"$K$deldiststring"mud_"$mud$bendiststring"mub_"$mub"popsize_"$initialPopsize"chromosomes_"$numberofchromosomes"L_"$chromosomesize"seed_"$seed"/"
fi

printf "directory path is %s \n" "$directory"

file1="popsnapshotfor_popsize_"$initialPopsize"_tskitstatus_"$tskitstatusstring"mub_"$mub".txt"

#checks if a previous snapshot of the simulation exist. Snapshots are saved as compressed files (.gz) to save space
prevsim=$([ -f $directory$file1".gz" ] && echo 1 || echo 0)

if [ $prevsim -eq 0 ]
then
	snapshot=0
elif [ $prevsim -eq 1 ]
then
	gzip -d $directory$file1".gz"
	snapshot=1
fi


SECONDS=0
echo "start of mutationload program"

# run mutationload program with arguments
./mutationload $timeSteps $initialPopsize $mud $chromosomesize $numberofchromosomes $bentodelratio $sb $bendist $typeofrun $slope $seed $K $fitnesstype $d0 $r $sdmin $tskitstatus $modularepis $elementsperl $snapshot $file1 $SdtoSbratio $deldist $rawdatafilesize

echo $SECONDS

echo "end of mutationload program"

#$snapshot $directory$file1

gzip $directory$file1