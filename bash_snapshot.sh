#!/bin/bash

timeSteps=10
initialPopsize=20000
mud=0.6
chromosomesize=50
numberofchromosomes=23
bentodelratio=0.0
sb="0.0"
#0 for point; 1 for exponential; and 2 for uniform
bendist=1
#0 for root; 1 for single
typeofrun=1
slope=0
seed=57
K=25000
#0 for relative; 1 for absolute
fitnesstype=1
d0=0.2
r=0.8
sdmin=0.01
#0 for runs without modular epistasis; 1 for runs with modular epistasis
modularepis=1
elementsperl=1


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

#mub is written as a formated double in mutation load program
mub=$(echo "$mud * $bentodelratio" | bc -l)
mub=$(printf "%.4f" $mub)

#creates 2 strings; directory refers to the folder where data for the specified parameters will be stored; file is the snapshot of the simulation at its end.
if [ $modularepis -eq 0 ]
then
	directory="datafor_"$fitnessstring"d0_"$d0"r_"$r"sdmin_"$sdmin"K_"$K"mud_"$mud$bendiststring"mub_"$mub"popsize_"$initialPopsize"chromosomes_"$numberofchromosomes"L_"$chromosomesize"totalt_"$timeSteps"seed_"$seed"/"
elif [ $modularepis -eq 1 ]
then
	directory="datafor_"$fitnessstring"elementsperlb_"$elementsperl"d0_"$d0"r_"$r"sdmin_"$sdmin"K_"$K"mud_"$mud$bendiststring"mub_"$mub"popsize_"$initialPopsize"chromosomes_"$numberofchromosomes"L_"$chromosomesize"totalt_"$timeSteps"seed_"$seed"/"
fi

file1='popsnapshotfor_Sb_'$sb'mub_'$mub".txt"

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

./mutationload $timeSteps $initialPopsize $mud $chromosomesize $numberofchromosomes $bentodelratio $sb $bendist $typeofrun $slope $seed $K $fitnesstype $d0 $r $sdmin $modularepis $elementsperl $snapshot $file1

#$snapshot $directory$file1

gzip $directory$file1


