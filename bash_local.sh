#!/bin/bash

#Mutationload variables
timeSteps=500000
initialPopsize=20000
mud=2.0
chromosomesize=50
numberofchromosomes=23
bentodelrateratio=0.0
sb="0.0"
#0 for point; 1 for exponential; and 2 for uniform
bendist=0
#0 for root; 1 for single
typeofrun=1
slope=0
seed=57
K=20000
#0 for relative; 1 for absolute
fitnesstype=1
r=0.98
i_init=400
s=0.01
tskit=0
#0 for relative; 1 for absolute
fitnesstype=0
r=0.98
i_init=400
s=0.01
#0 for runs without modular epistasis; 1 for runs with modular epistasis
modularepis=0
elementsperl=0
#0 for no tskit; 1 for tskit on; 2 for tskit on after burnin
tskitstatus=2
SdtoSbratio=0.029
#0 for Kim et al., 1 for exponential, 2 for point
deldist=1
#rawdata file size in datapoints
rawfilesize=10000
redinK=0


if [ $fitnesstype -eq 0 ]
then
	fitnessstring="relative_"
elif [ $fitnesstype -eq 1 ]
then
	fitnessstring="absolute_"
fi

if [ $tskit -eq 0 ]
then
    	tskitstring="OFF_"
elif [ $tskit -eq 1 ]
then
    	tskitstring="ON_"
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

if [ $deldist -eq 0 ]
then
	deldiststring="kim_"
elif [ $deldist -eq 1 ]
then
	deldiststring="exponential_"
elif [ $deldist -eq 2]
then
	deldiststring="point_"
fi

#mub is written as a formated double in mutation load program
mub=$(echo "$mud * $bentodelrateratio" | bc -l)
mub=$(printf "%.4f" $mub)

#creates 2 strings; directory refers to the folder where data for the specified parameters will be stored; file is the snapshot of the simulation at its end.
if [ $modularepis -eq 0 ]
then
	directory="datafor_"$fitnessstring"tskit_"$tskitstring"r_"$r"_iinit_"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"bendist_"$bendiststring"mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"_seed_"$seed"/"
elif [ $modularepis -eq 1 ]
then
	directory="datafor_"$fitnessstring"tskit_"$tskitstring"elementsperlb_"$elementsperl"_r_"$r"_iinit_"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"bendist_"$bendiststring"mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"_seed_"$seed"/"
fi

file1='popsnapshotfor_Sb_'$Sbname'mub_'$mub".txt"

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

./mutationload $timeSteps $initialPopsize $mud $chromosomesize $numberofchromosomes $bentodelrateratio $sb $bendist $typeofrun $slope $seed $K $fitnesstype $r $i_init $s $tskit $modularepis $elementsperl $snapshot $file1 $SdtoSbratio $deldist $rawfilesize $redinK

#$snapshot $directory$file1
echo $directory$file1

echo $SECONDS

echo "end of mutationload program"

gzip $directory$file1
