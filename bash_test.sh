#!/bin/bash

#SBATCH --output=job%A_%a.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=masel

#Mutationload variables
timeSteps=10
initialPopsize=2000
mud=2
chromosomesize=200
numberofchromosomes=23
bentodelratio=0
sb=1
#0 for point; 1 for exponential; and 2 for uniform
bendist=1
#0 for root sb; 1 for single; 2 for root Ncrit
typeofrun=1
slope=0
seed=24
K=20000
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
rawdatafilesize=10
#change in carrying capacity, type in the difference in popsize
redinmaxpopsize=0
#status of fixation calculation; 0 for OFF; 1 for ON
calcfixation=0


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
	deldiststring="kim_"
elif [ $deldist -eq 1 ]
then
	deldiststring="exponential_"
elif [ $deldist -eq 2]
then
	deldiststring="point_"
fi

#mub is written as a formated double in mutation load program
mub=$(echo "$mud * $bentodelratio" | bc -l)
mub=$(printf "%.4f" $mub)

#creates 2 strings; directory refers to the folder where data for the specified parameters will be stored; file is the snapshot of the simulation at its end.
if [ $modularepis -eq 0 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstring"r_"$r"_i_init"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"bendist_"$bendiststring"_mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"seed_"$seed"/"
elif [ $modularepis -eq 1 ]
then
	directory="datafor_"$fitnessstring"tskitstatus_"$tskitstring"elementsperlb_"$elementsperl"_r_"$r"_i_init"$i_init"_s_"$s"_K_"$K"_deldist_"$deldiststring"_bendist_"$bendiststring"_mub_"$mub"_chromnum_"$numberofchromosomes"_N0_"$initialPopsize"_mud_"$mud"_L_"$chromosomesize"seed_"$seed"/"
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
./mutationload $timeSteps $initialPopsize $mud $chromosomesize $numberofchromosomes $bentodelratio $sb $bendist $typeofrun $slope $seed $K $fitnesstype $r $i_init $s $tskitstatus $modularepis $elementsperl $snapshot $file1 $SdtoSbratio $deldist $rawdatafilesize $redinmaxpopsize $calcfixation

echo $SECONDS

echo "end of mutationload program"

#$snapshot $directory$file1

gzip $directory$file1
