#!/bin/bash
#SBATCH --job-name=mutationload
#SBATCH --output=job_%A_%a.out
#SBATCH --account=masel
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --mem=5gb                    
#SBATCH --time=07-00:00:00   
#SBATCH --partition=standard
#SBATCH --array=0-6
 
# SLURM Inherits your environment. cd $SLURM_SUBMIT_DIR not needed
pwd; hostname; date

#module load gsl
module restore default

cd ~/MutationLoad
# run Makefile
# make whereisrun=HPC

#array of values
arraytimeSteps=(100000 200000 300000 400000 500000 600000 700000)
arrayN=(1000 2000 3000 4000 5000 6000 7000)
#arraymud=(0.5 1.0 2.0 5.0 10.0)
#arrayratio=(0.002 0.002)

#Mutationload variables
timeSteps=${arraytimeSteps[${SLURM_ARRAY_TASK_ID}]}
initialPopsize=${arrayN[${SLURM_ARRAY_TASK_ID}]}
mud=2.0
chromosomesize=55
numberofchromosomes=22
bentodelratio=0.001
sb=0.005
#0 for point; 1 for exponential; and 2 for uniform
bendist=1
#0 for root sb; 1 for single; 2 for root Ncrit
typeofrun=1
slope=0
seed=24
K=10000
#0 for relative; 1 for absolute
fitnesstype=0
r=0.98
i_init=500
s=0.01
#0 for runs without modular epistasis; 1 for runs with modular epistasis
modularepis=0
elementsperl=0
#0 for no tskit; 1 for tskit on; 2 for tskit on after burnin
tskitstatus=1
SdtoSbratio=1.0
#0 for point; 1 for exponential
deldist=1
#rawdata file size in datapoints
rawdatafilesize=100
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
