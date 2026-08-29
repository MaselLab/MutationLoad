#!/bin/bash
#
# =============================================================================
#  MutationLoad - evolution of the mutation rate
#  General SLURM submission script for the University of Arizona HPC
#  (Puma / Ocelote / El Gato)
# =============================================================================
#
#  QUICK START
#     sbatch slurm_mutationload.sh                  # runs the whole array below
#     sbatch --array=0-4 slurm_mutationload.sh      # override the array range
#     squeue -u $USER                               # check on it
#     scancel <jobid>                               # kill it
#
#  WHAT THIS SCRIPT DOES
#     Runs one independent simulation per array task. Each task gets its own
#     output directory under $OUTROOT, so array tasks can never overwrite each
#     other's files even when they differ only in a parameter that does not
#     appear in the directory name.
#
#  HOW TO SWEEP A PARAMETER
#     Put the values you want in the sweep array below and set --array to match
#     (0 to N-1). The default sweeps the random seed, i.e. it produces
#     independent replicates of the same parameter set. Commented examples show
#     how to sweep the mutator parameters instead, or two parameters at once.
#
# =============================================================================
#  SBATCH DIRECTIVES  -  edit these first
# =============================================================================

#SBATCH --job-name=mutationload
#SBATCH --output=logs/mutationload_%A_%a.out
#SBATCH --error=logs/mutationload_%A_%a.err

#SBATCH --account=masel
#SBATCH --partition=standard
# For jobs that should not consume your group's monthly allocation, swap the
# partition for windfall (lower priority, pre-emptible). Check the current UA
# HPC documentation for the exact windfall settings before relying on it.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# The simulation is strictly single-threaded. Asking for more cores will not
# make it faster; it will only make the job wait longer in the queue.

#SBATCH --mem=8gb
# MEMORY SIZING. The dominant cost is the population itself:
#
#     bytes ~= popsize * numberofchromosomes * 2 * chromosomesize * B
#
#     B = 12  in global-mask mode      (8 for the fitness double + 4 for the state int)
#     B = 13  in inherited-mask mode   (+1 for the per-individual modifier mask byte)
#
# For the defaults below (popsize 20000, 23 chromosomes, 200 blocks) that is
# about 2.2 GB (global) or 2.4 GB (inherited). Tree-sequence recording adds more
# on top and grows with run length, so 8 GB is a safe starting point with
# tskitstatus 0 or 2. Raise it if you turn tskit on from generation 0.

#SBATCH --time=220:00:00
# Walltime. A 20000 x 20000 run takes a long time; check the standard
# partition's current maximum before increasing this.

#SBATCH --array=0-9
# One task per entry in the sweep array below. MUST be 0 to (number of entries - 1).

# Optional: email yourself when the job ends or dies.
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=YOUR_NETID@arizona.edu

# =============================================================================
#  ENVIRONMENT
# =============================================================================

set -u   # abort on an undefined variable rather than silently passing an empty
         # argument, which would shift every later positional argument

# At RUN time the job needs the GSL runtime library, not a compiler - the binary
# is already built by then. So load gsl here if your cluster provides it as a
# module, and do NOT load a compiler module: on this cluster "module load gcc"
# fails with "The following module(s) are unknown: gcc". Because there is no
# set -e, that failure does not stop the job, but it does clutter the .err file.
# Run "module avail gsl" to see what is available; leave this commented out if
# GSL is already on your default path (which is how the existing
# bash_hpc_array*.sh scripts in this repository ran).
# module load gsl

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
EXE="$SUBMITDIR/mutationload"

# Array index. Defaults to 0 so the script also works when submitted without
# --array, or when run directly for a quick test.
TASKID="${SLURM_ARRAY_TASK_ID:-0}"

mkdir -p "$SUBMITDIR/logs"

if [ ! -x "$EXE" ]; then
    echo "ERROR: $EXE not found or not executable. Build it first with 'make'." >&2
    exit 1
fi

# =============================================================================
#  THE SWEEP
# =============================================================================
# Default: 10 independent replicates of one parameter set, differing only in the
# random seed. Now that the GSL generator is seeded properly (see main.c), these
# really are independent.

seedarray=(101 102 103 104 105 106 107 108 109 110)
seed=${seedarray[$TASKID]}

# ---- Alternative sweeps: comment out the two lines above and use one of these.
#
# Sweep the mutator strength factor f (--array=0-5):
#   msfarray=(1.0 1.001 1.01 1.05 1.1 1.5)
#   mutator_strength_factor=${msfarray[$TASKID]}
#   seed=101
#
# Sweep the fraction of blocks carrying a modifier locus (--array=0-4):
#   mlfarray=(0.001 0.005 0.01 0.05 0.1)
#   modifier_locus_fraction=${mlfarray[$TASKID]}
#   seed=101
#
# Two parameters at once - 4 values of f x 5 seeds = 20 tasks (--array=0-19):
#   msfarray=(1.0 1.01 1.05 1.1)
#   seedarray=(101 102 103 104 105)
#   mutator_strength_factor=${msfarray[$(( TASKID / 5 ))]}
#   seed=${seedarray[$(( TASKID % 5 ))]}

# =============================================================================
#  SIMULATION PARAMETERS
# =============================================================================
# Anything set by the sweep above is guarded with ${VAR:=default} so the sweep
# wins and the default applies otherwise.

#--- general ---------------------------------------------------------------
timeSteps=20000
initialPopsize=20000
mud=2.1
chromosomesize=200
numberofchromosomes=23
bentodelratio=0
sb=1
#0 for point; 1 for exponential; 2 for uniform
bendist=1
#0 for root sb; 1 for single run; 2 for root Ncrit
typeofrun=1
#0 for no tskit; 1 for tskit from generation 0; 2 for tskit only after burn-in
tskitstatus=2
SdtoSbratio=0.029
#0 for Kim et al.; 1 for exponential; 2 for point
deldist=1
#0 for relative fitness; 1 for absolute (absolute runs are disabled in this build)
fitnesstype=0
: "${seed:=101}"

#--- absolute-fitness parameters (unused by relative runs, but they occupy
#    positional argument slots and must still be supplied) -------------------
slope=0
K=20000
r=0.98
i_init=400
s=0.01
rawdatafilesize=10
redinmaxpopsize=0
calcfixation=0

#--- modular epistasis: NOT SUPPORTED, both must stay 0 ---------------------
modularepis=0
elementsperl=0

#--- mutation-rate evolution ------------------------------------------------
# mu_deleterious = mud * f^n   and   mu_beneficial = mub * f^n,
# where n is the summed modifier state over the whole diploid genome.
# See sharedfunc_flag.h and general_bash_local.sh for the full explanation.

#f; 1.0 means mutator alleles have no effect on the mutation rate
: "${mutator_strength_factor:=1.0}"
#per-MODIFIER-LOCUS, per-gamete switch probability; 0 disables modifier evolution
: "${mutator_switch_rate:=0.0}"
#multiplier on the anti-mutator -> mutator rate; >1 favours mutators
: "${mutator_bias:=1.0}"
#p: fraction of linkage blocks carrying a modifier locus
: "${modifier_locus_fraction:=0.01}"
#0 = one global fixed mask; 1 = per-individual inherited mask
: "${modifier_mask_mode:=0}"
#0 = anti-mutator stored as 0 (n = mutator count)
#1 = anti-mutator stored as -1 (n = net sum, so mu can fall below mud)
: "${antimutator_encoding:=0}"
#q: fraction of modifier loci starting in the +1 state
: "${initial_mutator_fraction:=0.0}"

#--- optional per-individual tracking ---------------------------------------
# EXPENSIVE: each firing writes $initialPopsize rows. The population-level means
# and variances go to the raw data file regardless of this setting.
: "${trackindividuals:=0}"
: "${trackinterval:=100}"
: "${trackstartgen:=1}"

# =============================================================================
#  OUTPUT LOCATION
# =============================================================================
# Each array task gets its own directory. The program creates a further
# parameter-named subdirectory inside it and writes all of its output there.
#
# On UA HPC, /xdisk is the right place for large output. Point OUTROOT there if
# these runs will produce a lot of data, e.g.
#   OUTROOT=/xdisk/masel/$USER/mutationload
OUTROOT="${OUTROOT:-$SUBMITDIR/results}"
WORKDIR="$OUTROOT/job${SLURM_JOB_ID:-local}_task${TASKID}_seed${seed}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || { echo "ERROR: cannot cd to $WORKDIR" >&2; exit 1; }

# The snapshot argument pair. The snapshot/restart workflow belongs to the
# absolute-fitness runs and is not used here, so start fresh every time.
snapshot=0
file1="popsnapshotfor_popsize_${initialPopsize}_seed_${seed}.txt"

# =============================================================================
#  RUN
# =============================================================================

echo "==============================================================="
echo " job          : ${SLURM_JOB_ID:-local}  task ${TASKID}"
echo " host         : $(hostname)"
echo " started      : $(date)"
echo " workdir      : $WORKDIR"
echo " seed         : $seed"
echo " popsize      : $initialPopsize   timeSteps: $timeSteps   mud: $mud"
echo " modifier     : f=$mutator_strength_factor switch=$mutator_switch_rate bias=$mutator_bias"
echo "                p=$modifier_locus_fraction maskmode=$modifier_mask_mode enc=$antimutator_encoding q=$initial_mutator_fraction"
echo " tracking     : on=$trackindividuals interval=$trackinterval startgen=$trackstartgen"
echo "==============================================================="

SECONDS=0

# 36 positional arguments, in this exact order. Adding, removing or reordering
# any of them shifts everything after it - see AssignArgumentstoVar in main.c.
"$EXE" \
    "$timeSteps" "$initialPopsize" "$mud" "$chromosomesize" "$numberofchromosomes" \
    "$bentodelratio" "$sb" "$bendist" "$typeofrun" "$slope" \
    "$seed" "$K" "$fitnesstype" "$r" "$i_init" \
    "$s" "$tskitstatus" "$modularepis" "$elementsperl" "$snapshot" \
    "$file1" "$SdtoSbratio" "$deldist" "$rawdatafilesize" "$redinmaxpopsize" \
    "$calcfixation" "$mutator_strength_factor" "$mutator_switch_rate" "$mutator_bias" \
    "$modifier_locus_fraction" "$modifier_mask_mode" "$antimutator_encoding" \
    "$initial_mutator_fraction" "$trackindividuals" "$trackinterval" "$trackstartgen"

status=$?

echo "==============================================================="
echo " finished     : $(date)"
echo " elapsed      : ${SECONDS}s"
echo " exit status  : $status"
echo "==============================================================="

exit $status
