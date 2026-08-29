#ifndef SHAREFUNC_FLAG_H_INCLUDED
#define SHAREFUNC_FLAG_H_INCLUDED 1

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <err.h>
#include "dependencies/pcg_basic.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* =========================================================================
 * MODIFIER-LOCUS MODEL (mutation-rate evolution)
 * =========================================================================
 * Every linkage block carries a fitness component. A USER-CHOSEN FRACTION of
 * the linkage blocks additionally carries a "modifier locus". A modifier locus
 * is in one of two states; a block without a modifier locus is a "non-modifier"
 * block and is permanently inert.
 *
 *      state  +1                 -> mutator allele
 *      state   0  (or  -1)       -> anti-mutator allele   (see ANTIMUTATOR_* below)
 *      state   0                 -> non-modifier block, NEVER changes
 *
 * The realised mutation rates of an individual are
 *
 *      mu_deleterious = mu_d0 * f^n
 *      mu_beneficial  = mu_b0 * f^n
 *
 * with f = mutator_strength_factor and n = the sum of the modifier states over
 * the WHOLE DIPLOID genome (Individual.netModifierSum). Non-modifier blocks
 * contribute 0 to n by construction, so they never affect the mutation rate.
 *
 * Because the modifier state lives in the same per-block array as the fitness
 * effect and is copied with the SAME recombination breakpoints, a modifier
 * allele can never separate from the linkage block it sits on.
 * ========================================================================= */

/* --- How the "which blocks carry a modifier locus" mask is handled -------- */
/* Selected by the modifier_mask_mode command-line argument.                  */
#define MODIFIERMASK_GLOBAL     0  /* One mask, length = haploid genome length, drawn once at   */
                                   /* startup. Position j is a modifier locus in EVERY          */
                                   /* individual and on BOTH homologs. The mask can therefore   */
                                   /* never change and is not stored per individual.            */
#define MODIFIERMASK_INHERITED  1  /* Each founder individual gets its own independent mask,    */
                                   /* mirrored across that individual's two homologs. The mask  */
                                   /* is then inherited: it recombines with the fitness and     */
                                   /* state arrays, so which blocks are modifier loci evolves.  */

/* --- How the anti-mutator state is encoded -------------------------------- */
/* Selected by the antimutator_encoding command-line argument.                */
#define ANTIMUTATOR_AS_ZERO    0   /* anti-mutator stored as  0 -> contributes  0 to n.         */
                                   /* n is then simply the COUNT of mutator alleles and the     */
                                   /* mutation rate can only ever be >= mu0 (for f > 1).        */
#define ANTIMUTATOR_AS_MINUS1  1   /* anti-mutator stored as -1 -> contributes -1 to n.         */
                                   /* n is then the NET SUM (#mutators - #anti-mutators), so an */
                                   /* anti-mutator-heavy individual gets mu < mu0.              */

/* -------------------------------------------------------------------------
 * MutatorConfig
 * -------------------------------------------------------------------------
 * Everything the user can set about the modifier-locus / mutation-rate-evolution
 * model, bundled so the already very long simulation signatures do not grow by
 * another seven arguments. Filled once in main() from the command line and then
 * passed by value down to RunSimulationRel().
 * ------------------------------------------------------------------------- */
typedef struct{
    double strengthfactor;         /* f in mu = mu0 * f^n. 1.0 = mutator alleles have no effect. */
    double switchrate;             /* Per-MODIFIER-LOCUS, per-gamete probability of switching     */
                                   /* state. 0 disables modifier evolution entirely.             */
    double bias;                   /* Multiplier applied to the anti-mutator -> mutator rate      */
                                   /* relative to the mutator -> anti-mutator rate.              */
    double locusfraction;          /* p: fraction of linkage blocks that carry a modifier locus.  */
                                   /* Exactly round(p * haploid genome length) blocks are chosen  */
                                   /* uniformly at random, without replacement.                  */
    int    maskmode;               /* MODIFIERMASK_GLOBAL or MODIFIERMASK_INHERITED.             */
    int    antimutatorencoding;    /* ANTIMUTATOR_AS_ZERO or ANTIMUTATOR_AS_MINUS1.              */
    int    antimutatorstate;       /* Derived from antimutatorencoding: the integer actually      */
                                   /* stored for an anti-mutator allele, i.e. 0 or -1.           */
    double initialmutatorfraction; /* q: fraction of each haplotype's modifier loci that start in */
                                   /* the +1 (mutator) state at generation 0.                    */
} MutatorConfig;

/* -------------------------------------------------------------------------
 * TrackingConfig
 * -------------------------------------------------------------------------
 * Controls the OPTIONAL per-individual dump (fitness and mutation rates of every
 * individual). This is off by default because it writes popsize rows every time
 * it fires; with interval and startgen you can restrict it to just the window of
 * the run you actually want to plot.
 * ------------------------------------------------------------------------- */
typedef struct{
    int enabled;   /* 0 = off, 1 = on                                                   */
    int interval;  /* Dump every this many N-timesteps (generations). Must be >= 1.      */
    int startgen;  /* First generation (1-based, as printed in the raw data file) that   */
                   /* is eligible for dumping. Use this to skip the burn-in.            */
} TrackingConfig;

/* Global Individual struct definition */
typedef struct{
    double *fitnessArray; // Still need to change the name to WiArray
    int *mutatorArray;    // Per-block modifier STATE. +1 = mutator; 0 or -1 = anti-mutator
                          // (see ANTIMUTATOR_*); 0 = non-modifier block (always inert).
    char *modifierMask;   // Per-block flag: 1 = this block carries a modifier locus, 0 = it does not.
                          // ALLOCATED ONLY IN MODIFIERMASK_INHERITED MODE. In MODIFIERMASK_GLOBAL
                          // mode this stays NULL and the single shared global mask is consulted
                          // instead, which saves popsize * genomelength bytes.
    double fitness;      // Overall fitness (Wi)
    double mutationRate; // Actual DELETERIOUS mutation rate based on mutator loci (mu_d0 * f^n)
    double beneficialMutationRate; // Actual BENEFICIAL mutation rate based on mutator loci (mu_b0 * f^n)
    int mutatorCount;    // Number of blocks currently in the +1 (mutator) state.
    int modifierCount;   // Number of blocks that carry a modifier locus at all.
    int netModifierSum;  // n, i.e. the sum of mutatorArray over the whole diploid genome.
} Individual;

/* -------------------------------------------------------------------------
 * ModifierLocusIndex
 * -------------------------------------------------------------------------
 * Scratch index over ONE gamete, built for free by RecombineChromosomesIntoGamete
 * (which already walks every position of the gamete to copy it) and consumed by
 * ProduceMutatedGamete.
 *
 * It exists purely as a performance fix: the old code tested EVERY one of the
 * ~4600 gamete positions against the switch rate on every birth. With this index
 * the switching step only ever touches the modifier loci, and draws the NUMBER of
 * switches from a binomial instead of rolling a die per locus. See the comment
 * block above SwitchModifierLoci() in sharedfunc_flag.c.
 * ------------------------------------------------------------------------- */
typedef struct{
    int *mutatorpositions;      /* gamete positions whose state is +1 (mutator)            */
    int  nmutatorpositions;
    int *antimutatorpositions;  /* gamete positions that are modifier loci but not +1      */
    int  nantimutatorpositions;
} ModifierLocusIndex;

void MutateGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, double mutationeffectsize);

double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, Individual *wholepopulation, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer);

// Updated to accept mutator arrays, the inherited modifier mask, and mutator parameters
void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gameteFitness, int *parent1gameteMutators, char *parent1gameteMask, double *parent2gameteFitness, int *parent2gameteMutators, char *parent2gameteMask, int maxPopSize, int *pPopSize, int birthplace, Individual *wholepopulation, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer, const char *globalmodifiermask, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

// Individual helper functions
Individual createIndividual(double *fitnessArray, int *mutatorArray, char *modifierMask, int totalindividualgenomelength);
void UpdateIndividual(Individual *ind, int totalindividualgenomelength, const char *globalmodifiermask, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate);

// Updated Recombine and Mutate prototypes
void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gameteFitness, int *gameteMutators, char *gameteMask, const char *globalmodifiermask, ModifierLocusIndex *gameteModifierIndex, Individual *wholepopulation, int totalindividualgenomelength);

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double parent_specific_deleterious_rate, double parent_specific_beneficial_rate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gameteFitness, int *gameteMutators, ModifierLocusIndex *gameteModifierIndex, int antimutatorstate, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer);

#endif // SHAREFUNC_FLAG_H_INCLUDED
