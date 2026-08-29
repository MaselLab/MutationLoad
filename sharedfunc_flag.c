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
#include "relative_functions.h"
/* #include "absolute_functions.h"  - COMMENTED OUT (Tier 3): absolute_functions.c
 * is not in the Makefile and nothing in this file references any symbol it
 * declares. Including it only pulled in declarations for a translation unit that
 * is never linked. */
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* -------------------------------------------------------------------------
 * MAXMUTATIONSPERGAMETE  (fix for the fixed-size Sds[] buffer)
 * -------------------------------------------------------------------------
 * The deleterious effect sizes drawn for one gamete are held in a fixed-size
 * stack buffer so that no allocation happens on the hot path. The old buffer
 * was 30 entries, which was safe only because the mutation rate was constant.
 * Now that mu evolves upwards (mu = mu0 * f^n) an unlucky individual can draw
 * far more than 30 mutations, so the buffer is raised to 10000 entries
 * (10000 * sizeof(double) = 80 kB of stack, well inside the default 8 MB).
 *
 * If a draw ever exceeds this, the run ABORTS with a diagnostic rather than
 * silently truncating the mutations, because a truncated draw would bias the
 * realised mutation rate without any visible symptom. Hitting this limit means
 * mutator_strength_factor is too large for the chosen baseline rate.
 * ------------------------------------------------------------------------- */
#define MAXMUTATIONSPERGAMETE 10000

void MutateGamete(int tskitstatus, int isburninphaseover,  tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, double mutationeffectsize)
{
    /* NOTE (Tier 3): the isburninphaseover parameter is unused in this function.
     * In the original code it gated recording for absolute runs only; the relative
     * path now controls recording through the flag passed in the tskitstatus slot
     * (see item 7). The parameter is kept so this shared signature is unchanged. */
    tsk_id_t idofnewmutation;
    int mutatedsite = DetermineMutationSite(totalindividualgenomelength/2);
    if(isabsolute){
        gamete[mutatedsite] += (mutationeffectsize);
    }else{
        gamete[mutatedsite] += log(1 + mutationeffectsize);
    }
    char derivedstate[400];
    sprintf(derivedstate, "%.11f", mutationeffectsize);
    /* NOTE (item 7): callers in the relative path now pass a *recording-active*
     * flag in the tskitstatus slot rather than the raw tskitstatus, so that
     * tskitstatus == 2 ("record only after burn-in") is honoured here without
     * changing this function's contract. See PerformOneTimeStepRel(). */
    if (tskitstatus != 0){
        idofnewmutation = tsk_mutation_table_add_row(&treesequencetablecollection->mutations, wholepopulationsitesarray[mutatedsite], childnode, TSK_NULL, ((double) totaltimesteps - currenttimestep), derivedstate, 12, NULL, 0);
        check_tsk_error(idofnewmutation); 
    }
}

// PerformDeath remains mostly unchanged, just fixing array references passed as NULL in relative runs
double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, Individual *wholepopulation, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer)
{
    /* The isabsolute branch was an empty stub. main() aborts on absolute runs and
     * RunSimulationRel refuses them, so isabsolute is always false here and only
     * the relative body ever executed. Original structure preserved:
     *
     * if(isabsolute){
     *     // Absolute logic
     * }
     * else{
     */
    *psumofloads -= wholepopulation[victim].fitness;
    wholepopulation[victim].fitness = 0.0;
    /* } */
    Fen_set(wholepopulationselectiontree, maxPopSize, 0.0, victim);
    return 0.0;
}

/* NOTE (Tier 1): ismodular and elementsperlb are unused - modular epistasis is
 * not supported in this build and main() rejects it. The parameters are kept so
 * this shared signature is unchanged. */
void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gameteFitness, int *parent1gameteMutators, char *parent1gameteMask, double *parent2gameteFitness, int *parent2gameteMutators, char *parent2gameteMask, int maxPopSize, int *pPopSize, int birthplace, Individual *wholepopulation, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer, const char *globalmodifiermask, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate)
{
    int i;
    long double newwi;
    int halfgenome = totalindividualgenomelength/2;

    // Copy gametes into the individual at 'birthplace'
    for (i = 0; i < halfgenome; i++) {
        wholepopulation[birthplace].fitnessArray[i] = parent1gameteFitness[i];
        wholepopulation[birthplace].fitnessArray[halfgenome + i] = parent2gameteFitness[i];
        wholepopulation[birthplace].mutatorArray[i] = parent1gameteMutators[i];
        wholepopulation[birthplace].mutatorArray[halfgenome + i] = parent2gameteMutators[i];
    }

    /* Inherited-mask mode only: the "which blocks carry a modifier locus" mask
     * travels with the gamete, so it is copied here exactly like the state and
     * fitness arrays. In global-mask mode modifierMask is NULL everywhere and
     * the shared globalmodifiermask is consulted instead, so there is nothing
     * to copy and no per-individual memory is spent. */
    if (wholepopulation[birthplace].modifierMask != NULL && parent1gameteMask != NULL && parent2gameteMask != NULL) {
        for (i = 0; i < halfgenome; i++) {
            wholepopulation[birthplace].modifierMask[i] = parent1gameteMask[i];
            wholepopulation[birthplace].modifierMask[halfgenome + i] = parent2gameteMask[i];
        }
    }

    /* As in PerformDeath, the isabsolute branch was an empty stub and isabsolute
     * is always false in this build. Original structure preserved:
     *
     * if(isabsolute){
     *     // Absolute logic...
     * }
     * else{
     */
    // Re-calculate fitness since we just overwrote the arrays
    double currentlinkageblockssum = 0.0;
    for (i = 0; i < totalindividualgenomelength; i++) {
        currentlinkageblockssum += wholepopulation[birthplace].fitnessArray[i];
    }
    newwi = exp(currentlinkageblockssum);

    Fen_set(wholepopulationselectiontree, maxPopSize, newwi, birthplace);
    wholepopulation[birthplace].fitness = newwi;
    *psumofloads += newwi;
    /* } */
    
    // Update Cached Mutation Rate for new individual
    UpdateIndividual(&wholepopulation[birthplace], totalindividualgenomelength, globalmodifiermask, mutator_strength_factor, baseline_deleterious_rate, baseline_beneficial_rate);

    /* NOTE (item 7): as in MutateGamete, the relative path passes a
     * recording-active flag here, not the raw tskitstatus. */
    if (tskitstatus != 0){
        wholepopulationnodesarray[birthplace*2] = childnode1;
        wholepopulationnodesarray[birthplace*2 + 1] = childnode2; 
    }
}

Individual createIndividual(double *fitnessArray, int *mutatorArray, char *modifierMask, int totalindividualgenomelength){
    Individual ind;
    ind.fitnessArray = fitnessArray;
    ind.mutatorArray = mutatorArray;
    ind.modifierMask = modifierMask;   /* NULL in MODIFIERMASK_GLOBAL mode */
    ind.fitness = 0.0;
    ind.mutationRate = 0.0;
    ind.beneficialMutationRate = 0.0;
    ind.mutatorCount = 0;
    ind.modifierCount = 0;
    ind.netModifierSum = 0;
    return ind;
}

/* -------------------------------------------------------------------------
 * UpdateIndividual
 * -------------------------------------------------------------------------
 * Recomputes everything that is cached on an Individual: its fitness Wi, the
 * exponent n, the realised deleterious and beneficial mutation rates, and the
 * bookkeeping counts used for the per-generation output.
 *
 *   n  = sum of mutatorArray over the whole diploid genome
 *        (+1 per mutator allele; 0 or -1 per anti-mutator allele depending on
 *         the antimutator_encoding; 0 per non-modifier block)
 *
 *   mu_deleterious = mu_d0 * f^n        <- item 8: BOTH rates are now scaled by
 *   mu_beneficial  = mu_b0 * f^n           the same modifier equation, using the
 *                                          same f and the same n, differing only
 *                                          in the baseline rate.
 *
 * globalmodifiermask is the shared haploid-length mask in MODIFIERMASK_GLOBAL
 * mode and NULL in MODIFIERMASK_INHERITED mode (where ind->modifierMask is used).
 * It is only needed to count how many blocks carry a modifier locus at all,
 * which is reported but does not enter the rate calculation.
 * ------------------------------------------------------------------------- */
void UpdateIndividual(Individual *ind, int totalindividualgenomelength, const char *globalmodifiermask, double mutator_strength_factor, double baseline_deleterious_rate, double baseline_beneficial_rate){
    int i;
    long double currentlinkageblockssum = 0.0;
    int mutator_count = 0;
    int modifier_count = 0;
    int net_modifier_sum = 0;
    int halfgenome = totalindividualgenomelength / 2;

    /* One pass over the haploid length handles both homologs, so the global
     * mask (which is haploid-length and mirrored across homologs) can be read
     * with a single index and no modulo in the inner loop. */
    for (i = 0; i < halfgenome; i++){
        int stateA = ind->mutatorArray[i];
        int stateB = ind->mutatorArray[halfgenome + i];

        currentlinkageblockssum += ind->fitnessArray[i];
        currentlinkageblockssum += ind->fitnessArray[halfgenome + i];

        net_modifier_sum += stateA + stateB;
        if (stateA == 1) mutator_count++;
        if (stateB == 1) mutator_count++;

        if (globalmodifiermask != NULL) {
            /* MODIFIERMASK_GLOBAL: same mask on both homologs. */
            if (globalmodifiermask[i]) modifier_count += 2;
        } else if (ind->modifierMask != NULL) {
            /* MODIFIERMASK_INHERITED: each slot carries its own mask bit. */
            if (ind->modifierMask[i]) modifier_count++;
            if (ind->modifierMask[halfgenome + i]) modifier_count++;
        }
    }

    ind->fitness = exp(currentlinkageblockssum);
    ind->mutatorCount = mutator_count;
    ind->modifierCount = modifier_count;
    ind->netModifierSum = net_modifier_sum;
    // Formula: mu = mu0 * f^n
    ind->mutationRate = baseline_deleterious_rate * pow(mutator_strength_factor, (double) net_modifier_sum);
    ind->beneficialMutationRate = baseline_beneficial_rate * pow(mutator_strength_factor, (double) net_modifier_sum);
}

void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gameteFitness, int *gameteMutators, char *gameteMask, const char *globalmodifiermask, ModifierLocusIndex *gameteModifierIndex, Individual *wholepopulation, int totalindividualgenomelength)
{
    int recombinationsite, startchromosome, h, i, returnvaluefortskit;
    
    tsk_id_t parentnode1 = (tsk_id_t) 2*persontorecombine;
    tsk_id_t parentnode2 = (tsk_id_t) (2*persontorecombine + 1);
    int chromatid_len = totalindividualgenomelength / 2;

    /* The modifier index is rebuilt from scratch for every gamete. */
    if (gameteModifierIndex != NULL) {
        gameteModifierIndex->nmutatorpositions = 0;
        gameteModifierIndex->nantimutatorpositions = 0;
    }

    /* NOTE (item 7): the relative path passes a recording-active flag here. */
    if (tskitstatus != 0){
        *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
        check_tsk_error(*childnode);
    }

    for (h = 0; h < numberofchromosomes; h++) {
        startchromosome = pcg32_boundedrand(2); 
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); 
        
        if (tskitstatus != 0){
            if (startchromosome == 0){
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (double)(h*chromosomesize), (double)(h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
            }else{
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize), (h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
                returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                check_tsk_error(returnvaluefortskit);
            }
        }
        
        // Recombine Fitness AND Mutator Arrays
        if(!ismodular){
            /* The fitness effect, the modifier state and (in inherited-mask mode)
             * the modifier mask of a linkage block are all copied through the SAME
             * breakpoint, so a modifier allele can never be separated from its
             * linkage block. While we are already walking every position we also
             * build the modifier index used by the switching step - this is what
             * makes the switching step O(number of modifier loci) instead of
             * O(genome length). See ProduceMutatedGamete()/SwitchModifierLoci(). */
            for (i = 0; i < chromosomesize; i++) {
                int source_offset;
                int idx = h*chromosomesize + i;
                int ismodifier;

                if (i < recombinationsite) {
                    source_offset = (startchromosome == 0) ? 0 : chromatid_len;
                } else {
                    source_offset = (startchromosome == 0) ? chromatid_len : 0;
                }

                gameteFitness[idx]  = wholepopulation[persontorecombine].fitnessArray[source_offset + idx];
                gameteMutators[idx] = wholepopulation[persontorecombine].mutatorArray[source_offset + idx];

                if (globalmodifiermask != NULL) {
                    /* MODIFIERMASK_GLOBAL: mask is haploid-length and mirrored,
                     * so the copied slot has the same mask bit either way and
                     * nothing needs to be carried in the gamete. */
                    ismodifier = globalmodifiermask[idx];
                } else {
                    /* MODIFIERMASK_INHERITED: carry the mask bit with the block. */
                    ismodifier = wholepopulation[persontorecombine].modifierMask[source_offset + idx];
                    if (gameteMask != NULL) gameteMask[idx] = (char) ismodifier;
                }

                if (ismodifier && gameteModifierIndex != NULL) {
                    if (gameteMutators[idx] == 1) {
                        gameteModifierIndex->mutatorpositions[gameteModifierIndex->nmutatorpositions++] = idx;
                    } else {
                        gameteModifierIndex->antimutatorpositions[gameteModifierIndex->nantimutatorpositions++] = idx;
                    }
                }
            }
        } else {
             /* ---------------------------------------------------------------
              * MODULAR-EPISTASIS BRANCH - DISABLED.
              * ---------------------------------------------------------------
              * This project does not use modular epistasis (run with
              * modularepis = 0). The original block below was incomplete: it
              * only ever copied the FIRST half of each chromosome ("... similar
              * for second half" was never written) and its elementsperlb
              * indexing overran the gamete buffers. It is also not aware of the
              * modifier mask. Rather than leave code that would silently produce
              * wrong genotypes, the branch now aborts. The original lines are
              * preserved verbatim, commented out, immediately below.
              * ---------------------------------------------------------------
             for (i = 0; i < recombinationsite*elementsperlb; i++) {
                int source_offset = (startchromosome == 0) ? 0 : chromatid_len;
                int idx = h*chromosomesize*elementsperlb + i;
                gameteFitness[idx] = wholepopulation[persontorecombine].fitnessArray[source_offset + idx];
                // Assuming mutators align with elements per lb, or are just per block?
                // Based on context, mutators seem to be per linkage block.
                // If elementsperlb > 1, the mutator array indexing needs careful handling or mutators need to be per element.
                // Assuming 1-to-1 mapping for simplicity given the provided code context.
                gameteMutators[idx] = wholepopulation[persontorecombine].mutatorArray[source_offset + idx];
            }
            // ... similar for second half
              * --------------------------------------------------------------- */
            fprintf(stderr, "Error: modular epistasis (modularepis = 1) is not supported by the mutation-rate-evolution build. Run with modularepis = 0.\n");
            exit(1);
        }
    }
}

/* -------------------------------------------------------------------------
 * SwitchModifierLoci   (fix for the per-locus Bernoulli loop)
 * -------------------------------------------------------------------------
 * WHAT CHANGED AND WHY
 *
 * The old implementation looped over EVERY position of the gamete (~4600) and
 * drew one uniform random number per position to decide whether that position
 * switched state. For a 20000 x 20000 run that is roughly 3.7e12 random draws
 * spent almost entirely on positions that cannot switch at all.
 *
 * Two changes, neither of which alters the model:
 *
 *   1. Only modifier loci are considered. Non-modifier blocks can never change
 *      state, so testing them was pure waste. RecombineChromosomesIntoGamete
 *      already walks the gamete to copy it, so it builds the list of modifier
 *      positions (split by current state) for free.
 *
 *   2. Instead of one Bernoulli trial per eligible locus, the NUMBER of
 *      switches is drawn once from the exact Binomial distribution that those
 *      independent trials define, and then that many DISTINCT loci are chosen
 *      uniformly at random via a partial Fisher-Yates shuffle. This is
 *      distributionally identical to the per-locus loop but costs
 *      O(number of switches) instead of O(number of loci).
 *
 * Both directions use the rates the original code used:
 *      anti-mutator -> mutator :  mutator_switch_rate * mutator_bias
 *      mutator -> anti-mutator :  mutator_switch_rate
 *
 * Both draws are taken against the PRE-switch state, exactly as the original
 * per-locus loop did (each locus was evaluated once, against the state it had
 * on entry), so a locus can never be flipped twice in one call.
 *
 * antimutatorstate is 0 or -1 depending on the antimutator_encoding argument.
 * ------------------------------------------------------------------------- */
static void SwitchModifierLoci(int *gameteMutators, ModifierLocusIndex *gameteModifierIndex, int antimutatorstate, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma)
{
    unsigned int numberofswitches;
    unsigned int j;
    double uprate, downrate;

    if (gameteModifierIndex == NULL || mutator_switch_rate <= 0.0) return;

    uprate   = mutator_switch_rate * mutator_bias;
    downrate = mutator_switch_rate;
    if (uprate   > 1.0) uprate   = 1.0;
    if (downrate > 1.0) downrate = 1.0;

    /* anti-mutator -> mutator */
    if (gameteModifierIndex->nantimutatorpositions > 0 && uprate > 0.0) {
        numberofswitches = gsl_ran_binomial(randomnumbergeneratorforgamma, uprate, (unsigned int) gameteModifierIndex->nantimutatorpositions);
        for (j = 0; j < numberofswitches; j++) {
            /* Partial Fisher-Yates: swap a uniformly chosen not-yet-used entry
             * into slot j, guaranteeing j distinct positions after j steps. */
            int remaining = gameteModifierIndex->nantimutatorpositions - (int) j;
            int pick = (int) j + (int) pcg32_boundedrand((uint32_t) remaining);
            int chosen = gameteModifierIndex->antimutatorpositions[pick];
            gameteModifierIndex->antimutatorpositions[pick] = gameteModifierIndex->antimutatorpositions[j];
            gameteModifierIndex->antimutatorpositions[j] = chosen;
            gameteMutators[chosen] = 1;
        }
    }

    /* mutator -> anti-mutator */
    if (gameteModifierIndex->nmutatorpositions > 0 && downrate > 0.0) {
        numberofswitches = gsl_ran_binomial(randomnumbergeneratorforgamma, downrate, (unsigned int) gameteModifierIndex->nmutatorpositions);
        for (j = 0; j < numberofswitches; j++) {
            int remaining = gameteModifierIndex->nmutatorpositions - (int) j;
            int pick = (int) j + (int) pcg32_boundedrand((uint32_t) remaining);
            int chosen = gameteModifierIndex->mutatorpositions[pick];
            gameteModifierIndex->mutatorpositions[pick] = gameteModifierIndex->mutatorpositions[j];
            gameteModifierIndex->mutatorpositions[j] = chosen;
            gameteMutators[chosen] = antimutatorstate;
        }
    }
}

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double parent_specific_deleterious_rate, double parent_specific_beneficial_rate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gameteFitness, int *gameteMutators, ModifierLocusIndex *gameteModifierIndex, int antimutatorstate, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    /* See MAXMUTATIONSPERGAMETE at the top of this file: raised from 30 to 10000
     * because the mutation rate now evolves. */
    static double Sds[MAXMUTATIONSPERGAMETE];

    // 1. Fitness Mutations (Uses Parent's Mutation Rate)
    bool stayInWhileLoop = true;
    while (stayInWhileLoop) {
        stayInWhileLoop = false;
        numberofdeleteriousmutations = DetermineNumberOfMutations(parent_specific_deleterious_rate);

        /* Hard stop rather than a silent truncation - a truncated draw would
         * bias the realised mutation rate with no visible symptom. */
        if (numberofdeleteriousmutations > MAXMUTATIONSPERGAMETE) {
            fprintf(miscfilepointer, "\nFATAL: %d deleterious mutations drawn for one gamete of individual %d, which exceeds MAXMUTATIONSPERGAMETE (%d).\n", numberofdeleteriousmutations, parent, MAXMUTATIONSPERGAMETE);
            fprintf(miscfilepointer, "The parent's realised deleterious mutation rate was %g. Reduce mutator_strength_factor or the baseline mutation rate, or raise MAXMUTATIONSPERGAMETE in sharedfunc_flag.c.\n", parent_specific_deleterious_rate);
            fflush(miscfilepointer);
            fprintf(stderr, "FATAL: deleterious mutation count %d exceeds MAXMUTATIONSPERGAMETE (%d). See miscellaneous.txt.\n", numberofdeleteriousmutations, MAXMUTATIONSPERGAMETE);
            exit(1);
        }

        for (k = 0; k < numberofdeleteriousmutations; k++) {
            if (deleteriousdistribution == 0) {
                 Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646);
            } else if (deleteriousdistribution == 1) {
                 Sds[k] = gsl_ran_exponential(randomnumbergeneratorforgamma, Sd);
            } else {
                 Sds[k] = Sd;
            }
            if (!isabsolute && Sds[k] >= 1) {
                stayInWhileLoop = true;
                break;
            }
        }
    }

    for (k = 0; k < numberofdeleteriousmutations; k++) {
        double effect = isabsolute ? Sds[k] : -Sds[k];
        MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, effect);
    }
    
    /* ---------------------------------------------------------------------
     * BENEFICIAL MUTATIONS (item 8)
     * ---------------------------------------------------------------------
     * Two fixes here:
     *  (a) the number of beneficial mutations is now drawn from the PARENT'S
     *      realised beneficial rate mu_b0 * f^n, i.e. the same modifier
     *      equation as the deleterious rate with a different baseline, rather
     *      than from the unmodified population-wide beneficialmutationrate; and
     *  (b) beneficialdistribution is honoured again. The previous version had
     *      collapsed every distribution to a point effect of Sb ("Simplified
     *      for brevity"). The three branches below are restored from the
     *      pre-refactor implementation:
     *          0 -> point,       effect = Sb
     *          1 -> exponential, mean  = Sb
     *          2 -> uniform on [0, 2*Sb]
     * Effect sign convention is unchanged: negative under absolute fitness,
     * positive under relative fitness.
     * --------------------------------------------------------------------- */
    numberofbeneficialmutations = DetermineNumberOfMutations(parent_specific_beneficial_rate);

    if (beneficialdistribution == 0) {
        //point distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = Sb;
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, (isabsolute ? -generatedSb : generatedSb));
        }
    } else if (beneficialdistribution == 1) {
        //exponential distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, (isabsolute ? -generatedSb : generatedSb));
        }
    } else if (beneficialdistribution == 2) {
        //uniform distribution
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gameteFitness, (isabsolute ? -generatedSb : generatedSb));
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        fflush(miscfilepointer);
        exit(0);
    }

    /* 2. Modifier-locus switching (0/-1 <-> +1).
     *    Only modifier loci are eligible; non-modifier blocks stay at 0 forever.
     *    See the comment block above SwitchModifierLoci for what changed. */
    SwitchModifierLoci(gameteMutators, gameteModifierIndex, antimutatorstate, mutator_switch_rate, mutator_bias, randomnumbergeneratorforgamma);

    return true;
}
