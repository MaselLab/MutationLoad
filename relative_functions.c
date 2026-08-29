#include <stdio.h>
#include <stdlib.h>
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
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* =========================================================================
 * DrawModifierMask
 * =========================================================================
 * Chooses which linkage blocks carry a modifier locus.
 *
 * EXACTLY round(locusfraction * haploidgenomelength) blocks are chosen,
 * uniformly at random and without replacement, via a partial Fisher-Yates
 * shuffle. Using an exact count (rather than an independent coin flip per
 * block) means the NUMBER of modifier loci is fixed by the parameter, so runs
 * that differ only in random seed are directly comparable.
 *
 * The mask is haploid-length. How it is used depends on mutatorconfig.maskmode:
 *   MODIFIERMASK_GLOBAL    - one mask for the whole run, mirrored across both
 *                            homologs of every individual.
 *   MODIFIERMASK_INHERITED - this is called once per founder individual, and
 *                            the result is mirrored across that individual's
 *                            two homologs and thereafter inherited.
 * ========================================================================= */
void DrawModifierMask(char *mask, int haploidgenomelength, double locusfraction)
{
    int i, numbertochoose, pick, temp;
    int *positions;

    for (i = 0; i < haploidgenomelength; i++) mask[i] = 0;

    if (locusfraction <= 0.0) return;

    numbertochoose = (int) floor(locusfraction * (double) haploidgenomelength + 0.5);
    if (numbertochoose < 0) numbertochoose = 0;
    if (numbertochoose > haploidgenomelength) numbertochoose = haploidgenomelength;
    if (numbertochoose == 0) return;

    positions = malloc(sizeof(int) * haploidgenomelength);
    for (i = 0; i < haploidgenomelength; i++) positions[i] = i;

    for (i = 0; i < numbertochoose; i++) {
        pick = i + (int) pcg32_boundedrand((uint32_t)(haploidgenomelength - i));
        temp = positions[pick];
        positions[pick] = positions[i];
        positions[i] = temp;
        mask[positions[i]] = 1;
    }
    free(positions);
}

/* =========================================================================
 * SeedInitialMutatorStates
 * =========================================================================
 * Fills one HAPLOID state array from a haploid modifier mask.
 *
 *   - a block with no modifier locus is set to 0 and can never change;
 *   - exactly round(initialmutatorfraction * m) of the m modifier loci are set
 *     to +1 (mutator), chosen uniformly at random without replacement;
 *   - every other modifier locus is set to antimutatorstate (0 or -1, per the
 *     antimutator_encoding argument).
 * ========================================================================= */
void SeedInitialMutatorStates(int *stateshaploid, const char *mask, int haploidgenomelength, double initialmutatorfraction, int antimutatorstate)
{
    int i, nmodifier = 0, numbertochoose, pick, temp;
    int *modifierpositions;

    for (i = 0; i < haploidgenomelength; i++) {
        if (mask[i]) {
            stateshaploid[i] = antimutatorstate;
            nmodifier++;
        } else {
            stateshaploid[i] = 0;   /* non-modifier block: permanently inert */
        }
    }

    if (nmodifier == 0 || initialmutatorfraction <= 0.0) return;

    numbertochoose = (int) floor(initialmutatorfraction * (double) nmodifier + 0.5);
    if (numbertochoose < 0) numbertochoose = 0;
    if (numbertochoose > nmodifier) numbertochoose = nmodifier;
    if (numbertochoose == 0) return;

    modifierpositions = malloc(sizeof(int) * nmodifier);
    nmodifier = 0;
    for (i = 0; i < haploidgenomelength; i++) {
        if (mask[i]) modifierpositions[nmodifier++] = i;
    }

    for (i = 0; i < numbertochoose; i++) {
        pick = i + (int) pcg32_boundedrand((uint32_t)(nmodifier - i));
        temp = modifierpositions[pick];
        modifierpositions[pick] = modifierpositions[i];
        modifierpositions[i] = temp;
        stateshaploid[modifierpositions[i]] = 1;
    }
    free(modifierpositions);
}

/* =========================================================================
 * SeedTreeSequenceTables   (item 7)
 * =========================================================================
 * Puts a freshly initialised (or freshly cleared) table collection into the
 * state the simulation expects before it starts recording: sequence length set,
 * 2 x popsize sample nodes stamped at nodetime, and one site per haploid
 * linkage block with ancestral state 0.
 *
 * Factored out of InitializePopulationRel because tskitstatus == 2 needs to do
 * exactly this a second time, at the moment the burn-in ends.
 * ========================================================================= */
void SeedTreeSequenceTables(tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int haploidgenomelength, double nodetime)
{
    int i;
    treesequencetablecollection->sequence_length = (double) haploidgenomelength;
    for (i = 0; i < (2 * popsize); i++) {
        wholepopulationnodesarray[i] = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, nodetime, TSK_NULL, TSK_NULL, NULL, 0);
        check_tsk_error(wholepopulationnodesarray[i]);
    }
    for (i = 0; i < haploidgenomelength; i++) {
        wholepopulationsitesarray[i] = tsk_site_table_add_row(&treesequencetablecollection->sites, i, "0.00000000", 10, NULL, 0);
        check_tsk_error(wholepopulationsitesarray[i]);
    }
}

/* =========================================================================
 * WritePopulationModifierSummary   (item 5)
 * =========================================================================
 * Appends the mutation-rate-evolution columns to one row of the raw data file.
 * The caller has already written the leading columns and writes the newline.
 *
 * Columns appended, in order:
 *   Mean.deleterious.mutation.rate   mean over individuals of mu_d0 * f^n
 *   Mean.beneficial.mutation.rate    mean over individuals of mu_b0 * f^n
 *   Mean.net.modifier.sum            mean over individuals of n
 *   Mean.mutator.freq.perindividual  mean over individuals of (that individual's
 *                                    mutator alleles / its modifier slots)
 *   Var.mutator.freq.acrossindividuals   variance of the quantity above
 *   Mean.mutator.freq.perlocus       mean over LOCI of (mutator alleles at that
 *                                    locus / modifier slots at that locus,
 *                                    pooled over both homologs and all individuals)
 *   Var.mutator.freq.acrossloci      variance of the quantity above
 *
 * The two "per individual" columns and the two "per locus" columns answer
 * different questions: the first pair measures how much individuals differ from
 * each other in mutator load; the second is the classic site-frequency view and
 * measures how much the individual modifier loci differ from each other.
 *
 * locusmutatorcounts and locusmodifiercounts are caller-owned scratch buffers of
 * haploid genome length, passed in so they are not reallocated every generation.
 * The per-locus pass is O(popsize x haploid genome length); that is negligible
 * next to the O(popsize x genome length) recombination work already done every
 * single timestep, of which there are popsize per generation.
 * ========================================================================= */
void WritePopulationModifierSummary(FILE *rawdatafilepointer, Individual *wholepopulation, int popsize, int totalindividualgenomelength, const char *globalmodifiermask, int *locusmutatorcounts, int *locusmodifiercounts)
{
    int i, j;
    int halfgenome = totalindividualgenomelength / 2;
    double sumdelrate = 0.0, sumbenrate = 0.0, sumnet = 0.0;
    double sumfreq = 0.0, sumfreqsquared = 0.0;
    double meandelrate, meanbenrate, meannet, meanfreqind, varfreqind;
    double meanfreqlocus = 0.0, varfreqlocus = 0.0;
    int nlociwithmodifiers = 0;
    double locussum = 0.0, locussumsquared = 0.0;

    /* ---- per-individual pass ---- */
    for (i = 0; i < popsize; i++) {
        double freq;
        sumdelrate += wholepopulation[i].mutationRate;
        sumbenrate += wholepopulation[i].beneficialMutationRate;
        sumnet     += (double) wholepopulation[i].netModifierSum;
        freq = (wholepopulation[i].modifierCount > 0)
             ? ((double) wholepopulation[i].mutatorCount / (double) wholepopulation[i].modifierCount)
             : 0.0;
        sumfreq        += freq;
        sumfreqsquared += freq * freq;
    }
    meandelrate = sumdelrate / (double) popsize;
    meanbenrate = sumbenrate / (double) popsize;
    meannet     = sumnet     / (double) popsize;
    meanfreqind = sumfreq    / (double) popsize;
    varfreqind  = (sumfreqsquared / (double) popsize) - (meanfreqind * meanfreqind);
    if (varfreqind < 0.0) varfreqind = 0.0;   /* guard against round-off */

    /* ---- per-locus pass ---- */
    for (j = 0; j < halfgenome; j++) {
        locusmutatorcounts[j] = 0;
        locusmodifiercounts[j] = 0;
    }
    for (i = 0; i < popsize; i++) {
        for (j = 0; j < halfgenome; j++) {
            int ismodifierA, ismodifierB;
            if (globalmodifiermask != NULL) {
                ismodifierA = globalmodifiermask[j];
                ismodifierB = globalmodifiermask[j];
            } else {
                ismodifierA = wholepopulation[i].modifierMask[j];
                ismodifierB = wholepopulation[i].modifierMask[halfgenome + j];
            }
            if (ismodifierA) {
                locusmodifiercounts[j]++;
                if (wholepopulation[i].mutatorArray[j] == 1) locusmutatorcounts[j]++;
            }
            if (ismodifierB) {
                locusmodifiercounts[j]++;
                if (wholepopulation[i].mutatorArray[halfgenome + j] == 1) locusmutatorcounts[j]++;
            }
        }
    }
    for (j = 0; j < halfgenome; j++) {
        if (locusmodifiercounts[j] > 0) {
            double f = (double) locusmutatorcounts[j] / (double) locusmodifiercounts[j];
            locussum += f;
            locussumsquared += f * f;
            nlociwithmodifiers++;
        }
    }
    if (nlociwithmodifiers > 0) {
        meanfreqlocus = locussum / (double) nlociwithmodifiers;
        varfreqlocus  = (locussumsquared / (double) nlociwithmodifiers) - (meanfreqlocus * meanfreqlocus);
        if (varfreqlocus < 0.0) varfreqlocus = 0.0;
    }

    fprintf(rawdatafilepointer, ",%.12g,%.12g,%.6f,%.10f,%.10f,%.10f,%.10f",
            meandelrate, meanbenrate, meannet, meanfreqind, varfreqind, meanfreqlocus, varfreqlocus);
}

/* =========================================================================
 * WriteIndividualSnapshot   (item 5, optional detailed tracking)
 * =========================================================================
 * Dumps one row per individual. Controlled entirely by TrackingConfig, which is
 * off by default because at popsize = 20000 each firing writes 20000 rows.
 * ========================================================================= */
void WriteIndividualSnapshot(FILE *individualfilepointer, Individual *wholepopulation, int popsize, int generation)
{
    int k;
    for (k = 0; k < popsize; k++) {
        fprintf(individualfilepointer, "%d,%d,%.12g,%.12g,%.12g,%.12g,%d,%d,%d\n",
                generation,
                k + 1,
                wholepopulation[k].fitness,
                log(wholepopulation[k].fitness),
                wholepopulation[k].mutationRate,
                wholepopulation[k].beneficialMutationRate,
                wholepopulation[k].mutatorCount,
                wholepopulation[k].modifierCount,
                wholepopulation[k].netModifierSum);
    }
}

double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig)
{
    if(isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationRel within an absolute fitness program. \n");
        exit(0);
    }

    /* PerformOneTimeStepRel draws two DISTINCT parents, which is impossible with
     * a single individual. Refuse rather than spin forever in that loop. */
    if (popsize < 2) {
        fprintf(miscfilepointer, "\nError: popsize is %d. RunSimulationRel needs at least 2 individuals because a birth requires two distinct parents.\n", popsize);
        fflush(miscfilepointer);
        fprintf(stderr, "Error: popsize must be at least 2.\n");
        exit(1);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    FILE *nodefilepointer;
    FILE *edgefilepointer;
    FILE *sitefilepointer;
    FILE *mutationfilepointer;
    FILE *individualfilepointer = NULL;   /* optional per-individual tracking (item 5) */
    
    int i, j, k;
    
    /* The mutator/modifier parameters are now part of the raw data file name.
     * Without them, two runs that differ only in their mutator settings landed in
     * the same directory under the same file name and silently overwrote each
     * other, because MakeDirectoryName() does not know about them either. */
    char * rawdatafilename = (char *) malloc(400);
    strcpy(rawdatafilename, "rawdatafor");
    strcat(rawdatafilename, "Nxtimesteps"); strcat(rawdatafilename, Nxtimestepsname);
    strcat(rawdatafilename, "popsize"); strcat(rawdatafilename, popsizename);
    strcat(rawdatafilename, "mutrate"); strcat(rawdatafilename, delmutratename);
    strcat(rawdatafilename, "chromsize"); strcat(rawdatafilename, chromsizename);
    strcat(rawdatafilename, "chromnum"); strcat(rawdatafilename, chromnumname);
    strcat(rawdatafilename, "benmutrate"); strcat(rawdatafilename, mubname);
    strcat(rawdatafilename, "Sb"); strcat(rawdatafilename, Sbname);
    strcat(rawdatafilename, "msf"); strcat(rawdatafilename, mutator_strength_factorname);
    strcat(rawdatafilename, "msr"); strcat(rawdatafilename, mutator_switch_ratename);
    strcat(rawdatafilename, "mb"); strcat(rawdatafilename, mutator_biasname);
    /* modifier-locus settings, so runs differing only in those do not collide */
    {
        char modifiertag[120];
        snprintf(modifiertag, sizeof(modifiertag), "mlf%g_imf%g_mm%d_ae%d",
                 mutatorconfig.locusfraction, mutatorconfig.initialmutatorfraction,
                 mutatorconfig.maskmode, mutatorconfig.antimutatorencoding);
        strcat(rawdatafilename, modifiertag);
    }
    strcat(rawdatafilename, ".txt");

    rawdatafilepointer = fopen(rawdatafilename, "w");
    /* Header extended with the mutation-rate-evolution columns (item 5). */
    fprintf(rawdatafilepointer, "Nxtimesteps,Sum.of.wis,Variance.in.log.fitness,FractionSelectiveDeaths,FractionSelectiveDeaths_exponantiated,Mean.deleterious.mutation.rate,Mean.beneficial.mutation.rate,Mean.net.modifier.sum,Mean.mutator.freq.perindividual,Var.mutator.freq.acrossindividuals,Mean.mutator.freq.perlocus,Var.mutator.freq.acrossloci\n");
    
    char * summarydatafilename = (char *) malloc(200);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb"); strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub"); strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, "msf"); strcat(summarydatafilename, mutator_strength_factorname);
    strcat(summarydatafilename, "msr"); strcat(summarydatafilename, mutator_switch_ratename);
    strcat(summarydatafilename, "mb"); strcat(summarydatafilename, mutator_biasname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w");
    
    nodefilepointer = fopen("nodetable.txt", "w");
    edgefilepointer = fopen("edgetable.txt", "w");
    sitefilepointer = fopen("sitetable.txt", "w");
    mutationfilepointer = fopen("mutationtable.txt", "w");
    
    int totaltimesteps = Nxtimesteps * popsize;
    double currenttimestep = 0.0;
    
    Individual *wholepopulation;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = popsize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    int haploidgenomelength = totalindividualgenomelength / 2;
    
    wholepopulation = malloc(sizeof(Individual) * popsize);
    // Note: Initialization of internal arrays happens in InitializePopulationRel

    /* The shared modifier mask. In MODIFIERMASK_GLOBAL mode this is THE mask and
     * is passed down everywhere. In MODIFIERMASK_INHERITED mode every individual
     * carries its own copy instead, and every consumer is handed NULL so that it
     * knows to look at Individual.modifierMask. */
    char *globalmodifiermask = malloc(sizeof(char) * haploidgenomelength);
    const char *effectiveglobalmask = (mutatorconfig.maskmode == MODIFIERMASK_GLOBAL) ? globalmodifiermask : NULL;
    
    long double sumofwis;
    long double *psumofwis = &sumofwis;
    long double *wholepopulationwistree;
    wholepopulationwistree = malloc(sizeof(long double) * popsize);
    
    tsk_table_collection_t treesequencetablecollection;
    tsk_table_collection_t * tablepointer = &treesequencetablecollection;
    int returnvaluefortskfunctions = tsk_table_collection_init(&treesequencetablecollection, 0);
    check_tsk_error(returnvaluefortskfunctions);
    
    tsk_id_t *wholepopulationnodesarray;
    wholepopulationnodesarray = malloc(sizeof(tsk_id_t) * 2 * popsize);
    tsk_id_t wholepopulationsitesarray[totalindividualgenomelength / 2];

    /* ---------------------------------------------------------------------
     * sortedwisarray - COMMENTED OUT (item 10: allocated and freed but never
     * read or written anywhere in the program).
     * ---------------------------------------------------------------------
    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * popsize);
     * ------------------------------------------------------------------- */

    /* Scratch buffers for the per-locus mutator-frequency statistics, allocated
     * once instead of every generation. */
    int *locusmutatorcounts  = malloc(sizeof(int) * haploidgenomelength);
    int *locusmodifiercounts = malloc(sizeof(int) * haploidgenomelength);

    InitializePopulationRel(tskitstatus, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, wholepopulationwistree, wholepopulation, popsize, totalpopulationgenomelength, totaltimesteps, psumofwis, globalmodifiermask, mutatorconfig, miscfilepointer);
    
    // Set initial mutation rates for population based on their starting modifier states
    for(k = 0; k < popsize; k++) {
        UpdateIndividual(&wholepopulation[k], totalindividualgenomelength, effectiveglobalmask, mutatorconfig.strengthfactor, deleteriousmutationrate, beneficialmutationrate);
    }

    fprintf(miscfilepointer, "Modifier-locus setup: maskmode=%d (0=global,1=inherited), locusfraction=%g, modifier loci per haplotype=%d, antimutator state=%d, initial mutator fraction=%g, f=%g, switch rate=%g, bias=%g\n",
            mutatorconfig.maskmode, mutatorconfig.locusfraction, wholepopulation[0].modifierCount / 2,
            mutatorconfig.antimutatorstate, mutatorconfig.initialmutatorfraction,
            mutatorconfig.strengthfactor, mutatorconfig.switchrate, mutatorconfig.bias);
    fflush(miscfilepointer);

    /* Optional per-individual tracking file (item 5). */
    if (trackingconfig.enabled) {
        char * individualfilename = (char *) malloc(400);
        strcpy(individualfilename, "individualtrackingfor");
        strcat(individualfilename, "Sb"); strcat(individualfilename, Sbname);
        strcat(individualfilename, "mub"); strcat(individualfilename, mubname);
        strcat(individualfilename, "msf"); strcat(individualfilename, mutator_strength_factorname);
        strcat(individualfilename, "msr"); strcat(individualfilename, mutator_switch_ratename);
        strcat(individualfilename, "mb"); strcat(individualfilename, mutator_biasname);
        strcat(individualfilename, ".txt");
        individualfilepointer = fopen(individualfilename, "w");
        fprintf(individualfilepointer, "Generation,Individual,Wi,LogWi,DeleteriousMutationRate,BeneficialMutationRate,MutatorAlleleCount,ModifierLocusCount,NetModifierSum\n");
        free(individualfilename);
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * Nxtimesteps);
    
    // Parent gamete arrays - allocated once here
    double parent1gameteFitness[numberofchromosomes*chromosomesize], parent2gameteFitness[numberofchromosomes*chromosomesize];
    int parent1gameteMutators[numberofchromosomes*chromosomesize], parent2gameteMutators[numberofchromosomes*chromosomesize];

    /* Gamete-level modifier-mask buffers. Only needed in inherited-mask mode;
     * NULL otherwise so nothing is copied and nothing is spent. */
    char *parent1gameteMask = NULL, *parent2gameteMask = NULL;
    if (mutatorconfig.maskmode == MODIFIERMASK_INHERITED) {
        parent1gameteMask = malloc(sizeof(char) * haploidgenomelength);
        parent2gameteMask = malloc(sizeof(char) * haploidgenomelength);
    }

    /* Modifier-locus indices, rebuilt for each gamete during recombination and
     * consumed by the switching step. See SwitchModifierLoci in sharedfunc_flag.c. */
    ModifierLocusIndex parent1modifierindex, parent2modifierindex;
    parent1modifierindex.mutatorpositions     = malloc(sizeof(int) * haploidgenomelength);
    parent1modifierindex.antimutatorpositions = malloc(sizeof(int) * haploidgenomelength);
    parent1modifierindex.nmutatorpositions = 0;
    parent1modifierindex.nantimutatorpositions = 0;
    parent2modifierindex.mutatorpositions     = malloc(sizeof(int) * haploidgenomelength);
    parent2modifierindex.antimutatorpositions = malloc(sizeof(int) * haploidgenomelength);
    parent2modifierindex.nmutatorpositions = 0;
    parent2modifierindex.nantimutatorpositions = 0;
    
    size_t step = 1;
    double *last200Ntimestepsvariance;
    double *literallyjustlast200Ntimesteps;
    literallyjustlast200Ntimesteps = malloc(sizeof(double) * 200);
    last200Ntimestepsvariance = malloc(sizeof(double) * 200);
    for (k = 0; k < 200; k++) {
        literallyjustlast200Ntimesteps[k] = 0.0;
        last200Ntimestepsvariance[k] = 0.0;
    }
    double slopeofvariance;
    int isburninphaseover = 0;
    int didpopulationcrash = 0;
    int endofburninphase;
    int endofdelay = Nxtimesteps-1;
    int endofsimulation = Nxtimesteps-1;
    int Nxtimestepsafterburnin = 0;
    double arbitrarynumber;
    arbitrarynumber = (-1 * 0.007 / popsize);
    double slopeoflogfitness;    
    double varianceinlogfitness;   
    long double fitnessfittest;
    long double FractionSelectiveDeaths;
    long double FractionSelectiveDeaths_exponantiatebirthrates;

    /* ---------------------------------------------------------------------
     * TREE-SEQUENCE RECORDING STATE (item 7)
     * ---------------------------------------------------------------------
     * tskitstatus semantics, now honoured in the relative path:
     *      0 - never record
     *      1 - record from generation 0 (tables already seeded in
     *          InitializePopulationRel)
     *      2 - record ONLY after the burn-in phase has been called as over
     *
     * Previously the relative path tested only "tskitstatus != 0" / "> 0", so
     * mode 2 behaved exactly like mode 1 and recorded from the very first
     * timestep - which is also what the default bash script asks for.
     *
     * istskitrecording is the single flag that everything downstream keys off.
     * It is passed to the shared recording functions IN PLACE OF tskitstatus
     * (they only ever test it against 0), which keeps their contract unchanged
     * and leaves the absolute-fitness path untouched.
     * ------------------------------------------------------------------- */
    int istskitrecording = (tskitstatus == 1) ? 1 : 0;
    
    for (i = 0; i < Nxtimesteps; i++) {
        for (j = 0; j < popsize; j++) {
            currenttimestep += 1.0;            
            PerformOneTimeStepRel(istskitrecording, isabsolute, isburninphaseover, ismodular, elementsperlb, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, totaltimesteps, currenttimestep, wholepopulationwistree, wholepopulation, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gameteFitness, parent1gameteMutators, parent1gameteMask, &parent1modifierindex, parent2gameteFitness, parent2gameteMutators, parent2gameteMask, &parent2modifierindex, effectiveglobalmask, randomnumbergeneratorforgamma, miscfilepointer, mutatorconfig);  
        }
        
        varianceinlogfitness = CalculateVarianceInLogFitness(popsize, wholepopulation, *psumofwis);
        fitnessfittest = FindFittestWi(wholepopulation, popsize);
        FractionSelectiveDeaths = (fitnessfittest-(sumofwis/popsize))/fitnessfittest;
        FractionSelectiveDeaths_exponantiatebirthrates = (exp(fitnessfittest)-exp((sumofwis/popsize)))/exp(fitnessfittest);
        
        fprintf(rawdatafilepointer, "%d,%Lf,%.18f,%Lf,%Lf", i+1, *psumofwis, varianceinlogfitness, FractionSelectiveDeaths, FractionSelectiveDeaths_exponantiatebirthrates);
        WritePopulationModifierSummary(rawdatafilepointer, wholepopulation, popsize, totalindividualgenomelength, effectiveglobalmask, locusmutatorcounts, locusmodifiercounts);
        fprintf(rawdatafilepointer, "\n");
        fflush(rawdatafilepointer);

        /* Optional detailed per-individual dump (item 5). Fires only when
         * enabled, only from trackingconfig.startgen onwards, and then only
         * every trackingconfig.interval generations. */
        if (trackingconfig.enabled && individualfilepointer != NULL) {
            int generation = i + 1;
            if (generation >= trackingconfig.startgen &&
                ((generation - trackingconfig.startgen) % trackingconfig.interval) == 0) {
                WriteIndividualSnapshot(individualfilepointer, wholepopulation, popsize, generation);
                fflush(individualfilepointer);
            }
        }

        if (istskitrecording){
            if (i % 10 == 0) {
                returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                check_tsk_error(returnvaluefortskfunctions);
                returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*popsize), 0, NULL);
                check_tsk_error(returnvaluefortskfunctions);
                for (k = 0; k < (2*popsize); k++) {
                    wholepopulationnodesarray[k] = k;
                }   
            }
        }

        double c0, cov00, cov01, cov11, sumsq;
        if (isburninphaseover == 0) {
            UpdateLast200NTimeSteps(last200Ntimestepsvariance, varianceinlogfitness);
            UpdateLast200NTimeSteps(literallyjustlast200Ntimesteps, i+1);
            if (i > 199) {           
                slopeofvariance = 0.0;
                gsl_fit_linear(literallyjustlast200Ntimesteps, step, last200Ntimestepsvariance, step, 200, &c0, &slopeofvariance, &cov00, &cov01, &cov11, &sumsq);
                if (slopeofvariance < arbitrarynumber) {
                    endofburninphase = i;
                    endofdelay = endofburninphase + 500;
                    isburninphaseover = 1;
                    fprintf(miscfilepointer, "Burn-in phase called as ending in generation %d\n", i+1);
                    fprintf(summarydatafilepointer, "Burn-in phase called as ending in generation %d\n", i+1);

                    /* ---------------------------------------------------------
                     * item 7: tskitstatus == 2 starts recording HERE.
                     * ---------------------------------------------------------
                     * The tables are thrown away and re-seeded rather than
                     * merely un-gated. If they were only un-gated, the founder
                     * nodes added at generation 0 would still be sitting in the
                     * table stamped with generation-0 times, and every edge
                     * recorded from now on would hang off them, giving branch
                     * lengths that span the entire burn-in the simulation never
                     * actually recorded. Re-seeding makes the post-burn-in
                     * population the root of the recorded tree sequence, with
                     * node times stamped at the current timestep.
                     * ------------------------------------------------------- */
                    if (tskitstatus == 2) {
                        tsk_table_collection_free(&treesequencetablecollection);
                        returnvaluefortskfunctions = tsk_table_collection_init(&treesequencetablecollection, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        SeedTreeSequenceTables(&treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, haploidgenomelength, ((double) totaltimesteps - currenttimestep));
                        istskitrecording = 1;
                        fprintf(miscfilepointer, "tskitstatus == 2: tree sequence tables seeded at generation %d; recording starts now.\n", i+1);
                        fflush(miscfilepointer);
                    }
                }
            }
        }        
        
        if (typeofrun == 1) {
            if (i == 1999) {
                if (INDIVIDUALWIDATA == 1) {
                    fprintf(summarydatafilepointer, "Individual, Wi\n");
                    for (k = 0; k < popsize; k++) {
                        fprintf(summarydatafilepointer, "%d,%Lf\n", k+1, wholepopulation[k].fitness);
                    }
                }
            }
        }
        
        if (i > endofdelay) {
            logaveragefitnesseachNtimesteps[Nxtimestepsafterburnin] = log((double) *psumofwis / (double) popsize);
            Nxtimestepsafterburnin += 1;
        }
        
        long double currentfittestindividualswi = FindFittestWi(wholepopulation, popsize);
        if (currentfittestindividualswi < pow(10.0, -10.0)) {
            endofsimulation = i;
            i = Nxtimesteps;
            didpopulationcrash = 1;
        }
    }
    
    if(istskitrecording){
        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
        check_tsk_error(returnvaluefortskfunctions);
        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*popsize), 0, NULL);
        check_tsk_error(returnvaluefortskfunctions);
        for (k = 0; k < (2*popsize); k++) {
            wholepopulationnodesarray[k] = k;
        }
        fprintf(nodefilepointer, "is_sample time\n");
        for (k = 0; k < tablepointer->nodes.num_rows; k++) {
            if (k < (2*popsize)) fprintf(nodefilepointer, "1 %f\n", tablepointer->nodes.time[k]);
            else fprintf(nodefilepointer, "0 %f\n", tablepointer->nodes.time[k]);
        }
        fprintf(edgefilepointer, "left right parent child\n");
        for (k = 0; k < tablepointer->edges.num_rows; k++) {
            fprintf(edgefilepointer, "%f %f %d %d\n", tablepointer->edges.left[k], tablepointer->edges.right[k], tablepointer->edges.parent[k], tablepointer->edges.child[k]);
        }
        fprintf(sitefilepointer, "position ancestral_state\n");
        for (k = 0; k < tablepointer->sites.num_rows; k++) {
            fprintf(sitefilepointer, "%f 0.0\n", tablepointer->sites.position[k]);
        }
        fprintf(mutationfilepointer, "site node derived_state\n");
        for (k = 0; k < tablepointer->mutations.num_rows; k++) {
            fprintf(mutationfilepointer, "%d %d %.12s\n", tablepointer->mutations.site[k], tablepointer->mutations.node[k], (tablepointer->mutations.derived_state + k*12));
        }
    }

    if (didpopulationcrash == 0) endofsimulation = i;

    if (isburninphaseover == 1) {
        slopeoflogfitness = CalculateSlopeOfLogFitness(endofsimulation, endofdelay, logaveragefitnesseachNtimesteps);
        fprintf(summarydatafilepointer, "Slope of log(fitness) after the burn-in phase: %f\n", slopeoflogfitness);
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        fclose(nodefilepointer); fclose(edgefilepointer); fclose(sitefilepointer); fclose(mutationfilepointer);
        if (individualfilepointer != NULL) fclose(individualfilepointer);
        free(rawdatafilename); free(summarydatafilename); free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps); free(last200Ntimestepsvariance);
        
        for(k=0; k<popsize; k++){
            free(wholepopulation[k].fitnessArray);
            free(wholepopulation[k].mutatorArray);
            if (wholepopulation[k].modifierMask != NULL) free(wholepopulation[k].modifierMask);
        }
        free(wholepopulation);
        free(wholepopulationwistree);
        free(wholepopulationnodesarray);
        /* free(sortedwisarray);  - see the commented-out allocation above (item 10) */
        free(globalmodifiermask);
        free(locusmutatorcounts); free(locusmodifiercounts);
        if (parent1gameteMask != NULL) free(parent1gameteMask);
        if (parent2gameteMask != NULL) free(parent2gameteMask);
        free(parent1modifierindex.mutatorpositions); free(parent1modifierindex.antimutatorpositions);
        free(parent2modifierindex.mutatorpositions); free(parent2modifierindex.antimutatorpositions);
        tsk_table_collection_free(&treesequencetablecollection);
        return slopeoflogfitness;
    }

    if (isburninphaseover == 0) {
        fprintf(summarydatafilepointer, "End of burn-in phase not reached.");
        fclose(rawdatafilepointer); fclose(summarydatafilepointer);
        fclose(nodefilepointer); fclose(edgefilepointer); fclose(sitefilepointer); fclose(mutationfilepointer);
        if (individualfilepointer != NULL) fclose(individualfilepointer);
        free(rawdatafilename); free(summarydatafilename); free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps); free(last200Ntimestepsvariance);
        for(k=0; k<popsize; k++){
            free(wholepopulation[k].fitnessArray);
            free(wholepopulation[k].mutatorArray);
            if (wholepopulation[k].modifierMask != NULL) free(wholepopulation[k].modifierMask);
        }
        free(wholepopulation);
        free(wholepopulationwistree);
        free(wholepopulationnodesarray);
        /* free(sortedwisarray);  - see the commented-out allocation above (item 10) */
        free(globalmodifiermask);
        free(locusmutatorcounts); free(locusmodifiercounts);
        if (parent1gameteMask != NULL) free(parent1gameteMask);
        if (parent2gameteMask != NULL) free(parent2gameteMask);
        free(parent1modifierindex.mutatorpositions); free(parent1modifierindex.antimutatorpositions);
        free(parent2modifierindex.mutatorpositions); free(parent2modifierindex.antimutatorpositions);
        tsk_table_collection_free(&treesequencetablecollection);
        return -1.0;
    }

    return -1.0;   /* unreachable; silences -Wreturn-type */
}

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, Individual *wholepopulation, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gameteFitness, int *parent1gameteMutators, char *parent1gameteMask, ModifierLocusIndex *parent1modifierindex, double *parent2gameteFitness, int *parent2gameteMutators, char *parent2gameteMask, ModifierLocusIndex *parent2modifierindex, const char *globalmodifiermask, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, MutatorConfig mutatorconfig)
{
    /* NOTE (item 7): RunSimulationRel passes its istskitrecording flag in the
     * tskitstatus slot, so everything below records only when recording is
     * actually active. */
    int currentparent1, currentparent2, currentvictim;
    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    /* Guarded by the popsize >= 2 check in RunSimulationRel: with a single
     * individual there is no second distinct parent and this would spin forever. */
    while (currentparent1 == currentparent2) {
        if (popsize < 2) break;
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    }
    
    tsk_id_t childnode1, childnode2;
   
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, chromosomesize, numberofchromosomes, parent1gameteFitness, parent1gameteMutators, parent1gameteMask, globalmodifiermask, parent1modifierindex, wholepopulation, totalindividualgenomelength);
    
    /* Both rates come from the PARENT, because the mutations in this gamete
     * arise in the parent's germ line. Item 8: the beneficial rate is now scaled
     * by the same mu = mu0 * f^n equation as the deleterious rate. */
    double p1_delrate = wholepopulation[currentparent1].mutationRate;
    double p1_benrate = wholepopulation[currentparent1].beneficialMutationRate;
    
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, isabsolute, totalindividualgenomelength, p1_delrate, p1_benrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gameteFitness, parent1gameteMutators, parent1modifierindex, mutatorconfig.antimutatorstate, mutatorconfig.switchrate, mutatorconfig.bias, randomnumbergeneratorforgamma, miscfilepointer);
        
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, chromosomesize, numberofchromosomes, parent2gameteFitness, parent2gameteMutators, parent2gameteMask, globalmodifiermask, parent2modifierindex, wholepopulation, totalindividualgenomelength);
    
    double p2_delrate = wholepopulation[currentparent2].mutationRate;
    double p2_benrate = wholepopulation[currentparent2].beneficialMutationRate;
    
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, isabsolute, totalindividualgenomelength, p2_delrate, p2_benrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent2gameteFitness, parent2gameteMutators, parent2modifierindex, mutatorconfig.antimutatorstate, mutatorconfig.switchrate, mutatorconfig.bias, randomnumbergeneratorforgamma, miscfilepointer);
               
    /* TODO (bug 3, not yet fixed - awaiting sign-off): pPopSize is uninitialised
     * here and is passed straight into PerformDeath and PerformBirth. It is never
     * dereferenced on the relative-fitness code paths, so this is currently
     * harmless, but it is undefined behaviour and the compiler warns about it. */
    int *pPopSize; 
    
    PerformDeath(isabsolute, tskitstatus, isburninphaseover, popsize, pPopSize, currentvictim, deleteriousdistribution, wholepopulationwistree, wholepopulation, NULL, NULL, NULL, psumofwis, NULL, NULL, 0, 0, 0, 0, NULL, NULL, wholepopulationnodesarray, miscfilepointer);
    
    PerformBirth(tskitstatus, isburninphaseover, ismodular, elementsperlb, treesequencetablecollection, wholepopulationnodesarray, childnode1, childnode2, isabsolute, parent1gameteFitness, parent1gameteMutators, parent1gameteMask, parent2gameteFitness, parent2gameteMutators, parent2gameteMask, popsize, pPopSize, currentvictim, wholepopulation, totalindividualgenomelength, deleteriousdistribution, wholepopulationwistree, NULL, NULL, NULL, psumofwis, NULL, NULL, 0, 0, 0, 0, NULL, NULL, miscfilepointer, globalmodifiermask, mutatorconfig.strengthfactor, deleteriousmutationrate, beneficialmutationrate);
}

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, Individual *wholepopulation, int popsize, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis, char *globalmodifiermask, MutatorConfig mutatorconfig, FILE *miscfilepointer) 
{
    int i, j;
    double haploidgenomelength = (double) ((totalpopulationgenomelength / popsize) / 2);
    int genomelength = (totalpopulationgenomelength / popsize);
    int halfgenome = genomelength / 2;

    for (i = 0; i < popsize; i++){
        wholepopulation[i].fitnessArray = malloc(sizeof(double) * genomelength);
        wholepopulation[i].mutatorArray = malloc(sizeof(int) * genomelength);
        /* Per-individual modifier mask only exists in inherited-mask mode. In
         * global-mask mode this stays NULL, saving popsize * genomelength bytes
         * (about 184 MB at popsize 20000 with a 9200-block genome). */
        wholepopulation[i].modifierMask = (mutatorconfig.maskmode == MODIFIERMASK_INHERITED)
                                        ? malloc(sizeof(char) * genomelength)
                                        : NULL;
        wholepopulation[i].fitness = 1.0;
        wholepopulation[i].mutationRate = 0.0;
        wholepopulation[i].beneficialMutationRate = 0.0;
        wholepopulation[i].mutatorCount = 0;
        wholepopulation[i].modifierCount = 0;
        wholepopulation[i].netModifierSum = 0;
        wholepopulationwistree[i] = 1.0; 
    }

    for (i = 0; i < popsize; i++) {
        j = i + LSB(i+1);
        if (j < popsize) {
            wholepopulationwistree[j] += wholepopulationwistree[i];
        }
    }
    
    *psumofwis = (long double)popsize;
    
    for(j = 0; j < popsize; j++){
        for (i = 0; i < genomelength; i++){
            wholepopulation[j].fitnessArray[i] = 0.0;
            wholepopulation[j].mutatorArray[i] = 0; 
        }
    }

    /* =====================================================================
     * MODIFIER-LOCUS INITIALISATION
     * =====================================================================
     * Two modes, chosen by mutatorconfig.maskmode.
     *
     * MODIFIERMASK_GLOBAL
     *   One mask is drawn for the whole run and mirrored across both homologs
     *   of every individual. One set of starting mutator positions is drawn and
     *   applied identically to every haplotype, so the founding population is
     *   monomorphic at every modifier locus and all variation that appears later
     *   is generated by the simulation itself.
     *
     * MODIFIERMASK_INHERITED
     *   Each founder gets its OWN mask, mirrored across that individual's two
     *   homologs, so individuals differ in which blocks carry a modifier locus
     *   and recombination can genuinely reshuffle the mask over time. Since the
     *   modifier loci then sit at different places in different individuals,
     *   "the same starting positions everywhere" is not well defined; instead
     *   exactly round(q * m) of EACH individual's own m modifier loci start at
     *   +1, so every founder carries the same NUMBER of mutator alleles at
     *   different positions.
     * ===================================================================== */
    {
        char *haploidmask = malloc(sizeof(char) * halfgenome);
        int  *haploidstates = malloc(sizeof(int) * halfgenome);

        if (mutatorconfig.maskmode == MODIFIERMASK_GLOBAL) {
            DrawModifierMask(globalmodifiermask, halfgenome, mutatorconfig.locusfraction);
            SeedInitialMutatorStates(haploidstates, globalmodifiermask, halfgenome, mutatorconfig.initialmutatorfraction, mutatorconfig.antimutatorstate);
            for (j = 0; j < popsize; j++) {
                for (i = 0; i < halfgenome; i++) {
                    wholepopulation[j].mutatorArray[i] = haploidstates[i];
                    wholepopulation[j].mutatorArray[halfgenome + i] = haploidstates[i];
                }
            }
        } else {
            /* Global mask is unused in this mode; zero it so that nothing can
             * accidentally read a stale value out of it. */
            for (i = 0; i < halfgenome; i++) globalmodifiermask[i] = 0;

            for (j = 0; j < popsize; j++) {
                DrawModifierMask(haploidmask, halfgenome, mutatorconfig.locusfraction);
                SeedInitialMutatorStates(haploidstates, haploidmask, halfgenome, mutatorconfig.initialmutatorfraction, mutatorconfig.antimutatorstate);
                for (i = 0; i < halfgenome; i++) {
                    wholepopulation[j].modifierMask[i] = haploidmask[i];
                    wholepopulation[j].modifierMask[halfgenome + i] = haploidmask[i];
                    wholepopulation[j].mutatorArray[i] = haploidstates[i];
                    wholepopulation[j].mutatorArray[halfgenome + i] = haploidstates[i];
                }
            }
        }

        free(haploidmask);
        free(haploidstates);
    }
    
    /* item 7: only mode 1 seeds the tables now. Mode 2 seeds them when the
     * burn-in ends (see RunSimulationRel); mode 0 never does. */
    if (tskitstatus == 1){
        SeedTreeSequenceTables(treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, (int) haploidgenomelength, (double) totaltimesteps);
    }
}

// ChooseVictim, ChooseParentWithTree remain identical to original
int ChooseVictim(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);
    return randomindividual;
}

int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer)
{
    long double randomnumberofbirth;
    int newparent = 0;
    randomnumberofbirth = (ldexp(pcg32_random(), -32)) * sumofwis;
    int leftbound = 0, rightbound = popsize;
    if (leftbound > rightbound) {
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
        return -1;
    }
    newparent = (SearchTree(leftbound, rightbound, randomnumberofbirth, wholepopulationwistree));
    return newparent;
}

/* ---------------------------------------------------------------------------
 * CalculateWi - COMMENTED OUT (item 10: completely unused).
 * ---------------------------------------------------------------------------
 * Nothing calls this. Wi is recomputed inline in PerformBirth() and in
 * UpdateIndividual(), both of which read the Individual struct rather than a
 * pair of raw gamete arrays. Preserved here rather than deleted.
 *
double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength)
{
    double newwi = 0.0;
    long double currentlinkageblockssum = 0.0;
    int i;
    for (i = 0; i < (totalindividualgenomelength/2); i++) {
        currentlinkageblockssum += parent1gamete[i];
        currentlinkageblockssum += parent2gamete[i];
    }
    newwi = exp(currentlinkageblockssum);
    return newwi;
}
 * ------------------------------------------------------------------------- */
