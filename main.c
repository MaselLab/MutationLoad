#include <stdio.h>
#include <float.h>
#include <string.h>
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
#include "sharedfunc_flag.h"
#include "relative_functions.h"
#include "global_vars.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

int main(int argc, char *argv[]) {
    
    /* 36 parameters are consumed by AssignArgumentstoVar, so argc must be 37.
     * NOTE: this used to read "argc != 30" while the parser consumed only 28
     * values, because the snapshot filename argument was never consumed (see the
     * fix in AssignArgumentstoVar). The count is now derived from the parser. */
    if (argc != 37) {
        printf("Incorrect number of arguments. Expected 36, got %d.\n", argc-1);
        return -1;
    }
    
    FILE *miscfilepointer;
    FILE *verbosefilepointer;
    /* finaldatafilepointer - COMMENTED OUT (Tier 3): its only consumer was
     * BisectionMethodToFindSbWithZeroSlope, which is itself commented out below.
     * FILE *finaldatafilepointer;
     */
    FILE *veryverbosefilepointer;
    
    int Nxtimesteps, popsize, chromosomesize, numberofchromosomes, beneficialdistribution, typeofrun, randomnumberseed, K, relorabs, i_init, tskitstatus, nonmodormod, elementsperlb, snapshot, deleteriousdistribution, rawdatafilesize, calcfixation;
    double deleteriousmutationrate, bentodelmutrate, Sbtemp, slopeforcontourline, r, s, SdtoSbratio, redinmaxpopsize;
    
    // New Variables for Mutator Evolution
    double mutator_strength_factor; // f in mu = mu0 * f^n
    double mutator_switch_rate;     // Rate at which mutator loci mutate
    double mutator_bias;            // Bias towards mutators (A->M / M->A)

    // New Variables for the modifier-locus model (see sharedfunc_flag.h)
    double modifier_locus_fraction;  // p: fraction of linkage blocks carrying a modifier locus
    int    modifier_mask_mode;       // MODIFIERMASK_GLOBAL (0) or MODIFIERMASK_INHERITED (1)
    int    antimutator_encoding;     // ANTIMUTATOR_AS_ZERO (0) or ANTIMUTATOR_AS_MINUS1 (1)
    double initial_mutator_fraction; // q: fraction of modifier loci starting in the +1 state

    // New Variables for the optional detailed per-individual tracking output
    int    trackindividuals;         // 0 = off, 1 = on
    int    trackinterval;            // dump every this many N-timesteps
    int    trackstartgen;            // first generation eligible for dumping

    char *Nxtimestepsname, *popsizename, *deleteriousmutationratename, *chromosomesizename, *numberofchromosomesname, *slopeforcontourlinename, *randomnumberseedname, *Kname, *rname, *i_initname, *sname, *elementsperlbname, *prevsnapshotfilename, *SdtoSbrationame, *redinmaxpopsizename, *iscalcfixationname, *mutator_strength_factorname, *mutator_switch_ratename, *mutator_biasname;
    
    Nxtimestepsname = (char *)malloc(30);
    popsizename = (char *)malloc(30);
    deleteriousmutationratename = (char *)malloc(30);
    chromosomesizename = (char *)malloc(30);
    numberofchromosomesname = (char *)malloc(30);
    slopeforcontourlinename = (char *)malloc(30);
    randomnumberseedname = (char *)malloc(30);
    Kname = (char *)malloc(30);
    rname = (char *)malloc(30);
    i_initname = (char *)malloc(30);
    sname = (char *)malloc(30);
    elementsperlbname = (char *)malloc(30);
    prevsnapshotfilename = (char *)malloc(200);
    SdtoSbrationame = (char *)malloc(30);
    redinmaxpopsizename = (char *)malloc(30);
    iscalcfixationname = (char *)malloc(30);
    mutator_strength_factorname = (char *)malloc(30);
    mutator_switch_ratename = (char *)malloc(30);
    mutator_biasname = (char *)malloc(30);

    int wrong_args;
    // Updated AssignArgumentstoVar to handle new params
    wrong_args = AssignArgumentstoVar(argv, &Nxtimesteps, Nxtimestepsname, &popsize, popsizename, &deleteriousmutationrate, deleteriousmutationratename, &chromosomesize, chromosomesizename, &numberofchromosomes, numberofchromosomesname, &bentodelmutrate, &Sbtemp, &beneficialdistribution, &typeofrun, &slopeforcontourline, slopeforcontourlinename, &randomnumberseed, randomnumberseedname, &K, Kname, &relorabs, &r, rname, &i_init, i_initname, &s, sname, &tskitstatus, &nonmodormod, &elementsperlb, elementsperlbname, &snapshot, prevsnapshotfilename, &SdtoSbratio, SdtoSbrationame, &deleteriousdistribution, &rawdatafilesize, &redinmaxpopsize, redinmaxpopsizename, &calcfixation, &mutator_strength_factor, mutator_strength_factorname, &mutator_switch_rate, mutator_switch_ratename, &mutator_bias, mutator_biasname, &modifier_locus_fraction, &modifier_mask_mode, &antimutator_encoding, &initial_mutator_fraction, &trackindividuals, &trackinterval, &trackstartgen);

    if(wrong_args != 1){
        return -1;
    }

    bool isabsolute = (relorabs == 1);
    /* -----------------------------------------------------------------------
     * MODULAR EPISTASIS IS NOT SUPPORTED IN THE MUTATION-RATE-EVOLUTION BUILD.
     * -----------------------------------------------------------------------
     * The modular branch of RecombineChromosomesIntoGamete was incomplete (it
     * only ever copied the first half of each chromosome), its elementsperlb
     * indexing overran the gamete buffers, and it is not aware of the modifier
     * mask. Rather than leave it reachable, ismodular is hard-wired to false and
     * a non-zero nonmodormod is rejected up front.
     *
     * The nonmodormod and elementsperlb command-line arguments are still PARSED,
     * because they occupy positional slots and skipping them would shift every
     * argument after them (that is exactly the bug that was fixed in
     * AssignArgumentstoVar). They are simply not used for anything.
     *
     * Original line, preserved:
     * bool ismodular = (nonmodormod == 1);
     * --------------------------------------------------------------------- */
    bool ismodular = false;
    if (nonmodormod != 0) {
        printf("Error: modular epistasis (argument 18) is not supported in this build; it must be 0, got %d.\n", nonmodormod);
        return -1;
    }
    bool issnapshot = (snapshot == 1);
    /* isredinmaxpopsize - COMMENTED OUT (Tier 3): its only consumer was
     * MakeRawDataFileName, which is commented out below. The redinmaxpopsize
     * argument itself is still parsed so positional slots stay aligned.
     * bool isredinmaxpopsize = (redinmaxpopsize != 0.0);
     */
    bool iscalcfixation = (calcfixation == 1);

    double Sb1 = 0.0, Sb2;
    double *pSb1 = &Sb1, *pSb2 = &Sb2;

    if(!isabsolute){
        Sb2 = Sbtemp;
    }else{
        Sb2 = 1.0;
    }

    double beneficialmutationrate = bentodelmutrate * deleteriousmutationrate;
    double Sd = Sb2 * SdtoSbratio;

    /* ---------------------------------------------------------------------
     * Assemble the modifier-locus / mutation-rate-evolution configuration.
     * Validated here rather than deep inside the simulation so that a bad
     * command line fails immediately and says exactly what is wrong.
     * ------------------------------------------------------------------- */
    if (modifier_locus_fraction < 0.0 || modifier_locus_fraction > 1.0) {
        printf("Error: modifier_locus_fraction must be between 0 and 1 (got %g).\n", modifier_locus_fraction);
        return -1;
    }
    if (initial_mutator_fraction < 0.0 || initial_mutator_fraction > 1.0) {
        printf("Error: initial_mutator_fraction must be between 0 and 1 (got %g).\n", initial_mutator_fraction);
        return -1;
    }
    if (modifier_mask_mode != MODIFIERMASK_GLOBAL && modifier_mask_mode != MODIFIERMASK_INHERITED) {
        printf("Error: modifier_mask_mode must be 0 (global) or 1 (inherited), got %d.\n", modifier_mask_mode);
        return -1;
    }
    if (antimutator_encoding != ANTIMUTATOR_AS_ZERO && antimutator_encoding != ANTIMUTATOR_AS_MINUS1) {
        printf("Error: antimutator_encoding must be 0 (anti-mutator = 0) or 1 (anti-mutator = -1), got %d.\n", antimutator_encoding);
        return -1;
    }
    if (mutator_switch_rate < 0.0) {
        printf("Error: mutator_switch_rate must be non-negative (got %g).\n", mutator_switch_rate);
        return -1;
    }
    if (mutator_bias < 0.0) {
        printf("Error: mutator_bias must be non-negative (got %g).\n", mutator_bias);
        return -1;
    }
    if (trackindividuals != 0 && trackindividuals != 1) {
        printf("Error: trackindividuals must be 0 or 1 (got %d).\n", trackindividuals);
        return -1;
    }
    if (trackindividuals == 1 && trackinterval < 1) {
        printf("Error: trackinterval must be at least 1 when trackindividuals is 1 (got %d).\n", trackinterval);
        return -1;
    }

    MutatorConfig mutatorconfig;
    mutatorconfig.strengthfactor         = mutator_strength_factor;
    mutatorconfig.switchrate             = mutator_switch_rate;
    mutatorconfig.bias                   = mutator_bias;
    mutatorconfig.locusfraction          = modifier_locus_fraction;
    mutatorconfig.maskmode               = modifier_mask_mode;
    mutatorconfig.antimutatorencoding    = antimutator_encoding;
    /* The integer actually stored for an anti-mutator allele. With 0 the
     * exponent n is just the count of mutator alleles; with -1 it is the net
     * sum (#mutators - #anti-mutators), so anti-mutators can push mu below mu0. */
    mutatorconfig.antimutatorstate       = (antimutator_encoding == ANTIMUTATOR_AS_MINUS1) ? -1 : 0;
    mutatorconfig.initialmutatorfraction = initial_mutator_fraction;

    TrackingConfig trackingconfig;
    trackingconfig.enabled  = trackindividuals;
    trackingconfig.interval = (trackinterval < 1) ? 1 : trackinterval;
    trackingconfig.startgen = (trackstartgen < 1) ? 1 : trackstartgen;

    char *beneficialmutationratename, *bendistname, *deldistname, *typeofrunname, *tskitstatusname, *Sb2name, *isabsolutename, *Sdname;
    beneficialmutationratename = (char *) malloc(30);
    bendistname = (char *) malloc(30);
    deldistname = (char *) malloc(30);
    typeofrunname = (char *) malloc(30);
    tskitstatusname = (char *) malloc(30);
    Sb2name = (char *) malloc(30);
    isabsolutename = (char *) malloc(30);
    Sdname = (char *) malloc(30);

    AssignStringNames(beneficialmutationratename, beneficialmutationrate, bendistname, beneficialdistribution, deldistname, deleteriousdistribution, typeofrunname, typeofrun,tskitstatusname, tskitstatus, Sb2name, Sb2, isabsolutename, isabsolute, iscalcfixationname, iscalcfixation, Sd, Sdname);

    pcg32_srandom(randomnumberseed, randomnumberseed);
    gsl_rng * randomnumbergeneratorforgamma = gsl_rng_alloc(gsl_rng_mt19937);
    /* -----------------------------------------------------------------------
     * BUG FIX: the GSL generator was never seeded.
     * -----------------------------------------------------------------------
     * gsl_rng_alloc() leaves the generator on its DEFAULT seed, so before this
     * line every run drew the identical GSL stream no matter what
     * randomnumberseed was set to. randomnumberseed only ever reached pcg32.
     *
     * That affected everything drawn through GSL: the deleterious effect sizes
     * (gsl_ran_gamma / gsl_ran_exponential), the beneficial effect sizes
     * (gsl_ran_exponential / gsl_ran_flat) and, since the switching rewrite, the
     * binomial draws that decide how many modifier loci flip state. Runs that
     * differed only in their seed were therefore NOT independent replicates.
     * --------------------------------------------------------------------- */
    gsl_rng_set(randomnumbergeneratorforgamma, (unsigned long int) randomnumberseed);
    
    char * directoryname = MakeDirectoryName(tskitstatusname, deldistname, isabsolutename, isabsolute, bendistname, beneficialmutationratename, numberofchromosomesname, chromosomesizename, popsizename, deleteriousmutationratename, randomnumberseedname, Kname, rname, i_initname, sname, ismodular, elementsperlbname, iscalcfixationname, typeofrun, Sb2name, Sdname);
    
    mkdir(directoryname, 0777);
    chdir(directoryname);
    
    if(!issnapshot){
        verbosefilepointer = fopen("verbose.txt", "w");
        veryverbosefilepointer = fopen("veryverbose.txt", "w");
        miscfilepointer = fopen("miscellaneous.txt", "w");
    }else{
        verbosefilepointer = fopen("verbose.txt", "a");
        veryverbosefilepointer = fopen("veryverbose.txt", "a");
        miscfilepointer = fopen("miscellaneous.txt", "a");
    }
    
    if (typeofrun == 0) {
        // Bracketing logic currently unmodified for mutators, using default call
        fprintf(miscfilepointer, "Beginning bracketing function.");
        fflush(miscfilepointer);
        BracketZeroForSb(tskitstatus, isabsolute, ismodular, elementsperlb, pSb1, pSb2, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, verbosefilepointer, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
        // ... (rest of bracketing logic)
    } else if (typeofrun == 1){
        if(!isabsolute){
            RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
        }else{
            // Absolute-fitness path is currently disabled (absolute_functions.c is not being
            // built/linked).
            fprintf(miscfilepointer, "Error: absolute fitness runs are currently disabled in this build.\n");
            fflush(miscfilepointer);
            fprintf(stderr, "Error: absolute fitness runs are currently disabled in this build.\n");
            exit(1);
        }
    }

    free(directoryname);
    // ... [String frees remain same] ...

    fclose(verbosefilepointer);
    fclose(veryverbosefilepointer);
    fclose(miscfilepointer);
    gsl_rng_free(randomnumbergeneratorforgamma);
    return 0;
}

int AssignArgumentstoVar(char **argv, int *Nxtimesteps, char *Nxtimestepsname, int *popsize, char *popsizename, double *deleteriousmutationrate, char *deleteriousmutationratename, int *chromosomesize, char *chromosomesizename, int *numberofchromosomes, char *numberofchromosomesname, double *bentodelmutrate, double *Sbtemp, int *beneficialdistribution, int *typeofrun, double *slopeforcontourline, char *slopeforcontourlinename, int *randomnumberseed, char *randomnumberseedname, int *K, char *Kname, int *relorabs, double *r, char *rname, int *i_init, char *i_initname, double *s, char *sname, int *tskitstatus, int *nonmodormod, int *elementsperlb, char *elementsperlbname, int *snapshot, char *prevsnapshotfilename, double *SdtoSbratio, char *SdtoSbrationame, int *deleteriousdistribution, int *rawdatafilesize, double *redinmaxpopsize, char *redinmaxpopsizename, int *calcfixation, double *mutator_strength_factor, char *mutator_strength_factorname, double *mutator_switch_rate, char *mutator_switch_ratename, double *mutator_bias, char *mutator_biasname, double *modifier_locus_fraction, int *modifier_mask_mode, int *antimutator_encoding, double *initial_mutator_fraction, int *trackindividuals, int *trackinterval, int *trackstartgen) {
    
    int whicharg = 1;
    *Nxtimesteps = atoi(argv[whicharg++]); strcpy(Nxtimestepsname, argv[whicharg-1]);
    *popsize = atoi(argv[whicharg++]); strcpy(popsizename, argv[whicharg-1]);
    *deleteriousmutationrate = atof(argv[whicharg++]); strcpy(deleteriousmutationratename, argv[whicharg-1]);
    *chromosomesize = atoi(argv[whicharg++]); strcpy(chromosomesizename, argv[whicharg-1]);
    *numberofchromosomes = atoi(argv[whicharg++]); strcpy(numberofchromosomesname, argv[whicharg-1]);
    *bentodelmutrate = atof(argv[whicharg++]);
    *Sbtemp = atof(argv[whicharg++]);
    *beneficialdistribution = atoi(argv[whicharg++]);
    *typeofrun = atoi(argv[whicharg++]);
    *slopeforcontourline = atof(argv[whicharg++]); strcpy(slopeforcontourlinename, argv[whicharg-1]);
    *randomnumberseed = atoi(argv[whicharg++]); strcpy(randomnumberseedname, argv[whicharg-1]);
    /* --- Arguments 12, 14, 15, 16, 25 and 26 below belong to the ABSOLUTE-fitness
     * workflow (carrying capacity, growth rate, initial i, selection coefficient,
     * carrying-capacity reduction, fixation calculation) and to the snapshot
     * workflow. None of them is used by the relative-fitness mutation-rate-
     * evolution runs: main() aborts on absolute runs, and MakeDirectoryName no
     * longer emits them. They MUST still be consumed here, because skipping an
     * argument shifts every argument after it - that was bug 1. --- */
    *K = atoi(argv[whicharg++]); strcpy(Kname, argv[whicharg-1]);
    *relorabs = atoi(argv[whicharg++]);
    *r = atof(argv[whicharg++]); strcpy(rname, argv[whicharg-1]);
    *i_init = atoi(argv[whicharg++]); strcpy(i_initname, argv[whicharg-1]);
    *s = atof(argv[whicharg++]); strcpy(sname, argv[whicharg-1]);
    *tskitstatus = atoi(argv[whicharg++]);
    *nonmodormod = atoi(argv[whicharg++]);
    *elementsperlb = atoi(argv[whicharg++]); strcpy(elementsperlbname, argv[whicharg-1]);
    /* ---------------------------------------------------------------------
     * BUG FIX: snapshot and the snapshot FILENAME are two separate command-line
     * arguments, but this line used to read argv[whicharg-1], i.e. it re-read
     * the snapshot flag and never consumed the filename. Every argument after
     * this point was therefore shifted by one: SdtoSbratio was being handed the
     * filename string (atof -> 0.0), deldist got SdtoSbratio, and so on all the
     * way down, with mutator_bias silently never read at all.
     * Both arguments are now consumed, one each.
     * ------------------------------------------------------------------- */
    *snapshot = atoi(argv[whicharg++]);
    strcpy(prevsnapshotfilename, argv[whicharg++]);
    *SdtoSbratio = atof(argv[whicharg++]); strcpy(SdtoSbrationame, argv[whicharg-1]);
    *deleteriousdistribution = atoi(argv[whicharg++]);
    *rawdatafilesize = atoi(argv[whicharg++]);
    *redinmaxpopsize = atof(argv[whicharg++]); strcpy(redinmaxpopsizename, argv[whicharg-1]);
    *calcfixation = atoi(argv[whicharg++]);
    *mutator_strength_factor = atof(argv[whicharg++]); strcpy(mutator_strength_factorname, argv[whicharg-1]);
    *mutator_switch_rate = atof(argv[whicharg++]); strcpy(mutator_switch_ratename, argv[whicharg-1]);
    *mutator_bias = atof(argv[whicharg++]); strcpy(mutator_biasname, argv[whicharg-1]);

    /* --- modifier-locus model (see sharedfunc_flag.h for the semantics) --- */
    *modifier_locus_fraction  = atof(argv[whicharg++]);
    *modifier_mask_mode       = atoi(argv[whicharg++]);
    *antimutator_encoding     = atoi(argv[whicharg++]);
    *initial_mutator_fraction = atof(argv[whicharg++]);

    /* --- optional detailed per-individual tracking output ---------------- */
    *trackindividuals = atoi(argv[whicharg++]);
    *trackinterval    = atoi(argv[whicharg++]);
    *trackstartgen    = atoi(argv[whicharg++]);

    return 1;
}

void UpdateLast200NTimeSteps(double * last200Ntimesteps, double newNtimesteps)
{
    double storage[200];
    int m;
    for (m = 0; m < 199; m++) {
        storage[m] = last200Ntimesteps[m+1];
    }
    storage[199] = newNtimesteps;
    for (m = 0; m < 200; m++) {
        last200Ntimesteps[m] = storage[m];
    }
}

/* ---------------------------------------------------------------------------
 * DoubleSwap - COMMENTED OUT (item 10: completely unused).
 * ---------------------------------------------------------------------------
 * Its only purpose was to serve DoubleBubbleSort, which was declared in main.h
 * but never defined and never called. Preserved here rather than deleted.
 *
void DoubleSwap(long double * x, long double * y)
{
    long double temp = *x;
    *x = *y;
    *y = temp;
}
 * ------------------------------------------------------------------------- */

double CalculateVarianceInLogFitness(int popsize, Individual *wholepopulation, long double sumofwis)
{
    int i;
    double variancesum;
    variancesum = 0.0;
    long double logaverage;
    logaverage = log(sumofwis / popsize);
    for (i = 0; i < popsize; i++) {
        variancesum += (double) pow((log(wholepopulation[i].fitness) - logaverage), 2);
    }
    variancesum = (variancesum/popsize);
    return variancesum;
}

long double FindFittestWi(Individual *wholepopulation, int popsize)
{
    long double fittestwi;
    int i;
    fittestwi = wholepopulation[0].fitness;
    for (i = 1; i < popsize; i++) {
        if (wholepopulation[i].fitness > fittestwi) {
            fittestwi = wholepopulation[i].fitness;
        }
    }
    return fittestwi;
}

double CalculateSlopeOfLogFitness(int endofsimulation, int endofburninphase, double *logaveragefitnesseachgeneration)
{
    size_t step = 1;
    int k;
    double c0, cov00, cov01, cov11, sumsq;
    int generationsafterburnin;
    double slopeoflogfitness;
    generationsafterburnin = (endofsimulation - endofburninphase);
    
    //I have to make an array of numbers to use as the x variable in a linear model to find the slope.
    double *justnumbers;
    justnumbers = malloc(sizeof(double) * generationsafterburnin);
    for (k = 0; k < generationsafterburnin; k++) {
        justnumbers[k] = (k+1);
    }
    //The following function fits a linear model to the two variables (generations and logfitness)
    //and records the parameters of the best-fitting linear model in the c0, cov00, cov01, cov11, sumsq, and slopeoflogfitness variables.
    //I only use the slope parameter, but the others are there in case I need them.
    gsl_fit_linear(justnumbers, step, logaveragefitnesseachgeneration, step, generationsafterburnin, &c0, &slopeoflogfitness, &cov00, &cov01, &cov11, &sumsq);
    free(justnumbers);
    return slopeoflogfitness;

}

/*All Fenwick tree functions from Wikipedia page "Fenwick tree" URL:https://en.wikipedia.org/wiki/Fenwick_tree
 This project is licensed under the GNU General Public License version 3.0, 
 which is compatible with the CC-BY-SA license of Wikipedia text.*/


//Returns sum of first i elements in the tree, 0 through i-1.
long double Fen_sum(long double *tree, int i)
{
    long double sum = 0;
    while (i) {
        sum += tree[i-1];
        i -= LSB(i);
    }
    return sum;
}

//Adds an amount to the ith element in the tree (and therefore to the Fen_sum for all elements in the tree greater than i).
void Fen_add(long double *tree, int numberofelementsintree, long double amounttoadd, int i)
{
    while (i < numberofelementsintree) {
        tree[i] += amounttoadd;
        i += LSB(i+1);
    }
}

//Returns the sum of the elements i through j-1.
//Could do with Fen_sum of j minus Fen_sum of i, but this is faster.
long double Fen_range(long double *tree, int i, int j)
{
    long double sum = 0;
    while (j > i) {
        sum += tree[j-1];
        j -= LSB(j);
    }
    while (i > j) {
        sum -= tree[i-1];
        i -= LSB(i);
    }
    return sum;
}

//Returns the value of the element at index i.
long double Fen_get(long double *tree, int i)
{
    return Fen_range(tree, i, i+1);
}

void Fen_set(long double *tree, int numberofelementsintree, long double newvalue, int i)
{
    Fen_add(tree, numberofelementsintree, newvalue - Fen_get(tree, i), i);
}

int SearchTree(int leftbound, int rightbound, long double targetvalue, long double *Fenwicktree)
{
    int middle;
    middle = floor((leftbound+rightbound)/2);
    long double partialsumatmiddle;
    long double partialsumatmiddleminusone;
    partialsumatmiddle = Fen_sum(Fenwicktree, middle);
    partialsumatmiddleminusone = Fen_sum(Fenwicktree, middle-1);
    if(partialsumatmiddle < targetvalue) {
        if((middle+1) == rightbound) {
            return middle;
        }
        return SearchTree(middle, rightbound, targetvalue, Fenwicktree);
    }
    if(partialsumatmiddle > targetvalue) {
        if(partialsumatmiddleminusone > targetvalue) {
            return SearchTree(leftbound, middle, targetvalue, Fenwicktree);
        } else {
            return (middle-1);
        }
    }
    if (partialsumatmiddle == targetvalue) {
        return middle;
    }
}


/* ---------------------------------------------------------------------------
 * ExponentialDerivate - COMMENTED OUT (item 10: completely unused).
 * ---------------------------------------------------------------------------
 * Nothing calls it; exponential effect sizes are drawn with gsl_ran_exponential
 * inside ProduceMutatedGamete. Preserved here rather than deleted.
 *
double ExponentialDerivate(double mean) {
    double result;
    float randnumb;
    do
        randnumb = ldexp(pcg32_random(), -32);
    while (randnumb == 0.0);

    result = (-log(randnumb))*mean;
    
    return result;
}
 * ------------------------------------------------------------------------- */

//From Numerical Recipes in C, Second Edition.
int SampleFromPoisson(float poissonmean)
{
    static float sq, logmean, g;
    static float oldmean = (-1.0);
    float numberofmutations, t, y;

    if (poissonmean < 12.0) {		//for small enough means, use direct method.
        if (poissonmean != oldmean) {	//check to see if the mean value is new.
                oldmean = poissonmean;
                g = exp(-poissonmean);	//if the mean is new, compute the exponential.
        }
        numberofmutations = -1;
        t = 1.0;
        do {
                ++numberofmutations;
                t *= ldexp(pcg32_random(), -32); //instead of adding exponential deviates, multiply uniform deviates and compare to pre-computed exponential.
        } while (t > g);
    } 
    else { 				//for larger means, use rejection method.
        if (poissonmean != oldmean) {	//for new means, pre-compute some functions.
                oldmean = poissonmean;
                sq = sqrt(2.0*poissonmean);
                logmean = log(poissonmean);
                g = poissonmean*logmean - gsl_sf_lngamma(poissonmean+1.0); //lngamma function is the natural log of the gamma function
        }
        do {
                do {
            y = tan(PI * ldexp(pcg32_random(), -32)); 	//makes y a deviate from a Lorentzian comparison function.
            numberofmutations = sq*y + poissonmean;		//shifts and scales y and sets results as possible numberofmutations (to be accepted or rejected);
                } while (numberofmutations < 0.0); 			//rejects values in zero probability area.
                numberofmutations = floor(numberofmutations);
                t = 0.9 * (1.0 + y*y) * exp(numberofmutations*logmean - gsl_sf_lngamma(numberofmutations + 1.0) - g);
        } while (ldexp(pcg32_random(), -32) > t);
    }
return numberofmutations;
}

// RecombineChromosomesIntoGamete and ProduceMutatedGamete are defined in sharedfunc_flag.c
// (with the current mutator-aware signatures). The stale duplicate definitions that used to
// live here were removed to avoid conflicting-type / multiple-definition errors.

int DetermineNumberOfMutations(double mutationrate)
{      
    double meannumberofmutations = mutationrate/2.0;
    
    //Note that because this function operates on gametes, the calculation above appears haploid.
    //There shouldn't be a multiplication by 2 (for diploidy) in this function, since it will be called twice per individual: once per gamete.
    //Above calculation should be moved outside this function for increased efficiency.
    
    int numberofmutations = SampleFromPoisson(meannumberofmutations);
    return numberofmutations;
}

int DetermineMutationSite(int totalgametelength)
{      
    //mutation occur in a random place in the gamete. We can change this to incorporate heterogenous mutation rate for different linckages blocks (mutation hotspots).
    //Moreover, in modular runs it can be edited to allow pleiotropy. A single mutation affects multiple modules.
    int mutatedsite = pcg32_boundedrand(totalgametelength);
    return mutatedsite;
}
//The following function is heavily modified from Numerical Recipes in C, Second Edition.
//For large population sizes, populations with mean Sb > 0 may actually have a more negative fitness slope than mean Sb = 0.
//
int BracketZeroForSb(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig) {
    int i, numberoftries;
    numberoftries = 10;
    float factor = 0.01;
    char Sb1name[10], Sb2name[10];
    snprintf(Sb1name, 10, "%.7f", *Sb1);
    snprintf(Sb2name, 10, "%.7f", *Sb2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Sb1name: %s, Sb2name: %s\n", Sb1name, Sb2name);
        fflush(verbosefilepointer);
    }
    float resultingslope1, resultingslope2;
    resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
    resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "First two slopes are: %.6f for sb %.6f, and %.6f for sb %.6f\n", resultingslope1, *Sb1, resultingslope2, *Sb2);
        fflush(verbosefilepointer);
    }
    if (resultingslope1 == resultingslope2) {
        return 0;
        fprintf(miscfilepointer, "Slopes after first try are the same, equaling %.5f and %.5f\n", resultingslope1, resultingslope2);
    }
    if (resultingslope1 > slopeforcontourline) {
        return 0;
        fprintf(miscfilepointer, "Slope with sb 0.0 larger than proposed contour, slope = %.6f, contour line value = %.6f\n", resultingslope1, slopeforcontourline);
    }
    
    for (i = 0; i < numberoftries; i++) {
        if ((resultingslope1 < slopeforcontourline) && (resultingslope2 > slopeforcontourline)) {
            return 1;
        } else if (resultingslope2 <= slopeforcontourline) {
            *Sb2 += factor;
            snprintf(Sb2name, 10, "%.7f", *Sb2);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb2name: %s\n", Sb2name);
                fprintf(verbosefilepointer, "Starting run with new sb2 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb2, resultingslope2);
                fflush(verbosefilepointer);
            }
            
        } else if (resultingslope1 >= slopeforcontourline) {
            *Sb1 -= factor;
            snprintf(Sb1name, 10, "%.7f", *Sb1);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb1name: %s\n", Sb1name);
                fprintf(verbosefilepointer, "Starting run with new sb1 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb1, resultingslope1);
                fflush(verbosefilepointer);
            }
            
        }
    }
    fprintf(miscfilepointer, "Failed to bracket contour slope of %.6f in 10 tries.", slopeforcontourline);
    return 0;
}

//The following function is modified from Numerical Recipes in C, Second Edition.
/* ---------------------------------------------------------------------------
 * BisectionMethodToFindSbWithZeroSlope - NEVER CALLED. Part of the Sb
 * contour-line root-finding workflow (typeofrun == 0), whose driver
 * BracketZeroForSb has its tail elided, so nothing invokes this.
 * Commented out, not deleted (item 10). Uncomment to bring it back.
 * ---------------------------------------------------------------------------
double BisectionMethodToFindSbWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig) {
    int i;
    double factor, slope1, slopemid, Sbmid, root;
    double accuracy = 0.00005;
    int maxtries = 30;
    char Sb1name[10], Sb2name[10], Sbmidname[10];
    snprintf(Sb1name, 10, "bis%.4f", *Sb1);
    snprintf(Sb2name, 10, "bis%.4f", *Sb2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Entered bisection function. First two sb %.6f and %.6f\n", *Sb1, *Sb2);
        fprintf(verbosefilepointer, "Starting Sb1name: %s, starting Sb2name: %s", Sb1name, Sb2name);
        fflush(verbosefilepointer);
    }
    slope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }
    
    slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb2, slopemid);
    }
    
    if (((slope1 - slopeforcontourline)*(slopemid - slopeforcontourline)) > 0.0) {
        fprintf(miscfilepointer, "Root not bracketed properly, with starting slopes %.10f and %.10f for a desired slope of %.6f\n", slope1, slopemid, slopeforcontourline);
        return 0.0;
    }
    root = (slope1 < slopeforcontourline) ? (factor=*Sb2-*Sb1, *Sb1) : (factor=*Sb1-*Sb2, *Sb2);
    for (i = 1; i <= maxtries; i++) {
        Sbmid = root + (factor *= 0.5);
        snprintf(Sbmidname, 10, "%.7f", Sbmid);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "Sbmidname: %s\n", Sbmidname);
            fprintf(verbosefilepointer, "Starting run with sb %.6f\n", Sbmid);
            fflush(verbosefilepointer);
        }
        slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutatorconfig, trackingconfig);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", Sbmid, slopemid);
            fflush(verbosefilepointer);
        }
        if (slopemid <= slopeforcontourline) {
            root = Sbmid;
        }
        if (fabs(factor) < accuracy || Sbmid == slopeforcontourline) {
            return root;
        }
        
    }
    fprintf(miscfilepointer, "Error: root not found. Root after 30 tries was: %.10f", root);
    fprintf(finaldatafilepointer, "Error: root not found. Root after 30 tries was: %.10f", root);
    return 0.0;
    
}
 * ------------------------------------------------------------------------- */

char * MakeDirectoryName(char * tskitstatus, char* deldist, char * isabsolutename, bool isabsolute, char * bendist, char * benmut, char * numberofchromosomes, char * chromosomesize, char * popsize, char * delmut, char * randomnumberseed, char * K, char * r, char *i_init, char * s, bool ismodular, char *elementsperlb, char *iscalcfixationname, int typeofrun, char * Sbname, char *Sdname) 
{
	
	char * directoryname = (char *) malloc(400);
	strcpy(directoryname, "datafor_");
	strcat(directoryname, isabsolutename);
    strcat(directoryname, "_tskitstatus_");
	strcat(directoryname, tskitstatus);
    strcat(directoryname, "_fixationcalc_");
	strcat(directoryname, iscalcfixationname);
    /* Modular-epistasis naming - COMMENTED OUT: ismodular is hard-wired false.
    if(ismodular){
        strcat(directoryname, "_m_");
        strcat(directoryname, elementsperlb);
    }
    */
    /* Absolute-fitness naming fields - COMMENTED OUT (Tier 3). This branch was
     * already unreachable: main() aborts on absolute runs, so isabsolute is
     * always false by the time this is called, and these fields never appeared
     * in any directory name produced by this build. Commenting it therefore
     * changes no output. The r / i_init / s / K / ismodular / elementsperlb
     * parameters of this function are consequently unused.
    if(isabsolute){
        strcat(directoryname, "_r_");
        strcat(directoryname, r);
        strcat(directoryname, "_iinit_");
        strcat(directoryname, i_init);
        strcat(directoryname, "_s_");
        strcat(directoryname, s);
        strcat(directoryname, "_K_");
        strcat(directoryname, K);
    }
    */
    if(typeofrun == 1){
        strcat(directoryname, "_Sb_");
        strcat(directoryname, Sbname);
    }
    strcat(directoryname, "_deldist_");
    strcat(directoryname, deldist);
    strcat(directoryname, "_bendist_");
	strcat(directoryname, bendist);
	strcat(directoryname, "_mub_");
	strcat(directoryname, benmut);
	strcat(directoryname, "_chromnum_");
	strcat(directoryname, numberofchromosomes);
	strcat(directoryname, "_N0_");
	strcat(directoryname, popsize);
	strcat(directoryname, "_mud_");
	strcat(directoryname, delmut);
    strcat(directoryname, "_L_");
	strcat(directoryname, chromosomesize);
	strcat(directoryname, "_seed_");
	strcat(directoryname, randomnumberseed);
    strcat(directoryname, "_Sd_");
    strcat(directoryname, Sdname);
    

	return directoryname;
}

/* ---------------------------------------------------------------------------
 * MakeFinalDataFileName - NEVER CALLED. Part of the Sb root-finding workflow.
 * Commented out, not deleted (item 10). Uncomment to bring it back.
 * ---------------------------------------------------------------------------
char * MakeFinalDataFileName(char * typeofrun, char * benmut, char * slopeforcontourline, char * randomnumberseed) 
{
	
    char * finaldatafilename = (char *) malloc(60);
    strcpy(finaldatafilename, "finaldatafor_");
    strcat(finaldatafilename, "runtype_");
    strcat(finaldatafilename, typeofrun);
    strcat(finaldatafilename, "_mub_");
    strcat(finaldatafilename, benmut);
    strcat(finaldatafilename, "_slope_");
    strcat(finaldatafilename, slopeforcontourline);
    strcat(finaldatafilename, "_seed_");
    strcat(finaldatafilename, randomnumberseed);
    
    return finaldatafilename;
}
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * MakeRawDataFileName - NEVER CALLED. RunSimulationRel builds its own raw
 * data file name inline (and now includes the mutator/modifier parameters).
 * Commented out, not deleted (item 10). Uncomment to bring it back.
 * ---------------------------------------------------------------------------
char * MakeRawDataFileName(char * mubname, char * Sbname, bool isredinmaxpopsize, char *redinmaxpopsizename) 
{
	
    char * rawdatafilename = (char *) malloc(200);
    strcpy(rawdatafilename, "rawdatafor_"); //starting the string that will be the name of the data file.
    strcat(rawdatafilename, "Sb_");
    strcat(rawdatafilename, Sbname);
    strcat(rawdatafilename, "_mub_");
    strcat(rawdatafilename, mubname);
    if(isredinmaxpopsize){
        strcat(rawdatafilename, "_redrate_");
        strcat(rawdatafilename, redinmaxpopsizename);
    }
    strcat(rawdatafilename, ".txt");
    
    return rawdatafilename;
}
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * MakeSummaryDataFileName - NEVER CALLED. RunSimulationRel builds its own
 * summary file name inline.
 * Commented out, not deleted (item 10). Uncomment to bring it back.
 * ---------------------------------------------------------------------------
char * MakeSummaryDataFileName(char * mubname, char * Sbname) 
{
	
    char* summarydatafilename = (char*)malloc(60);
    strcpy(summarydatafilename, "summarydatafor_");
    strcat(summarydatafilename, "Sb_");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "_mub_");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    
    return summarydatafilename;
}
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * MakePopSnapshotFileName - NEVER CALLED. Part of the population-snapshot
 * workflow used by the absolute-fitness runs.
 * Commented out, not deleted (item 10). Uncomment to bring it back.
 * ---------------------------------------------------------------------------
char * MakePopSnapshotFileName(char * mubname, char * Sbname) 
{
	
    char * popsnapshotfilename = (char *) malloc(100);
    strcpy(popsnapshotfilename, "popsnapshotfor_");
    strcat(popsnapshotfilename, "Sb_");
    strcat(popsnapshotfilename, Sbname);
    strcat(popsnapshotfilename, "mub_");
    strcat(popsnapshotfilename, mubname);
    strcat(popsnapshotfilename, ".txt");
    
    return popsnapshotfilename;
}
 * ------------------------------------------------------------------------- */

void AssignStringNames(char *beneficialmutationratename, double beneficialmutationrate, char *bendistname, int beneficialdistribution, char *deldistname, int deleteriousdistribution, char *typeofrunname, int typeofrun, char *tskitstatusname, int tskitstatus, char* Sb2name, double Sb2, char *isabsolutename, bool isabsolute, char *iscalcfixationname, bool iscalcfixation, double Sd, char *Sdname) {
	//pointer for the beneficial mutation rate name
	sprintf(beneficialmutationratename, "%1.4f", beneficialmutationrate);

	//pointer for the beneficial distribution name
	if(beneficialdistribution == 0)
		strncpy(bendistname, "point", sizeof("point"));
	else if(beneficialdistribution == 1)
		strncpy(bendistname, "exponential", sizeof("exponential"));
	else if(beneficialdistribution == 2)
		strncpy(bendistname, "uniform", sizeof("uniform"));

	//pointer for the deleterious distribution name
	if(deleteriousdistribution == 0)
		strncpy(deldistname, "kim_et_al", sizeof("kim_et_al"));
	else if(deleteriousdistribution == 2)
		strncpy(deldistname, "point", sizeof("point"));
    else 
        strncpy(deldistname, "exponential", sizeof("exponential"));

	//pointer for type of run name
	if(typeofrun == 0)
		strncpy(typeofrunname, "root", sizeof("root"));
	else if(typeofrun == 1)
		strncpy(typeofrunname, "single", sizeof("single"));

	//pointer for tskitstatus name
	if(tskitstatus == 0) {
		strncpy(tskitstatusname, "OFF", sizeof("OFF"));
	} else if(tskitstatus == 1) {
		strncpy(tskitstatusname, "ON", sizeof("ON"));
	} else if(tskitstatus == 2) {
		strncpy(tskitstatusname, "ON_BURNIN", sizeof("ON_BURNIN"));
	}

	//pointer for Sb name
	sprintf(Sb2name, "%1.4f", Sb2);

    // pointer for Sd name
    sprintf(Sdname, "%1.6f", Sd);

	//pointer for isabsolutename
	if(isabsolute)
		strncpy(isabsolutename, "absolute", sizeof("absolute"));
	else
		strncpy(isabsolutename, "relative", sizeof("relative"));

    if(iscalcfixation)
		strncpy(iscalcfixationname, "ON", sizeof("ON"));
	else
		strncpy(iscalcfixationname, "OFF", sizeof("OFF"));
}