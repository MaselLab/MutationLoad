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
#include "absolute_functions.h"
#include "global_vars.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

int main(int argc, char *argv[]) {
    
    if (argc != 30) {
        printf("Incorrect number of arguments. Expected 30, got %d.\n", argc);
        return -1;
    }
    
    FILE *miscfilepointer;
    FILE *verbosefilepointer;
    FILE *finaldatafilepointer;
    FILE *veryverbosefilepointer;
    
    int Nxtimesteps, popsize, chromosomesize, numberofchromosomes, beneficialdistribution, typeofrun, randomnumberseed, K, relorabs, i_init, tskitstatus, nonmodormod, elementsperlb, snapshot, deleteriousdistribution, rawdatafilesize, calcfixation;
    double deleteriousmutationrate, bentodelmutrate, Sbtemp, slopeforcontourline, r, s, SdtoSbratio, redinmaxpopsize;
    
    // New Variables for Mutator Evolution
    double mutator_strength_factor; // f in mu = mu0 * f^n
    double mutator_switch_rate;     // Rate at which mutator loci mutate
    double mutator_bias;            // Bias towards mutators (A->M / M->A)

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
    wrong_args = AssignArgumentstoVar(argv, &Nxtimesteps, Nxtimestepsname, &popsize, popsizename, &deleteriousmutationrate, deleteriousmutationratename, &chromosomesize, chromosomesizename, &numberofchromosomes, numberofchromosomesname, &bentodelmutrate, &Sbtemp, &beneficialdistribution, &typeofrun, &slopeforcontourline, slopeforcontourlinename, &randomnumberseed, randomnumberseedname, &K, Kname, &relorabs, &r, rname, &i_init, i_initname, &s, sname, &tskitstatus, &nonmodormod, &elementsperlb, elementsperlbname, &snapshot, prevsnapshotfilename, &SdtoSbratio, SdtoSbrationame, &deleteriousdistribution, &rawdatafilesize, &redinmaxpopsize, redinmaxpopsizename, &calcfixation, &mutator_strength_factor, mutator_strength_factorname, &mutator_switch_rate, mutator_switch_ratename, &mutator_bias, mutator_biasname);

    if(wrong_args != 1){
        return -1;
    }

    bool isabsolute = (relorabs == 1);
    bool ismodular = (nonmodormod == 1);
    bool issnapshot = (snapshot == 1);
    bool isredinmaxpopsize = (redinmaxpopsize != 0.0);
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
        BracketZeroForSb(tskitstatus, isabsolute, ismodular, elementsperlb, pSb1, pSb2, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, verbosefilepointer, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        // ... (rest of bracketing logic)
    } else if (typeofrun == 1){
        if(!isabsolute){
            RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, mutator_strength_factor, mutator_switch_rate, mutator_bias);
        }else{
            // Absolute simulation call (unmodified here, ensure Absolute functions are updated if needed)
             RunSimulationAbs(issnapshot, prevsnapshotfilename, isredinmaxpopsize, redinmaxpopsizename, redinmaxpopsize, beneficialmutationratename, Sb2name, tskitstatus, ismodular, elementsperlb, isabsolute, Nxtimesteps, popsize, K, chromosomesize, numberofchromosomes, deleteriousmutationrate, Sd, deleteriousdistribution, beneficialmutationrate, Sb2, beneficialdistribution, r, i_init, s, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, iscalcfixation);
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

int AssignArgumentstoVar(char **argv, int *Nxtimesteps, char *Nxtimestepsname, int *popsize, char *popsizename, double *deleteriousmutationrate, char *deleteriousmutationratename, int *chromosomesize, char *chromosomesizename, int *numberofchromosomes, char *numberofchromosomesname, double *bentodelmutrate, double *Sbtemp, int *beneficialdistribution, int *typeofrun, double *slopeforcontourline, char *slopeforcontourlinename, int *randomnumberseed, char *randomnumberseedname, int *K, char *Kname, int *relorabs, double *r, char *rname, int *i_init, char *i_initname, double *s, char *sname, int *tskitstatus, int *nonmodormod, int *elementsperlb, char *elementsperlbname, int *snapshot, char *prevsnapshotfilename, double *SdtoSbratio, char *SdtoSbrationame, int *deleteriousdistribution, int *rawdatafilesize, double *redinmaxpopsize, char *redinmaxpopsizename, int *calcfixation, double *mutator_strength_factor, double *mutator_switch_rate, double *mutator_bias) {
    
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
    *K = atoi(argv[whicharg++]); strcpy(Kname, argv[whicharg-1]);
    *relorabs = atoi(argv[whicharg++]);
    *r = atof(argv[whicharg++]); strcpy(rname, argv[whicharg-1]);
    *i_init = atoi(argv[whicharg++]); strcpy(i_initname, argv[whicharg-1]);
    *s = atof(argv[whicharg++]); strcpy(sname, argv[whicharg-1]);
    *tskitstatus = atoi(argv[whicharg++]);
    *nonmodormod = atoi(argv[whicharg++]);
    *elementsperlb = atoi(argv[whicharg++]); strcpy(elementsperlbname, argv[whicharg-1]);
    *snapshot = atoi(argv[whicharg++]); strcpy(prevsnapshotfilename, argv[whicharg-1]);
    *SdtoSbratio = atof(argv[whicharg++]); strcpy(SdtoSbrationame, argv[whicharg-1]);
    *deleteriousdistribution = atoi(argv[whicharg++]);
    *rawdatafilesize = atoi(argv[whicharg++]);
    *redinmaxpopsize = atof(argv[whicharg++]); strcpy(redinmaxpopsizename, argv[whicharg-1]);
    *calcfixation = atoi(argv[whicharg++]);
    *mutator_strength_factor = atof(argv[whicharg++]); strcpy(mutator_strength_factorname, argv[whicharg-1]);
    *mutator_switch_rate = atof(argv[whicharg++]); strcpy(mutator_switch_ratename, argv[whicharg-1]);
    *mutator_bias = atof(argv[whicharg++]); strcpy(mutator_biasname, argv[whicharg-1]);

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

void DoubleSwap(long double * x, long double * y)
{
    long double temp = *x;
    *x = *y;
    *y = temp;
}

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


double ExponentialDerivate(double mean) {
    double result;
    float randnumb;
    do
        randnumb = ldexp(pcg32_random(), -32);
    while (randnumb == 0.0);

    result = (-log(randnumb))*mean;
    
    return result;
}

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

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, Individual *wholepopulation, int totalindividualgenomelength)
{
    int recombinationsite, startchromosome, startofindividual, h, i, returnvaluefortskit;
    startofindividual = persontorecombine * totalindividualgenomelength;
    
    tsk_id_t parentnode1 = (tsk_id_t) 2*persontorecombine;
    tsk_id_t parentnode2 = (tsk_id_t) (2*persontorecombine + 1);
        
    if (tskitstatus != 0){
        if (isabsolute){
            if (isburninphaseover != 0){
                double parenttime1 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode1]];
                double parenttime2 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode2]];    
    
                *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
                check_tsk_error(*childnode);
    
                double childtime = treesequencetablecollection->nodes.time[*childnode];
            }
        }else{
            double parenttime1 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode1]];
            double parenttime2 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode2]];    
    
            *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(*childnode);
    
            double childtime = treesequencetablecollection->nodes.time[*childnode];  
        }
    }

    for (h = 0; h < numberofchromosomes; h++) {
	    startchromosome = pcg32_boundedrand(2); //generates either a zero or a one to decide to start with chromosome 1 or 2.
	
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); //it doesn't make sense to do a recombination event before the first linkage block. Note that this will never break if the chromosome size is only one linkage block.
        
        //Tree sequence recording needs the recombination sites to add edges.
        if (tskitstatus != 0){
            if(isabsolute){
                if (isburninphaseover != 0){
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
            }else{
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
        }
        if(!ismodular){
            for (i = 0; i < recombinationsite; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize + i] = wholepopulation[persontorecombine].fitnessArray[h*chromosomesize + i];
                }
                else {
                    gamete[h*chromosomesize + i] = wholepopulation[persontorecombine].fitnessArray[totalindividualgenomelength/2 + (h*chromosomesize) + i];
                }
            }
            for (i = recombinationsite; i < chromosomesize; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize + i] = wholepopulation[persontorecombine].fitnessArray[totalindividualgenomelength/2 + (h*chromosomesize) + i];
                }
                else {
                    gamete[h*chromosomesize + i] = wholepopulation[persontorecombine].fitnessArray[startofindividual + (h*chromosomesize) + i];
                }
            }
        } else{
            for (i = 0; i < recombinationsite*elementsperlb; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulation[persontorecombine].fitnessArray[startofindividual + (h*chromosomesize*elementsperlb) + i];
                }
                else {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulation[persontorecombine].fitnessArray[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize*elementsperlb) + i];
                }
            }
            for (i = recombinationsite*elementsperlb; i < chromosomesize*elementsperlb; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulation[persontorecombine].fitnessArray[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize*elementsperlb) + i];
                }
                else {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulation[persontorecombine].fitnessArray[startofindividual + (h*chromosomesize*elementsperlb) + i];
                }
            }
        }
    }
}

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    double Sds[30];
    //Following lines stochastically generate a number of deleterious mutations drawn from a Poisson distribution with mean determined by the deleterious mutation rate
    //with effect sizes drawn from a gamma distribution with parameters taken from Kim et al 2017.
    
    // Note that deleteriousdistribution == 0 corresponds to Kim et al., == 1 corresponds to exponential,
    // and == 2 corresponds to point for th deleterious distribution. This loop also ensures that if we're
    // performing a relative run, we don't end up with any Sd's greater than or equal to 1.
    bool stayInWhileLoop = true;
    while (stayInWhileLoop) {
        // We update our looping flag to false so that the while loop will break, so long as we do not encounter
        // any lethal mutations within a relative run.
        stayInWhileLoop = false;
        numberofdeleteriousmutations = DetermineNumberOfMutations(deleteriousmutationrate);

        for (k = 0; k < numberofdeleteriousmutations; k++) {
            if (deleteriousdistribution == 0) {
                // Case for Kim et al.
                Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646); //Uses parameters for the gamma distribution of the selection coefficients of new mutations scaled to an inferred ancestral populations size. To produce the distribution of unscaled effect sizes, numbers drawn from this distribution must be divided by two times the ancestral population size for the population from which the distribution was derived (11,823 in this case). Data used to produce these fits were samples from 6503 individuals from the National Heart, Lung, and Blood Institute European-American dataset. Analysis of DFE from Kim et al. 2017.
            } else if (deleteriousdistribution == 1) {
                // Case for exponential
                Sds[k] = gsl_ran_exponential(randomnumbergeneratorforgamma, Sd);
            }  else if (deleteriousdistribution == 2) {
                // Case for point
                Sds[k] = Sd;
            }

            // We then check to see whether we encountered a lethal mutation, in which case we break
            // the inner for loop, and indicate that the outer while loop should not be broken. Note
            // that we need only check for this within a relative run
            if (!isabsolute && Sds[k] >= 1) {
                stayInWhileLoop = true;
                break;
            }
        }
    }

    

    //Adds the specified number of deleterious mutations to the gamete, recording the sites of each mutation for tree sequence recording.
    //Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the deleterious mutation effect is positive while for relative fitness the sign is negative
    for (k = 0; k < numberofdeleteriousmutations; k++) {
        if (isabsolute){
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, Sds[k]);
        }
        else{
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -Sds[k]);
        }
    }
    
    //Following lines stochastically generate a number of beneficial mutations drawn from a Poisson distribution with mean determined by the beneficial mutation rate.
    numberofbeneficialmutations = DetermineNumberOfMutations(beneficialmutationrate);
    //Adds the specified number of beneficial mutations, drawing Sb values from the specified distribution.
    //Sites of each mutation are added to the mutationsites array for tree sequence recording.
    ////Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the beneficial mutation effect is negative while for relative fitness the sign is positive
    //point distribution
    if (beneficialdistribution == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -Sb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, Sb);
            }
        }
    //exponential distribution
    } else if (beneficialdistribution == 1) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -generatedSb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, generatedSb);
            }
        }
    //uniform distribution
    } else if (beneficialdistribution == 2) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -generatedSb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, generatedSb);
            }
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        exit(0);
    }
    
    return true;
}

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
int BracketZeroForSb(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
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
    resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
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
            resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
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
            resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, mutator_switch_ratename, mutator_biasname, mutator_strength_factorname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
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
double BisectionMethodToFindSbWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
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
    slope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }
    
    slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
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
        slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
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

char * MakeDirectoryName(char * tskitstatus, char* deldist, char * isabsolutename, bool isabsolute, char * bendist, char * benmut, char * numberofchromosomes, char * chromosomesize, char * popsize, char * delmut, char * randomnumberseed, char * K, char * r, char *i_init, char * s, bool ismodular, char *elementsperlb, char *iscalcfixationname, int typeofrun, char * Sbname, char *Sdname) 
{
	
	char * directoryname = (char *) malloc(400);
	strcpy(directoryname, "datafor_");
	strcat(directoryname, isabsolutename);
    strcat(directoryname, "_tskitstatus_");
	strcat(directoryname, tskitstatus);
    strcat(directoryname, "_fixationcalc_");
	strcat(directoryname, iscalcfixationname);
    if(ismodular){
        strcat(directoryname, "_m_");
        strcat(directoryname, elementsperlb);
    }
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

