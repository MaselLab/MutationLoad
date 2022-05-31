/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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

int main(int argc, char *argv[]) {

    if (argc != 19) {
        printf("[Error]; Wrong number of arguments in program. It should be timeSteps; initialPopsize; genome-wide deleterious mutation rate; chromosomesize; numberofchromosomes; beneficial/deleterious mutation ratio; Sb; point, exponential or unifrom distribution; root or single; slope; seed; MaxPopSize; relative or absolute; d_0; epistasis r; sdmin; without modular epistasis or with it; elements per linkage block \n");
        return -1;
    }
    //declare the file pointers for the files used for printing across the program
    FILE *miscfilepointer;
    FILE *verbosefilepointer;
    FILE *finaldatafilepointer;
    FILE *veryverbosefilepointer;
    //number of time steps which crudely translate to generations. So far, 1000 is the value used
    int Nxtimesteps = atoi(argv[1]);
    //initial population size
    int popsize = atoi(argv[2]);
    //deleterious mutation rate; remember that this is the genome-wide mutation rate
    double deleteriousmutationrate = atof(argv[3]);
    //size of a single chromosome, i.e. number of blocks within a single chromosome. so far, we use 50 blocks
    int chromosomesize = atoi(argv[4]);
    //number of chromosomes; remember that this is total number of chromosomes, not ploidy -- all individuals will be diploid; to mimic a human genome, we use 23 chromosomes
    int numberofchromosomes = atoi(argv[5]); //
    //given that we don't have an exact idea what the beneficial mutation rate is in humans, we use instead a beneficial/deleterious mutation ratio.
    double bentodel = atof(argv[6]);
    
    char *delmutname = (char *) malloc(30);//create a pointer for the deleterious mutation rate
    sprintf(delmutname, "%1.1f", deleteriousmutationrate);
    //calculate the beneficial mutation rate from the inputted ratio
    double beneficialmutationrate = bentodel*deleteriousmutationrate;
    
    char *benmutname = (char *) malloc(30);//create a pointer for the beneficial mutation rate
    sprintf(benmutname, "%1.4f", beneficialmutationrate);
    //I have two parameters for Sb for the type of run that needs to have bracketed values of Sb.
    //In the case with just a single simulation being run, Sb2 here will be the value of Sb used.
    double Sb2 = atof(argv[7]);
    double *pSb2 = &Sb2;
    
    double Sb1 = 0.0;
    double *pSb1 = &Sb1;
    
    //0 for point; 1 for exponential; 2 for unifrom
    int beneficialdistribution = atoi(argv[8]);
    char *bendistname = (char *) malloc(30);
    if(beneficialdistribution == 0)
        strncpy(bendistname, "point", sizeof("point"));
    else if(beneficialdistribution == 1)
        strncpy(bendistname, "exponential", sizeof("exponential"));
    else if(beneficialdistribution == 2)
        strncpy(bendistname, "uniform", sizeof("uniform"));
    else
        strncpy(bendistname, "NA", sizeof("NA"));
    
    //0 is root; 1 is single
    int typeofrun = atoi(argv[9]);
    char *typeofrunname = (char *) malloc(30);
    if(typeofrun == 0)
        strncpy(typeofrunname, "root", sizeof("root"));
    else if(typeofrun == 1)
        strncpy(typeofrunname, "single", sizeof("single"));
    else
        strncpy(typeofrunname, "NA", sizeof("NA"));
    
    //the slope for the contour line; so far it is used in the relative scheme
    double slopeforcontourline = atof(argv[10]);
    // the see for the random number generators used in this program *don't forget to be consistent with the seed*
    int randomnumberseed = atoi(argv[11]);
    //the maximum population size K or carrying capacity. so far, we used 20000 as a value for K.
    int maxPopSize = atoi(argv[12]);
    if(maxPopSize < popsize){
        printf("[Error]maxPopSize is smaller than initial popsize \n");
        return -1;
    }
    
    //0 for relative simulation; 1 for absolute
    int relorabs = atoi(argv[13]);
    bool isabsolute;
    char *isabsolutename = (char *) malloc(30);
    if(relorabs == 0){
        isabsolute = false;
        strncpy(isabsolutename, "relative", sizeof("relative"));
    }
    else if(relorabs == 1){
        isabsolute = true;
        strncpy(isabsolutename, "absolute", sizeof("absolute"));
    }
    else{
        printf("[Error] 12th argument must be 0 or 1 \n");
        return -1;
    }
    //the initial death rate; so far, we test with two values: 0.5 and 0.2
    double d_0 = atof(argv[14]);
    if (d_0 >= 1 || d_0 <= 0){
        printf("[Error] initial death rate d0 should vary between 0 and 1");
        return -1;
    }   
    //the epistasis variable R; R varies between 0 and 1, 1 representing no epistasis and 0 is full epistasis
    double r = atof(argv[15]);
    if (r < 0 || r > 1){
        printf("[Error] epistasis degree variable R must be between 0 and 1");
        return -1;
    } 
    //the selection coefficient of a beneficial mutation at minimum death rate (i.e. max birth rate). so far we use the arbitrary value 0.01 
    double sdmin = atof(argv[16]);
    
    //Indicate whether epistasis is modular or not. 0 for non modular epistasis and 1 for modular epistasis. 
    int nonmodormod = atoi(argv[17]);
    bool ismodular;
    if(nonmodormod == 0){
        ismodular = false;
    }
    else if(nonmodormod == 1){
        ismodular = true;
    }
    else{
        printf("[Error] 17th argument (non modular epistasis or modular epistasis) must be 0 or 1 \n");
        return -1;
    }
    
    //Number of elements per linkage block
    int elementsperlb = atoi(argv[18]);
    
    pcg32_srandom(randomnumberseed, randomnumberseed); // seeds the random number generator.
    gsl_rng * randomnumbergeneratorforgamma = gsl_rng_alloc(gsl_rng_mt19937);
    //the gamma distribution function requires a gsl random number generator, which is set here.
    //it's a bit inelegant to have two different RNGs, which could be solved by using a different algorithm 
    //for choosing variates from a gamma distribution, instead of using the free one from gsl.
    
    char * directoryname = MakeDirectoryName(argv[1], argv[2], delmutname, argv[4], argv[5], benmutname, bendistname, argv[11], argv[12], isabsolute, isabsolutename, argv[14], argv[15], argv[16], ismodular, argv[18]);// this will create the directory name pointer using input parameter values
    
    mkdir(directoryname, 0777);//create the directory with the directory name pointer
    chdir(directoryname);//move into the created directory
    
    verbosefilepointer = fopen("verbose.txt", "w");	//opens the file to which to print verbose data.
    veryverbosefilepointer = fopen("veryverbose.txt", "w"); //opens the file to which to print very verbose data.
    miscfilepointer = fopen("miscellaneous.txt", "w"); //opens the file to which to print miscellaneous data.
    
   
    if (typeofrun == 0) {
        
        char * finaldatafilename = MakeFinalDataFileName(typeofrunname, benmutname, argv[10], argv[11]);//create the final data file name using the input parameter values
        finaldatafilepointer = fopen(finaldatafilename, "w");
        /*This type of run finds the Sb value for the given set of parameters
         that produces a population whose fitness stays almost exactly stable.
         It does this by finding values of Sb that lead to populations definitely
         increasing and definitely decreasing in fitness,
         and then searching between them until it finds a value of Sb that leads
         to a population with a long-term slope of fitness that is within an error term of zero.
         */
        double sbrequiredforzeroslopeoffitness;
        fprintf(miscfilepointer, "Beginning bracketing function.");
        fflush(miscfilepointer);
        BracketZeroForSb(isabsolute, pSb1, pSb2, argv[1], argv[2], delmutname, argv[4], argv[5], benmutname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, verbosefilepointer, miscfilepointer, veryverbosefilepointer);
        fprintf(miscfilepointer, "Finished bracketing function.");
        fflush(miscfilepointer);
        sbrequiredforzeroslopeoffitness = BisectionMethodToFindSbWithZeroSlope(isabsolute, pSb1, pSb2, argv[1], argv[2], delmutname, argv[4], argv[5], benmutname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, verbosefilepointer, finaldatafilepointer, veryverbosefilepointer);
        fprintf(finaldatafilepointer, "The value of Sb for which the slope of log fitness is zero with mub of %.10f is %.10f", beneficialmutationrate, sbrequiredforzeroslopeoffitness);
        
        free(finaldatafilename);
        fclose(finaldatafilepointer);
    
    } else if (typeofrun == 1) {
        if(isabsolute == 0){
            //This type of run just simulates a single population with the input parameters.
            RunSimulationRel(isabsolute, argv[1], argv[2], delmutname, argv[4], argv[5], benmutname, argv[7], typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
        }
        
        else{
            RunSimulationAbs(isabsolute, argv[1], argv[2], delmutname, argv[4], argv[5], benmutname, argv[7], typeofrun, Nxtimesteps, popsize, maxPopSize, d_0, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, randomnumbergeneratorforgamma, r, sdmin, miscfilepointer, veryverbosefilepointer);
        }
        
    } else{
        //One day maybe I'll have more types of runs.
        fprintf(miscfilepointer, "That type of run is not currently supported.");
        return -1;
    }
    
    free(directoryname);
    free(bendistname);
    free(typeofrunname);
    free(delmutname);
    free(benmutname);
    free(isabsolutename);
    
    fclose(verbosefilepointer);
    fclose(veryverbosefilepointer);
    fclose(miscfilepointer); 
    
    gsl_rng_free(randomnumbergeneratorforgamma);
    
    return 0;

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

//Inefficient algorithm (guaranteed to be order n^2). Improve algorithm if using more than once per simulation.
void DoubleBubbleSort(long double *arraytobesorted, int arraysize)
{
    int i, j;
    for (i = 0; i < arraysize-1; i++) {
        for (j = 0; j < arraysize-i-1; j++) {
            if (arraytobesorted[j] > arraytobesorted[j+1]) {
                DoubleSwap(&arraytobesorted[j], &arraytobesorted[j+1]);
            }
        }
    }
}

double CalculateVarianceInLogFitness(int popsize, long double *wholepopulationwisarray, long double sumofwis)
{
    int i;
    double variancesum;
    variancesum = 0.0;
    long double logaverage;
    logaverage = log(sumofwis / popsize);
    for (i = 0; i < popsize; i++) {
        variancesum += (double) pow((log(wholepopulationwisarray[i]) - logaverage), 2);
    }
    variancesum = (variancesum/popsize);
    return variancesum;
}

long double FindFittestWi(long double *wisarray, int popsize)
{
    long double fittestwi;
    int i;
    fittestwi = wisarray[0];
    for (i = 1; i < popsize; i++) {
        if (wisarray[i] > fittestwi) {
            fittestwi = wisarray[i];
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
void RecombineChromosomesIntoGamete(int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, double *wholepopulationgenomes, int totalindividualgenomelength)
{    
    int recombinationsite, startchromosome, startofindividual, h, i;
    startofindividual = persontorecombine * totalindividualgenomelength;
    
    
    for (h = 0; h < numberofchromosomes; h++) {
	startchromosome = pcg32_boundedrand(2); //generates either a zero or a one to decide to start with chromosome 1 or 2.
	
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); //it doesn't make sense to do a recombination event before the first linkage block. Note that this will never break if the chromosome size is only one linkage block.
        
        
        for (i = 0; i < recombinationsite; i++) {
            if (startchromosome == 0) {
                gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize) + i];
            }
            else {
                gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
            }
        }
        for (i = recombinationsite; i < chromosomesize; i++) {
            if (startchromosome == 0) {
                gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
            }
            else {
                gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize) + i];
            }
        }
    }
}

bool ProduceMutatedGamete(bool isabsolute, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double *gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    double Sds[30];
    
    //Following lines stochastically generate a number of deleterious mutations drawn from a Poisson distribution with mean determined by the deleterious mutation rate
    //with effect sizes drawn from a gamma distribution with parameters taken from Kim et al 2017.
    bool DontBreakWhileLoop = false;
    
    while(1){
        DontBreakWhileLoop = false;
        numberofdeleteriousmutations = DetermineNumberOfMutations(chromosomesize, numberofchromosomes, deleteriousmutationrate);
        
        for (k = 0; k < numberofdeleteriousmutations; k++) {
            Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646); //Uses parameters for the gamma distribution of the selection coefficients of new mutations scaled to an inferred ancestral populations size. To produce the distribution of unscaled effect sizes, numbers drawn from this distribution must be divided by two times the ancestral population size for the population from which the distribution was derived (11,823 in this case). Data used to produce these fits were samples from 6503 individuals from the National Heart, Lung, and Blood Institute European-American dataset. Analysis of DFE from Kim et al. 2017.
            //This gamma distribution can occasionally produce deleterious mutations with effect sizes larger than 1,
            //which would result in a gamete with fitness less than zero, which would break my algorithm.
            //The following if statement simply throws out any deleterious mutations with effect sizes larger than 1.
            if (Sds[k] >= 1) {
                //this mean a new lethal mutation appeared
                DontBreakWhileLoop = true;
                break;
            }
            //store the number of mutations excluded
        }
        //in absolute runs a lethal mutation means the new offspring 
        if(isabsolute && DontBreakWhileLoop){
            return false;
        }
        if (!DontBreakWhileLoop) 
            break;
    }
    
    //Adds the specified number of deleterious mutations to the gamete, recording the sites of each mutation for tree sequence recording.
    //Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the deleterious mutation effect is positive while for relative fitness the sign is negative
    for (k = 0; k < numberofdeleteriousmutations; k++) {
        if (isabsolute)
            MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, Sds[k]);
        else
            MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, -Sds[k]);
    }
    
    //Following lines stochastically generate a number of beneficial mutations drawn from a Poisson distribution with mean determined by the beneficial mutation rate.
    numberofbeneficialmutations = DetermineNumberOfMutations(chromosomesize, numberofchromosomes, beneficialmutationrate);
    //Adds the specified number of beneficial mutations, drawing Sb values from the specified distribution.
    //Sites of each mutation are added to the mutationsites array for tree sequence recording.
    ////Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the beneficial mutation effect is negative while for relative fitness the sign is positive
    //point distribution
    if (beneficialdistribution == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            if (isabsolute)
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, -Sb);
            else
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, Sb); 
        }
    //exponential distribution
    } else if (beneficialdistribution == 1) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            if (isabsolute)
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, -generatedSb);
            else
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, generatedSb);
        }
    //uniform distribution
    } else if (beneficialdistribution == 2) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            if (isabsolute)
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, -generatedSb);
            else
                MutateGamete(isabsolute, chromosomesize, numberofchromosomes, gamete, generatedSb);
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        exit(0);
    }
    
    return true;
}

int DetermineNumberOfMutations(int chromosomesize, int numberofchromosomes, double mutationrate)
{      
    double meannumberofmutations = mutationrate/2.0;
    
    //Note that because this function operates on gametes, the calculation above appears haploid.
    //There shouldn't be a multiplication by 2 (for diploidy) in this function, since it will be called twice per individual: once per gamete.
    //Above calculation should be moved outside this function for increased efficiency.
    
    int numberofmutations = SampleFromPoisson(meannumberofmutations);
    return numberofmutations;
}

//The following function is heavily modified from Numerical Recipes in C, Second Edition.
//For large population sizes, populations with mean Sb > 0 may actually have a more negative fitness slope than mean Sb = 0.
//
int BracketZeroForSb(bool isabsolute, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer) {
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
    resultingslope1 = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
    resultingslope2 = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
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
            resultingslope2 = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
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
            resultingslope1 = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb2, resultingslope2);
                fflush(verbosefilepointer);
            }
            
        }
    }
    fprintf(miscfilepointer, "Failed to bracket contour slope of %.6f in 10 tries.", slopeforcontourline);
    return 0;
}

//The following function is modified from Numerical Recipes in C, Second Edition.
double BisectionMethodToFindSbWithZeroSlope(bool isabsolute, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer) {
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
    slope1 = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }
    
    slopemid = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
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
        slopemid = RunSimulationRel(isabsolute, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer);
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
char * MakeDirectoryName(char * timeSteps, char * popsize, char * delmut, char * chromosomesize, char * numberofchromosomes, char * benmut, char * bendist, char * randomnumberseed, char * maxPopSize, bool isabsolute, char * isabsolutename, char * d_0, char * r, char * sdmin, bool ismodular, char *elementsperlb) 
{
	
	char * directoryname = (char *) malloc(200);
	strcpy(directoryname, "datafor_");
	strcat(directoryname, isabsolutename);
    strcat(directoryname, "_");
    if(ismodular){
        strcat(directoryname, "elementsperlb_");
        strcat(directoryname, elementsperlb);
    }
    if(isabsolute){
        strcat(directoryname, "d0_");
        strcat(directoryname, d_0);
        strcat(directoryname, "r_");
        strcat(directoryname, r);
        strcat(directoryname, "sdmin_");
        strcat(directoryname, sdmin);
        strcat(directoryname, "K_");
        strcat(directoryname, maxPopSize);
    }
    strcat(directoryname, "mud_");
	strcat(directoryname, delmut);
    strcat(directoryname, bendist);
    strcat(directoryname, "_");
    strcat(directoryname, "mub_");
	strcat(directoryname, benmut);
	strcat(directoryname, "popsize_");
	strcat(directoryname, popsize);
    strcat(directoryname, "chromosomes_");
	strcat(directoryname, numberofchromosomes);
    strcat(directoryname, "L_");
	strcat(directoryname, chromosomesize);
    strcat(directoryname, "totalt_");
	strcat(directoryname, timeSteps);
	strcat(directoryname, "seed_");
	strcat(directoryname, randomnumberseed);
	
	return directoryname;
}

char * MakeFinalDataFileName(char * typeofrun, char * benmut, char * slopeforcontourline, char * randomnumberseed) 
{
	
    char * finaldatafilename = (char *) malloc(60);
    strcpy(finaldatafilename, "finaldatafor_");
    strcat(finaldatafilename, "runtype_");
    strcat(finaldatafilename, typeofrun);
    strcat(finaldatafilename, "mub_");
    strcat(finaldatafilename, benmut);
    strcat(finaldatafilename, "slope_");
    strcat(finaldatafilename, slopeforcontourline);
    strcat(finaldatafilename, "seed_");
    strcat(finaldatafilename, randomnumberseed);
    
    return finaldatafilename;
}
