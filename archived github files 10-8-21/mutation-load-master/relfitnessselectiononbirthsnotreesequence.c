/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<stdio.h>
#include<float.h>
#include<string.h>
#include "pcg_basic.h"
#include<math.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<unistd.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_fit.h>
#include<err.h>

#define VERBOSE 0
#define VERYVERBOSE 1
#define MISCELLANEOUS 1
#define INDIVIDUALWIDATA 1
#define OVERLAPPINGGENERATIONS 1 //probably should be an input argument at some point.
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654

//bug checking files
FILE	*verbosefilepointer;
FILE	*veryverbosefilepointer;
FILE	*miscfilepointer;
//simulation data files
FILE	*rawdatafilepointer;
FILE    *summarydatafilepointer;
FILE    *finaldatafilepointer;

void UpdateLast200NTimeSteps(double *, double);
void DoubleSwap(long double *, long double *);
void DoubleBubbleSort(long double *, int);
double CalculateVarianceInLogFitness(int, long double *, long double);
long double FindFittestWi(long double *, int);
double CalculateSlopeOfLogFitness(int, int, double *);
long double Fen_sum(long double *, int);
void Fen_add(long double *, int, long double, int);
long double Fen_range(long double *, int, int);
long double Fen_get(long double *, int);
void Fen_set(long double *, int, long double, int);
void InitializePopulation(long double *, long double *, int , double *, int, int, long double *);
//int ChooseVictimWithoutTree(long double *, int , long double);
int SearchTree(int, int, long double, long double *);
int ChooseParentWithTree(long double *, int, long double);
int ChooseVictim(int);
void RecombineChromosomesIntoGamete(int, int, int, int , int , double *, double *, int, int *);
int SampleFromPoisson(float );
int DetermineNumberOfMutations(int, int, float);
void MutateGamete(int , int , double *, double);
double CalculateWi(double *, double *, int);
void ReplaceVictim(double *, double *, int, int, long double *, double *, int, long double *, long double *);
//void RecalculateSumOfDeathRates(long double *, int , int , int , struct individual *);
void ProduceMutatedRecombinedGamete(int, int, double *, int , int , int , int , double , double , double , char * , double *, gsl_rng *);
void PerformOneTimeStep(int, int, int, long double *, long double *, double *, long double *, int, int, int, double, double, double, char *, double *, double *, gsl_rng *);
double RunSimulation(char *, char *, char *, char *, char *, char *, char *, char *, int, int, int, int, double, double, double, char *, gsl_rng *, FILE *, FILE *, FILE *);
int BracketZeroForSb(double *, double *, char *, char *, char *, char *, char *, char *, char *, int , int, int, int, double, double, double, char *, gsl_rng *, FILE *, FILE *, FILE *);
double BisectionMethodToFindSbWithZeroSlope(double *, double *, char *, char *, char *, char *, char *, char *, char *, int , int, int, int, double, double, double, char *, gsl_rng *, FILE *, FILE *, FILE *);

void main(int argc, char *argv[]) {

    
    char * directoryname = (char *) malloc(200);
    strcpy(directoryname, "datafor");
    strcat(directoryname, argv[8]);
    strcat(directoryname, "mub");
    strcat(directoryname, argv[6]);
    strcat(directoryname, "chromosomes");
    strcat(directoryname, argv[5]);
    strcat(directoryname, "popsize");
    strcat(directoryname, argv[2]);
    strcat(directoryname, "mud");
    strcat(directoryname, argv[3]);
    strcat(directoryname, "seed");
    strcat(directoryname, argv[11]);
    int check1, check2;
    check1 = mkdir(directoryname, 0777);
    check2 = chdir(directoryname);

    int Nxtimesteps;
    Nxtimesteps = atoi(argv[1]);

    int popsize;
    popsize = atoi(argv[2]);

    double deleteriousmutationrate;
    deleteriousmutationrate = atof(argv[3]); //remember that this is the per-locus deleterious mutation rate, not the genome-wide mutation rate.
    
    int chromosomesize;
    chromosomesize = atoi(argv[4]);

    int numberofchromosomes;
    numberofchromosomes = atoi(argv[5]); //remember that this is total number of chromosomes, not ploidy -- all individuals will be diploid.

    double beneficialmutationrate;
    beneficialmutationrate = atof(argv[6]); //remember that this is the per-locus rate, not genome-wide.
    
    //I have two parameters for Sb for the type of run that needs to have bracketed values of Sb.
    //In the case with just a single simulation being run, Sb2 here will be the value of Sb used.
    
    double Sb2;
    Sb2 = atof(argv[7]);
    double *pSb2 = &Sb2;
    
    double Sb1;
    Sb1 = 0.0;
    double *pSb1 = &Sb1;
    
    char * beneficialdistribution = (char *) malloc(30);
    strcpy(beneficialdistribution, argv[8]);
    
    char * typeofrun = (char *) malloc(30);
    strcpy(typeofrun, argv[9]);
    
    double slopeforcontourline;
    slopeforcontourline = atof(argv[10]);
    
    int randomnumberseed;
    randomnumberseed = atoi(argv[11]);
    
    pcg32_srandom(randomnumberseed, randomnumberseed); // seeds the random number generator.
    gsl_rng * randomnumbergeneratorforgamma = gsl_rng_alloc(gsl_rng_mt19937);
    //the gamma distribution function requires a gsl random number generator, which is set here.
    //it's a bit inelegant to have two different RNGs, which could be solved by using a different algorithm 
    //for choosing variates from a gamma distribution, instead of using the free one from gsl.

    verbosefilepointer = fopen("verbose.txt", "w");	//opens the file to which to print verbose data.
    veryverbosefilepointer = fopen("veryverbose.txt", "w"); //opens the file to which to print very verbose data.
    miscfilepointer = fopen("miscellaneous.txt", "w"); //opens the file to which to print miscellaneous data.
    
    char * finaldatafilename = (char *) malloc(60);
    strcpy(finaldatafilename, "finaldatafor");
    strcat(finaldatafilename, "runtype");
    strcat(finaldatafilename, argv[9]);
    strcat(finaldatafilename, "mub");
    strcat(finaldatafilename, argv[6]);
    strcat(finaldatafilename, "slope");
    strcat(finaldatafilename, argv[10]);
    strcat(finaldatafilename, "seed");
    strcat(finaldatafilename, argv[11]);
    finaldatafilepointer = fopen(finaldatafilename, "w");
   
    if (strcmp(typeofrun, "root") == 0) {
        
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
        BracketZeroForSb(pSb1, pSb2, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
        fprintf(miscfilepointer, "Finished bracketing function.");
        fflush(miscfilepointer);
        sbrequiredforzeroslopeoffitness = BisectionMethodToFindSbWithZeroSlope(pSb1, pSb2, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
        fprintf(finaldatafilepointer, "The value of Sb for which the slope of log fitness is zero with mub of %.10f is %.10f", beneficialmutationrate, sbrequiredforzeroslopeoffitness);
    
    } else if (strcmp(typeofrun, "single") == 0) {
        
        //This type of run just simulates a single population with the input parameters.
        RunSimulation(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
        
    } else {
        
        //One day maybe I'll have more types of runs.
        fprintf(miscfilepointer, "That type of run is not currently supported.");
    }
    
    free(directoryname);
    free(finaldatafilename);
    free(beneficialdistribution);
    free(typeofrun);
    
    fclose(verbosefilepointer);
    fclose(veryverbosefilepointer);
    fclose(miscfilepointer); //closes data files
    //free(fitnessdistributiondatafilename);    
    gsl_rng_free(randomnumbergeneratorforgamma);
    fclose(finaldatafilepointer);

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

void InitializePopulation(long double *wholepopulationwistree, long double *wholepopulationwisarray, int populationsize, double *populationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis) {
    int i, j;
    
    double haploidgenomelength = (double) ((totalpopulationgenomelength / populationsize) / 2);
    
    for (i = 0; i < populationsize; i++) {
        wholepopulationwistree[i] = 1.0; //for relative fitness, all individuals start with probability of being chosen as a parent of 1/N
        wholepopulationwisarray[i] = 1.0;
    }
    //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
    for (i = 0; i < populationsize; i++) {
	j = i + LSB(i+1);
        if (j < populationsize) {
            wholepopulationwistree[j] += wholepopulationwistree[i];
        }
    }
    for (i = 0; i < totalpopulationgenomelength; i++) {
        populationgenomes[i] = 0.0;
    }
    *psumofwis = (long double) populationsize;
    
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

//The tree in the name is the Fenwick tree, which stores the fitnesses of individuals in the population.
//This function is where selection occurs -- individuals with higher-than-average fitness will be chosen more often as parents.
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis)
{
    long double randomnumberofbirth;
    int newparent = 0;
    
    randomnumberofbirth = (ldexp(pcg32_random(), -32)) * sumofwis;
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize;
    if (leftbound > rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    
    newparent = (SearchTree(leftbound, rightbound, randomnumberofbirth, wholepopulationwistree));
    return newparent;
}

//In this model, individuals die at random. There's no selection happening here.
int ChooseVictim(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);
    return randomindividual;
}

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(int totaltimesteps, int currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, double *populationgenomes, int totalindividualgenomelength, int * recombinationsites)
{
    int recombinationsite, startchromosome, startofindividual, h, i, returnvaluefortskit;
    startofindividual = persontorecombine * totalindividualgenomelength;
    
    
    for (h = 0; h < numberofchromosomes; h++) {
	startchromosome = pcg32_boundedrand(2); //generates either a zero or a one to decide to start with chromosome 1 or 2.
	
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); //it doesn't make sense to do a recombination event before the first linkage block. Note that this will never break if the chromosome size is only one linkage block.
        
        
	for (i = 0; i < recombinationsite; i++) {
            if (startchromosome == 0) {
		gamete[h*chromosomesize + i] = populationgenomes[startofindividual + (h*chromosomesize) + i];
            }
            else {
		gamete[h*chromosomesize + i] = populationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
            }
	}
        for (i = recombinationsite; i < chromosomesize; i++) {
            if (startchromosome == 0) {
                gamete[h*chromosomesize + i] = populationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
            }
            else {
                gamete[h*chromosomesize + i] = populationgenomes[startofindividual + (h*chromosomesize) + i];
            }
	}
    }
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
    } else { 				//for larger means, use rejection method.
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

int DetermineNumberOfMutations(int chromosomesize, int numberofchromosomes, float mutationrate)
{
    float meannumberofmutations = mutationrate * (float) chromosomesize * (float) numberofchromosomes;
    
    //Note that because this function operates on gametes, the calculation above appears haploid.
    //There shouldn't be a multiplication by 2 (for diploidy) in this function, since it will be called twice per individual: once per gamete.
    //Above calculation should be moved outside this function for increased efficiency.
    
    int numberofmutations = SampleFromPoisson(meannumberofmutations);
    return numberofmutations;
}

void MutateGamete(int chromosomesize, int numberofchromosomes, double *gamete, double mutationeffectsize)
{
    
    int randomchromosometomutate = pcg32_boundedrand(numberofchromosomes); //if we decide to include heterogenous rates of recombination/mutation, both of these will need to be replaced by a function that weights each linkage block's probability of mutating.
    int randomblocktomutate = pcg32_boundedrand(chromosomesize);
    int mutatedsite = randomchromosometomutate*chromosomesize + randomblocktomutate;
    gamete[mutatedsite] += log(1 + mutationeffectsize);
    char derivedstate[10];
    sprintf(derivedstate, "%.9f", mutationeffectsize);
    
}

//this function seems inefficient, but with recombination and mutation, I'm not sure there's a significantly easier way.
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

void ReplaceVictim(double *parent1gamete, double *parent2gamete, int currentpopsize, int currentvictim, long double *sumofwis, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationwistree, long double *wholepopulationwisarray)
{
    int i;
    double newwi;
	
    newwi = CalculateWi(parent1gamete, parent2gamete, totalindividualgenomelength);
	
    for (i = 0; i < (totalindividualgenomelength/2); i++) {
        wholepopulationgenomes[currentvictim*totalindividualgenomelength + i] = parent1gamete[i];
        wholepopulationgenomes[currentvictim*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
//It would be more efficient to build directly into victim slot, but right now it is possible for the victim to also be a parent, later add an if statement to more efficiently deal with more common case.
    }

    *sumofwis -= wholepopulationwisarray[currentvictim];
    Fen_set(wholepopulationwistree, currentpopsize, newwi, currentvictim);
    
    wholepopulationwisarray[currentvictim] = (long double) newwi;
    *sumofwis += (long double) newwi;

}

void ProduceMutatedRecombinedGamete(int totaltimesteps, int currenttimestep, double *wholepopulationgenomes, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, int currentparent, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char * beneficialdistribution, double *parentgamete, gsl_rng * randomnumbergeneratorforgamma)
{
    
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    double Sds[30];
    int recombinationsites[numberofchromosomes];
    
    //Following lines produce a gamete from parent 1 and add deleterious and beneficial mutations to the gamete.
    RecombineChromosomesIntoGamete(totaltimesteps, currenttimestep, currentparent, chromosomesize, numberofchromosomes, parentgamete, wholepopulationgenomes, totalindividualgenomelength, recombinationsites);
    
    //Following lines stochastically generate a number of deleterious mutations drawn from a Poisson distribution with mean determined by the deleterious mutation rate
    //with effect sizes drawn from a gamma distribution with parameters taken from Kim et al 2017.
    int DontBreakWhileLoop = 0;
    while (1) {
        DontBreakWhileLoop = 0;
        numberofdeleteriousmutations = DetermineNumberOfMutations(chromosomesize, numberofchromosomes, deleteriousmutationrate);
        for (k = 0; k < numberofdeleteriousmutations; k++) {
            Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646); //Uses parameters for the gamma distribution of the selection coefficients of new mutations scaled to an inferred ancestral populations size. To produce the distribution of unscaled effect sizes, numbers drawn from this distribution must be divided by two times the ancestral population size for the population from which the distribution was derived (11,823 in this case). Data used to produce these fits were samples from 6503 individuals from the National Heart, Lung, and Blood Institute European-American dataset. Analysis of DFE from Kim et al. 2017.
            //This gamma distribution can occasionally produce deleterious mutations with effect sizes larger than 1,
            //which would result in a gamete with fitness less than zero, which would break my algorithm.
            //The following if statement simply throws out any deleterious mutations with effect sizes larger than 1.
            if (Sds[k] >= 1) {
                DontBreakWhileLoop = 1;
                break;
            }     
        }
        if (DontBreakWhileLoop == 0) 
            break;
    }
    
    //Following lines stochastically generate a number of beneficial mutations drawn from a Poisson distribution with mean determined by the beneficial mutation rate.
    numberofbeneficialmutations = DetermineNumberOfMutations(chromosomesize, numberofchromosomes, beneficialmutationrate);
    
    //Adds the specified number of deleterious mutations to the gamete, recording the sites of each mutation for tree sequence recording.
    for (k = 0; k < numberofdeleteriousmutations; k++) {
        MutateGamete(chromosomesize, numberofchromosomes, parentgamete, -Sds[k]);
    }
    
    //Adds the specified number of beneficial mutations, drawing Sb values from the specified distribution.
    //Sites of each mutation are added to the mutationsites array for tree sequence recording.
    if (strcmp(beneficialdistribution, "point") == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, Sb);
        }
    } else if (strcmp(beneficialdistribution, "exponential") == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, generatedSb);
        }
    } else if (strcmp (beneficialdistribution, "uniform") == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, generatedSb);
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        for (k = 0; k < numberofbeneficialmutations; k++) {
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, Sb);
        }
    }
    
    
}

void PerformOneTimeStep(int popsize, int totaltimesteps, int currenttimestep, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char * beneficialdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma)
{

    int currentparent1, currentparent2, currentvictim;

    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis);
    while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis);
    }
    
   
    ProduceMutatedRecombinedGamete(totaltimesteps, currenttimestep, wholepopulationgenomes, chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent1, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, randomnumbergeneratorforgamma);
    ProduceMutatedRecombinedGamete(totaltimesteps, currenttimestep, wholepopulationgenomes, chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent2, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent2gamete, randomnumbergeneratorforgamma);
                
    ReplaceVictim(parent1gamete, parent2gamete, popsize, currentvictim, psumofwis, wholepopulationgenomes, totalindividualgenomelength, wholepopulationwistree, wholepopulationwisarray);
    
}

double RunSimulation(char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char * beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE * veryverbosefilepointer, FILE * verbosefilepointer, FILE * miscfilepointer)
{
    int i, j, k;
    
    char * rawdatafilename = (char *) malloc(200);
    strcpy(rawdatafilename, "rawdatafor"); //starting the string that will be the name of the data file.

    strcat(rawdatafilename, "Nxtimesteps"); //for adding values of generations to the data name.
    strcat(rawdatafilename, Nxtimestepsname);

    strcat(rawdatafilename, "popsize"); //for adding values of starting population sizes to the data name.
    strcat(rawdatafilename, popsizename);

    strcat(rawdatafilename, "mutrate"); //for adding values of mutation rate to the data name (remember that mutation rate is currently the per-locus rate, not per-genome).
    strcat(rawdatafilename, delmutratename);

    strcat(rawdatafilename, "chromsize"); //for adding values of chromosome size to the data name.
    strcat(rawdatafilename, chromsizename);

    strcat(rawdatafilename, "chromnum"); //for adding values of the number of chromosomes to the data name.
    strcat(rawdatafilename, chromnumname);
    
    strcat(rawdatafilename, "benmutrate"); //for adding values of the beneficial mutation rate to the data name.
    strcat(rawdatafilename, mubname);
    
    strcat(rawdatafilename, "Sb"); //for adding values of the beneficial mutation effect size to the data name.
    strcat(rawdatafilename, Sbname);
    
    strcat(rawdatafilename, ".txt");

    rawdatafilepointer = fopen(rawdatafilename, "w"); //opens the file to which to print data to be plotted.
    fprintf(rawdatafilepointer, "Nxtimesteps,Sum.of.wis,Variance.in.log.fitness\n");
    
    char * summarydatafilename = (char *) malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
    
    
    int totaltimesteps = Nxtimesteps * popsize;
    int currenttimestep = 0;
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = popsize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
    long double sumofwis;
    long double *psumofwis = &sumofwis;
    long double *wholepopulationwistree;
    wholepopulationwistree = malloc(sizeof(long double) * popsize);
       
    long double *wholepopulationwisarray;
    wholepopulationwisarray = malloc(sizeof(long double) * popsize);
    //The Fenwick tree does not store each individual's wi, but rather a collection of partial sums.
    //For debugging purposes and data that requires summarizing wis, storing the wis in an array is necessary.

    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * popsize);
    //In order to visualize the distribution of fitness in the population,
    //the individual fitnesses need to be sorted, which requires a separate array.

    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Entered simulation run.\n");
        fflush(veryverbosefilepointer);
    }
    
    InitializePopulation(wholepopulationwistree, wholepopulationwisarray, popsize, wholepopulationgenomes, totalpopulationgenomelength, totaltimesteps, psumofwis);
    /*Sets the initial population to have zeroes in all their linkage blocks,
    death rates equal to the baseline wi, and an identifier number.
    It also sums the wis and returns the sum.*/
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Population initialized.\n");
        fflush(veryverbosefilepointer);
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * Nxtimesteps);
    //In order to calculate the slope of degradation of fitness,
    //I need to store the average fitness each generation.
    
    long double currentfittestindividualswi;
    double parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
    
    //Following array is to store the variance in log(fitness) for twenty generations at a time for estimating the end of the burn-in phase.
    //The choice of twenty is completely arbitrary and could eventually be an input if it seems worth it.
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
    arbitrarynumber = (-1 * 0.007 / popsize); //using a number somewhere close to the mean of the DFE for deleterious mutations.
    double slopeoflogfitness;    
    double variancesum;
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
        fflush(veryverbosefilepointer);
    }
    
    //BEGIN THE SIMULATION FOR LOOP
    
    for (i = 0; i < Nxtimesteps; i++) {
        
        //Following code performs N rounds of paired births and deaths.
        for (j = 0; j < popsize; j++) {
            currenttimestep += 1;            
            PerformOneTimeStep(popsize, totaltimesteps, currenttimestep, wholepopulationwistree, wholepopulationwisarray, wholepopulationgenomes, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma);
            
        }
        
        //Following code calculates the variance in log(fitness) of the population after this generation of births and deaths.
        //May use an imprecise algorithm -- check before using as data.
        variancesum = CalculateVarianceInLogFitness(popsize, wholepopulationwisarray, *psumofwis);
        
        
        double c0, cov00, cov01, cov11, sumsq;
        if (isburninphaseover == 0) {
            UpdateLast200NTimeSteps(last200Ntimestepsvariance, variancesum);
            UpdateLast200NTimeSteps(literallyjustlast200Ntimesteps, i+1);
            if (i > 199) {           //to avoid calling the end of the burn-in phase at generation one
                                    //because of setting pre-simulation generations to zeroes
                                    //I just won't start looking for the end of the burn-in phase until 200 generations
                                    //This would be a mild problem if a simulation should end in 200 generations, but that shouldn't ever happen with the DFE I'm using.
                slopeofvariance = 0.0;
                gsl_fit_linear(literallyjustlast200Ntimesteps, step, last200Ntimestepsvariance, step, 200, &c0, &slopeofvariance, &cov00, &cov01, &cov11, &sumsq);
                if (slopeofvariance < arbitrarynumber) {
                    endofburninphase = i;
                    endofdelay = endofburninphase + 500;
                    isburninphaseover = 1;
                    fprintf(miscfilepointer, "Burn-in phase called as ending in generation %d\n", i+1);
                    fprintf(summarydatafilepointer, "Burn-in phase called as ending in generation %d\n", i+1);
                    if (VERBOSE == 1) {
                        fflush(miscfilepointer);
                        fflush(summarydatafilepointer);
                    }
                    
                }
            }
        }        
        
        //This is the main data output, currently the summed fitness and variance in log(fitness) in the population.
        fprintf(rawdatafilepointer, "%d,%Lf,%.18f\n", i+1, *psumofwis, variancesum);
        
        
        
        //This is to produce a histogram of the wis of the entire population from a single generation.
        //It's terrible and completely non-modular, but I just can't bring myself to add in two more user-input arguments.
        if (strcmp(typeofrun, "single") == 0) {
            if (i == 1999) {
                
                if (INDIVIDUALWIDATA == 1) {
                    if (VERYVERBOSE == 1) {
                        fprintf(veryverbosefilepointer, "Just before individual wi data lines.\n");
                        fflush(veryverbosefilepointer);
                    }                    
                    fprintf(summarydatafilepointer, "Individual, Wi\n");
                    for (k = 0; k < popsize; k++) {
                        fprintf(summarydatafilepointer, "%d,%Lf\n", k+1, wholepopulationwisarray[k]);
                    }
                }
            }
        }
        
        //If the burn-in phase has been called, wait 500 generations to start recording fitnesses.
        //This is to be sure that even when the beneficial rates/sizes are large, the only generations recorded will be from the uniformly sloping part of the simulation.
        //The average fitness from any generation after this delay period is recorded in the array of average fitnesses.
        if (i > endofdelay) {
            logaveragefitnesseachNtimesteps[Nxtimestepsafterburnin] = log((double) *psumofwis / (double) popsize);
            if (VERYVERBOSE == 1) {
                fprintf(veryverbosefilepointer, "log average fitness in generation %d, %d generations after burn-in, is: %f\n", i, Nxtimestepsafterburnin, logaveragefitnesseachNtimesteps[Nxtimestepsafterburnin]);
                fflush(veryverbosefilepointer);
            }
            Nxtimestepsafterburnin += 1;
        }
        
        //These lines ensure that the magnitude of fitness hasn't declined by too much.
        //At extremely small fitness values, floating-point math becomes imprecise.
        //These lines end the simulation if fitness declines below 10^-10, which should represent a completely degraded population.
        currentfittestindividualswi = FindFittestWi(wholepopulationwisarray, popsize);
        if (currentfittestindividualswi < pow(10.0, -10.0)) {
            fprintf(miscfilepointer, "\nFitness declined to less than 10^-10 during generation %d.", i+1);
            fprintf(summarydatafilepointer, "Fitness declined to catastrophic levels in generation %d.\n", i+1);
            endofsimulation = i;
            i = Nxtimesteps;
            didpopulationcrash = 1;
        }
    }
    
    //END OF SIMULATION FOR LOOP
   
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %.6f. Final population summed fitness was: %Lf\n", Sb, *psumofwis);
        fflush(veryverbosefilepointer);
    }
    if (didpopulationcrash == 0) {
        endofsimulation = i;
    }
    if (isburninphaseover == 1) {
        if (VERYVERBOSE == 1) {
            fprintf(veryverbosefilepointer, "Calculating slope of log fitness with the following parameters: endofsimulation = %d, endofdelay = %d, generationsafterburnin = %d\nLog fitness each generation: ", endofsimulation, endofdelay, Nxtimestepsafterburnin);
            for (j = 0; j < (endofsimulation - endofdelay); j++) {
                fprintf(veryverbosefilepointer, "%f ", logaveragefitnesseachNtimesteps[j]);
            }
            fprintf(veryverbosefilepointer, "\n");
            fflush(veryverbosefilepointer);
        }
        slopeoflogfitness = CalculateSlopeOfLogFitness(endofsimulation, endofdelay, logaveragefitnesseachNtimesteps);
        fprintf(summarydatafilepointer, "Slope of log(fitness) after the burn-in phase: %f\n", slopeoflogfitness);
        
        if (VERBOSE == 1) {
            fflush(summarydatafilepointer);
            fflush(rawdatafilepointer);
        }
        
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(sortedwisarray);
        
        return slopeoflogfitness;
    }
    if (isburninphaseover == 0) {
        fprintf(summarydatafilepointer, "End of burn-in phase not reached.");
        
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(sortedwisarray);
        
        return -1.0;
    }
}

//The following function is heavily modified from Numerical Recipes in C, Second Edition.
//For large population sizes, populations with mean Sb > 0 may actually have a more negative fitness slope than mean Sb = 0.
//
int BracketZeroForSb(double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, char * beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE * veryverbosefilepointer, FILE * verbosefilepointer, FILE * miscfilepointer) {
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
    resultingslope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
    resultingslope2 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
            resultingslope2 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
            resultingslope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
double BisectionMethodToFindSbWithZeroSlope(double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, char * beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE * veryverbosefilepointer, FILE * verbosefilepointer, FILE * miscfilepointer) {
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
    slope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }
    
    slopemid = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
        slopemid = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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