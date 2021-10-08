/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 //current bug is that the code is reading the data backwards.
 //tree goes with inverse sum of wis, Fen_set
 //data output is normal sum of wis, log fitness, calcuate variance, slop of fitness

#include<stdio.h>
#include<stdlib.h>
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

#define QUICK_DEATH 0
#define PERFORMBIRTHMARKER 0
#define PERFORMONETIMESTEPMARKERS 0
#define PERFORMDEATHMARKER 0
#define VERBOSE 0
#define VERYVERBOSE 0
#define RUNSIMULATIONMARKERS 1
#define MISCELLANEOUS 1
#define INDIVIDUALWIDATA 1
#define OVERLAPPINGGENERATIONS 1 //probably should be an input argument at some point.
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654
#define REDEFINECHROM 0

//bug checking files
FILE* verbosefilepointer;
FILE* veryverbosefilepointer;
FILE* miscfilepointer;
//simulation data files
FILE* rawdatafilepointer;
FILE* summarydatafilepointer;
FILE* finaldatafilepointer;

typedef struct individual {

    double deathRate;
    double wis;
    int index;
    double* genome;

}individual;

void arrayPositionFlip(int, int*, int);
double averageDeathRate(int*, struct individual*, int);
double calcAvgPopSize(double*);
void storePopsizeForGenerations(double*, int, int);
void UpdateLast200NTimeSteps(double*, double);
void UpdateLast200GenStepsAbs(int*, int);
void DoubleSwap(long double*, long double*);
void DoubleBubbleSort(long double*, int);
double CalculateVarianceInLogFitness(int, long double*, long double, struct individual*, int*);
long double FindFittestWi(int popsize, struct individual* popArray, int* indexArray);
double CalculateSlopeOfLogFitness(int, int, int, double*);
long double Fen_sum(long double*, int);
void Fen_add(long double*, int, long double, int);
long double Fen_range(long double*, int, int);
long double Fen_get(long double*, int);
void Fen_set(long double*, int, long double, int);
void InitializePopulation(long double*, int, int, long double*, long double*, struct individual*, int*, int*, int);
int SearchTree(int, int, long double, long double*, int);
int ChooseParent(int);
int ChooseVictimWithTree(long double*, int, long double, long double, int);
void RecombineChromosomesIntoGamete(int persontorecombine, int chromosomesize, int numberofchromosomes, double* gamete, int totalindividualgenomelength, int* recombinationsites, struct individual* popArray, int* indexArray, double currentTimeStep);
int SampleFromPoisson(float);
int DetermineNumberOfMutations(int, int, float);
void MutateGamete(int, int, double*, double);
double CalculateWi(double*, double*, int);
void ProduceMutatedRecombinedGamete(int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, int currentparent, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char* beneficialdistribution, double* parentgamete, gsl_rng* randomnumbergeneratorforgamma, struct individual* popArray, int* arrayOfIndexes, double currentTimeStep);
void PerformOneTimeStep(int* pPopSize, int currenttimestep, long double* wholepopulationwistree, double* wholepopulationgenomes, long double* psumofwis, long double* pInverseSumOfWis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char* beneficialdistribution, double* parent1gamete, double* parent2gamete, gsl_rng* randomnumbergeneratorforgamma, int birthBool, struct individual* popArray, int* freeIndexes, int* arrayOfIndexes, int eventNumber, int maxPopSize, double currentTimeStep, double averageDeathRate, long double* lastDeathRateGone);
double RunSimulation(char*, char*, char*, char*, char*, char*, char*, char*, int, int, int, int, int, double, double, double, char*, gsl_rng*, FILE*, FILE*, FILE*);
int BracketZeroForSb(double*, double*, char*, char*, char*, char*, char*, char*, char*, int, int, int, int, int, double, double, double, char*, gsl_rng*, FILE*, FILE*, FILE*);
double BisectionMethodToFindSbWithZeroSlope(double*, double*, char*, char*, char*, char*, char*, char*, char*, int, int, int, int, int, double, double, double, char*, gsl_rng*, FILE*, FILE*, FILE*);
double ExponentialDerivateOfUnitMeanOne(float idum);
double rateOfDeathsCalc(double, double);
double rateOfBirthsCalc(int, double, int, double);
int discoverEvent(double, double, double);
int monteCarloStep(int, double, double*, double, struct individual*, int, double);
double averageWiOfStructs(int*, int*, int);
void performDeath(int* pCurrentPopsize, int currentvictim, long double* sumofwis, long double* wholepopulationwistree, long double* pInverseSumOfWis, struct individual* popArray, int* arrayOfIndexes, int* arrayOfFreeIndexes, int MAX_POP_SIZE, double currentTimeStep, int eventNumber, int totalindividualgenomelength, long double* lastDeathRateGone);
void performBirth(double* parent1gamete, double* parent2gamete, int* pCurrentPopsize, int currentvictim, long double* pSumOfWis, int totalindividualgenomelength, long double* wholepopulationwistree, long double* pInverseSumOfWis, struct individual* popArray, int* freeIndexes, int* arrayOfIndexes, const int MAX_POP_SIZE, int eventNumber, double currentTimeStep, double populationAverage);
double sumOfWiOfStructs(int*, struct individual*, int);
void allocateMemoryForSizeOfGenome(int, int, struct individual**, struct individual*);

int main(int argc, char* argv[]) {

    if (argc != 13) {
        printf("[Error]; Wrong number of arguments in program. It should be timeSteps, initialPopsize, deleteriousmutationrate, chromosomesize, numberofchromosomes, beneficialmutationrate, Sb, beneficialdistribution, typeofrun, slope, seed, maxPopSize \n");
        return -1;
    }
    //develops directory name from arguments given in the scheduling script. Michael Osipov (8/14/2021)
    char* directoryname = (char*)malloc(200);
    strcpy(directoryname, "datafor");
    strcat(directoryname, argv[8]);
    strcat(directoryname, "mub");
    strcat(directoryname, argv[6]);
    strcat(directoryname, "chromosomes");
    strcat(directoryname, argv[5]);
    strcat(directoryname, "popsize");
    strcat(directoryname, argv[2]);
    strcat(directoryname, argv[3]);
    strcat(directoryname, "seed");
    strcat(directoryname, argv[11]);
    int check1, check2;
    check1 = mkdir(directoryname, 0777);
    check2 = chdir(directoryname);

    //Maximum number of timesteps to run the simulation. Michael Osipov (8/14/2021)
    int timeSteps;
    timeSteps = atoi(argv[1]);

    //Initial population size the simulation starts at. Michael Osipov (8/14/2021)
    int initialPopsize;
    initialPopsize = atoi(argv[2]);

    //Maximum population size or N in calculations. Michael Osipov (8/14/2021)
    int maxPopSize;
    maxPopSize = atoi(argv[12]);


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
    double* pSb2 = &Sb2;

    double Sb1;
    Sb1 = 0.0;
    double* pSb1 = &Sb1;

    char* beneficialdistribution = (char*)malloc(30);
    strcpy(beneficialdistribution, argv[8]);

    char* typeofrun = (char*)malloc(30);
    strcpy(typeofrun, argv[9]);

    double slopeforcontourline;
    slopeforcontourline = atof(argv[10]);

    //useful to change the seed if you wish to run the same code several times with seperate saved files.
    int randomnumberseed;
    randomnumberseed = atoi(argv[11]);

    pcg32_srandom(randomnumberseed, randomnumberseed); // seeds the random number generator.//set to one (6/16/2020)
    //pcg32_srandom(800, 800);
    gsl_rng* randomnumbergeneratorforgamma = gsl_rng_alloc(gsl_rng_mt19937);
    //the gamma distribution function requires a gsl random number generator, which is set here.
    //it's a bit inelegant to have two different RNGs, which could be solved by using a different algorithm 
    //for choosing variates from a gamma distribution, instead of using the free one from gsl.

    verbosefilepointer = fopen("verbose.txt", "w");	//opens the file to which to print verbose data.
    veryverbosefilepointer = fopen("veryverbose.txt", "w"); //opens the file to which to print very verbose data.
    miscfilepointer = fopen("miscellaneous.txt", "w"); //opens the file to which to print miscellaneous data.

    char* finaldatafilename = (char*)malloc(60);
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

    //root finding for values of Ub and Sb
    //this is not possible for a abs fitness run so its being ignored
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
        BracketZeroForSb(pSb1, pSb2, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], typeofrun, timeSteps, initialPopsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
        fprintf(miscfilepointer, "Finished bracketing function.");
        fflush(miscfilepointer);
        sbrequiredforzeroslopeoffitness = BisectionMethodToFindSbWithZeroSlope(pSb1, pSb2, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], typeofrun, timeSteps, initialPopsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
        fprintf(finaldatafilepointer, "The value of Sb for which the slope of log fitness is zero with mub of %.10f is %.10f", beneficialmutationrate, sbrequiredforzeroslopeoffitness);

    }
    else if (strcmp(typeofrun, "single") == 0) {

        //This type of run just simulates a single population with the input parameters.
        RunSimulation(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], typeofrun, timeSteps, initialPopsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);

    }
    else {

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

    return -1;
}

void arrayPositionFlip(int popSize, int* arrayOfIndexes, int currentVictim){

    /*
    * This function flips the integer value in two points withing the arrayOfIndexes 
    */

	int storageIntDead;
	int storageIntAlive;

	storageIntDead = arrayOfIndexes[currentVictim];
	storageIntAlive = arrayOfIndexes[popSize];

	arrayOfIndexes[popSize] = storageIntDead;
	arrayOfIndexes[currentVictim] = storageIntAlive;

}

double averageDeathRate(int* indexArray, struct individual* popArray, int popSize){

    /*
    * This function finds the arithmetic mean of the deathrates in the whole living population.
    */

    int i = 0;
    int currentIndex;

    double average = 0;
    double combinedDeathRate = 0;

    struct individual individualFrame;

    while (i < popSize) {

        currentIndex = indexArray[i];//goes through the sorted array of indexes to ensure there are no holes
        individualFrame = popArray[currentIndex];//currentindex values will always be locations in the poparray with living individuals
        combinedDeathRate = individualFrame.deathRate + combinedDeathRate;

        i++;

    }

    average = combinedDeathRate / popSize;

    return average;
}

//function calculates average popsize over a generatiom, generation is a inditerminate number used for testing
double calcAvgPopSize(double* popSizeArrayForAverage)
{
	int i;
	double averagePopSize;
	const double MAX_POP_SIZE = 10000;//remeber to have all functions inherit this variable this is just easier

	for(i = 0; i < MAX_POP_SIZE; i++){
		averagePopSize = averagePopSize + popSizeArrayForAverage[i];
	}

	averagePopSize = averagePopSize/MAX_POP_SIZE;//arithmetic mean

	return averagePopSize;

}

//function stores average population over a generation
void storePopsizeForGenerations(double* popSizeArrayForAverage, int popSize, int maxPopSize)
{
	int const MAX_POP_SIZE_MINUS_ONE = maxPopSize - 1;
	int i;
	double storage[maxPopSize];
	for(i = 0; i < MAX_POP_SIZE_MINUS_ONE; i++){
		storage[i] = popSizeArrayForAverage[i + 1];
	}
	storage[maxPopSize] = popSize;
	for(i = 0; i < maxPopSize; i++){
		popSizeArrayForAverage[i] = storage[i];
	}
}

void UpdateLast200NTimeSteps(double* last200NTimeSteps, double newNTimeSteps)
{
    double storage[200];
    int m;
    for (m = 0; m < 199; m++) {
        storage[m] = last200NTimeSteps[m + 1];
    }
    storage[199] = newNTimeSteps;
    for (m = 0; m < 200; m++) {
    	last200NTimeSteps[m] = storage[m];
    }
}

void UpdateLast200GenStepsAbs(int* last200Gen, int newAveragePopSize)
{
    int storage[200];
    int i;

    for(i = 0; i < 199; i++){
    	storage[i] = last200Gen[i + 1];
    }
    storage[199] = newAveragePopSize;
    for (i = 0; i < 200; i++) {
    	last200Gen[i] = storage[i];
    }

}

//swaps the addresses of two doubles within an array for DoubleBubbleSort
void DoubleSwap(long double* x, long double* y)
{
    long double temp = *x;
    *x = *y;
    *y = temp;
}

//Inefficient algorithm (guaranteed to be order n^2). Improve algorithm if using more than once per simulation.
void DoubleBubbleSort(long double* arrayToBeSorted, int arraySize)
{
    int i, j;
    for (i = 0; i < arraySize - 1; i++) {
        for (j = 0; j < arraySize - i - 1; j++) {
            if (arrayToBeSorted[j] > arrayToBeSorted[j + 1]) {
                DoubleSwap(&arrayToBeSorted[j], &arrayToBeSorted[j + 1]);
            }
        }
    }
}

//used originally as a bug testing function. If you think that an individual may have a much higher fitness than the average fitness value.
long double FindFittestWi(int popsize, struct individual* popArray, int* indexArray)
{
    long double fittestwi;
    int i;
    fittestwi = popArray[indexArray[0]].wis;
    for (i = 1; i < popsize; i++) {
        if (popArray[indexArray[i]].wis > fittestwi) {
            fittestwi = popArray[indexArray[i]].wis;
        }
    }

    return fittestwi;
}

double CalculateSlopeOfLogFitness(int endofsimulation, int endofburninphase, int numberOfTimeSteps, double* logaveragefitnesseachgeneration)
{
    size_t step = 1;
    int k;
    const int MAX_POP_SIZE = 10000;
    double c0, cov00, cov01, cov11, sumsq;
    int generationsafterburnin;
    double slopeoflogfitness;
    double slopeOfPopulationSize;//basic fix to checking burn in phase for a absolute population (8/3/2020)
    //generationsafterburnin = (endofsimulation - endofburninphase);//this equation was used for a relative fitness model
    generationsafterburnin = numberOfTimeSteps/MAX_POP_SIZE;//generations after burn in number of time steps

    //I have to make an array of numbers to use as the x variable in a linear model to find the slope.
    double* justnumbers;
    justnumbers = malloc(sizeof(double) * generationsafterburnin);
    //put an arbitrary number for generations
    for (k = 0; k < generationsafterburnin; k++) {
        justnumbers[k] = (k + 1);
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
long double Fen_sum(long double* tree, int i)
{
    long double sum = 0;
    while (i) {
        sum += tree[i - 1];
        i -= LSB(i);
    }
    return sum;
}

//modified to add inverse of fitness
//used only for future code
//will make code make no sense
//Adds an amount to the ith element in the tree (and therefore to the Fen_sum for all elements in the tree greater than i).

void Fen_add(long double* tree, int numberofelementsintree, long double amounttoadd, int i)
{

    while (i < numberofelementsintree) {
        tree[i] += amounttoadd;//changed from adding to subtraction 11/25/2019
        i += LSB(i + 1);

        fprintf(verbosefilepointer, " ");
    }
}

//Returns the sum of the elements i through j-1.
//Could do with Fen_sum of j minus Fen_sum of i, but this is faster.
long double Fen_range(long double* tree, int i, int j)
{
    long double sum = 0;
    while (j > i) {
        sum += tree[j - 1];
        j -= LSB(j);
    }
    while (i > j) {
        sum -= tree[i - 1];
        i -= LSB(i);
    }
    return sum;
}

//Returns the value of the element at index i.
long double Fen_get(long double* tree, int i)
{
    return Fen_range(tree, i, i + 1);
}

void Fen_set(long double* tree, int numberofelementsintree, long double newvalue, int i)
{
    Fen_add(tree, numberofelementsintree, newvalue - Fen_get(tree, i), i);
}

//taken from numerical recipes in C
//will be more official later put error here for reminder
//creating this as structs force something new to be done
double sumOfWiOfStructs(int* indexArray, struct individual* popArray, int popSize) {

    int i = 0;
    int currentIndex;

    double sum = 0;

    struct individual individualFrame;

    while (i < popSize) {

        currentIndex = indexArray[i];
        individualFrame = popArray[currentIndex];//redo this
        sum = individualFrame.wis + sum;

        i++;

    }

    return sum;
}

//finds the average fitness of all individuals in the population array
double averageWiOfStructs(int* indexArray, int* pPopArray, int popSize) {

    int i = 0;
    int currentIndex;

    double average = 0;

    struct individual individualFrame;

    while (i < popSize) {

        currentIndex = indexArray[i];
        individualFrame.index = pPopArray[currentIndex];
        average = individualFrame.wis + average;

        i++;

    }

    average = average / popSize;

    return average;
}

//add mean
//population of 10000 initial
double ExponentialDerivateOfUnitMeanOne(float idum) {

    float ran1(long* idum);

    float dum;



    do
        dum = ldexp(pcg32_random(), -32);//check if this is how pcg32 is called
    while (dum == 0.0);

    return -log(dum);

}

double rateOfDeathsCalc(double sumOfFitness, double D0) {

    int iterator = 0;

    double netFitness;
    double deathRate;

    iterator++;

    deathRate = D0 + sumOfFitness;

    return deathRate;

}

double rateOfBirthsCalc(int populationSize, double b, int MAX_POP_SIZE, double currentTimeStep) {

    double birthRate;
    double popSizeConvertedToDouble = populationSize * 1.0;

    birthRate = (b) * (popSizeConvertedToDouble) * (1 - (popSizeConvertedToDouble / MAX_POP_SIZE));

    return  birthRate;

}

int discoverEvent(double deathRate, double birthRate, double currentTimeStep) {

    int boolBirth;

    double combinedBirthDeathRate;
    double cutOffPoint;
    double randomNumber;

    randomNumber = ldexp(pcg32_random(), -32);

    //finds what percent of events at this time would be a death event to determine the cuttoff point.
    //if lands on the cutoff point death occurs, techincally skewing towards death but is very minor.
    combinedBirthDeathRate = deathRate + birthRate;
    cutOffPoint = deathRate/combinedBirthDeathRate;

    if (randomNumber > cutOffPoint) {
        boolBirth = 1;
    }
    else if (randomNumber <= cutOffPoint) {
        boolBirth = 0;
    }
    else {

    }

    return boolBirth;

}

//returns output of either birth or death
int monteCarloStep(int popSize, double sumWi, double* pTimeElapsed, double sumOfDeathRates, struct individual* popArray, int MAX_POP_SIZE, double currentTimeStep) {

    int boolVar;
    int randSeed = 1;

    double deathRate;
    double birthRate;
    double time;
    double mean;
    double deathRateAtPerfection;
    double randomNumber;

    double *pRandomNumber = &randomNumber;

    deathRateAtPerfection = popSize;

    const double b = 10;

    deathRate = sumOfDeathRates;

    //rate of births is calculated using equation used in lab write up
    birthRate = rateOfBirthsCalc(popSize, b, MAX_POP_SIZE, currentTimeStep);
	

    mean = ((1.0)/(deathRate + birthRate));
    time = ExponentialDerivateOfUnitMeanOne(randSeed);//draws a random number from a distribution with unit mean 1. This occurs because a even t is most likely to occur right after a previous event.
    time = (time) * (mean);

    *pTimeElapsed = time + *pTimeElapsed;

    //This is the actual monte carlo step
    boolVar = discoverEvent(deathRate, birthRate, currentTimeStep);

    return boolVar;

}

void allocateMemoryForSizeOfGenome(int maxPopSize, int genomeSize, struct individual** pointerTopopArray, struct individual *popArray) {

    int i = 0;
    int j = 0;

    const double ZEROVALUE = 0.0;

    while (i < maxPopSize) {

        pointerTopopArray[i]->genome = malloc(sizeof(double) * genomeSize);

        for (j = 0; j < genomeSize; j++) {
            pointerTopopArray[i]->genome[j] = ZEROVALUE;
            fprintf(veryverbosefilepointer, "%lf\n", popArray[i].genome[j]);
            fflush(veryverbosefilepointer);
        }

        fprintf(veryverbosefilepointer, "new individual \n");
        fflush(veryverbosefilepointer);

        i++;

    }

}

//discuss with joseph
//should a seperate fuction be created dealing with absolute fitness
void InitializePopulation(long double* wholepopulationwistree, int populationsize, int singleIndividualGenomeLength, long double* psumofwis, long double* pInverseSumOfWis, struct individual* popArray, int* indexArray, int* freeIndexes, int MAX_POP_SIZE)
{
    int i, j, k;

    double haploidgenomelength = (double)((singleIndividualGenomeLength) / 2);

    //right here is where part of the array is uninitialized
    for (i = 0; i < populationsize; i++) {

        popArray[i].wis = 1.0;

        popArray[i].deathRate = 1.0;

        wholepopulationwistree[i] = 1.0; //for relative fitness, all individuals start with probability of being chosen as a parent of 1/N this is a tree
        /*
        for (k = 0; k < singleIndividualGenomeLength; k++) {
            popArray[i].genome[k] = 0.0;
        }
        */

        popArray[i].index = i;
    }

    for (i = 0; i < MAX_POP_SIZE; i++) {

        if(i <= populationsize){
        	freeIndexes[i] = -1;
        	indexArray[i] = i;
        }
        else{
        	freeIndexes[i] = i;
        	indexArray[i] = -1;
        }

    }

    //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
    for (i = 0; i < populationsize; i++) {
        j = i + LSB(i + 1);
        if (j <= populationsize) {
            wholepopulationwistree[j] += wholepopulationwistree[i];
        }
    }

    *psumofwis = *psumofwis + (long double)populationsize;
    *pInverseSumOfWis = *pInverseSumOfWis + (long double)populationsize;
}

int SearchTree(int leftbound, int rightbound, long double targetvalue, long double* Fenwicktree, int eventNumber)
{
    //double sumOfInverseFitnesses;

    int middle;
    middle = floor((leftbound + rightbound) / 2);
    long double partialsumatmiddle;
    long double partialsumatmiddleminusone;
    partialsumatmiddle = Fen_sum(Fenwicktree, middle);
    partialsumatmiddleminusone = Fen_sum(Fenwicktree, middle - 1);

    if (partialsumatmiddle < targetvalue) {
        if ((middle + 1) == rightbound) {
            return middle;
        }
        return SearchTree(middle, rightbound, targetvalue, Fenwicktree, eventNumber);
    }
    if (partialsumatmiddle > targetvalue) {
        if (partialsumatmiddleminusone > targetvalue) {
            return SearchTree(leftbound, middle, targetvalue, Fenwicktree, eventNumber);
        }
        else {
            return (middle - 1);
        }
    }
    if (partialsumatmiddle == targetvalue) {
        return middle;
    }
}

//The tree in the name is the Fenwick tree, which stores the fitnesses of individuals in the population.
//This function is where selection occurs -- individuals with higher-than-average fitness will be chosen more often as parents.
//edit to make choose parent with tree a random selection
//changed name
int ChooseParent(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);

    return randomindividual;
}
//edit

//In this model, individuals die at random. There's no selection happening here.
//swapped choose parent with tree and choose victim
//returns number of vicim does not select from array
int ChooseVictimWithTree(long double* wholepopulationwistree, int popsize, long double sumofwis, long double inverseSumOfWis, int eventNumber)//using pinversesum intesting and is currently unchanged 11/18/2019
{
    long double randomnumberofdeath;
    int newVictim = 0;
    randomnumberofdeath = (long double)ldexp(pcg32_random(), -32) * (inverseSumOfWis);
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize;//change variable name
    if (leftbound >= rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    //the random death is causing a strange number
    newVictim = (SearchTree(leftbound, rightbound, randomnumberofdeath, wholepopulationwistree, eventNumber));//fixed possible error
    return newVictim;
}

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(int persontorecombine, int chromosomesize, int numberofchromosomes, double* gamete, int totalindividualgenomelength, int* recombinationsites, struct individual* popArray, int *indexArray, double currentTimeStep)
{
    int recombinationsite;
    int startchromosome;
    int startofindividual;
    int returnvaluefortskit;
    int indexIdentifier;
    int h, i;

    double recombinedGenome;
    double checkingData;
    double dataCheck = 100;

    indexIdentifier = indexArray[persontorecombine];

    if (REDEFINECHROM == 1) {
    	fprintf(veryverbosefilepointer, "\nFirst loop\n");
        fflush(veryverbosefilepointer);
    }

    for (h = 0; h < numberofchromosomes; h++) {
        startchromosome = pcg32_boundedrand(2);

        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);

        } while (recombinationsite == 0); //it doesn't make sense to do a recombination event before the first linkage block. Note that this will never break if the chromosome size is only one linkage block.
        

        for (i = 0; i < recombinationsite; i++) {
            if (startchromosome == 0) {

                gamete[h * chromosomesize + i] = popArray[indexIdentifier].genome[(h * chromosomesize) + i];

            }
            else {

                gamete[h * chromosomesize + i] = popArray[indexIdentifier].genome[totalindividualgenomelength / 2 + (h * chromosomesize) + i];

            }
        }

        for (i = recombinationsite; i < chromosomesize; i++) {
            if (startchromosome == 0) {//most likely will always be zero

                checkingData = popArray[indexIdentifier].genome[totalindividualgenomelength / 2 + (h * chromosomesize) + i];

                gamete[h * chromosomesize + i] = popArray[indexIdentifier].genome[totalindividualgenomelength / 2 + (h * chromosomesize) + i];

                if (checkingData > dataCheck) {
                    fprintf(veryverbosefilepointer, "data in flip position after event %lf\n\n\n", gamete[h * chromosomesize + i]);
                    fflush(veryverbosefilepointer);
                }

            }
            else {

                gamete[h * chromosomesize + i] = popArray[indexIdentifier].genome[(h * chromosomesize) + i];

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
    }
    else { 				//for larger means, use rejection method.
        if (poissonmean != oldmean) {	//for new means, pre-compute some functions.
            oldmean = poissonmean;
            sq = sqrt(2.0 * poissonmean);
            logmean = log(poissonmean);
            g = poissonmean * logmean - gsl_sf_lngamma(poissonmean + 1.0); //lngamma function is the natural log of the gamma function
        }
        do {
            do {
                y = tan(PI * ldexp(pcg32_random(), -32)); 	//makes y a deviate from a Lorentzian comparison function.
                numberofmutations = sq * y + poissonmean;		//shifts and scales y and sets results as possible numberofmutations (to be accepted or rejected);
            } while (numberofmutations < 0.0); 			//rejects values in zero probability area.
            numberofmutations = floor(numberofmutations);
            t = 0.9 * (1.0 + y * y) * exp(numberofmutations * logmean - gsl_sf_lngamma(numberofmutations + 1.0) - g);
        } while (ldexp(pcg32_random(), -32) > t);
    }

    return numberofmutations;
}

int DetermineNumberOfMutations(int chromosomesize, int numberofchromosomes, float mutationrate)
{
    float meannumberofmutations = mutationrate * (float)chromosomesize * (float)numberofchromosomes;

    //Note that because this function operates on gametes, the calculation above appears haploid.
    //There shouldn't be a multiplication by 2 (for diploidy) in this function, since it will be called twice per individual: once per gamete.
    //Above calculation should be moved outside this function for increased efficiency.

    int numberofmutations = SampleFromPoisson(meannumberofmutations);

    return numberofmutations;
}

void MutateGamete(int chromosomesize, int numberofchromosomes, double* gamete, double mutationeffectsize)
{

    int randomchromosometomutate = pcg32_boundedrand(numberofchromosomes); //if we decide to include heterogenous rates of recombination/mutation, both of these will need to be replaced by a function that weights each linkage block's probability of mutating.
    int randomblocktomutate = pcg32_boundedrand(chromosomesize);
    int mutatedsite = randomchromosometomutate * chromosomesize + randomblocktomutate;
    gamete[mutatedsite] += log(1 + mutationeffectsize);
    char derivedstate[100];
    sprintf(derivedstate, "%10.9f", mutationeffectsize);

}

//this function seems inefficient, but with recombination and mutation, I'm not sure there's a significantly easier way.
double CalculateWi(double* parent1gamete, double* parent2gamete, int totalindividualgenomelength)
{
    double newwi = 0.0;
    long double currentlinkageblockssum = 0.0;
    int i;

    for (i = 0; i < (totalindividualgenomelength / 2); i++) {
        currentlinkageblockssum += parent1gamete[i];
        currentlinkageblockssum += parent2gamete[i];
    }
    newwi = exp(currentlinkageblockssum);

    return newwi;
}

void performBirth(double* parent1gamete, double* parent2gamete, int* pCurrentPopsize, int currentvictim, long double* pSumOfWis, int totalindividualgenomelength, long double* wholepopulationwistree, long double* pInverseSumOfWis, struct individual* popArray, int* freeIndexes, int* arrayOfIndexes, const int MAX_POP_SIZE, int eventNumber, double currentTimeStep, double averageDeathRate)
{
    int i;
    int freeIndex;
    int indexOfStructure;
    int usedFreeIndex;
    double newWi;
    double newInverseWi;

    //parses through the free index array looking for a unused index
    for(i = 0; i < MAX_POP_SIZE; i++){
    	freeIndex = freeIndexes[i];

    	//case where the free index is not zero
    	if(freeIndex != -1){
            usedFreeIndex = i;
    		break;
    	}
    	else{

    	}
    }

    //error message
	if((VERBOSE == 1) && (freeIndex == 0)){
		fprintf(verbosefilepointer, "Perform birth takes from final index that is zero\n");
		fflush(verbosefilepointer);
	}

    //calculates the fitness for the new indivdual
    newWi = CalculateWi(parent1gamete, parent2gamete, totalindividualgenomelength);

    //Creates the new recombined genome for the born individual
	for (i = 0; i < (totalindividualgenomelength / 2); i++) {

		popArray[freeIndex].genome[i] = parent1gamete[i];
		popArray[freeIndex].genome[(totalindividualgenomelength / 2) + i] = parent2gamete[i];

	}

    newInverseWi = 1.0 / newWi;

    //Optional statement for quick population death
    if(QUICK_DEATH == 1){
    	*pSumOfWis = *pSumOfWis + newWi;
    	*pInverseSumOfWis += *pInverseSumOfWis + newInverseWi;
    }
    else{
    	*pSumOfWis = *pSumOfWis + newWi;
    	*pInverseSumOfWis = *pInverseSumOfWis + newInverseWi;
    }

    *pCurrentPopsize = *pCurrentPopsize + 1;

    if (*pCurrentPopsize > MAX_POP_SIZE) {
    	if(VERBOSE == 1){
    		fprintf(verbosefilepointer, "At event %d a birth occured over the population limit", eventNumber);
    		fflush(veryverbosefilepointer);
    	}
    }
    //fills out data for the newborn
    popArray[freeIndex].wis = newWi;
    popArray[freeIndex].deathRate = newInverseWi;
    popArray[freeIndex].index = freeIndex;

    if ((currentTimeStep > 0.672708) && (currentTimeStep < 0.683801)) {
        fprintf(veryverbosefilepointer, "death rate %lf \n", newInverseWi);
        fflush(veryverbosefilepointer);
    }

    //fills in data of indexes at newest point to ensure order
    arrayOfIndexes[*pCurrentPopsize] = freeIndex;

    freeIndexes[usedFreeIndex] = -1;

    Fen_add(wholepopulationwistree, MAX_POP_SIZE, newInverseWi, freeIndex);

}

void performDeath(int* pCurrentPopsize, int currentvictim, long double* sumofwis, long double* wholepopulationwistree, long double* pInverseSumOfWis, struct individual* popArray, int* arrayOfIndexes, int* arrayOfFreeIndexes, int MAX_POP_SIZE, double currentTimeStep, int eventNumber, int totalindividualgenomelength, long double *lastDeathRateGone) {
    const int ZERO = 0;
    const int NEGATIVE_ONE = -1;

    int i, j;
    int structFreed;
    int tempIndex;
    int openIndex;
    int sizeOfFreeIndexes;
    int newFreeIndex;
    int indexOfNewFreeIndex;
    double negativeInverseWi;

    //looking for a free index to fill out
    for (i = 0; i < MAX_POP_SIZE; i++) {
        indexOfNewFreeIndex = i;

        if (-1 == arrayOfFreeIndexes[indexOfNewFreeIndex]) {
            //if this reaches a free index with data corresponding data no error is sent
            break;
        }
        else {

        }
    }

    if ((VERBOSE == 1) && (arrayOfFreeIndexes[indexOfNewFreeIndex] != -1)) {
        fprintf(verbosefilepointer, "Perform death takes from final index that is already taken\n");
        fflush(verbosefilepointer);
    }

    //the victims index will be moved to the popsize a replaced
    openIndex = *pCurrentPopsize - 1;

    //flips the victims index with the last index of the array at popsize

    arrayPositionFlip(openIndex, arrayOfIndexes, currentvictim);

    //edits the values of sum of wis and inverse sum of wis
    *sumofwis = *sumofwis - popArray[arrayOfIndexes[openIndex]].wis;
    *pInverseSumOfWis = *pInverseSumOfWis - popArray[arrayOfIndexes[openIndex]].deathRate;
    *lastDeathRateGone = popArray[arrayOfIndexes[openIndex]].deathRate;

    //For the sum of the Fennwick tree this subtracts the value
    negativeInverseWi = (NEGATIVE_ONE) * (popArray[arrayOfIndexes[openIndex]].deathRate);//switched from fitness to death rate
    Fen_add(wholepopulationwistree, MAX_POP_SIZE, negativeInverseWi, currentvictim);

    //sends index within array of indexes to array of free indexes

    arrayOfFreeIndexes[indexOfNewFreeIndex] = arrayOfIndexes[openIndex];

    //As a death occured the popsize must decrease
    *pCurrentPopsize = *pCurrentPopsize - 1;
}

void ProduceMutatedRecombinedGamete(int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, int currentparent, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char* beneficialdistribution, double* parentgamete, gsl_rng* randomnumbergeneratorforgamma, struct individual* popArray, int* arrayOfIndexes, double currentTimeStep)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    double Sds[30];
    int recombinationsites[numberofchromosomes];

    //Following lines produce a gamete from parent 1 and add deleterious and beneficial mutations to the gamete.
    RecombineChromosomesIntoGamete(currentparent, chromosomesize, numberofchromosomes, parentgamete, totalindividualgenomelength, recombinationsites, popArray, arrayOfIndexes, currentTimeStep);

    //Following lines stochastically generate a number of deleterious mutations drawn from a Poisson distribution with mean determined by the deleterious mutation rate
    //with effect sizes drawn from a gamma distribution with parameters taken from Kim et al 2017.
    int DontBreakWhileLoop = 0;
    while (1) {
        DontBreakWhileLoop = 0;
        numberofdeleteriousmutations = DetermineNumberOfMutations(chromosomesize, numberofchromosomes, deleteriousmutationrate);//I have doughts a bug lies here.
        for (k = 0; k < numberofdeleteriousmutations; k++) {
            Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4) / 23646); //Uses parameters for the gamma distribution of the selection coefficients of new mutations scaled to an inferred ancestral populations size. To produce the distribution of unscaled effect sizes, numbers drawn from this distribution must be divided by two times the ancestral population size for the population from which the distribution was derived (11,823 in this case). Data used to produce these fits were samples from 6503 individuals from the National Heart, Lung, and Blood Institute European-American dataset. Analysis of DFE from Kim et al. 2017.
            if (REDEFINECHROM == 1) {
                fprintf(veryverbosefilepointer, "random numbertaken from generator %lf \n", Sds[k]);
                fflush(veryverbosefilepointer);
            }
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
    }
    else if (strcmp(beneficialdistribution, "exponential") == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, generatedSb);
        }
    }
    else if (strcmp(beneficialdistribution, "uniform") == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, generatedSb);
        }
    }
    else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        for (k = 0; k < numberofbeneficialmutations; k++) {
            MutateGamete(chromosomesize, numberofchromosomes, parentgamete, Sb);
        }
    }
}

void PerformOneTimeStep(int* pPopSize, int currenttimestep, long double* wholepopulationwistree, double* wholepopulationgenomes, long double* psumofwis, long double* pInverseSumOfWis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char* beneficialdistribution, double* parent1gamete, double* parent2gamete, gsl_rng* randomnumbergeneratorforgamma, int birthBool, struct individual* popArray, int* freeIndexes, int* arrayOfIndexes, int eventNumber, int maxPopSize, double currentTimeStep, double averageDeathRate, long double *plastDeathRateGone)
{

    if ((PERFORMONETIMESTEPMARKERS == 1) && (currentTimeStep >= 0.03)) {
        fprintf(veryverbosefilepointer, "pop size is currently %d\n", *pPopSize);
        fflush(veryverbosefilepointer);
    }

    int currentparent1, currentparent2, currentvictim;
    int popSize = *pPopSize;
    int i;

    const int BIRTH_OCCURS = 1;
    const int DEATH_OCCURS = 0;

    if (birthBool == BIRTH_OCCURS) {

        //choose parent returns a random int in the range of the popsize. This is not the position in the popArray but the position in the arrayOfIndexes
        currentparent1 = ChooseParent(popSize);//will corsepond with index of individual does not matter as the individual is random anyway putting an extra layer to the randomness should change nothing

        currentparent2 = ChooseParent(popSize);//need to select two parents

        if ((VERYVERBOSE == 1) && (currentTimeStep >= 0.03)) {
            fprintf(veryverbosefilepointer, "%d\n", arrayOfIndexes[currentparent2] );
            fprintf(veryverbosefilepointer, "Parent 2 choosen. \n");
            fflush(veryverbosefilepointer);
        }

        //Possiblity the same parent will be chosen, this is a simple work around that repeats the selection process on one until a new person is selected.
        while ( (currentparent1 == currentparent2) ) { //probably not ideal, since it'll never break with population sizes of zero or one.
            currentparent2 = ChooseParent(popSize);
        }

        if ((PERFORMONETIMESTEPMARKERS == 1)) {
            fprintf(veryverbosefilepointer, "Parental index 1 is %d\n", currentparent1);
            fprintf(veryverbosefilepointer, "Parental index 2 is %d\n", currentparent2);
            fprintf(veryverbosefilepointer, "current time step is %lf\n", currentTimeStep);
            fprintf(veryverbosefilepointer, "Passes choose parent final loop. \n");
            fflush(veryverbosefilepointer);
        }

        //this function mutates and recombindes the gametes.
        ProduceMutatedRecombinedGamete(chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent1, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, randomnumbergeneratorforgamma, popArray, arrayOfIndexes, currentTimeStep);
        ProduceMutatedRecombinedGamete(chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent2, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent2gamete, randomnumbergeneratorforgamma, popArray, arrayOfIndexes, currentTimeStep);

        fprintf(veryverbosefilepointer, "new individual /n/n");
        fflush(veryverbosefilepointer);

        if ((currentTimeStep > 0.672708) && (currentTimeStep < 0.683801)) {
            for (i = 0; i < totalindividualgenomelength; i++) {
                fprintf(veryverbosefilepointer, "%lf\n", popArray[currentparent1].genome[i]);
                fprintf(veryverbosefilepointer, "%lf\n", popArray[currentparent2].genome[i]);
                fflush(veryverbosefilepointer);
            }
        }

        if (PERFORMONETIMESTEPMARKERS == 1) {
            fprintf(veryverbosefilepointer, "Genes recombined. \n");
            fflush(veryverbosefilepointer);
        }

        //Replaces data of a dead individual with an alive one and moves data between freeIndexes and arrayOfIndexes
        performBirth(parent1gamete, parent2gamete, pPopSize, currentvictim, psumofwis, totalindividualgenomelength, wholepopulationwistree, pInverseSumOfWis, popArray, freeIndexes, arrayOfIndexes, maxPopSize, eventNumber, currentTimeStep, averageDeathRate);

    }
    else if (birthBool == DEATH_OCCURS) {

        //use of the fennwick tree to find the victim in the population
    	currentvictim = ChooseVictimWithTree(wholepopulationwistree, popSize, *psumofwis, *pInverseSumOfWis, eventNumber);

        //function moves index data from arrayOfIndexes to freeIndexes. No data is deleted or overwriten here as 
        performDeath(pPopSize, currentvictim, psumofwis, wholepopulationwistree, pInverseSumOfWis, popArray, arrayOfIndexes, freeIndexes, maxPopSize, currentTimeStep, eventNumber, totalindividualgenomelength, plastDeathRateGone);

    }
    else{
        //Putting an else function at the end of a if and else if section to ensure there is exit where nothing occurs and I can bug check for this. If you wish you could delete this.
    }

}

double RunSimulation(char* Nxtimestepsname, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* Sbname, char* typeofrun, int timeSteps, int initialPopSize, int maxPopSize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char* beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma, FILE* veryverbosefilepointer, FILE* verbosefilepointer, FILE* miscfilepointer)
{

    int i, j, k, w;

    char* rawdatafilename = (char*)malloc(sizeof(char) * 200);//editied slightly if everything blows up definitly this (11/25/2019)
    strcpy(rawdatafilename, "rawdatafor"); //starting the string that will be the name of the data file.

    strcat(rawdatafilename, "timeSteps"); //for adding values of generations to the data name.
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
    fprintf(rawdatafilepointer, "Time, Pop size, average death rate\n");
    fprintf(veryverbosefilepointer, "Wi, Free Index");
    fflush(veryverbosefilepointer);

    char* summarydatafilename = (char*)malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.

    int popSizeFiller;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    int birthBool;
    int popsize;
    int endofburninphase;
    int isburninphaseover = 0;
    int didpopulationcrash = 0;
    int endofdelay = timeSteps - 1;
    int endofsimulation = timeSteps - 1;
    int timeStepsAfterBurnin = 0;
    int EventsPreformed = 0;
    int avgPopsizeForOneRun = 0;
    int switchForChart = 0;

    double currenttimestep = 0;
    double numberOfTimeStepsBetweenEvents;
    double varianceOfPopulationForLastTenThousandTimeSteps;
    double c0, cov00, cov01, cov11, sumsq;
    double parent1gamete[numberofchromosomes * chromosomesize], parent2gamete[numberofchromosomes * chromosomesize];
    double slopeofvariance;
    double slopeoflogfitness;
    double variancesum;
    double variancePop;
    double averagePopsize;
    double arbitrarynumber = 1;
    double doubleTimeStep = timeSteps * 1.0;
    double previousSumOfWis;
    double changeOfSum;

    long double avgDeathRate = 0;
    long double sumofwis = 0;
    long double inverseSumOfWis = 0;
    long double currentfittestindividualswi;
    long double lastDeathRateGone = 0;

    int* arrayOfIndexes;
    int* arrayOfFreeIndexes;
    int* pPopSize;

    double* pCurrentTimeStep;
    double* popSizeArrayForAverage;
    double* wholepopulationgenomes;
    double* logaveragefitnesseachNtimesteps;
    double* last200Ntimestepsvariance;
    double* literallyjustlast200Ntimesteps;
    double* popVarianceArray;
    double* pNumberOfTimeStepsBetweenEvents;

    long double* psumofwis;
    long double* pInverseSumOfWis;
    long double* wholepopulationwistree;
    long double* wholepopulationwisarray;
    long double* plastDeathRateGone;

    struct individual emptyCase;

    struct individual* popArray;
    size_t step = 1;

    struct individual** pointerTopopArray = &popArray;

    totalpopulationgenomelength = maxPopSize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    popsize = initialPopSize;

    //wrapper for structure initilization
    emptyCase.deathRate = 0.0;
    emptyCase.genome = malloc(sizeof(double) * totalindividualgenomelength);
    for (j = 0; j < totalindividualgenomelength; j++) {
        emptyCase.genome[j] = 0.0;
    }
    emptyCase.index = 0;
    emptyCase.wis = 0.0;

    pCurrentTimeStep = &currenttimestep;
    pPopSize = &popsize;
    psumofwis = &sumofwis;
    pInverseSumOfWis = &inverseSumOfWis;
    pNumberOfTimeStepsBetweenEvents = &numberOfTimeStepsBetweenEvents;
    plastDeathRateGone = &lastDeathRateGone;

    arrayOfIndexes = (int*)malloc(sizeof(int) * maxPopSize);//Array responsible for ordering the popArray, contains all living individuals
    arrayOfFreeIndexes = (int*)malloc(sizeof(int) * maxPopSize);//Array to keep track of holes within the array of indexes
    wholepopulationgenomes = (double*)malloc(sizeof(double) * totalpopulationgenomelength);
    popArray = (struct individual*)(malloc(sizeof(struct individual) * maxPopSize));//Array of data for every individual in the population
    popSizeArrayForAverage = (double*)malloc(sizeof(double) * maxPopSize);
    wholepopulationwistree = (long double*)malloc(sizeof(long double) * maxPopSize);
    logaveragefitnesseachNtimesteps = (double*)malloc(sizeof(double) * timeSteps);
    literallyjustlast200Ntimesteps = (double*)malloc(sizeof(double) * 200);//use malloc
    last200Ntimestepsvariance = (double*)malloc(sizeof(double) * 200);//use malloc
    popVarianceArray = (double*)malloc(sizeof(double) * maxPopSize);

    for (k = 0; k < 200; k++) {
        literallyjustlast200Ntimesteps[k] = 0.0;
        last200Ntimestepsvariance[k] = 0.0;
    }

    fprintf(veryverbosefilepointer, "Entered simulation run.\n");
    fflush(veryverbosefilepointer);

    i = 0;
    j = 0;

    //assignment of wrapper to all indexes to popArray
    for (i = 0; i < maxPopSize; i++) {

        popArray[i] = emptyCase;

    }

    //assignment of data to popArray within genome
    for (i = 0; i < maxPopSize; i++) {

        popArray[i].genome = malloc(sizeof(double) * totalindividualgenomelength);

        for (j = 0; j < totalindividualgenomelength; j++) {
            popArray[i].genome[j] = 0;
        }

    }

    fprintf(veryverbosefilepointer, "Memory allocated for genomes. \n");
    fflush(veryverbosefilepointer);

    //assignment of data to popArray for index, wis, and deathrate
    InitializePopulation(wholepopulationwistree, initialPopSize, totalindividualgenomelength, psumofwis, pInverseSumOfWis, popArray, arrayOfIndexes, arrayOfFreeIndexes, maxPopSize);//sets up all data within the population for a run. As this initializes data I think it should be a seperate funtion.

    fprintf(veryverbosefilepointer, "Population initialized.\n");
    fflush(veryverbosefilepointer);

    if (RUNSIMULATIONMARKERS == 1) {
        fprintf(verbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
        fflush(verbosefilepointer);
    }

    //BEGIN THE SIMULATION FOR LOOP

    if (RUNSIMULATIONMARKERS == 1) {
        fprintf(veryverbosefilepointer, "Enters simulation run \n");
        fflush(veryverbosefilepointer);
    }
    //doubleTimeStep is a placeholder name, this is a maximum number of time steps you wish for the run to take. 
    while (currenttimestep <= doubleTimeStep) {

    	if(popsize < 3){
        		fprintf(summarydatafilepointer, "Population died during run at time step %d", EventsPreformed);
        		fflush(summarydatafilepointer);

        		break;
        }

        birthBool = monteCarloStep(popsize, sumofwis, pCurrentTimeStep, inverseSumOfWis, popArray, maxPopSize, currenttimestep);//This is the monte carlo step. This decides if a birth or a death event takes place by returning a 0 or 1

        //function that is responsible for all calls that would occur during a birth and death simulation.
		PerformOneTimeStep(pPopSize, currenttimestep, wholepopulationwistree, wholepopulationgenomes, psumofwis, pInverseSumOfWis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma, birthBool, popArray, arrayOfFreeIndexes, arrayOfIndexes, EventsPreformed, maxPopSize, currenttimestep, avgDeathRate, plastDeathRateGone);

		avgDeathRate = inverseSumOfWis/popsize;

        previousSumOfWis = inverseSumOfWis;

        if (RUNSIMULATIONMARKERS == 1) {

            if (avgDeathRate < 0) {
                fprintf(veryverbosefilepointer, "average death rate become negtive\n");
                fprintf(veryverbosefilepointer, "%Lf\n", inverseSumOfWis);
                fprintf(veryverbosefilepointer, "Occurs in time step %lf\n", currenttimestep);
                fprintf(veryverbosefilepointer, "Death rate is equal to %Lf\n", avgDeathRate);
                fflush(veryverbosefilepointer);

                fprintf(summarydatafilepointer, "\n\nValues of death Rate\n");
                for (k = 0; k < popsize; k++) {
                    //fprintf(veryverbosefilepointer, "%d,%lf\n", k + 1, popArray[arrayOfIndexes[k]].deathRate);
                    //fflush(veryverbosefilepointer);
                }

                exit(0);
            }

            fflush(veryverbosefilepointer);
            
        }

		if(((EventsPreformed % 1000) == 0) || ((currenttimestep > 2040) && (currenttimestep < 2060))){

			fprintf(rawdatafilepointer, "%lf,", currenttimestep);
			fprintf(rawdatafilepointer, "%d,", popsize);
			fprintf(rawdatafilepointer, "%Lf\n", avgDeathRate);
			fflush(rawdatafilepointer);

		}

        //This is to produce a histogram of the wis of the entire population from a single generation.
        //It's terrible and completely non-modular, but I just can't bring myself to add in two more user-input arguments.

        if (strcmp(typeofrun, "single") == 0) {

            if ( (currenttimestep >= 27.0) && (switchForChart == 0) ) {

            	fprintf(veryverbosefilepointer, "Just before individual wi data lines.\n");
            	fflush(veryverbosefilepointer);

                if (INDIVIDUALWIDATA == 1) {
                    //fprintf(summarydatafilepointer, "Individual, death rates\n");
                	fprintf(summarydatafilepointer, "Values in array of indexes\n");
                    for (k = 0; k < maxPopSize; k++) {
                    	fprintf(summarydatafilepointer, "%d,%d\n", k + 1, arrayOfIndexes[k]);
                        fflush(summarydatafilepointer);
                    }

                    fprintf(summarydatafilepointer, "\n\nValues in array of free indexes\n");
                    for (k = 0; k < maxPopSize; k++) {
                    	fprintf(summarydatafilepointer, "%d,%d\n", k + 1, arrayOfFreeIndexes[k]);
                        fflush(summarydatafilepointer);
                    }

                    fprintf(summarydatafilepointer, "\n\nValues of death Rate\n");
                    for (k = 0; k < popsize; k++) {
                        fprintf(summarydatafilepointer, "%d,%lf\n", k + 1, popArray[arrayOfIndexes[k]].deathRate);
                        fflush(summarydatafilepointer);
                    }
                }

                switchForChart = 1;

            }

        }

        //If the burn-in phase has been called, wait 500 generations to start recording fitnesses.
        //This is to be sure that even when the beneficial rates/sizes are large, the only generations recorded will be from the uniformly sloping part of the simulation.
        //The average fitness from any generation after this delay period is recorded in the array of average fitnesses.

        if (i > endofdelay) {
            logaveragefitnesseachNtimesteps[timeStepsAfterBurnin] = log((double)*psumofwis / (double)popsize);//This refers to a completed time step not the time steps relative to the algorithm.

            if (VERYVERBOSE == 1) {
                fprintf(veryverbosefilepointer, "log average fitness in event number %d, %d timesteps after burn-in, is: %f\n", EventsPreformed, timeStepsAfterBurnin, logaveragefitnesseachNtimesteps[timeStepsAfterBurnin]);
                fflush(veryverbosefilepointer);
            }
            timeStepsAfterBurnin += 1;
        }

        /*
        if (RUNSIMULATIONMARKERS == 1) {
            fprintf(veryverbosefilepointer, "Past end of delay \n");
            fflush(veryverbosefilepointer);
        }
        */

        //These lines ensure that the magnitude of fitness hasn't declined by too much.
        //At extremely small fitness values, floating-point math becomes imprecise.
        //These lines end the simulation if fitness declines below 10^-10, which should represent a completely degraded population.
        currentfittestindividualswi = FindFittestWi(popsize, popArray, arrayOfIndexes);
        if (currentfittestindividualswi < pow(10.0, -10.0)) {

            fprintf(miscfilepointer, "\nFitness declined to less than 10^-10 during generation %d.", EventsPreformed + 1);
            fprintf(summarydatafilepointer, "Fitness declined to catastrophic levels in generation %d.\n", EventsPreformed + 1);

            endofsimulation = EventsPreformed;
            EventsPreformed = timeSteps;
            didpopulationcrash = EventsPreformed;

            break;

        }
        /*
        if (RUNSIMULATIONMARKERS == 1) {
            fprintf(veryverbosefilepointer, "Current fittest past \n");
            fflush(veryverbosefilepointer);
        }
        */
        EventsPreformed++;

        storePopsizeForGenerations(popSizeArrayForAverage, popsize, maxPopSize);

        /*
        if (RUNSIMULATIONMARKERS == 1) {
            fprintf(veryverbosefilepointer, "Store popsize for generations completed \n");
            fflush(veryverbosefilepointer);
        }
        */

        if(EventsPreformed % maxPopSize == 0){

        	averagePopsize = calcAvgPopSize(popSizeArrayForAverage);
        	i++;

        	//moved for absolute fitness as there should be a check at the end of every generation also the code should just work a bit better
            if (isburninphaseover == 0) {
            	UpdateLast200NTimeSteps(last200Ntimestepsvariance, averagePopsize);
                UpdateLast200NTimeSteps(literallyjustlast200Ntimesteps, i + 1);

                /*
                if (RUNSIMULATIONMARKERS == 1) {
                    fprintf(veryverbosefilepointer, "Array update performed \n");
                    fflush(veryverbosefilepointer);
                }
                */

                //to avoid calling the end of the burn-in phase at generation one
                //because of setting pre-simulation generations to zeroes
                //I just won't start looking for the end of the burn-in phase until 200 generations
                //This would be a mild problem if a simulation should end in 200 generations, but that shouldn't ever happen with the DFE I'm using.
                if ((EventsPreformed/maxPopSize) > 200) {
                    slopeofvariance = 0.0;
                    gsl_fit_linear(literallyjustlast200Ntimesteps, step, popSizeArrayForAverage, step, 200, &c0, &slopeofvariance, &cov00, &cov01, &cov11, &sumsq);
                    if (slopeofvariance < arbitrarynumber) {
                        endofburninphase = EventsPreformed;
                        endofdelay = endofburninphase + 500;
                        isburninphaseover = 1;
                        fprintf(miscfilepointer, "Burn-in phase called as ending in generation %d\n", EventsPreformed + 1);
                        fprintf(summarydatafilepointer, "Burn-in phase called as ending in generation %d\n", EventsPreformed + 1);
                        if (VERBOSE == 1) {
                            fflush(miscfilepointer);
                            fflush(summarydatafilepointer);
                        }

                    }
                }
            }

        }

    }

    //END OF SIMULATION FOR LOOP


    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %.6f. Final population summed fitness was: %Lf\n", Sb, *psumofwis);
        fflush(veryverbosefilepointer);
    }
    if (didpopulationcrash == 0) {
        endofsimulation = EventsPreformed;
    }
    if (isburninphaseover == 1) {
    	/*
        if (VERYVERBOSE == 1) {
            fprintf(veryverbosefilepointer, "Calculating slope of log fitness with the following parameters: endofsimulation = %d, endofdelay = %d, generationsafterburnin = %d\nLog fitness each generation: ", endofsimulation, endofdelay, timeStepsAfterBurnin);
            for (j = 0; j < (endofsimulation - endofdelay); j++) {
                fprintf(veryverbosefilepointer, "%f ", logaveragefitnesseachNtimesteps[j]);
            }
            fprintf(veryverbosefilepointer, "\n");
            fflush(veryverbosefilepointer);
        }
        */

        slopeoflogfitness = CalculateSlopeOfLogFitness(endofsimulation, endofdelay, currenttimestep, logaveragefitnesseachNtimesteps);
        fprintf(summarydatafilepointer, "Slope of log(fitness) after the burn-in phase: %f\n", slopeoflogfitness);

        /*
        if (VERBOSE == 1) {
            fflush(summarydatafilepointer);
            fflush(rawdatafilepointer);
        }
        */

        fclose(rawdatafilepointer);
        fclose(summarydatafilepointer);
        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(arrayOfFreeIndexes);
        free(arrayOfIndexes);
        for (i = 0; i < maxPopSize; i++) {
            free(popArray[i].genome);
        }
        free(popArray);

        return slopeoflogfitness;//should return average fitness
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
        free(arrayOfFreeIndexes);
        free(arrayOfIndexes);
        for (i = 0; i < maxPopSize; i++) {
            free(popArray[i].genome);
        }
        free(popArray);

        return -1.0;
    }
}

//The following function is heavily modified from Numerical Recipes in C, Second Edition.
//For large population sizes, populations with mean Sb > 0 may actually have a more negative fitness slope than mean Sb = 0.
//
int BracketZeroForSb(double* Sb1, double* Sb2, char* Nxtimestepsname, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* typeofrun, int timeSteps, int initialPopSize, int maxPopSize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, char* beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma, FILE* veryverbosefilepointer, FILE* verbosefilepointer, FILE* miscfilepointer) {
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
    resultingslope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, timeSteps, initialPopSize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
    resultingslope2 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, timeSteps, initialPopSize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
        }
        else if (resultingslope2 <= slopeforcontourline) {
            *Sb2 += factor;
            snprintf(Sb2name, 10, "%.7f", *Sb2);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb2name: %s\n", Sb2name);
                fprintf(verbosefilepointer, "Starting run with new sb2 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope2 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, timeSteps, initialPopSize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb2, resultingslope2);
                fflush(verbosefilepointer);
            }

        }
        else if (resultingslope1 >= slopeforcontourline) {
            *Sb1 -= factor;
            snprintf(Sb1name, 10, "%.7f", *Sb1);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb1name: %s\n", Sb1name);
                fprintf(verbosefilepointer, "Starting run with new sb1 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, timeSteps, initialPopSize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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
double BisectionMethodToFindSbWithZeroSlope(double* Sb1, double* Sb2, char* Nxtimestepsname, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* typeofrun, int timeSteps, int popsize, int maxPopSize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, char* beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma, FILE* veryverbosefilepointer, FILE* verbosefilepointer, FILE* miscfilepointer) {
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
    slope1 = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, timeSteps, popsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }

    slopemid = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, timeSteps, popsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb2, slopemid);
    }

    if (((slope1 - slopeforcontourline) * (slopemid - slopeforcontourline)) > 0.0) {
        fprintf(miscfilepointer, "Root not bracketed properly, with starting slopes %.10f and %.10f for a desired slope of %.6f\n", slope1, slopemid, slopeforcontourline);
        return 0.0;
    }
    root = (slope1 < slopeforcontourline) ? (factor = *Sb2 - *Sb1, *Sb1) : (factor = *Sb1 - *Sb2, *Sb2);
    for (i = 1; i <= maxtries; i++) {
        Sbmid = root + (factor *= 0.5);
        snprintf(Sbmidname, 10, "%.7f", Sbmid);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "Sbmidname: %s\n", Sbmidname);
            fprintf(verbosefilepointer, "Starting run with sb %.6f\n", Sbmid);
            fflush(verbosefilepointer);
        }
        slopemid = RunSimulation(Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, typeofrun, timeSteps, popsize, maxPopSize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, randomnumbergeneratorforgamma, veryverbosefilepointer, verbosefilepointer, miscfilepointer);
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

