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
#include "absolute_functions.h"
#include "sharedfunc_flag.h"
#include "global_vars.h"
#include "main.h"

double RunSimulationAbs(bool isabsolute, char* maxTimeStepsname, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* Sbname, int typeofrun, int maxTimeSteps, int initialPopSize, int maxPopSize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma)
{

    if(!isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationAbs within a relative fitness program. \n");
        exit(0);
    }
    
    int i, j, k, w;

    char* rawdatafilename = (char*)malloc(sizeof(char) * 200);//editied slightly if everything blows up definitly this (11/25/2019)
    strcpy(rawdatafilename, "rawdatafor"); //starting the string that will be the name of the data file.

    strcat(rawdatafilename, "maxTimeSteps"); //for adding values of generations to the data name.
    strcat(rawdatafilename, maxTimeStepsname);

    strcat(rawdatafilename, "initialPopSize"); //for adding values of starting population sizes to the data name.
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

    char* summarydatafilename = (char*)malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.

    int popsize = initialPopSize;
    int *pPopSize;
    pPopSize = &popsize;
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = maxPopSize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
    long double inversesumofloads;
    long double *pInversesumofloads = &inversesumofloads;
    long double *wholepopulationloadstree;
    wholepopulationloadstree = malloc(sizeof(long double) * maxPopSize);
       
    long double *wholepopulationdeathratesarray;
    wholepopulationdeathratesarray = malloc(sizeof(long double) * maxPopSize);
    //The Fenwick tree does not store each individual's wi, but rather a collection of partial sums.
    //For debugging purposes and data that requires summarizing wis, storing the wis in an array is necessary.

    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * maxPopSize);
    //In order to visualize the distribution of fitness in the population,
    //the individual fitnesses need to be sorted, which requires a separate array.
    
    bool *wholepopulationisfree;
    wholepopulationisfree = malloc(sizeof(bool) * maxPopSize);
    
    int *wholepopulationindex;
    wholepopulationindex = malloc(sizeof(int) * maxPopSize);
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Entered simulation run.\n");
    }
    
    //psumofloads is only used in relative fitness scenarios, here it is just initialized
    long double *psumofloads, *wholepopulationwisarray;
    
    //assignment of data to popArray for index, wis, and deathrate
    InitializePopulation(isabsolute, wholepopulationloadstree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, initialPopSize, maxPopSize, wholepopulationgenomes, totalpopulationgenomelength, psumofloads, pInversesumofloads);//sets up all data within the population for a run. As this initializes data I think it should be a seperate funtion.

    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Population initialized.\n");
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * maxTimeSteps);
    //In order to calculate the slope of degradation of fitness,
    //I need to store the average fitness each generation.
    long double currentfittestindividualswi;
    
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
    int endofdelay = maxTimeSteps-1;
    int endofsimulation = maxTimeSteps-1;
    int Nxtimestepsafterburnin = 0;
    double arbitrarynumber = (-1.0 * 0.007 / maxPopSize); //using a number somewhere close to the mean of the DFE for deleterious mutations.
    double slopeoflogfitness;    
    double variancesum;
    
    bool birthBool;
    double currenttime = 0;
    double *pCurrentTime;
    pCurrentTime = &currenttime;
    
    double parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
    }
    
    
    //BEGIN THE SIMULATION FOR LOOP
    
    for (i = 0; i < maxTimeSteps; i++) {
        
        if(popsize < 3){
        		fprintf(summarydatafilepointer, "Population died during run at time step %d", i);
                break;
        }
        
        if(popsize >= (maxPopSize-1)){
        		fprintf(summarydatafilepointer, "Population achieved its maximum population size %d", i);
//                 printf("entro \n");
                break;
        }
        
        birthBool = monteCarloStep(popsize, pCurrentTime, inversesumofloads, maxPopSize);//This is the monte carlo step. This decides if a birth or a death event takes place by returning a 0 or 1
        
        PerformOneEventAbs(isabsolute, birthBool, maxPopSize, pPopSize, wholepopulationgenomes, wholepopulationloadstree, wholepopulationdeathratesarray, wholepopulationisfree, wholepopulationindex, pInversesumofloads, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma);
        
//         if(*pCurrentTime < 0.1)
//             ¸
//         printf("%d %d \n", birthBool, popsize);
//         
//     for(i = 0; i < maxPopSize; i++)
//         printf("%2d  ", (i+1));
//     printf("\n");
//     
//     for(i = 0; i < maxPopSize; i++)
//         printf("%2d  ", wholepopulationisfree[i]);
//     printf("\n \n");
//     
//     
//     for(i = 0; i < maxPopSize; i++)
//         printf("%2d  ", wholepopulationindex[i]);
//     printf("\n \n");
        
//         for(j = 0; j < maxPopSize; j++)
//             printf("%d ", wholepopulationisfree[j]);
//         printf("\n");
         
//          printf("%f\n", *pCurrentTime);
        
    }
    
//     printf("salio \n");
    
    //END OF SIMULATION FOR LOOP
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %f \n", Sb);
        fprintf(veryverbosefilepointer, "Time elapsed: %f\n", currenttime);
        fprintf(veryverbosefilepointer, "Final population size was: %d\n", popsize);
    }
    
    
    fclose(rawdatafilepointer);
    fclose(summarydatafilepointer);
    free(rawdatafilename);
    free(summarydatafilename);
    
    free(logaveragefitnesseachNtimesteps);
    free(literallyjustlast200Ntimesteps);
    free(last200Ntimestepsvariance);
    
    free(wholepopulationgenomes);
    free(wholepopulationloadstree);
    free(wholepopulationdeathratesarray);
    free(sortedwisarray);
    free(wholepopulationisfree);
    free(wholepopulationindex);

//     return slopeoflogfitness;//should return average fitness
    return 0;
       
}


//returns output of either birth or death
bool monteCarloStep(int popsize, double *pCurrentTime, double sumOfDeathRates, int maxPopSize) {

    double deathRate;
    double birthRate;
    double time;
    double mean;

    const double b_0 = 10;

    deathRate = sumOfDeathRates;

    //rate of births is calculated using equation used in lab write up
    birthRate = rateOfBirthsCalc(popsize, maxPopSize, b_0);

    mean = ((1.0)/(deathRate + birthRate));
    
    time = ExponentialDerivate(mean);//draws a random number from a exponential distribution with an specified mean. This occurs because a even t is most likely to occur right after a previous event.

    *pCurrentTime += time;

    //This is the actual monte carlo step
    return discoverEvent(deathRate, birthRate);

}

double rateOfBirthsCalc(int popsize, int maxPopSize, double b_0) {

    double birthRate;
    birthRate = (b_0) * (double)popsize * (1 - ((double)popsize/(double)maxPopSize));
    return  birthRate;

}

bool discoverEvent(double deathRate, double birthRate) {

    bool boolBirth;
    
    double cutOffPoint, randomNumber;

    randomNumber = ldexp(pcg32_random(), -32);

    //if lands on the cutoff point death occurs, techincally skewing towards death but is very minor.
    cutOffPoint = deathRate/(deathRate + birthRate);

    if (randomNumber > cutOffPoint)
        boolBirth = 1;
    else
        boolBirth = 0;
    
    return boolBirth;

}


void PerformOneEventAbs(bool isabsolute, bool birthBool, int maxPopSize, int *pPopSize, double *wholepopulationgenomes, long double *wholepopulationloadstree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, long double *pInversesumofloads, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution,  double* parent1gamete, double* parent2gamete, gsl_rng* randomnumbergeneratorforgamma)
{
    if(isabsolute == 0){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs within a non absolute fitness simulation. \n");
        exit(0);
    }
    
    int randparent1, randparent2, currentparent1, currentparent2, currentvictim;
    int i;

    const int BIRTH_OCCURS = 1;
    const int DEATH_OCCURS = 0;
        
    int victim;
    
    if(*pPopSize < 2){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs with a population size lower than 2. \n");
        exit(0);
    }
    
    if(*pPopSize >= maxPopSize && birthBool == BIRTH_OCCURS){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs with a population size greater than maxPopSize. \n");
        exit(0);
    }
    
    //next pointers are only used in relative fitness scenarios, here they are just initialized
    long double *psumofloads, *wholepopulationwisarray;
    
    if (birthBool == BIRTH_OCCURS){
        randparent1 = ChooseParent(*pPopSize);
        currentparent1 = wholepopulationindex[randparent1];
        do{
            randparent2 = ChooseParent(*pPopSize);
        } while(randparent2 == randparent1);
        currentparent2 = wholepopulationindex[randparent2];
        
        if(wholepopulationisfree[currentparent1] || wholepopulationisfree[currentparent2]){
            fprintf(miscfilepointer, "\n Trying to give birth to a new individual using a parent that is not present in PerformOneEventAbs. \n");
            exit(0);
        }
        

//         printf("%f \n", deleteriousmutationrate);
        
        ProduceMutatedRecombinedGamete(wholepopulationgenomes, chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent1, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, randomnumbergeneratorforgamma);
        ProduceMutatedRecombinedGamete(wholepopulationgenomes, chromosomesize, numberofchromosomes, totalindividualgenomelength, currentparent2, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent2gamete, randomnumbergeneratorforgamma);
        
        
        PerformBirth(isabsolute, parent1gamete, parent2gamete, maxPopSize, pPopSize, victim, wholepopulationgenomes, totalindividualgenomelength, wholepopulationloadstree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, pInversesumofloads);
    }
    
    
    
    else{
//         use of the fennwick tree to find the victim in the population, note that since fenwick tree stores all the population (non occupied spaces have a fitness of 0) there is not need to use the wholepopulationindex, fenwick search already gives you the space that the selected individual occupies.
    	victim = ChooseVictimWithTree(wholepopulationloadstree, *pPopSize, maxPopSize, *pInversesumofloads);
        
//         printf("%d\n", victim);
//         
//         for(i = 0; i < maxPopSize; i++)
//             printf("%2d  ", wholepopulationisfree[i]);
//         printf("\n");
        
        PerformDeath(isabsolute, maxPopSize, pPopSize, victim, wholepopulationloadstree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, pInversesumofloads);
    }


}

int ChooseParent(int populationsize)
{
    return pcg32_boundedrand(populationsize);
}

int ChooseVictimWithTree(long double* wholepopulationloadstree, int popsize, int maxPopSize, long double inversesumofloads)
{
    long double randomnumberofdeath;
    int newVictim = 0;
    randomnumberofdeath = (long double)ldexp(pcg32_random(), -32) * (inversesumofloads);
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = maxPopSize;
    
    if (leftbound >= rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    
    
    newVictim = SearchTree(leftbound, rightbound, randomnumberofdeath, wholepopulationloadstree);
    
    return newVictim;
}

int findinindex(int *wholepopulationindex, int which, int tam){
    int i;
    bool placefound = false;
    
    for(i = 0; i < tam; i++){
        if(wholepopulationindex[i] == which){
            placefound = true;
            return i;
        }
    }
    
    if(!placefound){
        fprintf(miscfilepointer, "\n Array of indexes is corrupted in findinindex. \n");
        exit(0);
    }
}

void indexArrayFlipDeath(int *wholepopulationindex, int placeinindex, int popsize){

    /*
    * This function flips the integer value from the last in the array of  wholepopulationindex to the place where an individual died. This way wholepopulationindex is mantained without unoccupied spaces
    */
	int indexlast = wholepopulationindex[(popsize-1)];
    int indexvictim = wholepopulationindex[placeinindex];
    
    wholepopulationindex[placeinindex] = indexlast;
    wholepopulationindex[(popsize-1)] = indexvictim;
}
