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
#include "sharedfunc_flag.h"
#include "main.h"


double RunSimulationRel(bool isabsolute, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer)
{
    if(isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationRel within an absolute fitness program. \n");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    
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
    
    InitializePopulationRel(wholepopulationwistree, wholepopulationwisarray, popsize, wholepopulationgenomes, totalpopulationgenomelength, psumofwis);
    
    printf("hasta aqui bien \n");
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
            PerformOneTimeStepRel(popsize, wholepopulationwistree, wholepopulationwisarray, wholepopulationgenomes, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);
            
        }
        
        //esta linea
        printf(" %Lf \n", sumofwis);
        
        
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
        if (typeofrun == 1) {
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

void PerformOneTimeStepRel(int popsize, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    bool isabsolute = 0;
    
    int currentparent1, currentparent2, currentvictim;

    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    }
    
   
    RecombineChromosomesIntoGamete(currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(isabsolute, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, randomnumbergeneratorforgamma, miscfilepointer);
        
    RecombineChromosomesIntoGamete(currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(isabsolute, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);
               
    //next variables are only used for absolute simulations, however since InitializePopulation is a shared function I must initialize them here, even if I don't use them again
    long double *pInverseSumOfWis;
    bool *wholepopulationisfree;
    int *wholepopulationindex;
    long double *wholepopulationdeathratesarray;
    int *pPopSize;
    double d_0;
    
    PerformDeath(isabsolute, popsize, pPopSize, currentvictim, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, miscfilepointer);
    
    PerformBirth(isabsolute, parent1gamete, parent2gamete, popsize, pPopSize, currentvictim, wholepopulationgenomes, totalindividualgenomelength, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, d_0, miscfilepointer);
    
}

void InitializePopulationRel(long double *wholepopulationwistree, long double *wholepopulationwisarray, int popsize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double * psumofwis) 
{
    int i, j;
        
    for (i = 0; i < popsize; i++){
        wholepopulationwistree[i] = 1.0; //all individuals start with load 1 (probability of being chosen to produce an offspring or to die of 1/N). 
        wholepopulationwisarray[i] = 1.0;
    }
    //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
    for (i = 0; i < popsize; i++) {
        j = i + LSB(i+1);
        if (j < popsize) {
            wholepopulationwistree[j] += wholepopulationwistree[i];
        }
    }
    
    *psumofwis = (long double)popsize;
    
    for (i = 0; i < totalpopulationgenomelength; i++){
        wholepopulationgenomes[i] = 0.0;
    }
    
}

//In this model, individuals die at random. There's no selection happening here.
int ChooseVictim(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);
    return randomindividual;
}

//The tree in the name is the Fenwick tree, which stores the fitnesses of individuals in the population.
//This function is where selection occurs -- individuals with higher-than-average fitness will be chosen more often as parents.
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer)
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
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
        return -1;
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    
    newparent = (SearchTree(leftbound, rightbound, randomnumberofbirth, wholepopulationwistree));
    return newparent;
}

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
