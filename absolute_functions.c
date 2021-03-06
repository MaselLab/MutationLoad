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
#include "main.h"

double RunSimulationAbs(bool issnapshot, char *prevsnapshotfilename, bool isabsolute, char* maxTimename, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* Sbname, int typeofrun, int maxTime, int initialPopSize, int maxPopSize, double d_0, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma, double r, double sdmin,  FILE *miscfilepointer, FILE *veryverbosefilepointer)
{

    if(!isabsolute){
        fprintf(miscfilepointer, "\nError: Trying to use RunSimulationAbs within a relative fitness program.");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    
    int i, j, k, w;

    char* rawdatafilename = (char*)malloc(sizeof(char) * 200);//edited slightly if everything blows up definitely this (11/25/2019)
    strcpy(rawdatafilename, "rawdatafor_"); //starting the string that will be the name of the data file.

    strcat(rawdatafilename, "Sb_"); //for adding values of the beneficial mutation effect size to the data name.
    strcat(rawdatafilename, Sbname);
    strcat(rawdatafilename, "mub_");
    strcat(rawdatafilename, mubname);
    strcat(rawdatafilename, ".txt");

    char* summarydatafilename = (char*)malloc(100);
    strcpy(summarydatafilename, "summarydatafor_");
    strcat(summarydatafilename, "Sb_");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub_");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    
    int popsize;
    double t;
    int *pPopSize = &popsize;
    double *pCurrenttime = &t;
    
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = maxPopSize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
    long double sumofdeathrates;
    long double sumofdeathratessquared;
    long double *psumofdeathrates;
    psumofdeathrates = &sumofdeathrates;
    long double *psumofdeathratessquared;
    psumofdeathratessquared = &sumofdeathratessquared;
    long double *wholepopulationselectiontree;
    wholepopulationselectiontree = malloc(sizeof(long double) * maxPopSize);
       
    long double *wholepopulationdeathratesarray;
    wholepopulationdeathratesarray = malloc(sizeof(long double) * maxPopSize);
    //The Fenwick tree does not store each individual's wi, but rather a collection of partial sums.
    //For debugging purposes and data that requires summarizing wis, storing the wis in an array is necessary.
    bool *wholepopulationisfree;
    wholepopulationisfree = malloc(sizeof(bool) * maxPopSize);
    
    int *wholepopulationindex;
    wholepopulationindex = malloc(sizeof(int) * maxPopSize);
    
    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * maxPopSize);
    //In order to visualize the distribution of fitness in the population,
    //the individual fitnesses need to be sorted, which requires a separate array.
    
    //These parameters specify the intervals to print rawdata from the simulations
    int printeach = 100;
    int printtime = 0;
    
    if(!issnapshot){
        rawdatafilepointer = fopen(rawdatafilename, "w"); //opens the file to which to print data to be plotted.
        fprintf(rawdatafilepointer, "Time\tPop_size\tMean_death_rate\tVariance_death_rate\tMean_birth_rate\n");
        
        summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
        popsize = initialPopSize;
        t = 0.0;
        //assignment of data to popArray for index, wis, and deathrate
        InitializePopulationAbs(wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, initialPopSize, maxPopSize, wholepopulationgenomes, totalpopulationgenomelength, psumofdeathrates, psumofdeathratessquared, d_0);//sets up all data within the population for a run. As this initializes data I think it should be a separate function.
    }else{
        rawdatafilepointer = fopen(rawdatafilename, "a");
        summarydatafilepointer = fopen(summarydatafilename, "a"); //opens the file to which to print summary data.
        FILE *prevsnapshotfilepointer;
        prevsnapshotfilepointer = fopen(prevsnapshotfilename, "r");
        if(prevsnapshotfilepointer == NULL){
            fprintf(miscfilepointer, "\nError: There was an error oppening the snapshot file");
            exit(0);
        }
        
        InitializeWithSnapshotAbs(wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationindex,wholepopulationisfree, maxPopSize, wholepopulationgenomes, totalpopulationgenomelength, psumofdeathrates, psumofdeathratessquared, pPopSize, pCurrenttime, prevsnapshotfilepointer, miscfilepointer);
        
        if (VERYVERBOSE == 1) {
            fprintf(veryverbosefilepointer, "Started population with snapshot file.\n");
        }
        //the previous simulation time had to be added to the max time of the simulation
        maxTime += t;
        printtime += t + printeach;
        
        fclose(prevsnapshotfilepointer);
    }
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Entered simulation run.\n");
    }
    
    //variables used to define death rates
    double const b_0 = 1.0;
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Population initialized.\n");
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * maxTime);
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
    int endofdelay = maxTime-1;
    int endofsimulation = maxTime-1;
    int Nxtimestepsafterburnin = 0;
    double arbitrarynumber = (-1.0 * 0.007 / maxPopSize); //using a number somewhere close to the mean of the DFE for deleterious mutations.
    double slopeoflogfitness;    
    double variancesum;
    
    bool birthhappens;
    
    double parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
    
    double birthrate;
    
    double *arrayofbirthrates;
    arrayofbirthrates = malloc(sizeof(double)*maxPopSize);
    
    calcRateofBirths(arrayofbirthrates, maxPopSize, b_0);
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
    }
    
    
    //BEGIN THE SIMULATION FOR LOOP
    
    while (t < maxTime) {
        
        if(popsize < 3){
        		fprintf(summarydatafilepointer, "Population died during run at time %f", t);
                break;
        }
        
        if(popsize >= (maxPopSize-1)){
        		fprintf(summarydatafilepointer, "Population achieved its maximum population size at time %f", t);
                break;
        }
        
        birthhappens = monteCarloStep(arrayofbirthrates, popsize, sumofdeathrates, pCurrenttime);//This is the monte carlo step. This decides if a birth or a death event takes place by returning a 0 or 1
        
        PerformOneEventAbs(isabsolute, birthhappens, maxPopSize, pPopSize, wholepopulationgenomes, wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationisfree, wholepopulationindex, psumofdeathrates,psumofdeathratessquared, d_0, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma, r, sdmin, miscfilepointer);
        
        if(t > printtime){
            birthrate = arrayofbirthrates[popsize];
            fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\n", t, popsize, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize));
            fflush(rawdatafilepointer);
            printtime += printeach;
        }        
        
    }  
    //END OF SIMULATION FOR LOOP
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %f \n", Sb);
        fprintf(veryverbosefilepointer, "Time elapsed: %f\n", t);
        fprintf(veryverbosefilepointer, "Final population size was: %d\n\n", popsize);
    }
    
    
    FILE *popsnapshotfilepointer;
    
    char* popsnapshotfilename = (char*)malloc(100);
    strcpy(popsnapshotfilename, "popsnapshotfor_");
    strcat(popsnapshotfilename, "Sb_");
    strcat(popsnapshotfilename, Sbname);
    strcat(popsnapshotfilename, "mub_");
    strcat(popsnapshotfilename, mubname);
    strcat(popsnapshotfilename, ".txt");
    popsnapshotfilepointer = fopen(popsnapshotfilename, "w"); //opens the file to which to print summary data.
    
    WritePopSnapshot(wholepopulationgenomes, totalpopulationgenomelength, sumofdeathrates, sumofdeathratessquared, wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationisfree, wholepopulationindex, maxPopSize, popsize, t, popsnapshotfilepointer);
    
    fclose(rawdatafilepointer);
    fclose(summarydatafilepointer);
    fclose(popsnapshotfilepointer);
    free(rawdatafilename);
    free(summarydatafilename);
    free(popsnapshotfilename);
    
    free(logaveragefitnesseachNtimesteps);
    free(literallyjustlast200Ntimesteps);
    free(last200Ntimestepsvariance);
    
    free(wholepopulationgenomes);
    free(wholepopulationselectiontree);
    free(wholepopulationdeathratesarray);
    free(sortedwisarray);
    free(wholepopulationisfree);
    free(wholepopulationindex);

//     return slopeoflogfitness;//should return average fitness
    return 0;
       
}


void calcRateofBirths(double *arrayofbirthrates, int maxPopSize, double b_0) {
    int i;
    for(i = 0; i < maxPopSize; i++){
        arrayofbirthrates[i] = (b_0) * (double)i * (1 - ((double)i/(double)maxPopSize));
    }
}

//returns output of either birth or death
bool monteCarloStep(double *arrayofbirthrates, int popsize, double sumofdeathrates, double *pCurrenttime) {

    double deathRate;
    double birthRate;
    double timestep;

    deathRate = sumofdeathrates;

    //rate of births is calculated using equation used in lab write up
    birthRate = arrayofbirthrates[popsize];
    
    timestep = 1/(deathRate + birthRate);//To make time steps dynamical, we directly use the inverse of the sum of the rates as value for a time step.

    *pCurrenttime += timestep;

    //This is the actual monte carlo step
    return discoverEvent(deathRate, birthRate);

}

bool discoverEvent(double deathRate, double birthRate) {

    bool boolBirth;
    
    double cutOffPoint, randomNumber;

    randomNumber = ldexp(pcg32_random(), -32);

    //if lands on the cutoff point death occurs, technically skewing towards death but is very minor.
    cutOffPoint = deathRate/(deathRate + birthRate);

    if (randomNumber > cutOffPoint)
        boolBirth = 1;
    else
        boolBirth = 0;
    
    return boolBirth;

}


bool PerformOneEventAbs(bool isabsolute, bool birthhappens, int maxPopSize, int *pPopSize, double *wholepopulationgenomes, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, long double *psumofdeathrates, long double *psumofdeathratessquared, double d_0, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution,  double* parent1gamete, double* parent2gamete, gsl_rng* randomnumbergeneratorforgamma, double r, double sdmin, FILE *miscfilepointer)
{
    if(isabsolute == 0){
        fprintf(miscfilepointer, "\nError: Trying to use PerformOneEventAbs within a non absolute fitness simulation.");
        exit(0);
    }
    
    int randparent1, randparent2, currentparent1, currentparent2, currentvictim;
    int i;
    bool nolethalmut = true;
        
    int victim;
    
    if(*pPopSize < 2){
        fprintf(miscfilepointer, "\nError: Trying to use PerformOneEventAbs with a population size lower than 2.");
        exit(0);
    }
    
    if(*pPopSize >= maxPopSize && birthhappens){
        fprintf(miscfilepointer, "\nError:  Trying to use PerformOneEventAbs with a population size greater than maxPopSize.");
        exit(0);
    }
    
    //next pointers are only used in relative fitness scenarios. Here, they are just initialized
    long double *psumofloads, *wholepopulationwisarray;
    
    if (birthhappens){
        randparent1 = ChooseParent(*pPopSize);
        currentparent1 = wholepopulationindex[randparent1];
        do{
            randparent2 = ChooseParent(*pPopSize);
        } while(randparent2 == randparent1);
        currentparent2 = wholepopulationindex[randparent2];
        
        if(wholepopulationisfree[currentparent1] || wholepopulationisfree[currentparent2]){
            fprintf(miscfilepointer, "\nError:  Trying to give birth to a new individual using a parent that is not present in PerformOneEventAbs.");
            exit(0);
        }
        

        RecombineChromosomesIntoGamete(currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
        nolethalmut = ProduceMutatedGamete(isabsolute, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent1gamete, randomnumbergeneratorforgamma, miscfilepointer);
        if(!nolethalmut)
            return false;
        
        RecombineChromosomesIntoGamete(currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
        nolethalmut = ProduceMutatedGamete(isabsolute, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);
        if(!nolethalmut)
            return false;
        
        PerformBirth(isabsolute, parent1gamete, parent2gamete, maxPopSize, pPopSize, victim, wholepopulationgenomes, totalindividualgenomelength, wholepopulationselectiontree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, psumofdeathrates, psumofdeathratessquared, d_0, r, sdmin, miscfilepointer);
    }
    
    
    
    else{
//         use of the fennwick tree to find the victim in the population, note that since fenwick tree stores all the population (non occupied spaces have a fitness of 0) there is not need to use the wholepopulationindex, fenwick search already gives you the space that the selected individual occupies.
    	victim = ChooseVictimWithTree(wholepopulationselectiontree, *pPopSize, maxPopSize, *psumofdeathrates, miscfilepointer);
        

        PerformDeath(isabsolute, maxPopSize, pPopSize, victim, wholepopulationselectiontree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, psumofdeathrates, psumofdeathratessquared, miscfilepointer);
    }
    
    return true;

}


void InitializePopulationAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, double d_0)
{    
    int i, j;

    for (i = 0; i < initialPopSize; i++) {
        wholepopulationselectiontree[i] = d_0; //all individuals start with death rate d_0. Currently d_0 is initialized in global_vars.h might be a better idea to initialize this variable in the OneStep function
        
        wholepopulationdeathratesarray[i] = d_0;
        
        wholepopulationindex[i] = i;
        
        wholepopulationisfree[i] = false;
    }
    
    for (i = initialPopSize; i < maxPopSize; i++) {
        wholepopulationselectiontree[i] = 0.0;
        
        wholepopulationdeathratesarray[i] = 0.0;
                                
        wholepopulationindex[i] = i;
        
        wholepopulationisfree[i] = true;
    }
    //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page. Absolute trees needs to contemplate all popSize, including non occupied spaces. 
    for (i = 0; i < maxPopSize; i++) {
        j = i + LSB(i+1);
        if (j < maxPopSize) {
            wholepopulationselectiontree[j] += wholepopulationselectiontree[i];
        }
    }
    
    *psumofdeathrates = (long double) initialPopSize * d_0;
    
    *psumofdeathratessquared = (long double) initialPopSize * pow(d_0, 2);
    
    for (i = 0; i < totalpopulationgenomelength; i++){
        wholepopulationgenomes[i];
    }
    
}

//used to initializa a population using a file that stores all the important variables
void InitializeWithSnapshotAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, int *pPopSize, double *pCurrenttime, FILE *prevsnapshotfilepointer, FILE *miscfilepointer)
{
    int i, j;

    char strtemp[200];
    
    //gets the genomes of the whole population from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    if(strcmp(strtemp, "Genomes:")){
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers does not incluide Genomes: check that the snapshot file is in the correct parameters folder");
        exit(0);
    }
    for (i = 0; i < totalpopulationgenomelength; i++){
        fscanf(prevsnapshotfilepointer, "%lf", &wholepopulationgenomes[i]);
    }
    
    //gets the sum of death rates from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    if(strcmp(strtemp, "Sum_of_death_rates:")){
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers does not incluide Sum_of_death_rates: check that the snapshot file is in the correct parameters folder or that there are no errors on the file");
        exit(0);
    }
    fscanf(prevsnapshotfilepointer, "%Lf", psumofdeathrates);
    
    //gets the sum of death rates squared from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    fscanf(prevsnapshotfilepointer, "%Lf", psumofdeathratessquared);
    
    //gets the selection tree from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    for (i = 0; i < maxPopSize; i++){
        fscanf(prevsnapshotfilepointer, "%Lf", &wholepopulationselectiontree[i]);
    }
    
    //gets the selection tree from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    for (i = 0; i < maxPopSize; i++){
        fscanf(prevsnapshotfilepointer, "%Lf", &wholepopulationdeathratesarray[i]);
    }
    
    //since fscanf does not read bools, first store the value in an int
    int tempspace;
    //gets the free spaces in the population from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    for (i = 0; i < maxPopSize; i++){
        fscanf(prevsnapshotfilepointer, "%d", &tempspace);
        wholepopulationisfree[i] = tempspace;
    }
    
    //gets the index array from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    for (i = 0; i < maxPopSize; i++){
        fscanf(prevsnapshotfilepointer, "%d", &wholepopulationindex[i]);
    }
    
    //gets the popsize from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    fscanf(prevsnapshotfilepointer, "%d", pPopSize);
    
    //gets the time from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    if(strcmp(strtemp, "Time:")){
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers does not incluide Time: check that the snapshot file is in the correct parameters folder or that there are no errors on the file");
        exit(0);
    }
    fscanf(prevsnapshotfilepointer, "%lf", pCurrenttime);
   
}

int ChooseParent(int populationsize)
{
    return pcg32_boundedrand(populationsize);
}

int ChooseVictimWithTree(long double* wholepopulationselectiontree, int popsize, int maxPopSize, long double inversesumofloads, FILE *miscfilepointer)
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
    
    
    newVictim = SearchTree(leftbound, rightbound, randomnumberofdeath, wholepopulationselectiontree);
    
    return newVictim;
}

int findinindex(int *wholepopulationindex, int which, int tam, FILE *miscfilepointer){
    int i;
    bool placefound = false;
    
    for(i = 0; i < tam; i++){
        if(wholepopulationindex[i] == which){
            placefound = true;
            return i;
        }
    }
    
    if(!placefound){
        fprintf(miscfilepointer, "\nError: Array of indexes is corrupted in findinindex. \n");
        exit(0);
    }
}

void indexArrayFlipDeath(int *wholepopulationindex, int placeinindex, int popsize){

    /*
    This function flips the integer value from the last in the array of  wholepopulationindex to the place where an individual died. This way wholepopulationindex is mantained without unoccupied spaces
    */
	int indexlast = wholepopulationindex[(popsize-1)];
    int indexvictim = wholepopulationindex[placeinindex];
    
    wholepopulationindex[placeinindex] = indexlast;
    wholepopulationindex[(popsize-1)] = indexvictim;
}

double CalculateDeathRate(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, double d_0, double r, double sdmin)
{
    double inddeathrate = 0.0;
    double currentlinkageblocksload = 0.0;
    int i;

    for (i = 0; i < (totalindividualgenomelength/2); i++) {
        currentlinkageblocksload += parent1gamete[i];
        currentlinkageblocksload += parent2gamete[i];
    }
    inddeathrate = d_0 + sdmin * (1 - pow(r, -currentlinkageblocksload))/log(r);
    return inddeathrate;
}

void WritePopSnapshot(double *wholepopulationgenomes, int totalpopulationgenomelength, long double sumofdeathrates, long double sumofdeathratessquared, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, int maxPopSize, int popsize, double t, FILE *popsnapshotfilepointer)
{
    int i;
    
    fprintf(popsnapshotfilepointer, "Genomes:\n");
    for (i = 0; i < totalpopulationgenomelength; i++)
        fprintf(popsnapshotfilepointer, "%lf\t", wholepopulationgenomes[i]);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Sum_of_death_rates:\n");
    fprintf(popsnapshotfilepointer, "%Lf", sumofdeathrates);
    fprintf(popsnapshotfilepointer, "\n");
    fprintf(popsnapshotfilepointer, "Sum_of_death_rates_squared:\n");
    fprintf(popsnapshotfilepointer, "%Lf", sumofdeathratessquared);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Selection_tree:\n");
    for (i = 0; i < maxPopSize; i++)
        fprintf(popsnapshotfilepointer, "%Lf\t", wholepopulationselectiontree[i]);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Deaths_array:\n");
    for (i = 0; i < maxPopSize; i++)
        fprintf(popsnapshotfilepointer, "%Lf\t", wholepopulationdeathratesarray[i]);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Pop_free_spaces:\n");
    for (i = 0; i < maxPopSize; i++)
        fprintf(popsnapshotfilepointer, "%i\t", wholepopulationisfree[i]);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Pop_index:\n");
    for (i = 0; i < maxPopSize; i++)
        fprintf(popsnapshotfilepointer, "%i\t", wholepopulationindex[i]);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Popsize:\n");
    fprintf(popsnapshotfilepointer, "%i", popsize);
    fprintf(popsnapshotfilepointer, "\n");
    
    fprintf(popsnapshotfilepointer, "Time:\n");
    fprintf(popsnapshotfilepointer, "%lf", t);
    fprintf(popsnapshotfilepointer, "\n");
}
