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
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, double mutator_strength_factor, double mutator_switch_rate, double mutator_bias)
{
    if(isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationRel within an absolute fitness program. \n");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    FILE *nodefilepointer;
    FILE *edgefilepointer;
    FILE *sitefilepointer;
    FILE *mutationfilepointer;
    
    int i, j, k;
    
    char * rawdatafilename = (char *) malloc(200);
    strcpy(rawdatafilename, "rawdatafor");
    strcat(rawdatafilename, "Nxtimesteps"); strcat(rawdatafilename, Nxtimestepsname);
    strcat(rawdatafilename, "popsize"); strcat(rawdatafilename, popsizename);
    strcat(rawdatafilename, "mutrate"); strcat(rawdatafilename, delmutratename);
    strcat(rawdatafilename, "chromsize"); strcat(rawdatafilename, chromsizename);
    strcat(rawdatafilename, "chromnum"); strcat(rawdatafilename, chromnumname);
    strcat(rawdatafilename, "benmutrate"); strcat(rawdatafilename, mubname);
    strcat(rawdatafilename, "Sb"); strcat(rawdatafilename, Sbname);
    strcat(rawdatafilename, ".txt");

    rawdatafilepointer = fopen(rawdatafilename, "w");
    fprintf(rawdatafilepointer, "Nxtimesteps,Sum.of.wis,Variance.in.log.fitness,FractionSelectiveDeaths,FractionSelectiveDeaths_exponantiated\n");
    
    char * summarydatafilename = (char *) malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb"); strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub"); strcat(summarydatafilename, mubname);
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
    
    wholepopulation = malloc(sizeof(Individual) * popsize);
    // Note: Initialization of internal arrays happens in InitializePopulationRel
    
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

    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * popsize);

    InitializePopulationRel(tskitstatus, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, wholepopulationwistree, wholepopulation, popsize, totalpopulationgenomelength, totaltimesteps, psumofwis);
    
    // Set initial mutation rates for population based on 0 mutator load
    for(k = 0; k < popsize; k++) {
        UpdateIndividual(&wholepopulation[k], totalindividualgenomelength, mutator_strength_factor, deleteriousmutationrate);
    }
    
    double *logaveragefitnesseachNtimesteps;
    logaveragefitnesseachNtimesteps = malloc(sizeof(double) * Nxtimesteps);
    
    // Parent gamete arrays - allocated once here
    double parent1gameteFitness[numberofchromosomes*chromosomesize], parent2gameteFitness[numberofchromosomes*chromosomesize];
    int parent1gameteMutators[numberofchromosomes*chromosomesize], parent2gameteMutators[numberofchromosomes*chromosomesize];
    
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
    
    for (i = 0; i < Nxtimesteps; i++) {
        for (j = 0; j < popsize; j++) {
            currenttimestep += 1.0;            
            PerformOneTimeStepRel(tskitstatus, isabsolute, isburninphaseover, ismodular, elementsperlb, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, totaltimesteps, currenttimestep, wholepopulationwistree, wholepopulation, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gameteFitness, parent1gameteMutators, parent2gameteFitness, parent2gameteMutators, randomnumbergeneratorforgamma, miscfilepointer, mutator_strength_factor, mutator_switch_rate, mutator_bias);  
        }
        
        varianceinlogfitness = CalculateVarianceInLogFitness(popsize, wholepopulation, *psumofwis);
        fitnessfittest = FindFittestWi(wholepopulation, popsize);
        FractionSelectiveDeaths = (fitnessfittest-(sumofwis/popsize))/fitnessfittest;
        FractionSelectiveDeaths_exponantiatebirthrates = (exp(fitnessfittest)-exp((sumofwis/popsize)))/exp(fitnessfittest);
        
        fprintf(rawdatafilepointer, "%d,%Lf,%.18f,%Lf,%Lf\n", i+1, *psumofwis, varianceinlogfitness, FractionSelectiveDeaths, FractionSelectiveDeaths_exponantiatebirthrates);
        fflush(rawdatafilepointer);

        if (tskitstatus > 0){
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
    
    if(tskitstatus > 0){
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
        free(rawdatafilename); free(summarydatafilename); free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps); free(last200Ntimestepsvariance);
        
        for(k=0; k<popsize; k++){
            free(wholepopulation[k].fitnessArray);
            free(wholepopulation[k].mutatorArray);
        }
        free(wholepopulation);
        free(wholepopulationwistree);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        tsk_table_collection_free(&treesequencetablecollection);
        return slopeoflogfitness;
    }

    if (isburninphaseover == 0) {
        fprintf(summarydatafilepointer, "End of burn-in phase not reached.");
        fclose(rawdatafilepointer); fclose(summarydatafilepointer);
        fclose(nodefilepointer); fclose(edgefilepointer); fclose(sitefilepointer); fclose(mutationfilepointer);
        free(rawdatafilename); free(summarydatafilename); free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps); free(last200Ntimestepsvariance);
        for(k=0; k<popsize; k++){
            free(wholepopulation[k].fitnessArray);
            free(wholepopulation[k].mutatorArray);
        }
        free(wholepopulation);
        free(wholepopulationwistree);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        tsk_table_collection_free(&treesequencetablecollection);
        return -1.0;
    }
}

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, Individual *wholepopulation, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gameteFitness, int *parent1gameteMutators, double *parent2gameteFitness, int *parent2gameteMutators, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, double mutator_strength_factor, double mutator_switch_rate, double mutator_bias)
{
    int currentparent1, currentparent2, currentvictim;
    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    while (currentparent1 == currentparent2) {
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    }
    
    tsk_id_t childnode1, childnode2;
   
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, chromosomesize, numberofchromosomes, parent1gameteFitness, parent1gameteMutators, wholepopulation, totalindividualgenomelength);
    
    double p1_rate = wholepopulation[currentparent1].mutationRate; 
    
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, isabsolute, totalindividualgenomelength, p1_rate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gameteFitness, parent1gameteMutators, mutator_switch_rate, mutator_bias, randomnumbergeneratorforgamma, miscfilepointer);
        
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, chromosomesize, numberofchromosomes, parent2gameteFitness, parent2gameteMutators, wholepopulation, totalindividualgenomelength);
    
    double p2_rate = wholepopulation[currentparent2].mutationRate;
    
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, isabsolute, totalindividualgenomelength, p2_rate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent2gameteFitness, parent2gameteMutators, mutator_switch_rate, mutator_bias, randomnumbergeneratorforgamma, miscfilepointer);
               
    int *pPopSize; 
    
    PerformDeath(isabsolute, tskitstatus, isburninphaseover, popsize, pPopSize, currentvictim, deleteriousdistribution, wholepopulationwistree, wholepopulation, NULL, NULL, NULL, psumofwis, NULL, NULL, 0, 0, 0, 0, NULL, NULL, wholepopulationnodesarray, miscfilepointer);
    
    PerformBirth(tskitstatus, isburninphaseover, ismodular, elementsperlb, treesequencetablecollection, wholepopulationnodesarray, childnode1, childnode2, isabsolute, parent1gameteFitness, parent1gameteMutators, parent2gameteFitness, parent2gameteMutators, popsize, pPopSize, currentvictim, wholepopulation, totalindividualgenomelength, deleteriousdistribution, wholepopulationwistree, NULL, NULL, NULL, psumofwis, NULL, NULL, 0, 0, 0, 0, NULL, NULL, miscfilepointer, mutator_strength_factor, deleteriousmutationrate);
}

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, Individual *wholepopulation, int popsize, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis) 
{
    int i, j;
    double haploidgenomelength = (double) ((totalpopulationgenomelength / popsize) / 2);
    int genomelength = (totalpopulationgenomelength / popsize);

    for (i = 0; i < popsize; i++){
        wholepopulation[i].fitnessArray = malloc(sizeof(double) * genomelength);
        wholepopulation[i].mutatorArray = malloc(sizeof(int) * genomelength);
        wholepopulation[i].fitness = 1.0;
        wholepopulation[i].mutationRate = 0.0;
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
    
    if (tskitstatus > 0){
        treesequencetablecollection->sequence_length = haploidgenomelength;
        for (i = 0; i < (2 * popsize); i++) {
            wholepopulationnodesarray[i] = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, totaltimesteps, TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(wholepopulationnodesarray[i]);
        }
        for (i = 0; i < haploidgenomelength; i++) {
            wholepopulationsitesarray[i] = tsk_site_table_add_row(&treesequencetablecollection->sites, i, "0.00000000", 10, NULL, 0);
            check_tsk_error(wholepopulationsitesarray[i]);
        }
    }
}

// ChooseVictim, ChooseParentWithTree, CalculateWi remain identical to original
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