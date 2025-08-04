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

double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize)
{
    if(isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationRel within an absolute fitness program. \n");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    //tree sequence data files
    FILE *nodefilepointer;
    FILE *edgefilepointer;
    FILE *sitefilepointer;
    FILE *mutationfilepointer;
    
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
    fprintf(rawdatafilepointer, "Nxtimesteps,Sum.of.wis,Variance.in.log.fitness,FractionSelectiveDeaths,FractionSelectiveDeaths_exponantiated\n");
    
    char * summarydatafilename = (char *) malloc(100);
    strcpy(summarydatafilename, "summarydatafor");
    strcat(summarydatafilename, "Sb");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "mub");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
    
    nodefilepointer = fopen("nodetable.txt", "w");
    edgefilepointer = fopen("edgetable.txt", "w");
    sitefilepointer = fopen("sitetable.txt", "w");
    mutationfilepointer = fopen("mutationtable.txt", "w");
    
    int totaltimesteps = Nxtimesteps * popsize;
    double currenttimestep = 0.0;
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
    
    //Following lines are from the tskit library.
    //Initializes the tables that make up the tree sequence recording.
    tsk_table_collection_t treesequencetablecollection;
    tsk_table_collection_t * tablepointer = &treesequencetablecollection;
    int returnvaluefortskfunctions = tsk_table_collection_init(&treesequencetablecollection, 0);
    check_tsk_error(returnvaluefortskfunctions);
    
    tsk_id_t *wholepopulationnodesarray;
    wholepopulationnodesarray = malloc(sizeof(tsk_id_t) * 2 * popsize);
    //The extant nodes need to have explicit identification in order to add edges between parents and children nodes.
    //Each node is only a single set of chromosomes, so the 2 here assumes diploidy.
    
    tsk_id_t wholepopulationsitesarray[totalindividualgenomelength / 2];
    //The number of sites is the number of linkage blocks in a single set of chromosomes (haploid), and won't change over the course of the simulation.
       
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
    
    InitializePopulationRel(tskitstatus, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, wholepopulationwistree, wholepopulationwisarray, popsize, wholepopulationgenomes, totalpopulationgenomelength, totaltimesteps, psumofwis);
    
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
    double varianceinlogfitness;   
    long double fitnessfittest;
    long double FractionSelectiveDeaths;
    long double FractionSelectiveDeaths_exponantiatebirthrates;
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
        fflush(veryverbosefilepointer);
    }
    
    //BEGIN THE SIMULATION FOR LOOP
    
    for (i = 0; i < Nxtimesteps; i++) {
        
        //Following code performs N rounds of paired births and deaths.
        for (j = 0; j < popsize; j++) {
            currenttimestep += 1.0;            
            PerformOneTimeStepRel(tskitstatus, isabsolute, isburninphaseover, ismodular, elementsperlb, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, popsize, totaltimesteps, currenttimestep,wholepopulationwistree, wholepopulationwisarray, wholepopulationgenomes, psumofwis, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gamete, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);  
        }
        
        //Following code calculates the variance in log(fitness) of the population after this generation of births and deaths.
        //May use an imprecise algorithm -- check before using as data.
        varianceinlogfitness = CalculateVarianceInLogFitness(popsize, wholepopulationwisarray, *psumofwis);
        fitnessfittest = FindFittestWi(wholepopulationwisarray, popsize);
        FractionSelectiveDeaths = (fitnessfittest-(sumofwis/popsize))/fitnessfittest;
	FractionSelectiveDeaths_exponantiatebirthrates = (exp(fitnessfittest)-exp((sumofwis/popsize)))/exp(fitnessfittest);
        
        //This is the main data output, currently the summed fitness and variance in log(fitness) in the population.
        fprintf(rawdatafilepointer, "%d,%Lf,%.18f,%Lf,%Lf\n", i+1, *psumofwis, varianceinlogfitness, FractionSelectiveDeaths, FractionSelectiveDeaths_exponantiatebirthrates);
        fflush(rawdatafilepointer);

        //fprintf(rawdatafilepointer, "%d \n", i+1);
        //fflush(rawdatafilepointer);

        //Following code periodically simplifies the tree sequence tables, which requires the tables to be sorted each time.
        //The tskit API calls out sorting each time as inefficient, but they haven't yet uploaded an example of how to do it differently. The API currently suggests that they would use the C++ API, which I don't want to deal with.
        //Simplify interval is currently set to 5. Should be either an input parameter or global variable at some point, probably.
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
   //I don't guarantee that the simplify interval matches up with the number of Nxtimesteps
    //Tree sequence recording requires that tables are sorted on the back end,
    //so I sort once again here at the end to ensure that all tables are sorted before they're read to file.
    //This might be inefficient, I'm not sure.
    if(tskitstatus > 0){
        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
        check_tsk_error(returnvaluefortskfunctions);

        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*popsize), 0, NULL);
        check_tsk_error(returnvaluefortskfunctions);
    
        for (k = 0; k < (2*popsize); k++) {
            wholepopulationnodesarray[k] = k;
        }
    
    //Printing out the node table in a way readable by python on the back end.
        fprintf(nodefilepointer, "is_sample time\n");
        for (k = 0; k < tablepointer->nodes.num_rows; k++) {
            if (k < (2*popsize)) {
                fprintf(nodefilepointer, "1 %f\n", tablepointer->nodes.time[k]);
            } else {
                fprintf(nodefilepointer, "0 %f\n", tablepointer->nodes.time[k]);
            }
        }
    
    //Printing out the edge table in a way readable by python on the back end.
        fprintf(edgefilepointer, "left right parent child\n");
        for (k = 0; k < tablepointer->edges.num_rows; k++) {
            fprintf(edgefilepointer, "%f %f %d %d\n", tablepointer->edges.left[k], tablepointer->edges.right[k], tablepointer->edges.parent[k], tablepointer->edges.child[k]);
        }
    
    //Printing out the site table in a way readable by python.
        fprintf(sitefilepointer, "position ancestral_state\n");
        for (k = 0; k < tablepointer->sites.num_rows; k++) {
            fprintf(sitefilepointer, "%f 0.0\n", tablepointer->sites.position[k]);
        }
    
    //Printing out the mutation table in a way readable by python.
        fprintf(mutationfilepointer, "site node derived_state\n");
        for (k = 0; k < tablepointer->mutations.num_rows; k++) {
            fprintf(mutationfilepointer, "%d %d %.12s\n", tablepointer->mutations.site[k], tablepointer->mutations.node[k], (tablepointer->mutations.derived_state + k*12));
        }
    }

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
        fclose(nodefilepointer);
        fclose(edgefilepointer);
        fclose(sitefilepointer);
        fclose(mutationfilepointer);

        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        
        tsk_table_collection_free(&treesequencetablecollection);
        
        return slopeoflogfitness;
    }

    if (isburninphaseover == 0) {
        fprintf(summarydatafilepointer, "End of burn-in phase not reached.");
        
        fclose(rawdatafilepointer); 
        fclose(summarydatafilepointer);
        fclose(nodefilepointer);
        fclose(edgefilepointer);
        fclose(sitefilepointer);
        fclose(mutationfilepointer);

        free(rawdatafilename);
        free(summarydatafilename);
        free(logaveragefitnesseachNtimesteps);
        free(literallyjustlast200Ntimesteps);
        free(last200Ntimestepsvariance);
        free(wholepopulationgenomes);
        free(wholepopulationwistree);
        free(wholepopulationwisarray);
        free(wholepopulationnodesarray);
        free(sortedwisarray);
        
        tsk_table_collection_free(&treesequencetablecollection);
        
        return -1.0;
    }
}

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int currentparent1, currentparent2, currentvictim;

    currentvictim = ChooseVictim(popsize);
    currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
        currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, *psumofwis, miscfilepointer);
    }
    
    tsk_id_t childnode1, childnode2;
   
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode1, totaltimesteps, currenttimestep, currentparent1, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gamete, randomnumbergeneratorforgamma, miscfilepointer);
        
    RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
    ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode2, totaltimesteps, currenttimestep, currentparent2, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);
               
    //next variables are only used for absolute simulations, however since InitializePopulation is a shared function I must initialize them here, even if I don't use them again
    long double *pInverseSumOfWis;
    long double *pInverseSumOfWissquared;
    long double *psumofload;
    long double *psumofloadsquared;
    bool *wholepopulationisfree;
    int *wholepopulationindex;
    long double *wholepopulationdeathratesarray;
    int *pPopSize;
	double b_0, r, s;
	int i_init;
    
    PerformDeath(isabsolute, tskitstatus, isburninphaseover, popsize, pPopSize, currentvictim, deleteriousdistribution, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, pInverseSumOfWissquared, b_0, r, i_init, s, psumofload, psumofloadsquared, wholepopulationnodesarray, miscfilepointer);
    
    PerformBirth(tskitstatus, isburninphaseover, ismodular, elementsperlb, treesequencetablecollection, wholepopulationnodesarray, childnode1, childnode2, isabsolute, parent1gamete, parent2gamete, popsize, pPopSize, currentvictim, wholepopulationgenomes, totalindividualgenomelength, deleteriousdistribution, wholepopulationwistree, wholepopulationwisarray, wholepopulationdeathratesarray,wholepopulationindex, wholepopulationisfree, psumofwis, pInverseSumOfWis, pInverseSumOfWissquared, b_0, r, i_init, s, psumofload, psumofloadsquared, miscfilepointer);
    
}

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, long double *wholepopulationwisarray, int popsize, double *wholepopulationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis) 
{
    int i, j;
    
    double haploidgenomelength = (double) ((totalpopulationgenomelength / popsize) / 2);
    
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
    //The following lines initialize the node table for tree sequence recording.
    //Note that nodes here are single sets of chromosomes, so the 2x popsize here assumes diploidy.
    if (tskitstatus > 0){
	treesequencetablecollection->sequence_length = haploidgenomelength;

        for (i = 0; i < (2 * popsize); i++) {
            wholepopulationnodesarray[i] = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, totaltimesteps, TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(wholepopulationnodesarray[i]);
        }
    
    //The following lines add a site to the tree sequence recording site table corresponding to each linkage block, with ancestral state of 0.
        for (i = 0; i < haploidgenomelength; i++) {
            wholepopulationsitesarray[i] = tsk_site_table_add_row(&treesequencetablecollection->sites, i, "0.00000000", 10, NULL, 0);
            check_tsk_error(wholepopulationsitesarray[i]);
        }
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
