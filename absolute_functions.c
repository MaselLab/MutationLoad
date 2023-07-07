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
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

double RunSimulationAbs(bool issnapshot, char *prevsnapshotfilename, bool isredinmaxpopsize, char *redinmaxpopsizename, double redinmaxpopsize, char* mubname, char* Sbname, int tskitstatus, bool ismodular, int elementsperlb, bool isabsolute, int maxTime, int initialPopSize, int K, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double Sd, int deleteriousdistribution, double beneficialmutationrate, double Sb, int beneficialdistribution, double r, int i_init, double s, gsl_rng* randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, bool iscalcfixation)
{

    if(!isabsolute){
        fprintf(miscfilepointer, "\n Trying to use RunSimulationAbs within a relative fitness program. \n");
        exit(0);
    }
    
    FILE *rawdatafilepointer;
    FILE *summarydatafilepointer;
    //tree sequence data files
    FILE *nodefilepointer;
    FILE *edgefilepointer;
    FILE *sitefilepointer;
    FILE *mutationfilepointer;
    FILE *tskitdatafilepointer;
    int i, j, k, w;

    char* rawdatafilename =  MakeRawDataFileName(mubname, Sbname, isredinmaxpopsize, redinmaxpopsizename);

    char* summarydatafilename = MakeSummaryDataFileName(mubname, Sbname);
    
    summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
    
    nodefilepointer = fopen("nodetable.txt", "w");
    edgefilepointer = fopen("edgetable.txt", "w");
    sitefilepointer = fopen("sitetable.txt", "w");
    mutationfilepointer = fopen("mutationtable.txt", "w");
    
    if (iscalcfixation)
    {
        tskitdatafilepointer = fopen("tskitdata.txt", "w"); //opens the file to which to print fixed mutaiton data.
        fprintf(tskitdatafilepointer, "Mutation\tTime_app\tTime_fix\n");
    }
    
    int kappa;
    if(r == 1.0){
        kappa = K/(s*i_init);
    }
    else{
        if(deleteriousdistribution == 0){
            kappa = (1-r)*K/(s);
        }
        else if(deleteriousdistribution == 1){
            kappa = -(log(r)*K)/(s);
        }
    }
    int maxPopSize = kappa;
	//variables used to define birth rate
	double const b_0 = 1.0;
	double birthrate;

    double *arrayofbirthrates;
	arrayofbirthrates = malloc(sizeof(double)*maxPopSize);

    calcRateofBirths(arrayofbirthrates, maxPopSize, kappa, b_0);

    int simplifyat = 0;
    int genafterburnin = 50000;
    int simplifyeach = 2000;
    int printeach = maxTime/rawdatafilesize;
    int printtime = 0;
    double redtime = 55000;
    double redtimeeach = redinmaxpopsize;
    if (tskitstatus == 2)
    {
        printtime = genafterburnin;
        if (iscalcfixation)
        {
            simplifyat = genafterburnin + 3000;
            simplifyeach = 3000;
        }else{
            simplifyat = genafterburnin;
        }
    }

    int updatesumofdeathrateseach = 1000;
    int updatesumtime = 0;
    
    int totaltimesteps = maxTime;
    int popsize = initialPopSize;
    int *pPopSize;
    pPopSize = &popsize;
    
    double t = 0.0;
    double *pCurrenttime = &t;
    
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    double *parent1gamete, *parent2gamete;

    if(!ismodular){
        totalpopulationgenomelength = maxPopSize * numberofchromosomes * 2 * chromosomesize;
        totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
        wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
        parent1gamete = malloc(sizeof(double) * numberofchromosomes*chromosomesize); 
        parent2gamete = malloc(sizeof(double) * numberofchromosomes*chromosomesize); 
    }
    else{
        //in modular runs the length of chromosome increases by a factor of elementsperlb (elements per linkage block)
        totalpopulationgenomelength = maxPopSize * numberofchromosomes * 2 * chromosomesize * elementsperlb;
        totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize * elementsperlb;
        wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
        parent1gamete = malloc(sizeof(double) * numberofchromosomes*chromosomesize*elementsperlb); 
        parent2gamete = malloc(sizeof(double) * numberofchromosomes*chromosomesize*elementsperlb); 
    }

    long double sumofdeathrates;
    long double sumofdeathratessquared;
    long double sumofload;
    long double sumofloadsquared;
    long double *psumofdeathrates;
    psumofdeathrates = &sumofdeathrates;
    long double *psumofdeathratessquared;
    psumofdeathratessquared = &sumofdeathratessquared;
    long double *psumofload;
    psumofload = &sumofload;
    long double *psumofloadsquared;
    psumofloadsquared = &sumofloadsquared;
    long double *wholepopulationselectiontree;
    wholepopulationselectiontree = malloc(sizeof(long double) * maxPopSize);
    
    //Following lines are from the tskit library.
    //Initializes the tables that make up the tree sequence recording.
    tsk_table_collection_t treesequencetablecollection;
    tsk_table_collection_t * tablepointer = &treesequencetablecollection;
    int returnvaluefortskfunctions = tsk_table_collection_init(&treesequencetablecollection, 0);
    check_tsk_error(returnvaluefortskfunctions);
    
    tsk_id_t *wholepopulationnodesarray;
    if (tskitstatus != 0){
        wholepopulationnodesarray = malloc(sizeof(tsk_id_t) * 2 * maxPopSize);
    }
    //The extant nodes need to have explicit identification in order to add edges between parents and children nodes.
    //Each node is only a single set of chromosomes, so the 2 here assumes diploidy.
    
    tsk_id_t wholepopulationsitesarray[totalindividualgenomelength / 2];
    //The number of sites is the number of linkage blocks in a single set of chromosomes (haploid), and won't change over the course of the simulation.
    
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

    //assignment of data to popArray for index, wis, and deathrate
    if(!issnapshot){
        rawdatafilepointer = fopen(rawdatafilename, "w"); //opens the file to which to print data to be plotted.
        fprintf(rawdatafilepointer, "Time\tPop_size\tMean_death_rate\tVar_death_rate\tMean_birth_rate\tMean_load\tVar_load\n");
        
        summarydatafilepointer = fopen(summarydatafilename, "w"); //opens the file to which to print summary data.
        popsize = initialPopSize;
        t = 0.0;
        //assignment of data to popArray for index, wis, and deathrate
        InitializePopulationAbs(tskitstatus, &treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, initialPopSize, maxPopSize, totaltimesteps, deleteriousdistribution, wholepopulationgenomes, totalpopulationgenomelength, psumofdeathrates, psumofdeathratessquared, psumofload, psumofloadsquared, b_0, r, i_init, s);//sets up all data within the population for a run. As this initializes data I think it should be a separate function.
    }else{
        rawdatafilepointer = fopen(rawdatafilename, "a");
        summarydatafilepointer = fopen(summarydatafilename, "a"); //opens the file to which to print summary data.
        FILE *prevsnapshotfilepointer;
        prevsnapshotfilepointer = fopen(prevsnapshotfilename, "r");
        if(prevsnapshotfilepointer == NULL){
            fprintf(miscfilepointer, "\nError: There was an error oppening the snapshot file");
            exit(0);
        }
        
        InitializeWithSnapshotAbs(wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, maxPopSize, wholepopulationgenomes, totalpopulationgenomelength, psumofdeathrates, psumofdeathratessquared, psumofload, psumofloadsquared, pPopSize, pCurrenttime, prevsnapshotfilepointer, miscfilepointer);
        
        if (VERYVERBOSE == 1) {
            fprintf(veryverbosefilepointer, "Started population with snapshot file.\n");
        }
        //the previous simulation time had to be added to the max time of the simulation
        maxTime += t;
        printtime += t + printeach;
        updatesumtime += t + updatesumofdeathrateseach;
	    redtime += t + redinmaxpopsize;

        fclose(prevsnapshotfilepointer);
    }

    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Population initialized.\n");
	    fflush(veryverbosefilepointer);
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
    if (tskitstatus == 1){
        isburninphaseover = 1;
    }
    
    int didpopulationcrash = 0;
    int endofburninphase;
    int endofdelay = maxTime-1;
    int endofsimulation = maxTime-1;
    int Nxtimestepsafterburnin = 0;
    double arbitrarynumber = (-1.0 * 0.007 / maxPopSize); //using a number somewhere close to the mean of the DFE for deleterious mutations.
    double slopeoflogfitness;    
    double variancesum;
    bool birthhappens;
    
    //variables used for limiting printing of popsize based on lines
    int N = 0;
    int linesize = 200;
	int index = 0;
	int upperlimit = 0;
    int lowerlimit = maxPopSize;

    int *arrayofpopsizes;
	arrayofpopsizes = malloc(sizeof(int)*maxPopSize);

    for(i = 0; i < maxPopSize; i++){
        arrayofpopsizes[i] = N;
        N += linesize;
    }
    
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Variables initialized, preparing to begin simulation.\n");
	    fflush(veryverbosefilepointer);
    }
    //BEGIN THE SIMULATION FOR LOOP
    
    while (t < maxTime) {
        if(popsize < 3){
        	fprintf(summarydatafilepointer, "Population died during run at time %f", t);
            birthrate = arrayofbirthrates[popsize];
			fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\n", t, popsize, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize), (sumofload/(double)popsize), ((sumofloadsquared/(double)popsize) - (long double) pow((sumofloadsquared/(double)popsize),2)));
            break;
        }
        
        if(popsize >= (maxPopSize-1)){
        	fprintf(summarydatafilepointer, "Population achieved its maximum population size at time %f", t); 
            birthrate = arrayofbirthrates[popsize];
			fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\n", t, popsize, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize), (sumofload/(double)popsize), ((sumofloadsquared/(double)popsize) - (long double) pow((sumofloadsquared/(double)popsize),2)));
            break;
        }
        
        if(tskitstatus == 2 && t >= genafterburnin){
            isburninphaseover = 1;
        }

        if (isredinmaxpopsize && t >= redtime){
		    K = --K;
        	if(r == 1.0){
        		kappa = K/(s*i_init);
    		} else{
			    if(deleteriousdistribution == 0){
            		kappa = (1-r)*K/(s);
        		} else if(deleteriousdistribution == 1){
            		kappa = -(log(r)*K)/(s);
        		}
		}
            calcRateofBirths(arrayofbirthrates, maxPopSize, kappa, b_0);
		    redtime += redinmaxpopsize;
        }
        
        birthhappens = monteCarloStep(arrayofbirthrates, popsize, pCurrenttime, sumofdeathrates);//This is the monte carlo step. This decides if a birth or a death event takes place by returning a 0 or 1
        
        PerformOneEventAbs(tskitstatus, isburninphaseover, ismodular, elementsperlb, &treesequencetablecollection,  wholepopulationnodesarray, wholepopulationsitesarray, pPopSize, pCurrenttime, wholepopulationgenomes, wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationisfree, wholepopulationindex, psumofdeathrates, psumofdeathratessquared, psumofload, psumofloadsquared, parent1gamete, parent2gamete, totaltimesteps, isabsolute, birthhappens, maxPopSize, chromosomesize, numberofchromosomes, totalindividualgenomelength, deleteriousmutationrate, Sd, deleteriousdistribution, beneficialmutationrate, Sb, beneficialdistribution, b_0, r, i_init, s, randomnumbergeneratorforgamma, miscfilepointer);

        if (tskitstatus == 2){
            if (isburninphaseover != 0){
                //This SECTION covers the search for fixed mutations and their time of fixation
                if (iscalcfixation)
                {
                    if (t > simplifyat)
                    {
                        //sort the table
                        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence table collection sorted. \n");
                        fflush(miscfilepointer);
                        //reindex the table
                        returnvaluefortskfunctions = tsk_table_collection_build_index(&treesequencetablecollection, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence table collection index built. \n");
                        fflush(miscfilepointer);
                        tsk_table_collection_print_state(&treesequencetablecollection, miscfilepointer);
                        //change sample status of nodes alive in the population
                        for (k = 0; k < (2*maxPopSize); k++) {
                            tablepointer->nodes.flags[wholepopulationnodesarray[k]] = TSK_NODE_IS_SAMPLE;
                        }
                        //create tree sequence object
                        tsk_treeseq_t treesequence;
                        returnvaluefortskfunctions = tsk_treeseq_init(&treesequence, &treesequencetablecollection, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence object initialized. \n");
                        fflush(miscfilepointer);
                        fprintf(miscfilepointer, "\n sample number %lld. \n", (long long) tsk_treeseq_get_num_samples(&treesequence));
                        //create tree object
                        tsk_tree_t tree;
                        int ret;
                        //initialize a tree object
                        ret = tsk_tree_init(&tree, &treesequence, 0);
                        check_tsk_error(ret);
                        fprintf(miscfilepointer, "\n tree object initialized. \n");
                        fflush(miscfilepointer);

                        for (ret= tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree))
                        {
                            fprintf(miscfilepointer, "\n analyzing tree \n");
                            fflush(miscfilepointer);
                            fprintf(miscfilepointer, "\n (within tree) number of roots of tree %lld is %lld. \n", (long long) tree.index,  (long long) tsk_tree_get_num_roots(&tree));
                            fflush(miscfilepointer);
                            //operation on a tree is to find the MRCA of the sample and then all the nodes between the MRCA and the root
                            FindFixedMutations(&tree, tskitdatafilepointer, miscfilepointer);
                            fprintf(miscfilepointer, "\n analyzing next tree \n");
                            fflush(miscfilepointer);
                        }
                        check_tsk_error(ret);
                        fprintf(miscfilepointer, "\n tree sequence analysis completed. \n");
                        fflush(miscfilepointer);
                        tsk_tree_free(&tree);
                        tsk_treeseq_free(&treesequence);
                        
                        //sort the table
                        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        //simplify after finishing the fixation part
                        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*maxPopSize), 0, NULL);
                        check_tsk_error(returnvaluefortskfunctions);
 
                        for (k = 0; k < (2*maxPopSize); k++) {
                            wholepopulationnodesarray[k] = k;
                        }

                        simplifyat += simplifyeach;
                    }
                }
                //printing section
                if (t > printtime){
                    birthrate = arrayofbirthrates[popsize];
                    if (t > updatesumtime){
                        sumofdeathrates = Fen_sum(wholepopulationselectiontree, maxPopSize);
                        updatesumtime += updatesumofdeathrateseach;
                    }
                    if (popsize > upperlimit || popsize < lowerlimit){
                        index = SearchforPopsize(arrayofpopsizes, maxPopSize, popsize);
                        upperlimit = arrayofpopsizes[index];
                        lowerlimit = arrayofpopsizes[index-1];
                
                        fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\t%Lf\t%Lf\n", t, (upperlimit+lowerlimit)/2, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize), (sumofload/(double)popsize));
                    }
                    if((int)t % 1000 == 0){
                        fflush(rawdatafilepointer);
                    }
                }
            }
        } else if (tskitstatus == 1){
            if (isburninphaseover != 0){
                //This SECTION covers the search for fixed mutations and their time of fixation
                if (iscalcfixation)
                {
                    if (t > simplifyat)
                    {
                        //sort the table
                        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence table collection sorted. \n");
                        fflush(miscfilepointer);
                        //reindex the table
                        returnvaluefortskfunctions = tsk_table_collection_build_index(&treesequencetablecollection, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence table collection index built. \n");
                        fflush(miscfilepointer);
                        tsk_table_collection_print_state(&treesequencetablecollection, miscfilepointer);
                        //change sample status of nodes alive in the population
                        for (k = 0; k < (2*maxPopSize); k++) {
                            tablepointer->nodes.flags[wholepopulationnodesarray[k]] = TSK_NODE_IS_SAMPLE;
                        }
                        //create tree sequence object
                        tsk_treeseq_t treesequence;
                        returnvaluefortskfunctions = tsk_treeseq_init(&treesequence, &treesequencetablecollection, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        fprintf(miscfilepointer, "\n tree sequence object initialized. \n");
                        fflush(miscfilepointer);
                        fprintf(miscfilepointer, "\n sample number %lld. \n", (long long) tsk_treeseq_get_num_samples(&treesequence));
                        //create tree object
                        tsk_tree_t tree;
                        int ret;
                        //initialize a tree object
                        ret = tsk_tree_init(&tree, &treesequence, 0);
                        check_tsk_error(ret);
                        fprintf(miscfilepointer, "\n tree object initialized. \n");
                        fflush(miscfilepointer);

                        for (ret= tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree))
                        {
                            fprintf(miscfilepointer, "\n analyzing tree \n");
                            fflush(miscfilepointer);
                            fprintf(miscfilepointer, "\n (within tree) number of roots of tree %lld is %lld. \n", (long long) tree.index,  (long long) tsk_tree_get_num_roots(&tree));
                            fflush(miscfilepointer);
                            //operation on a tree is to find the MRCA of the sample and then all the nodes between the MRCA and the root
                            FindFixedMutations(&tree, tskitdatafilepointer, miscfilepointer);
                            fprintf(miscfilepointer, "\n moving to next tree \n");
                            fflush(miscfilepointer);
                        }
                        fprintf(miscfilepointer, "\n tree sequence analysis completed. \n");
                        fflush(miscfilepointer);
                        tsk_tree_free(&tree);
                        tsk_treeseq_free(&treesequence);

                        //sort the table
                        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                        check_tsk_error(returnvaluefortskfunctions);
                        //simplify after finishing the fixation part
                        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*maxPopSize), 0, NULL);
                        check_tsk_error(returnvaluefortskfunctions);

                        for (k = 0; k < (2*maxPopSize); k++) {
                            wholepopulationnodesarray[k] = k;
                        }

                        simplifyat += simplifyeach;
                    }
                } else{
                    if (t > simplifyat) {
                        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
                        check_tsk_error(returnvaluefortskfunctions);
            
                        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*maxPopSize), 0, NULL);
                        check_tsk_error(returnvaluefortskfunctions);

                        for (k = 0; k < (2*maxPopSize); k++) {
                            wholepopulationnodesarray[k] = k;
                        }
                        simplifyat += simplifyeach;
                    }
                }
            }
            //printing section
            if(t > printtime){
                birthrate = arrayofbirthrates[popsize];
                if (t > updatesumtime){
                    sumofdeathrates = Fen_sum(wholepopulationselectiontree, maxPopSize);
                    updatesumtime += updatesumofdeathrateseach;
                }
                fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\t%Lf\t%Lf\n", t, popsize, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize), (sumofload/(double)popsize));
                if((int)t % 1000 == 0){
                    fflush(rawdatafilepointer);
                }
                printtime += printeach;
            }
        } else{
            if(t > printtime){
                birthrate = arrayofbirthrates[popsize];
                if (t > updatesumtime){
                    sumofdeathrates = Fen_sum(wholepopulationselectiontree, maxPopSize);
                    updatesumtime += updatesumofdeathrateseach;
                }
                fprintf(rawdatafilepointer, "%f\t%d\t%Lf\t%Lf\t%f\t%Lf\t%Lf\n", t, popsize, (sumofdeathrates/(double)popsize), ((sumofdeathratessquared/(double)popsize) - (long double) pow((sumofdeathrates/(double)popsize),2)), (birthrate/(double)popsize), (sumofload/(double)popsize));
                if((int)t % 1000 == 0){
                    fflush(rawdatafilepointer);
                }
                printtime += printeach;
            }
        }
    }
    
    //END OF SIMULATION FOR LOOP
    //Tree sequence recording requires that tables are sorted on the back end,
    //so I sort once again here at the end to ensure that all tables are sorted before they're read to file.
    //This might be inefficient, I'm not sure.
    if(tskitstatus != 0){
        printf("Sort at final generation %lld: (%lld nodes %lld edges)",
            (long long) t,
            (long long) tablepointer->nodes.num_rows,
            (long long) tablepointer->edges.num_rows);
        returnvaluefortskfunctions = tsk_table_collection_sort(&treesequencetablecollection, NULL, 0);
        check_tsk_error(returnvaluefortskfunctions);

        returnvaluefortskfunctions = tsk_table_collection_simplify(&treesequencetablecollection, wholepopulationnodesarray, (2*maxPopSize), 0, NULL);
        check_tsk_error(returnvaluefortskfunctions);
    
        printf(" -> (%lld nodes %lld edges)\n",
                (long long) tablepointer->nodes.num_rows,
                (long long) tablepointer->edges.num_rows);
    
        for (k = 0; k < (2*maxPopSize); k++) {
            wholepopulationnodesarray[k] = k;
        }
    
    //Printing out the node table in a way readable by python on the back end.
        fprintf(nodefilepointer, "is_sample time\n");
        for (k = 0; k < tablepointer->nodes.num_rows; k++) {
            if (k < (2*popsize)) {
                fprintf(nodefilepointer, "1 %f\n", tablepointer->nodes.time[k]);
                fflush(nodefilepointer);
            } else {
                fprintf(nodefilepointer, "0 %f\n", tablepointer->nodes.time[k]);
                fflush(nodefilepointer);
            }
        }
    
    //Printing out the edge table in a way readable by python on the back end.
        fprintf(edgefilepointer, "left right parent child\n");
        for (k = 0; k < tablepointer->edges.num_rows; k++) {
            fprintf(edgefilepointer, "%f %f %d %d\n", tablepointer->edges.left[k], tablepointer->edges.right[k], tablepointer->edges.parent[k], tablepointer->edges.child[k]);
            fflush(edgefilepointer);
        }
    
    //Printing out the site table in a way readable by python.
        fprintf(sitefilepointer, "position ancestral_state\n");
        for (k = 0; k < tablepointer->sites.num_rows; k++) {
            fprintf(sitefilepointer, "%f 0.0\n", tablepointer->sites.position[k]);
            fflush(sitefilepointer);
        }
    
    //Printing out the mutation table in a way readable by python.
        fprintf(mutationfilepointer, "site node time derived_state\n");
        for (k = 0; k < tablepointer->mutations.num_rows; k++) {
            fprintf(mutationfilepointer, "%d %d %f %.12s\n", tablepointer->mutations.site[k], tablepointer->mutations.node[k], tablepointer->mutations.time[k], (tablepointer->mutations.derived_state + k*12));
            fflush(mutationfilepointer);
        }
    }
    if (VERYVERBOSE == 1) {
        fprintf(veryverbosefilepointer, "Finished simulation with mean sb %f \n", Sb);
        fprintf(veryverbosefilepointer, "Time elapsed: %f\n", t);
        fprintf(veryverbosefilepointer, "Final population size was: %d\n", popsize);
    }
    
    FILE *popsnapshotfilepointer;
    
    char* popsnapshotfilename = MakePopSnapshotFileName(mubname, Sbname);
    popsnapshotfilepointer = fopen(popsnapshotfilename, "w"); //opens the file to which to print summary data.
    
    WritePopSnapshot(wholepopulationgenomes, totalpopulationgenomelength, sumofdeathrates, sumofdeathratessquared, sumofload, sumofloadsquared, wholepopulationselectiontree, wholepopulationdeathratesarray, wholepopulationisfree, wholepopulationindex, maxPopSize, popsize, t, popsnapshotfilepointer);
    
    fclose(rawdatafilepointer);
    fclose(summarydatafilepointer);
    fclose(nodefilepointer);
    fclose(edgefilepointer);
    fclose(sitefilepointer);
    fclose(mutationfilepointer);
    fclose(popsnapshotfilepointer);

    if (iscalcfixation)
    {
        fclose(tskitdatafilepointer);
    }
    
    
    free(popsnapshotfilename);
    free(rawdatafilename);
    free(summarydatafilename);
    
	free(arrayofbirthrates);
    free(arrayofpopsizes);
    free(logaveragefitnesseachNtimesteps);
    free(literallyjustlast200Ntimesteps);
    free(last200Ntimestepsvariance);
    
    free(wholepopulationgenomes);
    free(parent1gamete);
    free(parent2gamete);
    free(wholepopulationselectiontree);
    free(wholepopulationdeathratesarray);
    free(sortedwisarray);
    free(wholepopulationisfree);
    free(wholepopulationindex);
    if (tskitstatus != 0){
        free(wholepopulationnodesarray);
    }
    
    tsk_table_collection_free(&treesequencetablecollection);
    
    return 0;     
}

//calculates every possible birth rate for every possible value of population size within range [0-MaxPopSize]
void calcRateofBirths(double *arrayofbirthrates, int maxPopSize, int kappa, double b_0) {
    int i;
    for(i = 0; i < maxPopSize; i++){
        arrayofbirthrates[i] = (b_0) * (double)i * (1 - ((double)i/(double)kappa));
    }
}

//searches for where the current popsize is an array of offset lines to determine upper and lower limit. This function serves the module of limiting excessive printing.
int SearchforPopsize(int *a, int n, int popsize) {
    int i = 0;

    while (i < n && a[i] < popsize) i++;

    return i;
}

//returns output of either birth or death
bool monteCarloStep(double *arrayofbirthrates, int popsize, double *pCurrentTime, double sumofdeathrates) {

    double deathRate;
    double birthRate;
    double timestep;

    deathRate = sumofdeathrates;

    //rate of births is calculated using equation used in lab write up
    birthRate = arrayofbirthrates[popsize];
    
    timestep = 1/(deathRate + birthRate);//To make time steps dynamical, we directly use the inverse of the sum of the rates as value for a time step.

    *pCurrentTime += timestep;

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


bool PerformOneEventAbs(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int *pPopSize, double * pCurrenttime, double *wholepopulationgenomes, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, double* parent1gamete, double* parent2gamete, int totaltimesteps, bool isabsolute, bool birthhappens, int maxPopSize, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double Sd, int deleteriousdistribution, double beneficialmutationrate, double Sb, int beneficialdistribution,  double b_0, double r, int i_init, double s, gsl_rng* randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    if(isabsolute == 0){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs within a non absolute fitness simulation. \n");
        exit(0);
    }
    
    int randparent1, randparent2, currentparent1, currentparent2, currentvictim;
    int i;
    bool nolethalmut = true;
        
    int victim, birthplace;
    
    if(*pPopSize < 2){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs with a population size lower than 2. \n");
        exit(0);
    }
    
    if(*pPopSize >= maxPopSize && birthhappens){
        fprintf(miscfilepointer, "\n Trying to use PerformOneEventAbs with a population size greater than maxPopSize. \n");
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
            fprintf(miscfilepointer, "\n Trying to give birth to a new individual using a parent that is not present in PerformOneEventAbs. \n");
            exit(0);
        }
        
        tsk_id_t childnode1, childnode2;
        
        double currenttime;

        RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode1, totaltimesteps, currenttime, currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
        nolethalmut = ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode1, totaltimesteps, currenttime, currentparent1, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent1gamete, randomnumbergeneratorforgamma, miscfilepointer);
        if(!nolethalmut)
            return false;
        
        RecombineChromosomesIntoGamete(isabsolute, tskitstatus, ismodular, elementsperlb, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, &childnode2, totaltimesteps, currenttime, currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
        nolethalmut = ProduceMutatedGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationnodesarray, wholepopulationsitesarray, &childnode2, totaltimesteps, currenttime, currentparent2, isabsolute, totalindividualgenomelength, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, parent2gamete, randomnumbergeneratorforgamma, miscfilepointer);
        if(!nolethalmut)
            return false;
        
        PerformBirth(tskitstatus, isburninphaseover, ismodular, elementsperlb, treesequencetablecollection, wholepopulationnodesarray, childnode1, childnode2, isabsolute, parent1gamete, parent2gamete, maxPopSize, pPopSize, birthplace, wholepopulationgenomes, totalindividualgenomelength, deleteriousdistribution, wholepopulationselectiontree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, psumofdeathrates, psumofdeathratessquared, b_0, r, i_init, s, psumofload, psumofloadsquared, miscfilepointer);
            }
     
    else{
//         use of the fennwick tree to find the victim in the population, note that since fenwick tree stores all the population (non occupied spaces have a fitness of 0) there is not need to use the wholepopulationindex, fenwick search already gives you the space that the selected individual occupies.
    	victim = ChooseVictimWithTree(wholepopulationselectiontree, *pPopSize, maxPopSize, *psumofdeathrates, miscfilepointer);
        
//         printf("%d\n", victim);
//         
//         for(i = 0; i < maxPopSize; i++)
//             printf("%2d  ", wholepopulationisfree[i]);
//         printf("\n");
        
        PerformDeath(isabsolute, tskitstatus, isburninphaseover, maxPopSize, pPopSize, victim, deleteriousdistribution, wholepopulationselectiontree, wholepopulationwisarray, wholepopulationdeathratesarray, wholepopulationindex, wholepopulationisfree, psumofloads, psumofdeathrates,psumofdeathratessquared, b_0, r, i_init, s, psumofload, psumofloadsquared, wholepopulationnodesarray, miscfilepointer);
    }
    
    return true;

}


void InitializePopulationAbs(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, int totaltimesteps, int deleteriousdistribution, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, double b_0, double r, int i_init, double s)
{    
    int i, j;

    double haploidgenomelength = (double) ((totalpopulationgenomelength / maxPopSize) / 2);
    
    double starting_load;

    if(r == 1.0){
        starting_load = b_0 - s*i_init;
    }
    else{
        if(deleteriousdistribution == 0){
            starting_load = b_0 - s*(1 - pow(r, i_init))/(1-r);
        }
        else if(deleteriousdistribution == 1){
            starting_load = b_0 - s*(pow(r, (i_init -1)) -1)/(log(r));
        }
    }
    for (i = 0; i < initialPopSize; i++) {
        wholepopulationselectiontree[i] = starting_load; //all individuals start with death rate d_0. Currently d_0 is set through the input parameters
        
        wholepopulationdeathratesarray[i] = starting_load;
        
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
    
    *psumofload = (long double) initialPopSize * i_init;

    *psumofloadsquared = (long double) initialPopSize * pow(i_init, 2); 

    *psumofdeathrates = (long double) initialPopSize * starting_load;
    
    *psumofdeathratessquared = (long double) initialPopSize * pow(starting_load, 2);
    
    for (i = 0; i < totalpopulationgenomelength; i++){
        wholepopulationgenomes[i] = 0.0;
    }
    if (tskitstatus != 0){
        treesequencetablecollection->sequence_length = haploidgenomelength;
    
    //The following lines initialize the node table for tree sequence recording.
    //Note that nodes here are single sets of chromosomes, so the 2x popsize here assumes diploidy.
        for (i = 0; i < (2 * maxPopSize); i++) {
            wholepopulationnodesarray[i] = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, totaltimesteps, TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(wholepopulationnodesarray[i]);
        }
    
    //The following lines add a site to the tree sequence recording site table corresponding to each linkage block, with ancestral state of 0.
        for (i = 0; i < haploidgenomelength; i++) {
            wholepopulationsitesarray[i] = tsk_site_table_add_row(&treesequencetablecollection->sites, i, "0.0000000000", 12, NULL, 0);
            check_tsk_error(wholepopulationsitesarray[i]);
        }
    }
}

//used to initialize a population using a file that stores all the important variables
void InitializeWithSnapshotAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, int *pPopSize, double *pCurrenttime, FILE *prevsnapshotfilepointer, FILE *miscfilepointer)
{
    int i, j;

    char strtemp[200];
    
    //gets the genomes of the whole population from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    if(strcmp(strtemp, "Genomes:")){
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers do not include Genomes: check that the snapshot file is in the correct parameters folder");
        exit(0);
    }
    for (i = 0; i < totalpopulationgenomelength; i++){
        fscanf(prevsnapshotfilepointer, "%lf", &wholepopulationgenomes[i]);
    }
    
    //gets the sum of death rates from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    if(strcmp(strtemp, "Sum_of_death_rates:")){
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers do not include Sum_of_death_rates: check that the snapshot file is in the correct parameters folder or that there are no errors in the file");
        exit(0);
    }
    fscanf(prevsnapshotfilepointer, "%Lf", psumofdeathrates);
    
    //gets the sum of death rates squared from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    fscanf(prevsnapshotfilepointer, "%Lf", psumofdeathratessquared);

    //gets the sum of load from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    fscanf(prevsnapshotfilepointer, "%Lf", psumofload);
    
    //gets the sum of load sqaured from the snapshot file
    fscanf(prevsnapshotfilepointer, "%s" , strtemp);
    fscanf(prevsnapshotfilepointer, "%Lf", psumofloadsquared);
    
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
        fprintf(miscfilepointer, "\nError: Corrupted snapshot file, headers do not include Time: check that the snapshot file is in the correct parameters folder or that there are no errors in the file");
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
        fprintf(miscfilepointer, "\n Array of indexes is corrupted in findinindex. \n");
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

double CalculateDeathRate(bool ismodular, int elementsperlb, double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, int deleteriousdistribution,double b_0, double r, int i_init, double s)
{
    double inddeathrate = 0.0;
    double currentlinkageblocksload = 0.0;
    double currentlinkageblocksloadmod[elementsperlb];
    double sumofloadmod = 0;
    int totallinkageblocks;
    int i, m;

    if (!ismodular){
        //since the load is calculated per gamete, the total number of linkages blocks is half the total genome length of an ind
        totallinkageblocks = totalindividualgenomelength/2;
        for (i = 0; i < totallinkageblocks; i++) {
            currentlinkageblocksload += parent1gamete[i];
            currentlinkageblocksload += parent2gamete[i];
        }

	//load is calculated as the number of mutations per individual. An heterozygous should have 1/2 mutation, and an homozygous 1 mutation. Thus we divide the sum of the load of the two sister chromosomes by 2
	currentlinkageblocksload = currentlinkageblocksload/2; 
        if(r == 1.0){
            inddeathrate = b_0 - s*(i_init - currentlinkageblocksload);
        }   
        else{
            if(deleteriousdistribution == 0){
                inddeathrate = b_0 - s*(1 - pow(r, (i_init -currentlinkageblocksload)))/(1-r);
            }
            else if(deleteriousdistribution == 1)
                inddeathrate = b_0 - s*(pow(r, (i_init -currentlinkageblocksload -1)) -1)/(log(r));
            }
    }
    else{
        printf("Error no code for CalculateDeathRate function with modularity yet");
        exit(0);       
    }
    return inddeathrate;
}

void FindFixedMutations(const tsk_tree_t *tree, FILE *tskitdatafilepointer, FILE *miscfilepointer)
{
    const tsk_id_t *parent = tree->parent;
    //find sample IDs
    const tsk_id_t *samples = tsk_treeseq_get_samples(tree->tree_sequence);
    tsk_size_t num_samples = tsk_treeseq_get_num_samples(tree->tree_sequence);
    fprintf(miscfilepointer, "\n (within tree) number of roots of tree is %d. \n", tsk_treeseq_get_num_samples(tree->tree_sequence));
    fflush(miscfilepointer);
    //find MRCA of sample
    tsk_id_t mrca;
    tsk_id_t u, v; 
    tsk_size_t j;
    int ret, numnodes;
    int numroots = tsk_tree_get_num_roots(tree);
    fprintf(miscfilepointer, "\n (within tree) number of roots of tree is %d. \n", tsk_tree_get_num_roots(tree));
    fflush(miscfilepointer);
    
    mrca = samples[0];
    fprintf(miscfilepointer, "\n (within tree) MRCA loop start. \n");
    fflush(miscfilepointer);
    //Loop through the sample and find the pairwise MRCA until finding the MRCA of all the sample
    
    for (j = 1; j < num_samples; j++) {   
        ret = tsk_tree_get_mrca(tree, mrca, samples[j], &mrca);
        check_tsk_error(ret);
        if (mrca == TSK_NULL || mrca == samples[j])
            break;
    }

    fprintf(miscfilepointer, "\n (within tree) MRCA %lld found. \n", mrca);
    fflush(miscfilepointer);
    //Find all nodes from MRCA to root
    tsk_id_t *nodes = malloc(tsk_tree_get_size_bound(tree) * sizeof(*nodes));
    if (nodes == NULL) {
        errx(EXIT_FAILURE, "Out of memory");
    }
    
    ret = tsk_tree_get_depth(tree, mrca, &numnodes);
    check_tsk_error(ret);
    fprintf(miscfilepointer, "\n (within tree) Depth calculated. \n");
    fflush(miscfilepointer);

    for (int i = 0; i < numnodes+1; i++)
    { 
        nodes[i] = mrca;
        mrca = parent[mrca];
    }
    fprintf(miscfilepointer, "\n (within tree) all nodes between MRCA and root found. \n");
    fflush(miscfilepointer);
    //Find and print mutations and time of appearance occuring at each node
    const tsk_id_t *mutation_nodes = tree->tree_sequence->tables->mutations.node;
    const tsk_id_t *sites = tree->tree_sequence->tables->mutations.site;
    if (mutation_nodes == NULL) {
        errx(EXIT_FAILURE, "Out of memory 2");
    }

    double fixationtime;
    fprintf(miscfilepointer, "\n (within tree) starting fixation time calculation. \n");
    fflush(miscfilepointer);
    for(j = 0; j < numnodes+1; j++)
    {
        for(int k = 0; k < tree->tree_sequence->tables->mutations.num_rows; k++)
        {
            while (sites[k] == tree->index)
            {
                if(nodes[j] == mutation_nodes[k])
                {
                    fixationtime = FindFixationTime(tree, tree->tree_sequence->tables->mutations.time[k]);
                    fprintf(tskitdatafilepointer, "%.12s\t%f\t%f\n", (tree->tree_sequence->tables->mutations.derived_state + k*12), tree->tree_sequence->tables->mutations.time[k], fixationtime);
                    fflush(tskitdatafilepointer);
                }
            } 
        }
    }
    fprintf(miscfilepointer, "\n (within tree) fixation time calculation finished. \n");
    fflush(miscfilepointer);
    free(nodes);
}

double FindFixationTime(const tsk_tree_t *tree, double timeappearance)
{
    int ret;
    tsk_size_t num_nodes, k;
    const double *node_time = tree->tree_sequence->tables->nodes.time;
    tsk_node_t node;
    tsk_id_t leaf;
    tsk_edge_t edge;
    double ret_tbl;
    double x;

    for (k = 0; k < tree->tree_sequence->tables->nodes.num_rows; k++) {
        while (node_time[k] >= timeappearance)
        {
            if (node_time[k] == timeappearance)
            {
                ret = tsk_treeseq_get_node(tree->tree_sequence, k, &node);
                check_tsk_error(ret);

                FindDeepestNode(tree, node.id, &leaf);
                ret_tbl = node_time[k] - node_time[leaf];
            
                if(ret_tbl > x)
                    x = ret_tbl;
            }
        }  
    }

    for (k = tree->tree_sequence->tables->edges.num_rows; k > 0 ; k--)
    {
        ret = tsk_treeseq_get_edge(tree->tree_sequence, k, &edge);
        check_tsk_error(ret);

        if (node_time[edge.parent] > timeappearance)
        {
            FindDeepestNode(tree, node.id, &leaf);
                ret_tbl = node_time[k] - node_time[leaf];
            
                if(ret_tbl > x)
                    x = ret_tbl;
        }
    }
    return x;
}

// Recursive function to find the deepest node in the tree
void FindDeepestNodeHelper(const tsk_tree_t *tree, tsk_id_t node, int level, int* deepestLevel, tsk_id_t* deepestNode) {
    if (node == TSK_NULL) {
        return;
    }
    if (level > *deepestLevel) {
        *deepestLevel = level;
        *deepestNode = node;
    }
    FindDeepestNodeHelper(tree, tree->left_child[node], level + 1, deepestLevel, deepestNode);
    FindDeepestNodeHelper(tree, tree->right_child[node], level + 1, deepestLevel, deepestNode);
}

// Function to find the deepest node in the tree
void FindDeepestNode(const tsk_tree_t *tree, tsk_id_t root, tsk_id_t* leaf) {
    int deepestLevel = -1;
    tsk_id_t deepestNode = TSK_NULL;
    FindDeepestNodeHelper(tree, root, 0, &deepestLevel, &deepestNode);
    *leaf = deepestNode;
}

void WritePopSnapshot(double *wholepopulationgenomes, int totalpopulationgenomelength, long double sumofdeathrates, long double sumofdeathratessquared, long double sumofload, long double sumofloadsquared, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, int maxPopSize, int popsize, double t, FILE *popsnapshotfilepointer)
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
    fprintf(popsnapshotfilepointer, "Sum_of_load:\n");
    fprintf(popsnapshotfilepointer, "%Lf", sumofload);
    fprintf(popsnapshotfilepointer, "\n");
    fprintf(popsnapshotfilepointer, "Sum_of_load_squared:\n");
    fprintf(popsnapshotfilepointer, "%Lf", sumofloadsquared);
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
