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
#include "absolute_functions.h"
#include "sharedfunc_flag.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>

void MutateGamete(int tskitstatus, int isburninphaseover,  tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double * pCurrenttime, bool isabsolute, int totalindividualgenomelength, double *gamete, double mutationeffectsize)
{
    tsk_id_t idofnewmutation;
    
    int mutatedsite = DetermineMutationSite(totalindividualgenomelength/2);
    if(isabsolute)
        gamete[mutatedsite] += (mutationeffectsize);
    else
        gamete[mutatedsite] += log(1 + mutationeffectsize);
    char derivedstate[400];
    sprintf(derivedstate, "%.11f", mutationeffectsize);
    
    if (tskitstatus != 0){
        if (isburninphaseover != 0){
            idofnewmutation = tsk_mutation_table_add_row(&treesequencetablecollection->mutations, wholepopulationsitesarray[mutatedsite], childnode, TSK_NULL, ((double) totaltimesteps - *pCurrenttime), derivedstate, 12, NULL, 0);
            check_tsk_error(idofnewmutation);
        }
    }
}

double PerformDeath(bool isabsolute, int maxPopSize, int *pPopSize, int victim, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, FILE *miscfilepointer)
{
    int placeinindex;
    if(isabsolute){
        if(wholepopulationisfree[victim]){
            fprintf(miscfilepointer, "\n Array of free indexes is corrupted in PerformDeath. \n");
            exit(0);
        }
        
        *psumofdeathrates -= wholepopulationdeathratesarray[victim];
        *psumofdeathratessquared -= pow(wholepopulationdeathratesarray[victim], 2);
        
        wholepopulationisfree[victim] = true;
        wholepopulationdeathratesarray[victim] = 0.0;
        
        //Joseph way of doing the tree just with popsize is better, once this is working change it
        placeinindex = findinindex(wholepopulationindex, victim, *pPopSize, miscfilepointer);
        indexArrayFlipDeath(wholepopulationindex, placeinindex, *pPopSize);
        
        
        *pPopSize -= 1;
    }
    
    else{
        *psumofloads -= wholepopulationwisarray[victim];
        wholepopulationwisarray[victim] = 0.0;
    }
    
    Fen_set(wholepopulationselectiontree, maxPopSize, 0.0, victim);
}


void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, FILE *miscfilepointer)
{
    int i;
    
    long double newwi;
    
    long double inddeathrate;
    
//     printf("%Lf %Lf \n", newwi, newInverse);
    
    bool placefound = false;
    
    if(isabsolute){
        if(*pPopSize >= maxPopSize){
            fprintf(miscfilepointer, "\n Trying to use PerformBirth with a population size greater than maxPopSize. \n");
            exit(0);
        }
        
        birthplace = wholepopulationindex[*pPopSize];
        
        if(!wholepopulationisfree[birthplace]){
            printf("wholepopulationisfree is corrupted in PerformBirth \n");
            exit(0);
        }
        

        inddeathrate = (long double) CalculateDeathRate(ismodular, elementsperlb, parent1gamete, parent2gamete, totalindividualgenomelength, deleteriousdistribution, b_0, r, i_init, s);

        
        if(inddeathrate < 0.0){
            fprintf(miscfilepointer, "\n The individual death rate of a new born is less than 0.0 (he is even more than immortal). \n");
            exit(0);
        }
    }
    else{
        newwi = (long double) CalculateWi(parent1gamete, parent2gamete, totalindividualgenomelength);
    }
    
    for (i = 0; i < (totalindividualgenomelength/2); i++) {
        wholepopulationgenomes[birthplace*totalindividualgenomelength + i] = parent1gamete[i];
        wholepopulationgenomes[birthplace*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
    }
        
    if(isabsolute){
        Fen_set(wholepopulationselectiontree, maxPopSize, inddeathrate, birthplace);
        wholepopulationdeathratesarray[birthplace] = inddeathrate;
        wholepopulationisfree[birthplace] = false;
        
        *psumofdeathrates += inddeathrate;
        *psumofdeathratessquared += pow(inddeathrate, 2);        
        *pPopSize += 1;
    }
    else{
        Fen_set(wholepopulationselectiontree, maxPopSize, newwi, birthplace);
        wholepopulationwisarray[birthplace] = newwi;
        
        *psumofloads += newwi;
    }
    if (tskitstatus != 0){
        if (isburninphaseover !=0){
            wholepopulationnodesarray[birthplace*2] = childnode1;
            wholepopulationnodesarray[birthplace*2 + 1] = childnode2;
        }      
    }
}
