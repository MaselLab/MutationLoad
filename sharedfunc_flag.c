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

void MutateGamete(bool isabsolute, int chromosomesize, int numberofchromosomes, double *gamete, double mutationeffectsize)
{
    int randomchromosometomutate = pcg32_boundedrand(numberofchromosomes); //if we decide to include heterogenous rates of recombination/mutation, both of these will need to be replaced by a function that weights each linkage block's probability of mutating.
    int randomblocktomutate = pcg32_boundedrand(chromosomesize);
    int mutatedsite = randomchromosometomutate*chromosomesize + randomblocktomutate;
    if(isabsolute)
        gamete[mutatedsite] += (mutationeffectsize);
    else
        gamete[mutatedsite] += log(1 + mutationeffectsize);

}

double PerformDeath(bool isabsolute, int maxPopSize, int *pPopSize, int victim, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, FILE *miscfilepointer)
{
    int placeinindex;
    if(isabsolute){
        if(wholepopulationisfree[victim]){
            fprintf(miscfilepointer, "\n Array of free indexes is corrupted in PerformDeath. \n");
            exit(0);
        }
        
        *psumofdeathrates -= wholepopulationdeathratesarray[victim];
        
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

void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, double d_0, double r, double sdmin, FILE *miscfilepointer)
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
        
        inddeathrate = CalculateDeathRate(parent1gamete, parent2gamete, totalindividualgenomelength, d_0, r, sdmin);
        
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
                
        *pPopSize += 1;
    }
    else{
        Fen_set(wholepopulationselectiontree, maxPopSize, newwi, birthplace);
        wholepopulationwisarray[birthplace] = newwi;
        
        *psumofloads += newwi;
    }
    
}
