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
#include "global_vars.h"
#include "main.h"


void InitializePopulation(bool isabsolute, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double * psumofloads, long double* pInversesumofloads) {
    
    int i, j;
    
    for (i = 0; i < initialPopSize; i++) {
        wholepopulationloadstree[i] = 1.0; //all individuals start with load 1 (probability of being chosen to produce an offspring or to die of 1/N). wholepopulationloadstree is mult((1+s_m)(1+s_f)) for relative fitness simulations and mult(1/(1+s_m)(1+s_f)) for absolute fitness simulations.
    }
    
    if(isabsolute){
        for (i = 0; i < initialPopSize; i++) {
            wholepopulationdeathratessarray[i] = 1.0;
            
            wholepopulationindex[i] = i;
            
            wholepopulationisfree[i] = false;
        }
        
        for (i = initialPopSize; i < maxPopSize; i++) {
            wholepopulationloadstree[i] = 0.0;
            
            wholepopulationdeathratessarray[i] = 0.0;
                                    
            wholepopulationindex[i] = i;
            
            wholepopulationisfree[i] = true;
        }
        //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page. Absolute trees needs to contemplate all popSize, including non occupied spaces. 
        for (i = 0; i < maxPopSize; i++) {
            j = i + LSB(i+1);
            if (j < maxPopSize) {
                wholepopulationloadstree[j] += wholepopulationloadstree[i];
            }
        }
        
        *pInversesumofloads = (long double)initialPopSize;
    }
    
    else{
        for (i = 0; i < initialPopSize; i++)
            wholepopulationwisarray[i] = 1.0;
        //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
        for (i = 0; i < initialPopSize; i++) {
            j = i + LSB(i+1);
            if (j < initialPopSize) {
                wholepopulationloadstree[j] += wholepopulationloadstree[i];
            }
        }
        
        *psumofloads = (long double)initialPopSize;
    }
    
    for (i = 0; i < totalpopulationgenomelength; i++){
        wholepopulationgenomes[i] = 0.0;
    }
    
}


double PerformDeath(bool isabsolute, int maxPopSize, int *pPopSize, int victim, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *pInversesumofloads)
{
    int placeinindex;
    if(isabsolute){
        if(wholepopulationisfree[victim]){
            fprintf(miscfilepointer, "\n Array of free indexes is corrupted in PerformDeath. \n");
            exit(0);
        }
        
        *pInversesumofloads -= wholepopulationdeathratessarray[victim];
        
        wholepopulationisfree[victim] = true;
        wholepopulationdeathratessarray[victim] = 0.0;
        
        //Joseph way of doing the tree just with popsize is better, once this is working change it
        placeinindex = findinindex(wholepopulationindex, victim, *pPopSize);
        indexArrayFlipDeath(wholepopulationindex, placeinindex, *pPopSize);
        
        
        *pPopSize -= 1;
    }
    
    else{
        *psumofloads -= wholepopulationwisarray[victim];
        wholepopulationwisarray[victim] = 0.0;
    }
    
    Fen_set(wholepopulationloadstree, maxPopSize, 0.0, victim);
}

void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *pInversesumofloads)
{
    int i;
    
    long double newload = (long double) CalculateLoad(parent1gamete, parent2gamete, totalindividualgenomelength);
    
    long double newInverse = 1.0 / (long double) newload;
    
//     printf("%Lf %Lf \n", newload, newInverse);
    
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
    }
    
    for (i = 0; i < (totalindividualgenomelength/2); i++) {
            wholepopulationgenomes[birthplace*totalindividualgenomelength + i] = parent1gamete[i];
            wholepopulationgenomes[birthplace*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
    }
        
    if(isabsolute){
        Fen_set(wholepopulationloadstree, maxPopSize, newInverse, birthplace);
        wholepopulationdeathratessarray[birthplace] = newInverse;
        wholepopulationisfree[birthplace] = false;
        
        *pInversesumofloads += newInverse;
                
        *pPopSize += 1;
    }
    else{
        Fen_set(wholepopulationloadstree, maxPopSize, newload, birthplace);
        wholepopulationwisarray[birthplace] = newload;
        
        *psumofloads += newload;
    }
    
}
