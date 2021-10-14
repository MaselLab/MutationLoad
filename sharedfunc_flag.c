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
#include "global_vars.h"
#include "main.h"


void InitializePopulation(bool isabsolute, long double *wholepopulationwistree, long double *wholepopulationwisarray, int populationsize, double *populationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis) {
    if(!isabsolute){
        int i, j;
        
        double haploidgenomelength = (double) ((totalpopulationgenomelength / populationsize) / 2);
        
        for (i = 0; i < populationsize; i++) {
            wholepopulationwistree[i] = 1.0; //for relative fitness, all individuals start with probability of being chosen as a parent of 1/N
            wholepopulationwisarray[i] = 1.0;
        }
        //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
        for (i = 0; i < populationsize; i++) {
        j = i + LSB(i+1);
            if (j < populationsize) {
                wholepopulationwistree[j] += wholepopulationwistree[i];
            }
        }
        for (i = 0; i < totalpopulationgenomelength; i++) {
            populationgenomes[i] = 0.0;
        }
        *psumofwis = (long double) populationsize;
    //     printf("%LF \n", *psumofwis);
    }
}

double PerformDeath(bool isabsolute, int currentpopsize, int currentvictim, long double *sumofwis, long double *wholepopulationwistree, long double *wholepopulationwisarray)
{
    if(!isabsolute){
        *sumofwis -= wholepopulationwisarray[currentvictim];
        Fen_set(wholepopulationwistree, currentpopsize, 0, currentvictim);
        wholepopulationwisarray[currentvictim] = (long double) 0.0;
    }
}

void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int currentpopsize, int currentvictim, long double *sumofwis, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationwistree, long double *wholepopulationwisarray)
{
    if(!isabsolute){
        int i;
        double newwi;
        
        newwi = CalculateWiRel(parent1gamete, parent2gamete, totalindividualgenomelength);
        
        for (i = 0; i < (totalindividualgenomelength/2); i++) {
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + i] = parent1gamete[i];
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
    //It would be more efficient to build directly into victim slot, but right now it is possible for the victim to also be a parent, later add an if statement to more efficiently deal with more common case.
        }

        Fen_set(wholepopulationwistree, currentpopsize, newwi, currentvictim);
        
        wholepopulationwisarray[currentvictim] = (long double) newwi;
        *sumofwis += (long double) newwi;
    }
}
