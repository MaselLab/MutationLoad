#ifndef SHAREFUNC_FLAG_H_INCLUDED
#define SHAREFUNC_FLAG_H_INCLUDED 1

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
#include "global_vars.h"
#include "main.h"


void InitializePopulation(bool isabsolute, long double *wholepopulationwistree, long double *wholepopulationwisarray, int populationsize, double *populationgenomes, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis);

double PerformDeath(bool isabsolute, int currentpopsize, int currentvictim, long double *sumofwis, long double *wholepopulationwistree, long double *wholepopulationwisarray);
void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int currentpopsize, int currentvictim, long double *sumofwis, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationwistree, long double *wholepopulationwisarray);



#endif // SHAREFUNC_FLAG_H_INCLUDED

