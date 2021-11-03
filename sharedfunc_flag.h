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
#include "absolute_functions.h"
#include "global_vars.h"
#include "main.h"


void InitializePopulation(bool isabsolute, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofloads, long double* pInversesumofloads);

double PerformDeath(bool isabsolute, int maxPopSize, int *pPopSize, int victim, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *pInversesumofloads);

void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationloadstree, long double *wholepopulationwisarray, long double *wholepopulationdeathratessarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *pInversesumofloads);



#endif // SHAREFUNC_FLAG_H_INCLUDED

