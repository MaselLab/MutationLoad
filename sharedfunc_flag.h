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
#include "main.h"

void MutateGamete(bool isabsolute, int chromosomesize, int numberofchromosomes, double *gamete, double mutationeffectsize);

double PerformDeath(bool isabsolute, int maxPopSize, int *pPopSize, int victim, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, FILE *miscfilepointer);

void PerformBirth(bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double d_0, double r, double sdmin, FILE *miscfilepointer);



#endif // SHAREFUNC_FLAG_H_INCLUDED

