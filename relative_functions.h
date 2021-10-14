#ifndef RELATIVE_FUNCTIONS_H_INCLUDED
#define RELATIVE_FUNCTIONS_H_INCLUDED 1

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
#include "sharedfunc_flag.h"
#include "global_vars.h"
#include "main.h"


double RunSimulationRel(char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma);

void PerformOneTimeStepRel(int popsize, int totaltimesteps, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma);

int ChooseVictim(int populationsize);
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis);

double CalculateWiRel(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);

#endif // RELATIVE_FUNCTIONS_H_INCLUDED
