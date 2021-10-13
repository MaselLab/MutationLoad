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


//We might want to delete the file arguments because they are global variables UH
double RunSimulation(char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char * beneficialdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE * veryverbosefilepointer, FILE * verbosefilepointer, FILE * miscfilepointer);
//We might want to remove currenttimestep because this argument is not used at all
void PerformOneTimeStep(int popsize, int totaltimesteps, int currenttimestep, long double *wholepopulationwistree, long double *wholepopulationwisarray, double *wholepopulationgenomes, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, char * beneficialdistribution, double *parent1gamete, double *parent2gamete, gsl_rng * randomnumbergeneratorforgamma);

int ChooseVictim(int populationsize);
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis);

void ReplaceVictim(double *parent1gamete, double *parent2gamete, int currentpopsize, int currentvictim, long double *sumofwis, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationwistree, long double *wholepopulationwisarray);

double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);

#endif // RELATIVE_FUNCTIONS_H_INCLUDED

