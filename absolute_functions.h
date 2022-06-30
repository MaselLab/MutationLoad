#ifndef ABSOLUTE_FUNCTIONS_H_INCLUDED
#define ABSOLUTE_FUNCTIONS_H_INCLUDED 1

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
#include "main.h"

double RunSimulationAbs(bool issnapshot, char *prevsnapshotfilename, bool isabsolute, char* maxTimeStepsname, char* popsizename, char* delmutratename, char* chromsizename, char* chromnumname, char* mubname, char* Sbname, int typeofrun, int maxTimeSteps, int initialPopSize, int maxPopSize, double d_0, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, gsl_rng* randomnumbergeneratorforgamma, double r, double sdmin, FILE *miscfilepointer, FILE *veryverbosefilepointer);

void calcRateofBirths(double *arrayofbirthrates, int maxPopSize, double b_0);
bool monteCarloStep(double *arrayofbirthrates, int popsize, double sumofdeathrates, double *pCurrenttime);
bool discoverEvent(double deathRate, double birthRate);

bool PerformOneEventAbs(bool isabsolute, bool birthhappens, int maxPopSize, int *pPopSize, double *wholepopulationgenomes, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, long double *psumofdeathrates, long double *psumofdeathratessquared, double d_0, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution,  double* parent1gamete, double* parent2gamete, gsl_rng* randomnumbergeneratorforgamma, double r, double sdmin, FILE *miscfilepointer);


void InitializePopulationAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double* psumofdeathrates, long double *psumofdeathratessquared, double d_0);

void InitializeWithSnapshotAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, int *pPopSize, double *pCurrenttime, FILE *prevsnapshotfilepointer, FILE *miscfilepointer);

int ChooseParent(int populationsize);
int ChooseVictimWithTree(long double *wholepopulationselectiontree, int popsiez, int maxPopSize, long double inversesumofloads, FILE *miscfilepointer);

int findinindex(int *wholepopulationindex, int which, int tam, FILE *miscfilepointer);
void indexArrayFlipDeath(int *wholepopulationindex, int placeinindex, int popsize);

double CalculateDeathRate(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, double d_0, double r, double sdmin);

void WritePopSnapshot(double *wholepopulationgenomes, int totalpopulationgenomelength, long double sumofdeathrates, long double sumofdeathratessquared, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, int maxPopSize, int popsize, double t, FILE *popsnapshotfilepointer);

#endif // ABSOLUTE_FUNCTIONS_H_INCLUDED
