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
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

void MutateGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, double mutationeffectsize);

double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer);

void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gamete, double *parent2gamete, int maxPopSize, int *pPopSize, int birthplace, double *wholepopulationgenomes, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationwisarray, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer);



#endif // SHAREFUNC_FLAG_H_INCLUDED

