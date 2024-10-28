#ifndef ABSOLUTE_FUNCTIONS_H_INCLUDED
#define ABSOLUTE_FUNCTIONS_H_INCLUDED 1

#include <stdio.h>
#include <stdlib.h>
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
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

double RunSimulationAbs(bool issnapshot, char *prevsnapshotfilename, bool isredinmaxpopsize, char *redinmaxpopsizename, double redinmaxpopsize, char* mubname, char* Sbname, int tskitstatus, bool ismodular, int elementsperlb, bool isabsolute, int maxTime, int initialPopSize, int K, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double Sd, int deleteriousdistribution, double beneficialmutationrate, double Sb, int beneficialdistribution, double r, int i_init, double s, gsl_rng* randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, bool iscalcfixation);

void calcRateofBirths(double *arrayofbirthrates, int maxPopSize, int kappa, double b_0);
bool monteCarloStep(double *arrayofbirthrates, int popsize, double *pCurrentTime, double sumofdeathrates);
bool discoverEvent(double deathRate, double birthRate);


bool PerformOneEventAbs(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int *pPopSize, double * pCurrenttime, double *wholepopulationgenomes, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, double* parent1gamete, double* parent2gamete, int totaltimesteps, bool isabsolute, bool birthhappens, int maxPopSize, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double Sd, int deleteriousdistribution, double beneficialmutationrate, double Sb, int beneficialdistribution,  double b_0, double r, int i_init, double s, gsl_rng* randomnumbergeneratorforgamma, FILE *miscfilepointer);


void InitializePopulationAbs(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int initialPopSize, int maxPopSize, int totaltimesteps, int deleteriousdistribution, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, double b_0, double r, int i_init, double s);

void InitializeWithSnapshotAbs(long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, int maxPopSize, double *wholepopulationgenomes, int totalpopulationgenomelength, long double *psumofdeathrates, long double *psumofdeathratessquared, long double *psumofload, long double *psumofloadsquared, int *pPopSize, double *pCurrenttime, FILE *prevsnapshotfilepointer, FILE *miscfilepointer);

int ChooseParent(int populationsize);
int ChooseVictimWithTree(long double *wholepopulationselectiontree, int popsiez, int maxPopSize, long double inversesumofloads, FILE *miscfilepointer);

int findinindex(int *wholepopulationindex, int which, int tam, FILE *miscfilepointer);
void indexArrayFlipDeath(int *wholepopulationindex, int placeinindex, int popsize);

double CalculateDeathRate(bool ismodular, int elementsperlb, double *parent1gamete, double *parent2gamete, int totalindividualgenomelength, int deleteriousdistribution, double b_0, double r, int i_init, double s);

int SearchforPopsize(int *a, int n, int popsize);

void FindFixedMutations(const tsk_tree_t *tree, FILE *tskitdatafilepointer, FILE *miscfilepointer);
double FindFixationTime(const tsk_tree_t *tree, double timeappearance);
void FindDeepestNodeHelper(const tsk_tree_t *tree, tsk_id_t node, int level, int* deepestLevel, tsk_id_t* deepestNode);
void FindDeepestNode(const tsk_tree_t *tree, tsk_id_t root, tsk_id_t* leaf);

void WritePopSnapshot(double *wholepopulationgenomes, int totalpopulationgenomelength, long double sumofdeathrates, long double sumofdeathratessquared, long double sumofload, long double sumofloadsquared, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, bool *wholepopulationisfree, int *wholepopulationindex, int maxPopSize, int popsize, double t, FILE *popsnapshotfilepointer);

#endif // ABSOLUTE_FUNCTIONS_H_INCLUDED
