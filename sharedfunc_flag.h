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
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

// Global Individual struct definition
typedef struct{
    double *fitnessArray; // Still need to change the name to WiArray
    int *mutatorArray;    // 0 = Anti-mutator, 1 = Mutator
    double fitness;      // Overall fitness (Wi)
    double mutationRate; // Actual mutation rate based on mutator loci
} Individual;

void MutateGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationsitesarray, tsk_id_t childnode, int totaltimesteps, double currenttimestep, bool isabsolute, int totalindividualgenomelength, double *gamete, double mutationeffectsize);

double PerformDeath(bool isabsolute, int tskitstatus, int isburninphaseover, int maxPopSize, int *pPopSize, int victim, int deleteriousdistribution, long double *wholepopulationselectiontree, Individual *wholepopulation, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, tsk_id_t * wholepopulationnodesarray, FILE *miscfilepointer);

// Updated to accept mutator arrays and mutator parameters
void PerformBirth(int tskitstatus, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t childnode1, tsk_id_t childnode2, bool isabsolute, double *parent1gameteFitness, int *parent1gameteMutators, double *parent2gameteFitness, int *parent2gameteMutators, int maxPopSize, int *pPopSize, int birthplace, Individual *wholepopulation, int totalindividualgenomelength, int deleteriousdistribution, long double *wholepopulationselectiontree, long double *wholepopulationdeathratesarray, int *wholepopulationindex, bool *wholepopulationisfree, long double *psumofloads, long double *psumofdeathrates, long double *psumofdeathratessquared, double b_0, double r,  int i_init, double s, long double *psumofload, long double *psumofloadsquared, FILE *miscfilepointer, double mutator_strength_factor, double baseline_mutation_rate);

// Individual helper functions
Individual createIndividual(double *fitnessArray, int *mutatorArray, int totalindividualgenomelength);
void UpdateIndividual(Individual *ind, int totalindividualgenomelength, double mutator_strength_factor, double baseline_mutation_rate);

// Updated Recombine and Mutate prototypes
void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gameteFitness, int *gameteMutators, Individual *wholepopulation, int totalindividualgenomelength);

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double parent_specific_mutation_rate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gameteFitness, int *gameteMutators, double mutator_switch_rate, double mutator_bias, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer);

#endif // SHAREFUNC_FLAG_H_INCLUDED