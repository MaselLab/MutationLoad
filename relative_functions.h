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
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

/* The three loose mutator doubles that used to sit at the end of these
 * signatures are now carried in MutatorConfig (see sharedfunc_flag.h), together
 * with the new modifier-locus settings. TrackingConfig carries the optional
 * per-individual dump settings. Both are passed by value. */
double RunSimulationRel(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * Sbname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig);

void PerformOneTimeStepRel(int tskitstatus, bool isabsolute, int isburninphaseover, bool ismodular, int elementsperlb, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int totaltimesteps, double currenttimestep, long double *wholepopulationwistree, Individual *wholepopulation, long double * psumofwis, int chromosomesize, int numberofchromosomes, int totalindividualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *parent1gameteFitness, int *parent1gameteMutators, char *parent1gameteMask, ModifierLocusIndex *parent1modifierindex, double *parent2gameteFitness, int *parent2gameteMutators, char *parent2gameteMask, ModifierLocusIndex *parent2modifierindex, const char *globalmodifiermask, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, MutatorConfig mutatorconfig);

void InitializePopulationRel(int tskitstatus, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, long double *wholepopulationwistree, Individual *wholepopulation, int popsize, int totalpopulationgenomelength, int totaltimesteps, long double * psumofwis, char *globalmodifiermask, MutatorConfig mutatorconfig, FILE *miscfilepointer);

int ChooseVictim(int populationsize);
int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis, FILE *miscfilepointer);

/* --- modifier-locus helpers (mutation-rate evolution) --------------------- */
void DrawModifierMask(char *mask, int haploidgenomelength, double locusfraction);
void SeedInitialMutatorStates(int *stateshaploid, const char *mask, int haploidgenomelength, double initialmutatorfraction, int antimutatorstate);

/* --- tree-sequence bootstrap (item 7: honour tskitstatus == 2) ------------ */
void SeedTreeSequenceTables(tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, int popsize, int haploidgenomelength, double nodetime);

/* --- output helpers (item 5) ---------------------------------------------- */
void WritePopulationModifierSummary(FILE *rawdatafilepointer, Individual *wholepopulation, int popsize, int totalindividualgenomelength, const char *globalmodifiermask, int *locusmutatorcounts, int *locusmodifiercounts);
void WriteIndividualSnapshot(FILE *individualfilepointer, Individual *wholepopulation, int popsize, int generation);

/* ---------------------------------------------------------------------------
 * CalculateWi - COMMENTED OUT (item 10: completely unused).
 * ---------------------------------------------------------------------------
 * Superseded by the inline Wi recomputation inside PerformBirth() and by
 * UpdateIndividual(), both of which work on the Individual struct rather than on
 * a pair of raw gamete arrays. Kept here, commented, rather than deleted.
 *
 * double CalculateWi(double *parent1gamete, double *parent2gamete, int totalindividualgenomelength);
 * ------------------------------------------------------------------------- */

#endif // RELATIVE_FUNCTIONS_H_INCLUDED
