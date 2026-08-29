#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED 1

#include <stdio.h>
#include <float.h>
#include <string.h>
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
#include "relative_functions.h"
#include "global_vars.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

//fitness
void UpdateLast200NTimeSteps(double *last200Ntimesteps, double newNtimesteps);
/* DoubleSwap / DoubleBubbleSort - COMMENTED OUT (item 10: completely unused).
 * DoubleBubbleSort was never even defined, only declared; DoubleSwap exists only
 * to serve it. Nothing in the program sorts fitnesses. Kept, commented, rather
 * than deleted.
 * void DoubleSwap(long double *x, long double *y);
 * void DoubleBubbleSort(long double *arraytobesorted, int arraysize);
 */
double CalculateVarianceInLogFitness(int popsize, Individual *wholepopulation, long double sumofwis);
long double FindFittestWi(Individual *wholepopulation, int popsize);
double CalculateSlopeOfLogFitness(int endofsimulation, int endofburninphase, double *logaveragefitnesseachgeneration);

//Fenwick trees
long double Fen_sum(long double *tree, int i);
void Fen_add(long double *tree, int numberofelementsintree, long double amounttoadd, int i);
long double Fen_range(long double *tree, int i, int j);
long double Fen_get(long double *tree, int i);
void Fen_set(long double *tree, int numberofelementsintree, long double newvalue, int i);
int SearchTree(int leftbound, int rightbound, long double targetvalue, long double *Fenwicktree);

//Distribution
/* ExponentialDerivate - COMMENTED OUT (item 10: completely unused). Effect sizes
 * are drawn with gsl_ran_exponential instead. Kept, commented, rather than deleted.
 * double ExponentialDerivate(double mean);
 */
int SampleFromPoisson(float poissonmean);

//Gametes.
//UH change to two recombination sites.
//RecombineChromosomesIntoGamete and ProduceMutatedGamete are declared in sharedfunc_flag.h (included above),
//which now carries the mutator-aware signatures; do not redeclare them here.
int DetermineNumberOfMutations(double mutationrate);
int DetermineMutationSite(int totalgametelength);

//root simulations. UH might move them to shared with flags and remove files pointers (they are declared globally)
int BracketZeroForSb(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig);
/* BisectionMethodToFindSbWithZeroSlope - COMMENTED OUT (item 10: never called).
 * double BisectionMethodToFindSbWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, char * mutator_switch_ratename, char * mutator_biasname, char * mutator_strength_factorname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize, MutatorConfig mutatorconfig, TrackingConfig trackingconfig);
 */
//make directory name and final data file name function. GOAL is to increase modularity in program
char * MakeDirectoryName(char * tskitstatus, char* deldist, char * isabsolutename, bool isabsolute, char * bendist, char * benmut, char * numberofchromosomes, char * chromosomesize, char * popsize, char * delmut, char * randomnumberseed, char * K, char * r, char *i_init, char * s, bool ismodular, char *elementsperlb, char *iscalcfixationname, int typeofrun, char * Sbname, char *Sdname);
/* MakeFinalDataFileName - COMMENTED OUT (item 10: never called).
 * char * MakeFinalDataFileName(char * typeofrun, char * benmut, char * slopeforcontourline, char * randomnumberseed);
 */
/* MakeRawDataFileName - COMMENTED OUT (item 10: never called).
 * char * MakeRawDataFileName(char * mubname, char * Sbname, bool isredinmaxpopsize, char *redinmaxpopsizename);
 */
/* MakeSummaryDataFileName - COMMENTED OUT (item 10: never called).
 * char * MakeSummaryDataFileName(char * mubname, char * Sbname);
 */
/* MakePopSnapshotFileName - COMMENTED OUT (item 10: never called).
 * char * MakePopSnapshotFileName(char * mubname, char * Sbname);
 */
//functions to assign arguments passed in the command line to their proper variables
int AssignArgumentstoVar(char **argv, int *Nxtimesteps, char *Nxtimestepsname, int *popsize, char *popsizename, double *deleteriousmutationrate, char *deleteriousmutationratename, int *chromosomesize, char *chromosomesizename, int *numberofchromosomes, char *numberofchromosomesname, double *bentodelmutrate, double *Sbtemp, int *beneficialdistribution, int *typeofrun, double *slopeforcontourline, char *slopeforcontourlinename, int *randomnumberseed, char *randomnumberseedname, int *K, char *Kname, int *relorabs, double *r, char *rname, int *i_init, char *i_initname, double *s, char *sname, int *tskitstatus, int *nonmodormod, int *elementsperlb, char *elementsperlbname, int *snapshot, char *prevsnapshotfilename, double *SdtoSbratio, char *SdtoSbrationame, int *deleteriousdistribution, int *rawdatafilesize, double *redinmaxpopsize, char *redinmaxpopsizename, int *calcfixation, double *mutator_strength_factor, char *mutator_strength_factorname, double *mutator_switch_rate, char *mutator_switch_ratename, double *mutator_bias, char *mutator_biasname, double *modifier_locus_fraction, int *modifier_mask_mode, int *antimutator_encoding, double *initial_mutator_fraction, int *trackindividuals, int *trackinterval, int *trackstartgen);
void AssignStringNames(char *beneficialmutationratename, double beneficialmutationrate, char *bendistname, int beneficialdistribution, char *deldistname, int deleteriousdistribution, char *typeofrunname, int typeofrun, char *tskitstatusname, int tskitstatus, char* Sb2name, double Sb2, char *isabsolutename, bool isabsolute, char *iscalcfixationname, bool iscalcfixation, double Sd, char *Sdname);

#endif // MAIN_H_INCLUDED



