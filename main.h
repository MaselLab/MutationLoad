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

//fitness
void UpdateLast200NTimeSteps(double *last200Ntimesteps, double newNtimesteps);
void DoubleSwap(long double *x, long double *y);
void DoubleBubbleSort(long double *arraytobesorted, int arraysize);
double CalculateVarianceInLogFitness(int popsize, long double *wholepopulationwisarray, long double sumofwis);
long double FindFittestWi(long double *wisarray, int popsize);
double CalculateSlopeOfLogFitness(int endofsimulation, int endofburninphase, double *logaveragefitnesseachgeneration);

//Fenwick trees
long double Fen_sum_squared(long double *tree, int i);
long double Fen_sum(long double *tree, int i);
void Fen_add(long double *tree, int numberofelementsintree, long double amounttoadd, int i);
long double Fen_range(long double *tree, int i, int j);
long double Fen_get(long double *tree, int i);
void Fen_set(long double *tree, int numberofelementsintree, long double newvalue, int i);
int SearchTree(int leftbound, int rightbound, long double targetvalue, long double *Fenwicktree);

//Distribution
double ExponentialDerivate(double mean);
int SampleFromPoisson(float poissonmean);

//Gametes.
//UH change to two recombination sites.
void RecombineChromosomesIntoGamete(int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double * pCurrenttime, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, double *wholepopulationgenomes, int totalindividualgenomelength);
bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover,tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double * pCurrenttime, int parent, bool isabsolute, int individualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer);
int DetermineNumberOfMutations(double mutationrate);
int DetermineMutationSite(int totalgametelength);

//root simulations.
int BracketZeroForSb(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize);
double BisectionMethodToFindSbWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize);
//make directory name and final data file name function. GOAL is to increase modularity in program
char * MakeDirectoryName(char * tskitstatus, char* deldist, char * isabsolutename, bool isabsolute, char * bendist, char * benmut, char * numberofchromosomes, char * chromosomesize, char * popsize, char * delmut, char * randomnumberseed, char * K, char * r, char *i_init, char * s, bool ismodular, char *elementsperlb);
char * MakeFinalDataFileName(char * typeofrun, char * benmut, char * slopeforcontourline, char * randomnumberseed);
char * MakeRawDataFileName(char * mubname, char * Sbname);
char * MakeSummaryDataFileName(char * mubname, char * Sbname);
char * MakePopSnapshotFileName(char * mubname, char * Sbname);
//functions to assign arguments passed in the command line to their proper variables
int AssignArgumentstoVar(char **argv, int *Nxtimesteps, char *Nxtimestepsname, int *popsize, char *popsizename, double *deleteriousmutationrate, char *deleteriousmutationratename, int *chromosomesize, char *chromosomesizename, int *numberofchromosomes, char *numberofchromosomesname, double *bentodelmutrate, double *Sbtemp, int *beneficialdistribution, int *typeofrun, double *slopeforcontourline, char *slopeforcontourlinename, int *randomnumberseed, char *randomnumberseedname, int *K, char *Kname, int *relorabs, double *r, char *rname, int *i_init, char *i_initname, double *s, char *sname, int *tskitstatus, int *nonmodormod, int *elementsperlb, char *elementsperlbname, int *snapshot, char *prevsnapshotfilename, double *SdtoSbratio, char *SdtoSbrationame, int *deleteriousdistribution, int *rawdatafilesize, int *redinmaxpopsize);
void AssignStringNames(char *beneficialmutationratename, double beneficialmutationrate, char *bendistname, int beneficialdistribution, char *deldistname, int deleteriousdistribution, char *typeofrunname, int typeofrun, char *tskitstatusname, int tskitstatus, char* Sb2name, double Sb2, char *isabsolutename, bool isabsolute);

#endif // MAIN_H_INCLUDED


