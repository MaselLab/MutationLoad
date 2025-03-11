/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include <stdlib.h>
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
#include "absolute_functions.h"
#include "global_vars.h"
#include "main.h"
#include <tskit.h>
#include <tskit/tables.h>
#include <kastore.h>
#include <tskit/core.h>
#include <tskit/trees.h>

int main(int argc, char *argv[]) {

    if (argc != 27) {
        printf("[Error]; Wrong number of arguments in program. usage: timeSteps; initialPopsize; genome-wide deleterious mutation rate; chromosome size; number of chromosomes; beneficial/deleterious mutation ratio; Sb; beneficial distribution; type of run; slope; seed; MaxPopSize; relative or absolute; epistasis r; i_init; deltad_0; tskit status; without modular epistasis or with it; elements per linkage block; snapshot; snapshot file name; deleterious/beneficial s ratio; point or exponential deleterious distribution; size of raw data file; reduction in MaxPopSize; Fixation Calculation status \n");
        return -1;
    }
    //declare the file pointers for the files used for printing across the program
    FILE *miscfilepointer;
    FILE *verbosefilepointer;
    FILE *finaldatafilepointer;
    FILE *veryverbosefilepointer;
    
    int Nxtimesteps, popsize, chromosomesize, numberofchromosomes, beneficialdistribution, typeofrun, randomnumberseed, K, relorabs, i_init, tskitstatus, nonmodormod, elementsperlb, snapshot, deleteriousdistribution, rawdatafilesize, calcfixation;
	double deleteriousmutationrate, bentodelmutrate, Sbtemp, slopeforcontourline, r, s, SdtoSbratio, redinmaxpopsize;

	char *Nxtimestepsname, *popsizename, *deleteriousmutationratename, *chromosomesizename, *numberofchromosomesname, *slopeforcontourlinename, *randomnumberseedname, *Kname, *rname, *i_initname, *sname, *elementsperlbname, *prevsnapshotfilename, *SdtoSbrationame, *redinmaxpopsizename, *iscalcfixationname;

    Nxtimestepsname = (char *)malloc(30);
	popsizename = (char *)malloc(30);
	deleteriousmutationratename = (char *)malloc(30);
	chromosomesizename = (char *)malloc(30);
	numberofchromosomesname = (char *)malloc(30);
	slopeforcontourlinename = (char *)malloc(30);
	randomnumberseedname = (char *)malloc(30);
	Kname = (char *)malloc(30);
	rname = (char *)malloc(30);
	i_initname = (char *)malloc(30);
	sname = (char *)malloc(30);
	elementsperlbname = (char *)malloc(30);
	prevsnapshotfilename = (char *)malloc(200);
	SdtoSbrationame = (char *)malloc(30);
    redinmaxpopsizename = (char *)malloc(30);
    iscalcfixationname = (char *)malloc(30);

    int wrong_args;
	wrong_args = AssignArgumentstoVar(argv, &Nxtimesteps, Nxtimestepsname, &popsize, popsizename, &deleteriousmutationrate, deleteriousmutationratename, &chromosomesize, chromosomesizename, &numberofchromosomes, numberofchromosomesname, &bentodelmutrate, &Sbtemp, &beneficialdistribution, &typeofrun, &slopeforcontourline, slopeforcontourlinename, &randomnumberseed, randomnumberseedname, &K, Kname, &relorabs, &r, rname, &i_init, i_initname, &s, sname, &tskitstatus, &nonmodormod, &elementsperlb, elementsperlbname, &snapshot, prevsnapshotfilename, &SdtoSbratio, SdtoSbrationame, &deleteriousdistribution, &rawdatafilesize, &redinmaxpopsize, redinmaxpopsizename, &calcfixation);
    //main run error check
    if(wrong_args == -1){
		return -1;
    }

	bool isabsolute;
    //checking absolute vs relative fitness model choice
	if(relorabs == 0){
		isabsolute = false;
    }else if(relorabs == 1){
		isabsolute = true;
    }
    //checking if modular chromosome model or not
	bool ismodular;
	if(nonmodormod == 0){
		ismodular = false;
    }else if(nonmodormod == 1){
		ismodular = true;
    }
    //checking if snapshot of previous run is being used or fresh run
	bool issnapshot;
	if(snapshot == 0){
		issnapshot = false;
    }else if(snapshot == 1){
		issnapshot = true;
    }
    //checking if reduction in max popsize is being run or not
    bool isredinmaxpopsize;
    if(redinmaxpopsize == 0.0){
        isredinmaxpopsize = false;
    } else{
        isredinmaxpopsize = true;
    }
    //checking if fixation time calculation module using tree sequences is being used or not
    bool iscalcfixation;
	if(calcfixation == 0){
		iscalcfixation = false;
    }else if(calcfixation == 1){
		iscalcfixation = true;
    }
    //I have two parameters for Sb for the type of run that needs to have bracketed values of Sb.
	//In the case with just a single simulation being run, Sb2 here will be the value of Sb used.
	double Sb1;
	double *pSb1 = &Sb1;

	double Sb2;
	double *pSb2 = &Sb2;

	Sb1 = 0.0;
	if(!isabsolute){
		Sb2 = Sbtemp;
    }else{
		//For absolute runs mean sb is set to 1. Beneficial DFE has a mean 1.0
		Sb2 = 1.0;
	}

    //I have two parameters for N for the type of run that needs to have bracketed values of N.
	int N1;
    N1 = (int) popsize*0.25;//we choose the lower N to be 25% lower than the input N
	int *pN1 = &N1;

	int N2;
    N2 = popsize;
	int *pN2 = &N2;

	//calculate the beneficial mutation rate from the inputted ratio
	double beneficialmutationrate = bentodelmutrate*deleteriousmutationrate;
	double Sd = Sb2*SdtoSbratio;

	char *beneficialmutationratename, *bendistname, *deldistname, *typeofrunname, *tskitstatusname, *Sb2name, *isabsolutename;
	beneficialmutationratename = (char *) malloc(30);
	bendistname = (char *) malloc(30);
	deldistname = (char *) malloc(30);
	typeofrunname = (char *) malloc(30);
	tskitstatusname = (char *) malloc(30);
	Sb2name = (char *) malloc(30);
	isabsolutename = (char *) malloc(30);

	AssignStringNames(beneficialmutationratename, beneficialmutationrate, bendistname, beneficialdistribution, deldistname, deleteriousdistribution, typeofrunname, typeofrun,tskitstatusname, tskitstatus, Sb2name, Sb2, isabsolutename, isabsolute, iscalcfixationname, iscalcfixation);

    pcg32_srandom(randomnumberseed, randomnumberseed); // seeds the random number generator.
    gsl_rng * randomnumbergeneratorforgamma = gsl_rng_alloc(gsl_rng_mt19937);
    //the gamma distribution function requires a gsl random number generator, which is set here.
    //it's a bit inelegant to have two different RNGs, which could be solved by using a different algorithm 
    //for choosing variates from a gamma distribution, instead of using the free one from gsl.
    
    char * directoryname = MakeDirectoryName(tskitstatusname, deldistname, isabsolutename, isabsolute, bendistname, beneficialmutationratename, numberofchromosomesname, chromosomesizename, popsizename, deleteriousmutationratename, randomnumberseedname, Kname, rname, i_initname, sname, ismodular, elementsperlbname, iscalcfixationname, typeofrun, Sb2name);// this will create the directory name pointer using input parameter values
    
    mkdir(directoryname, 0777);//create the directory with the directory name pointer
    chdir(directoryname);//move into the created directory
    
    //open files to which to print: verbose, very verbose and misc data. According to snapshot status, files are opened as write (w) or append (a)
    if(!issnapshot){
        verbosefilepointer = fopen("verbose.txt", "w");
        veryverbosefilepointer = fopen("veryverbose.txt", "w");
        miscfilepointer = fopen("miscellaneous.txt", "w");
    }else{
        verbosefilepointer = fopen("verbose.txt", "a");
        veryverbosefilepointer = fopen("veryverbose.txt", "a");
        miscfilepointer = fopen("miscellaneous.txt", "a");
    }
    
    
   //START OF RUNS
    if (typeofrun == 0) {
        fprintf(miscfilepointer, "Running root finding for Sb");
        //create and open the final data file name using the input parameter values
        char * finaldatafilename = MakeFinalDataFileName(typeofrunname, beneficialmutationratename, slopeforcontourlinename, randomnumberseedname);
        finaldatafilepointer = fopen(finaldatafilename, "w");
        /*This type of run finds the Sb value for the given set of parameters
         that produces a population whose fitness stays almost exactly stable.
         It does this by finding values of Sb that lead to populations definitely
         increasing and definitely decreasing in fitness,
         and then searching between them until it finds a value of Sb that leads
         to a population with a long-term slope of fitness that is within an error term of zero.
         */
        double sbrequiredforzeroslopeoffitness;
        
        fprintf(miscfilepointer, "Beginning bracketing function.");
        
        fflush(miscfilepointer);
        
        BracketZeroForSb(tskitstatus, isabsolute, ismodular, elementsperlb, pSb1, pSb2, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, verbosefilepointer, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        
        fprintf(miscfilepointer, "Finished bracketing function.");
        
        fflush(miscfilepointer);
        sbrequiredforzeroslopeoffitness = BisectionMethodToFindSbWithZeroSlope(tskitstatus, isabsolute, ismodular, elementsperlb, pSb1, pSb2, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, verbosefilepointer, finaldatafilepointer, veryverbosefilepointer, rawdatafilesize);
        
        fprintf(finaldatafilepointer, "The value of Sb for which the slope of log fitness is zero with mub of %.10f is %.10f", beneficialmutationrate, sbrequiredforzeroslopeoffitness);

        free(finaldatafilename);
		fclose(finaldatafilepointer);
    
    } else if (typeofrun == 1){
        double slope;
        fprintf(miscfilepointer, "Running single run simulation");
        if(!isabsolute){
            //This type of run just simulates a single population with the input parameters.
            slope = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        }else{
            RunSimulationAbs(issnapshot, prevsnapshotfilename, isredinmaxpopsize, redinmaxpopsizename, redinmaxpopsize, beneficialmutationratename, Sb2name, tskitstatus, ismodular, elementsperlb, isabsolute, Nxtimesteps, popsize, K, chromosomesize, numberofchromosomes, deleteriousmutationrate, Sd, deleteriousdistribution, beneficialmutationrate, Sb2, beneficialdistribution, r, i_init, s, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize, iscalcfixation);
        }
        fprintf(miscfilepointer, "slope of run is %.10f", slope);
        
    } else if (typeofrun == 2){
        fprintf(miscfilepointer, "Running root finding for N crit\n");
        //create and open the final data file name using the input parameter values
        char * finaldatafilename = MakeFinalDataFileName(typeofrunname, beneficialmutationratename, slopeforcontourlinename, randomnumberseedname);
        finaldatafilepointer = fopen(finaldatafilename, "w");
        /*This type of run finds the Ncrit value for the given set of parameters
         that produces a population whose fitness stays almost exactly stable.
         It does this by finding values of Ncrit that lead to populations definitely
         increasing and definitely decreasing in fitness,
         and then searching between them until it finds a value of Ncrit that leads
         to a population with a long-term slope of fitness that is within an error term of zero.
         */
        int wrong_bounds;
        double Ncritrequiredforzeroslopeoffitness;
        
        fprintf(miscfilepointer, "Beginning bracketing function.\n");
        
        fflush(miscfilepointer);
        
        wrong_bounds = BracketZeroForNcrit(tskitstatus, isabsolute, ismodular, elementsperlb, pN1, pN2, Nxtimestepsname, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, typeofrun, Nxtimesteps, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, Sb2name, Sb2, deleteriousdistribution, randomnumbergeneratorforgamma, verbosefilepointer, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        
        if(wrong_bounds == -1){
            fprintf(miscfilepointer, "Bracketing function failed.\n");
		    return -1;
        }
        
        fprintf(miscfilepointer, "Finished bracketing function.\n");
        
        fflush(miscfilepointer);
        fprintf(miscfilepointer, "Beginning Secant function.\n");

        fflush(miscfilepointer);
        Ncritrequiredforzeroslopeoffitness = SecantMethodToFindNcritWithZeroSlope(tskitstatus, isabsolute, ismodular, elementsperlb, pN1, pN2, Nxtimestepsname, deleteriousmutationratename, chromosomesizename, numberofchromosomesname, beneficialmutationratename, typeofrun, Nxtimesteps, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, slopeforcontourline, beneficialdistribution, Sd, Sb2name, Sb2, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, verbosefilepointer, finaldatafilepointer, veryverbosefilepointer, rawdatafilesize);
        
        fprintf(finaldatafilepointer, "The value of Ncrit for which the slope of log fitness is zero with mub of %.10f, sb of %.10f, Ud of %.4f is %.2f", beneficialmutationrate, Sb2, deleteriousmutationrate, Ncritrequiredforzeroslopeoffitness);
        fflush(finaldatafilepointer);

        free(finaldatafilename);
		fclose(finaldatafilepointer);
    } else{
        //One day maybe I'll have more types of runs.
        fprintf(miscfilepointer, "That type of run is not currently supported.");
        return -1;
    }
    
    free(directoryname);

	free(Nxtimestepsname);
	free(popsizename);
	free(deleteriousmutationratename);
	free(chromosomesizename);
	free(numberofchromosomesname);
	free(slopeforcontourlinename);
	free(randomnumberseedname);
	free(Kname);
	free(rname);
	free(i_initname);
	free(sname);
	free(elementsperlbname);
	free(prevsnapshotfilename);
	free(SdtoSbrationame);
    free(iscalcfixationname);
	free(beneficialmutationratename);
	free(bendistname);
	free(deldistname);
	free(typeofrunname);
	free(tskitstatusname);
	free(Sb2name);
	free(isabsolutename);

    fclose(verbosefilepointer);
    fclose(veryverbosefilepointer);
    fclose(miscfilepointer);
    
    gsl_rng_free(randomnumbergeneratorforgamma);
    
    return 0;
}

void UpdateLast200NTimeSteps(double * last200Ntimesteps, double newNtimesteps)
{
    double storage[200];
    int m;
    for (m = 0; m < 199; m++) {
        storage[m] = last200Ntimesteps[m+1];
    }
    storage[199] = newNtimesteps;
    for (m = 0; m < 200; m++) {
        last200Ntimesteps[m] = storage[m];
    }
}

void DoubleSwap(long double * x, long double * y)
{
    long double temp = *x;
    *x = *y;
    *y = temp;
}

//Inefficient algorithm (guaranteed to be order n^2). Improve algorithm if using more than once per simulation.
void DoubleBubbleSort(long double *arraytobesorted, int arraysize)
{
    int i, j;
    for (i = 0; i < arraysize-1; i++) {
        for (j = 0; j < arraysize-i-1; j++) {
            if (arraytobesorted[j] > arraytobesorted[j+1]) {
                DoubleSwap(&arraytobesorted[j], &arraytobesorted[j+1]);
            }
        }
    }
}

double CalculateVarianceInLogFitness(int popsize, long double *wholepopulationwisarray, long double sumofwis)
{
    int i;
    double variancesum;
    variancesum = 0.0;
    long double logaverage;
    logaverage = log(sumofwis / popsize);
    for (i = 0; i < popsize; i++) {
        variancesum += (double) pow((log(wholepopulationwisarray[i]) - logaverage), 2);
    }
    variancesum = (variancesum/(popsize-1)); //note that in the article (Mawass et al. 2024) the Bessel correction is missing here
    //However, given that we use the variance measure only to decide when to end the burn-in period, the effect of its absence is very minimal in practice over long runs
    return variancesum;
}

long double FindFittestWi(long double *wisarray, int popsize)
{
    long double fittestwi;
    int i;
    fittestwi = wisarray[0];
    for (i = 1; i < popsize; i++) {
        if (wisarray[i] > fittestwi) {
            fittestwi = wisarray[i];
        }
    }
    return fittestwi;
}

double CalculateSlopeOfLogFitness(int endofsimulation, int endofburninphase, double *logaveragefitnesseachgeneration)
{
    size_t step = 1;
    int k;
    double c0, cov00, cov01, cov11, sumsq;
    int generationsafterburnin;
    double slopeoflogfitness;
    generationsafterburnin = (endofsimulation - endofburninphase);
    
    //I have to make an array of numbers to use as the x variable in a linear model to find the slope.
    double *justnumbers;
    justnumbers = malloc(sizeof(double) * generationsafterburnin);
    for (k = 0; k < generationsafterburnin; k++) {
        justnumbers[k] = (k+1);
    }
    //The following function fits a linear model to the two variables (generations and logfitness)
    //and records the parameters of the best-fitting linear model in the c0, cov00, cov01, cov11, sumsq, and slopeoflogfitness variables.
    //I only use the slope parameter, but the others are there in case I need them.
    gsl_fit_linear(justnumbers, step, logaveragefitnesseachgeneration, step, generationsafterburnin, &c0, &slopeoflogfitness, &cov00, &cov01, &cov11, &sumsq);
    free(justnumbers);
    return slopeoflogfitness;

}

/*All Fenwick tree functions from Wikipedia page "Fenwick tree" URL:https://en.wikipedia.org/wiki/Fenwick_tree
 This project is licensed under the GNU General Public License version 3.0, 
 which is compatible with the CC-BY-SA license of Wikipedia text.*/


//Returns sum of first i elements in the tree, 0 through i-1.
long double Fen_sum(long double *tree, int i)
{
    long double sum = 0;
    while (i) {
        sum += tree[i-1];
        i -= LSB(i);
    }
    return sum;
}

//Adds an amount to the ith element in the tree (and therefore to the Fen_sum for all elements in the tree greater than i).
void Fen_add(long double *tree, int numberofelementsintree, long double amounttoadd, int i)
{
    while (i < numberofelementsintree) {
        tree[i] += amounttoadd;
        i += LSB(i+1);
    }
}

//Returns the sum of the elements i through j-1.
//Could do with Fen_sum of j minus Fen_sum of i, but this is faster.
long double Fen_range(long double *tree, int i, int j)
{
    long double sum = 0;
    while (j > i) {
        sum += tree[j-1];
        j -= LSB(j);
    }
    while (i > j) {
        sum -= tree[i-1];
        i -= LSB(i);
    }
    return sum;
}

//Returns the value of the element at index i.
long double Fen_get(long double *tree, int i)
{
    return Fen_range(tree, i, i+1);
}

void Fen_set(long double *tree, int numberofelementsintree, long double newvalue, int i)
{
    Fen_add(tree, numberofelementsintree, newvalue - Fen_get(tree, i), i);
}

int SearchTree(int leftbound, int rightbound, long double targetvalue, long double *Fenwicktree)
{
    int middle;
    middle = floor((leftbound+rightbound)/2);
    long double partialsumatmiddle;
    long double partialsumatmiddleminusone;
    partialsumatmiddle = Fen_sum(Fenwicktree, middle);
    partialsumatmiddleminusone = Fen_sum(Fenwicktree, middle-1);
    if(partialsumatmiddle < targetvalue) {
        if((middle+1) == rightbound) {
            return middle;
        }
        return SearchTree(middle, rightbound, targetvalue, Fenwicktree);
    }
    if(partialsumatmiddle > targetvalue) {
        if(partialsumatmiddleminusone > targetvalue) {
            return SearchTree(leftbound, middle, targetvalue, Fenwicktree);
        } else {
            return (middle-1);
        }
    }
    if (partialsumatmiddle == targetvalue) {
        return middle;
    }
}


double ExponentialDerivate(double mean) {
    double result;
    float randnumb;
    do
        randnumb = ldexp(pcg32_random(), -32);
    while (randnumb == 0.0);

    result = (-log(randnumb))*mean;
    
    return result;
}

//From Numerical Recipes in C, Second Edition.
int SampleFromPoisson(float poissonmean)
{
    static float sq, logmean, g;
    static float oldmean = (-1.0);
    float numberofmutations, t, y;

    if (poissonmean < 12.0) {		//for small enough means, use direct method.
        if (poissonmean != oldmean) {	//check to see if the mean value is new.
                oldmean = poissonmean;
                g = exp(-poissonmean);	//if the mean is new, compute the exponential.
        }
        numberofmutations = -1;
        t = 1.0;
        do {
                ++numberofmutations;
                t *= ldexp(pcg32_random(), -32); //instead of adding exponential deviates, multiply uniform deviates and compare to pre-computed exponential.
        } while (t > g);
    } 
    else { 				//for larger means, use rejection method.
        if (poissonmean != oldmean) {	//for new means, pre-compute some functions.
                oldmean = poissonmean;
                sq = sqrt(2.0*poissonmean);
                logmean = log(poissonmean);
                g = poissonmean*logmean - gsl_sf_lngamma(poissonmean+1.0); //lngamma function is the natural log of the gamma function
        }
        do {
                do {
            y = tan(PI * ldexp(pcg32_random(), -32)); 	//makes y a deviate from a Lorentzian comparison function.
            numberofmutations = sq*y + poissonmean;		//shifts and scales y and sets results as possible numberofmutations (to be accepted or rejected);
                } while (numberofmutations < 0.0); 			//rejects values in zero probability area.
                numberofmutations = floor(numberofmutations);
                t = 0.9 * (1.0 + y*y) * exp(numberofmutations*logmean - gsl_sf_lngamma(numberofmutations + 1.0) - g);
        } while (ldexp(pcg32_random(), -32) > t);
    }
return numberofmutations;
}

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(bool isabsolute, int tskitstatus, bool ismodular, int elementsperlb, int isburninphaseover, tsk_table_collection_t * treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, double *wholepopulationgenomes, int totalindividualgenomelength)
{
    int recombinationsite, startchromosome, startofindividual, h, i, returnvaluefortskit;
    startofindividual = persontorecombine * totalindividualgenomelength;
    
    tsk_id_t parentnode1 = (tsk_id_t) 2*persontorecombine;
    tsk_id_t parentnode2 = (tsk_id_t) (2*persontorecombine + 1);
        
    if (tskitstatus != 0){
        if (isabsolute){
            if (isburninphaseover != 0){
                double parenttime1 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode1]];
                double parenttime2 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode2]];    
    
                *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
                check_tsk_error(*childnode);
    
                double childtime = treesequencetablecollection->nodes.time[*childnode];
            }
        }else{
            double parenttime1 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode1]];
            double parenttime2 = treesequencetablecollection->nodes.time[wholepopulationnodesarray[parentnode2]];    
    
            *childnode = tsk_node_table_add_row(&treesequencetablecollection->nodes, 0, ((double) totaltimesteps - currenttimestep), TSK_NULL, TSK_NULL, NULL, 0);
            check_tsk_error(*childnode);
        }
    }

    for (h = 0; h < numberofchromosomes; h++) {
	    startchromosome = pcg32_boundedrand(2); //generates either a zero or a one to decide to start with chromosome 1 or 2.
	
        do {
            recombinationsite = pcg32_boundedrand(chromosomesize);
        } while (recombinationsite == 0); //it doesn't make sense to do a recombination event before the first linkage block. Note that this will never break if the chromosome size is only one linkage block.
        
        //Tree sequence recording needs the recombination sites to add edges.
        if (tskitstatus != 0){
            if(isabsolute){
                if (isburninphaseover != 0){
                    if (startchromosome == 0){
                        returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (double)(h*chromosomesize), (double)(h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                        check_tsk_error(returnvaluefortskit);
                        returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                        check_tsk_error(returnvaluefortskit);
                    }else{
                        returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize), (h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                        check_tsk_error(returnvaluefortskit);
                        returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                        check_tsk_error(returnvaluefortskit);
                    }
                }
            } else{
                if (startchromosome == 0){
                    returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (double)(h*chromosomesize), (double)(h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                    check_tsk_error(returnvaluefortskit);
                    returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                    check_tsk_error(returnvaluefortskit);
                }else{
                    returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize), (h*chromosomesize + recombinationsite), wholepopulationnodesarray[parentnode2], *childnode, NULL, 0);
                    check_tsk_error(returnvaluefortskit);
                    returnvaluefortskit = tsk_edge_table_add_row(&treesequencetablecollection->edges, (h*chromosomesize + recombinationsite), ((h+1)*chromosomesize), wholepopulationnodesarray[parentnode1], *childnode, NULL, 0);
                    check_tsk_error(returnvaluefortskit);
                }
            }
        }
        if(!ismodular){
            for (i = 0; i < recombinationsite; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize) + i];
                }
                else {
                    gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize) + i];
                }
            }
            for (i = recombinationsite; i < chromosomesize; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize) + i];
                }
                else {
                    gamete[h*chromosomesize + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize) + i];
                }
            }
        } else{
            for (i = 0; i < recombinationsite*elementsperlb; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize*elementsperlb) + i];
                }
                else {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize*elementsperlb) + i];
                }
            }
            for (i = recombinationsite*elementsperlb; i < chromosomesize*elementsperlb; i++) {
                if (startchromosome == 0) {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulationgenomes[startofindividual + totalindividualgenomelength/2 + (h*chromosomesize*elementsperlb) + i];
                }
                else {
                    gamete[h*chromosomesize*elementsperlb + i] = wholepopulationgenomes[startofindividual + (h*chromosomesize*elementsperlb) + i];
                }
            }
        }
    }
}

bool ProduceMutatedGamete(int tskitstatus, int isburninphaseover, tsk_table_collection_t *treesequencetablecollection, tsk_id_t * wholepopulationnodesarray, tsk_id_t * wholepopulationsitesarray, tsk_id_t * childnode, int totaltimesteps, double currenttimestep, int parent, bool isabsolute, int individualgenomelength, double deleteriousmutationrate, double beneficialmutationrate, double Sb, int beneficialdistribution, double Sd, int deleteriousdistribution, double *gamete, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer)
{
    int k, numberofbeneficialmutations, numberofdeleteriousmutations;
    double generatedSb;
    double Sds[30];
    //Following lines stochastically generate a number of deleterious mutations drawn from a Poisson distribution with mean determined by the deleterious mutation rate
    //with effect sizes drawn from a gamma distribution with parameters taken from Kim et al 2017.
    bool DontBreakWhileLoop = false;
    
    while(1){
        DontBreakWhileLoop = false;
        numberofdeleteriousmutations = DetermineNumberOfMutations(deleteriousmutationrate);
        if(!isabsolute){
            for (k = 0; k < numberofdeleteriousmutations; k++) {
                Sds[k] = (gsl_ran_gamma(randomnumbergeneratorforgamma, 0.169, 1327.4)/23646); //Uses parameters for the gamma distribution of the selection coefficients of new mutations scaled to an inferred ancestral populations size. To produce the distribution of unscaled effect sizes, numbers drawn from this distribution must be divided by two times the ancestral population size for the population from which the distribution was derived (11,823 in this case). Data used to produce these fits were samples from 6503 individuals from the National Heart, Lung, and Blood Institute European-American dataset. Analysis of DFE from Kim et al. 2017.
            //This gamma distribution can occasionally produce deleterious mutations with effect sizes larger than 1,
            //which would result in a gamete with fitness less than zero, which would break my algorithm.
            //The following if statement simply throws out any deleterious mutations with effect sizes larger than 1.
                if (Sds[k] >= 1) {
                //this mean a new lethal mutation appeared
                    DontBreakWhileLoop = true;
                    break;
                }
            }
        }else{
            for (k = 0; k < numberofdeleteriousmutations; k++) {
                if(deleteriousdistribution == 0){
                    Sds[k] = Sd;
                }
                else if(deleteriousdistribution == 1)
                    Sds[k] = gsl_ran_exponential(randomnumbergeneratorforgamma, Sd);
                }
            }
            if (!DontBreakWhileLoop) 
                break;
    }

    //Adds the specified number of deleterious mutations to the gamete, recording the sites of each mutation for tree sequence recording.
    //Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the deleterious mutation effect is positive while for relative fitness the sign is negative
    for (k = 0; k < numberofdeleteriousmutations; k++) {
        if (isabsolute){
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, Sds[k]);
        }
        else{
            MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -Sds[k]);
        }
    }
    
    //Following lines stochastically generate a number of beneficial mutations drawn from a Poisson distribution with mean determined by the beneficial mutation rate.
    numberofbeneficialmutations = DetermineNumberOfMutations(beneficialmutationrate);
    //Adds the specified number of beneficial mutations, drawing Sb values from the specified distribution.
    //Sites of each mutation are added to the mutationsites array for tree sequence recording.
    ////Mutation effect sign depends on fitness scheme, for absolute fitness the sign of the beneficial mutation effect is negative while for relative fitness the sign is positive
    //point distribution
    if (beneficialdistribution == 0) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -Sb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, Sb);
            }
        }
    //exponential distribution
    } else if (beneficialdistribution == 1) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            generatedSb = gsl_ran_exponential(randomnumbergeneratorforgamma, Sb);
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -generatedSb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, generatedSb);
            }
        }
    //uniform distribution
    } else if (beneficialdistribution == 2) {
        for (k = 0; k < numberofbeneficialmutations; k++) {
            double upperlimitforuniform = (2 * Sb);
            generatedSb = gsl_ran_flat(randomnumbergeneratorforgamma, 0, upperlimitforuniform);
            if (isabsolute){
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, -generatedSb);
            }
            else{
                MutateGamete(tskitstatus, isburninphaseover, treesequencetablecollection, wholepopulationsitesarray, *childnode, totaltimesteps, currenttimestep, isabsolute, individualgenomelength, gamete, generatedSb);
            }
        }
    } else {
        fprintf(miscfilepointer, "Error: type of distribution for beneficial effect sizes not recognized.");
        exit(0);
    }
    
    return true;
}

int DetermineNumberOfMutations(double mutationrate)
{      
    double meannumberofmutations = mutationrate/2.0;
    
    //Note that because this function operates on gametes, the calculation above appears haploid.
    //There shouldn't be a multiplication by 2 (for diploidy) in this function, since it will be called twice per individual: once per gamete.
    //Above calculation should be moved outside this function for increased efficiency.
    
    int numberofmutations = SampleFromPoisson(meannumberofmutations);
    return numberofmutations;
}

int DetermineMutationSite(int totalgametelength)
{      
    //mutation occur in a random place in the gamete. We can change this to incorporate heterogenous mutation rate for different linckages blocks (mutation hotspots).
    //Moreover, in modular runs it can be edited to allow pleiotropy. A single mutation affects multiple modules.
    int mutatedsite = pcg32_boundedrand(totalgametelength);
    return mutatedsite;
}
//The following function is heavily modified from Numerical Recipes in C, Second Edition.
//For large population sizes, populations with mean Sb > 0 may actually have a more negative fitness slope than mean Sb = 0.
//
int BracketZeroForSb(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double *Sb1, double *Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
    int i, numberoftries;
    numberoftries = 10;
    float factor = 0.01;
    char Sb1name[10], Sb2name[10];
    snprintf(Sb1name, 10, "%.7f", *Sb1);
    snprintf(Sb2name, 10, "%.7f", *Sb2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Sb1name: %s, Sb2name: %s\n", Sb1name, Sb2name);
        fflush(verbosefilepointer);
    }
    float resultingslope1, resultingslope2;
    resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "First two slopes are: %.6f for sb %.6f, and %.6f for sb %.6f\n", resultingslope1, *Sb1, resultingslope2, *Sb2);
        fflush(verbosefilepointer);
    }
    if (resultingslope1 == resultingslope2) {
        fprintf(miscfilepointer, "Slopes after first try are the same, equaling %.5f and %.5f\n", resultingslope1, resultingslope2);
        return 0;
    }
    if (resultingslope1 > slopeforcontourline) {
        fprintf(miscfilepointer, "Slope with sb 0.0 larger than proposed contour, slope = %.6f, contour line value = %.6f\n", resultingslope1, slopeforcontourline);
        return 0;
    }
    
    for (i = 0; i < numberoftries; i++) {
        if ((resultingslope1 < slopeforcontourline) && (resultingslope2 > slopeforcontourline)) {
            return 1;
        } else if (resultingslope2 <= slopeforcontourline) {
            *Sb2 += factor;
            snprintf(Sb2name, 10, "%.7f", *Sb2);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb2name: %s\n", Sb2name);
                fprintf(verbosefilepointer, "Starting run with new sb2 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb2, resultingslope2);
                fflush(verbosefilepointer);
            }
            
        } else if (resultingslope1 >= slopeforcontourline) {
            *Sb1 -= factor;
            snprintf(Sb1name, 10, "%.7f", *Sb1);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "New Sb1name: %s\n", Sb1name);
                fprintf(verbosefilepointer, "Starting run with new sb1 = %.6f\n", *Sb2);
                fflush(verbosefilepointer);
            }
            resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for sb %.6f = %.6f\n", *Sb2, resultingslope2);
                fflush(verbosefilepointer);
            }
            
        }
    }
    fprintf(miscfilepointer, "Failed to bracket contour slope of %.6f in 10 tries.", slopeforcontourline);
    return 0;
}

//The following function is modified from Numerical Recipes in C, Second Edition.
double BisectionMethodToFindSbWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, double * Sb1, double * Sb2, char * Nxtimestepsname, char * popsizename, char * delmutratename, char * chromsizename, char * chromnumname, char * mubname, int typeofrun, int Nxtimesteps, int popsize, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, int deleteriousdistribution, gsl_rng * randomnumbergeneratorforgamma, FILE *miscfilepointer, FILE *verbosefilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
    int i;
    double factor, slope1, slopemid, Sbmid, root;
    double accuracy = 0.00005;
    int maxtries = 30;
    char Sb1name[10], Sb2name[10], Sbmidname[10];
    snprintf(Sb1name, 10, "bis%.4f", *Sb1);
    snprintf(Sb2name, 10, "bis%.4f", *Sb2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Entered bisection function. First two sb %.6f and %.6f\n", *Sb1, *Sb2);
        fprintf(verbosefilepointer, "Starting Sb1name: %s, starting Sb2name: %s", Sb1name, Sb2name);
        fflush(verbosefilepointer);
    }
    slope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb1name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb1, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb1, slope1);
        fflush(verbosefilepointer);
    }
    
    slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sb2name, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, *Sb2, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", *Sb2, slopemid);
    }
    
    if (((slope1 - slopeforcontourline)*(slopemid - slopeforcontourline)) > 0.0) {
        fprintf(miscfilepointer, "Root not bracketed properly, with starting slopes %.10f and %.10f for a desired slope of %.6f\n", slope1, slopemid, slopeforcontourline);
        return 0.0;
    }
    root = (slope1 < slopeforcontourline) ? (factor=*Sb2-*Sb1, *Sb1) : (factor=*Sb1-*Sb2, *Sb2);
    for (i = 1; i <= maxtries; i++) {
        Sbmid = root + (factor *= 0.5);
        snprintf(Sbmidname, 10, "%.7f", Sbmid);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "Sbmidname: %s\n", Sbmidname);
            fprintf(verbosefilepointer, "Starting run with sb %.6f\n", Sbmid);
            fflush(verbosefilepointer);
        }
        slopemid = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, popsizename, delmutratename, chromsizename, chromnumname, mubname, Sbmidname, typeofrun, Nxtimesteps, popsize, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sbmid, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "Finished run with sb %.6f, resulting in a slope of %.6f\n", Sbmid, slopemid);
            fflush(verbosefilepointer);
        }
        if (slopemid <= slopeforcontourline) {
            root = Sbmid;
        }
        if (fabs(factor) < accuracy || Sbmid == slopeforcontourline) {
            return root;
        }
        
    }
    fprintf(miscfilepointer, "Error: root not found. Root after 30 tries was: %.10f", root);
    fprintf(finaldatafilepointer, "Error: root not found. Root after 30 tries was: %.10f", root);
    return 0.0;
    
}

//This is the second bracketing function to find Ncrit which is the population size where the net fitness flux of 0 for a given set of parameters.
int BracketZeroForNcrit(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, int *N1, int *N2, char *Nxtimestepsname, char *delmutratename, char *chromsizename, char *chromnumname, char *mubname, int typeofrun, int Nxtimesteps, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, char *Sbname, double Sb, int deleteriousdistribution, gsl_rng *randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
    //setting up variables
    int i, numberoftriesforlower, numberoftriesforupper, generations;
    double resultingslope1, resultingslope2;
    int gen = 10;//scale number of generations based on tested N
    //number of attempts of finding bounds
    numberoftriesforlower = 3;
    numberoftriesforupper = 3;
    //factor to increase or decrease the next attempted N from the initial unsuccessful guesses
    int factor = 2;

    char N1name[6], N2name[6];
    snprintf(N1name, 6, "%d", *N1);
    snprintf(N2name, 6, "%d", *N2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "N1name: %s, N2name: %s\n", N1name, N2name);
        fflush(verbosefilepointer);
    }

    generations = *N1 * gen;
    fprintf(verbosefilepointer, "generations for N1 %d: %d\n", *N1, generations);
    fflush(verbosefilepointer);
    //testing lower bound N1
    resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N1name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, *N1, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    
    generations = *N2 * gen;
    fprintf(verbosefilepointer, "generations for N2 %d: %d\n", *N2, generations);
    fflush(verbosefilepointer);
    //testing higher bound N2
    resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N2name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, *N2, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);

    //different print statments based on results of initial guesses
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "First two slopes are: %.6f for N %d, and %.6f for N %d\n", resultingslope1, *N1, resultingslope2, *N2);
        fflush(verbosefilepointer);
    }
    if (resultingslope1 == resultingslope2) {
        fprintf(miscfilepointer, "Slopes after first try are the same, equaling %.6f and %.6f\n", resultingslope1, resultingslope2);
        return -1;
    }
    if (resultingslope1 > slopeforcontourline) {
        fprintf(miscfilepointer, "Slope with inital lower bound N1 %d larger than proposed contour, slope = %.6f, contour line value = %.6f. Lower bound N should have slope < contour line\n",*N1, resultingslope1, slopeforcontourline);
    }
    if (resultingslope2 < slopeforcontourline) {
        fprintf(miscfilepointer, "Slope with inital upper bound N2 %d lower than proposed contour, slope = %.6f, contour line value = %.6f. Upper bound N should have slope > contour line\n",*N2, resultingslope2, slopeforcontourline);
    }
    
    if ((resultingslope1 < slopeforcontourline) && (resultingslope2 > slopeforcontourline)) {
        return 0;//if guesses are successful, everything is good to go 
    } else if (resultingslope2 <= slopeforcontourline) {
        //this block of code will attempt to find a new upper bound
        for (i = 0; i < numberoftriesforupper; i++) {
            //if initial upper bound guess is unsuccessful
            *N1 = *N2;//if the first guess of an upper bound is wrong, then we make it the new lower bound
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Previous upper bound N2 is now new lower bound N1 : %s\n", N2name);
                fflush(verbosefilepointer);
            }
            *N2 = *N2 * factor;//calculate new upper limit guess by multiplying old guess by preset factor
            generations = *N2 * gen;//set number of generations
            snprintf(N2name, 6, "%d", *N2);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Testing new upper bound N2: %s\n", N2name);
                fprintf(verbosefilepointer, "Starting run with new N2 = %d\n", *N2);
                fflush(verbosefilepointer);
            }
            //calculating slope of new tested guess of upper limit
            resultingslope2 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N2name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, *N2, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for new upper bound N2 %d = %.6f\n", *N2, resultingslope2);
                fflush(verbosefilepointer);
            }
            if (resultingslope2 > slopeforcontourline)
            {
                //if a new guess of upper limit was found
                fprintf(verbosefilepointer, "new upper bound N2 found %d\n", *N2);
                fflush(verbosefilepointer);
                break;
            } else if (i == numberoftriesforupper && resultingslope2 <= slopeforcontourline){
                //if no new guess of upper limit is found within the max number of tries
                fprintf(verbosefilepointer, "No new upper bound N2 found after %d attempts\n", numberoftriesforupper);
                fflush(verbosefilepointer);
                return -1;
            }    
        }          
    } else if (resultingslope1 >= slopeforcontourline) {
        //this block of code will attempt to find a new lower bound
        for (i = 0; i < numberoftriesforlower; i++) {
            //if initial lower bound guess is unsuccessful
            *N2 = *N1;//if the first guess of a lower bound is wrong, then we make it the new upper bound
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Previous lower bound N1 is now new upper bound N2 : %s\n", N1name);
                fflush(verbosefilepointer);
            }
            *N1 = *N1 / factor;//calculate new lower limit guess by dividing old guess by preset factor
            generations = *N1 * gen;//set number of generations
            snprintf(N1name, 6, "%d", *N1);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Testing new lower bound N1: %s\n", N1name);
                fprintf(verbosefilepointer, "Starting run with new N1 = %d\n", *N1);
                fflush(verbosefilepointer);
            }
            //calculating slope of new tested guess of upper limit
            resultingslope1 = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N1name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, Nxtimesteps, *N1, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
            if (VERBOSE == 1) {
                fprintf(verbosefilepointer, "Slope for new lower bound N1 %d = %.6f\n", *N1, resultingslope1);
                fflush(verbosefilepointer);
            }
            if (resultingslope1 < slopeforcontourline)
            {
                //if a new guess of upper limit was found
                fprintf(verbosefilepointer, "new lower bound N1 found %d\n", *N1);
                fflush(verbosefilepointer);
                break;
            } else if (i == numberoftriesforlower && resultingslope1 >= slopeforcontourline){
                //if no new guess of upper limit is found within the max number of tries
                fprintf(verbosefilepointer, "No new lower bound N1 found after %d attempts\n", numberoftriesforlower);
                fflush(verbosefilepointer);
                return -1;
            } 
        }  
    }
}

//The following function is modified from Numerical Recipes in C, Second Edition.
double SecantMethodToFindNcritWithZeroSlope(int tskitstatus, bool isabsolute, bool ismodular, int elementsperlb, int *N1, int *N2, char *Nxtimestepsname, char *delmutratename, char * chromsizename, char *chromnumname, char *mubname, int typeofrun, int Nxtimesteps, int chromosomesize, int numberofchromosomes, double deleteriousmutationrate, double beneficialmutationrate, double slopeforcontourline, int beneficialdistribution, double Sd, char *Sbname, double Sb, int deleteriousdistribution, gsl_rng *randomnumbergeneratorforgamma, FILE *verbosefilepointer, FILE *miscfilepointer, FILE *finaldatafilepointer, FILE *veryverbosefilepointer, int rawdatafilesize) {
    int i, Nl, root, generations;
    int gen = 10;
    double *N = NULL;
    double *v = NULL;   
    // Variable to store the number of elements in the array
    int num_elements = 0;
    int filtered_count = 0;
    double slope, slopel, slopetemp;
    //We choose a threshold of 25% change in slope. This choice is somewhat arbitrary
    double percent_threshold = 25.0;
    //Number of maximum tries to find Ncrit 
    int maxtries = 30;
    char N1name[6], N2name[6], rootname[6];
    snprintf(N1name, 6, "%d", *N1);
    snprintf(N2name, 6, "%d", *N2);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Entered secant function. First two N %d and %d\n", *N1, *N2);
        fprintf(verbosefilepointer, "Starting N1name: %s, starting N2name: %s\n", N1name, N2name);
        fflush(verbosefilepointer);
    }
    //number of generations is a function of the attempted census pop size
    generations = *N1 * gen;
    //calculate the slope of fitness change as a function of the lower pop size
    slopel = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N1name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, *N1, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with N1 %d, resulting in a slope of %.6f\n", *N1, slopel);
        fflush(verbosefilepointer);
    }

    // Dynamically reallocate memory for the array to store the new number
    num_elements++;
    N = (double *)realloc(N, num_elements * sizeof(double));
    v = (double *)realloc(v, num_elements * sizeof(double));

    // Check if memory allocation was successful
    if (v == NULL || N == NULL) {
        printf("Memory allocation failed\n");
        exit(1); // Exit the program due to memory allocation failure
    }

    // Store the input number in the array
    N[num_elements - 1] = (double) *N1;
    v[num_elements - 1] = slopel;
    
    generations = *N2 * gen;
    //calculate the slope of fitness change as a function of the higher pop size
    slope = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N2name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, *N2, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
    if (VERBOSE == 1) {
        fprintf(verbosefilepointer, "Finished run with N2 %d, resulting in a slope of %.6f\n", *N2, slope);
    }
    
    num_elements++;
    N = (double *)realloc(N, num_elements * sizeof(double));
    v = (double *)realloc(v, num_elements * sizeof(double));

    // Check if memory allocation was successful
    if (v == NULL || N == NULL) {
        printf("Memory allocation failed\n");
        exit(1); // Exit the program due to memory allocation failure
    }

    // Store the input number in the array
    N[num_elements - 1] = (double) *N2;
    v[num_elements - 1] = slope;
    //check if the bracketing was not successful. This is just a double check out of caution
    if (((slopel - slopeforcontourline)*(slope - slopeforcontourline)) > 0.0) {
        fprintf(miscfilepointer, "Root not bracketed properly, with starting slopes %.10f and %.10f for a desired slope of %.6f\n", slopel, slope, slopeforcontourline);
        return 0.0;
    }

    //Secant calculation here
    if (abs(slopel) < abs(slope)){
        root = *N1;
        Nl = *N2;
        slopetemp = slopel;
        slopel = slope;
        slope = slopetemp;
    } else{
        Nl = *N1;
        root = *N2;
    }
    
    for (i = 1; i <= maxtries; i++) {
        int dx = (Nl-root)*slope/(slope-slopel);
        Nl=root;
        slopel=slope;
        root += dx;
        snprintf(rootname, 6, "%d", root);
        if (VERBOSE == 1) {
            fprintf(verbosefilepointer, "rootname: %s\n", rootname);
            fprintf(verbosefilepointer, "Starting run with N %d\n", root);
            fflush(verbosefilepointer);
        }

        if (root == Nl) {
            fprintf(miscfilepointer, "Ncrit found based on secant method converging on same N: %d\n", root);
            break;
        }
        generations = root * gen;

        slope = RunSimulationRel(tskitstatus, isabsolute, ismodular, elementsperlb, Nxtimestepsname, N2name, delmutratename, chromsizename, chromnumname, mubname, Sbname, typeofrun, generations, root, chromosomesize, numberofchromosomes, deleteriousmutationrate, beneficialmutationrate, Sb, beneficialdistribution, Sd, deleteriousdistribution, randomnumbergeneratorforgamma, miscfilepointer, veryverbosefilepointer, rawdatafilesize);
        
        fprintf(verbosefilepointer, "Finished run with N %d, resulting in a slope of %.6f\n", root, slope);
        fflush(verbosefilepointer);
        
        num_elements++;
        N = (double *)realloc(N, num_elements * sizeof(double));
        v = (double *)realloc(v, num_elements * sizeof(double));

        // Check if memory allocation was successful
        if (v == NULL || N == NULL) {
            printf("Memory allocation failed\n");
            exit(1); // Exit the program due to memory allocation failure
        }

        // Store the input number in the array
        N[num_elements - 1] = (double) root;
        v[num_elements - 1] = slope;
        
        double relative_change = fabs((slope - slopel) / slopel) * 100.0;

        if (relative_change < percent_threshold) {
            fprintf(miscfilepointer, "Ncrit found based on secant method due to relative change being less than %.1f: %d\n", percent_threshold, root);
            break;
        }
    }
    // Create arrays to store filtered N and v
    double *filtered_N = malloc(num_elements * sizeof(double));
    double *filtered_v = malloc(num_elements * sizeof(double));

    if (filtered_N == NULL || filtered_v == NULL) {
        printf("Memory allocation failed\n");
        exit(1); // Exit the program due to memory allocation failure
    }

    // Iterate through recorded data
    for (int i = 0; i < num_elements; ++i) {
        // Check if N is within a factor of 3 from the root
        if (N[i] >= root / 3 && N[i] <= root * 3) {
            // Store the pair in filtered arrays
            filtered_N[filtered_count] = N[i];
            filtered_v[filtered_count] = v[i];
            filtered_count++;
        }
    }
    // Free memory for the original N and v arrays
    free(N);
    free(v);

    // Resize the arrays to the actual number of valid elements
    filtered_N = realloc(filtered_N, filtered_count * sizeof(double));
    filtered_v = realloc(filtered_v, filtered_count * sizeof(double));

    if (filtered_N == NULL || filtered_v == NULL) {
        printf("Memory reallocation failed\n");
        exit(1); // Exit the program due to memory reallocation failure
    }
    // defining variables related to linear model fit
    size_t step = 1;
    double cov00, cov01, cov11, sumsq, Ncrit;
    double b, m;
    //The following function fits a linear model to the two variables (popsize and slope of log fitness)
    //and records the parameters of the best-fitting linear model in the intercept, cov00, cov01, cov11, sumsq, and slope variables.
    //I use the slope and intercept parameters,  to find the root as -intercept/slope.
    gsl_fit_linear(filtered_N, step, filtered_v, step, filtered_count, &b, &m, &cov00, &cov01, &cov11, &sumsq);
    fprintf(miscfilepointer, "the slope of the linear curve is %.10f and the intercept is %.10f\n", b, m);
    fflush(miscfilepointer);
    Ncrit = -b/m;
    fprintf(miscfilepointer, "Ncrit found based on linear curve fitting method: %.2f", Ncrit);
    free(filtered_N);
    free(filtered_v);
    return Ncrit;
}
//The next few functions serve to aggregate arguments to build the name of the directory or raw files names for the particular simulation run
char * MakeDirectoryName(char * tskitstatus, char* deldist, char * isabsolutename, bool isabsolute, char * bendist, char * benmut, char * numberofchromosomes, char * chromosomesize, char * popsize, char * delmut, char * randomnumberseed, char * K, char * r, char *i_init, char * s, bool ismodular, char *elementsperlb, char *iscalcfixationname, int typeofrun, char * Sbname) 
{
	
	char * directoryname = (char *) malloc(400);
	strcpy(directoryname, "datafor_");
	strcat(directoryname, isabsolutename);
    strcat(directoryname, "_tskitstatus_");
	strcat(directoryname, tskitstatus);
    strcat(directoryname, "_fixationcalc_");
	strcat(directoryname, iscalcfixationname);
    if(ismodular){
        strcat(directoryname, "_m_");
        strcat(directoryname, elementsperlb);
    }
    if(isabsolute){
        strcat(directoryname, "_r_");
        strcat(directoryname, r);
        strcat(directoryname, "_iinit_");
        strcat(directoryname, i_init);
        strcat(directoryname, "_s_");
        strcat(directoryname, s);
        strcat(directoryname, "_K_");
        strcat(directoryname, K);
    }
    if(typeofrun == 1){
        strcat(directoryname, "_Sb_");
        strcat(directoryname, Sbname);
    }
    if(typeofrun == 2){
        strcat(directoryname, "_Ncrit");
        strcat(directoryname, "_Sb_");
        strcat(directoryname, Sbname);
    }
    strcat(directoryname, "_deldist_");
    strcat(directoryname, deldist);
    strcat(directoryname, "_bendist_");
	strcat(directoryname, bendist);
	strcat(directoryname, "_mub_");
	strcat(directoryname, benmut);
	strcat(directoryname, "_chromnum_");
	strcat(directoryname, numberofchromosomes);
	strcat(directoryname, "_N0_");
	strcat(directoryname, popsize);
	strcat(directoryname, "_mud_");
	strcat(directoryname, delmut);
    strcat(directoryname, "_L_");
	strcat(directoryname, chromosomesize);
	strcat(directoryname, "_seed_");
	strcat(directoryname, randomnumberseed);

	return directoryname;
}

char * MakeFinalDataFileName(char * typeofrun, char * benmut, char * slopeforcontourline, char * randomnumberseed) 
{
	
    char * finaldatafilename = (char *) malloc(60);
    strcpy(finaldatafilename, "finaldatafor_");
    strcat(finaldatafilename, "runtype_");
    strcat(finaldatafilename, typeofrun);
    strcat(finaldatafilename, "_mub_");
    strcat(finaldatafilename, benmut);
    strcat(finaldatafilename, "_slope_");
    strcat(finaldatafilename, slopeforcontourline);
    strcat(finaldatafilename, "_seed_");
    strcat(finaldatafilename, randomnumberseed);
    
    return finaldatafilename;
}

char * MakeRawDataFileName(char * mubname, char * Sbname, bool isredinmaxpopsize, char *redinmaxpopsizename) 
{
	
    char * rawdatafilename = (char *) malloc(200);
    strcpy(rawdatafilename, "rawdatafor_"); //starting the string that will be the name of the data file.
    strcat(rawdatafilename, "Sb_");
    strcat(rawdatafilename, Sbname);
    strcat(rawdatafilename, "_mub_");
    strcat(rawdatafilename, mubname);
    if(isredinmaxpopsize){
        strcat(rawdatafilename, "_redrate_");
        strcat(rawdatafilename, redinmaxpopsizename);
    }
    strcat(rawdatafilename, ".txt");
    
    return rawdatafilename;
}
char * MakeSummaryDataFileName(char * mubname, char * Sbname) 
{
	
    char* summarydatafilename = (char*)malloc(60);
    strcpy(summarydatafilename, "summarydatafor_");
    strcat(summarydatafilename, "Sb_");
    strcat(summarydatafilename, Sbname);
    strcat(summarydatafilename, "_mub_");
    strcat(summarydatafilename, mubname);
    strcat(summarydatafilename, ".txt");
    
    return summarydatafilename;
}

char * MakePopSnapshotFileName(char * mubname, char * Sbname) 
{
	
    char * popsnapshotfilename = (char *) malloc(100);
    strcpy(popsnapshotfilename, "popsnapshotfor_");
    strcat(popsnapshotfilename, "Sb_");
    strcat(popsnapshotfilename, Sbname);
    strcat(popsnapshotfilename, "mub_");
    strcat(popsnapshotfilename, mubname);
    strcat(popsnapshotfilename, ".txt");
    
    return popsnapshotfilename;
}

int AssignArgumentstoVar(char **argv, int *Nxtimesteps, char *Nxtimestepsname, int *popsize, char *popsizename, double *deleteriousmutationrate, char *deleteriousmutationratename, int *chromosomesize, char *chromosomesizename, int *numberofchromosomes, char *numberofchromosomesname, double *bentodelmutrate, double *Sbtemp, int *beneficialdistribution, int *typeofrun, double *slopeforcontourline, char *slopeforcontourlinename, int *randomnumberseed, char *randomnumberseedname, int *K, char *Kname, int *relorabs, double *r, char *rname, int *i_init, char *i_initname, double *s, char *sname, int *tskitstatus, int *nonmodormod, int *elementsperlb, char *elementsperlbname, int *snapshot, char *prevsnapshotfilename, double *SdtoSbratio, char *SdtoSbrationame, int *deleteriousdistribution, int *rawdatafilesize, double *redinmaxpopsize, char *redinmaxpopsizename, int *calcfixation) {
	//whicharg is used for simplification in future modifications at the order of the arguments in the code
	int whicharg;
	whicharg=1;
	//number of time steps that the simulation last.
	*Nxtimesteps = atoi(argv[whicharg]);
	strcpy(Nxtimestepsname, argv[whicharg]);
	whicharg=2;
	//initial population size
	*popsize = atoi(argv[whicharg]);
	strcpy(popsizename, argv[whicharg]);
	whicharg=3;
	//deleterious mutation rate; remember that this is the genome-wide mutation rate
	*deleteriousmutationrate = atof(argv[whicharg]);
	strcpy(deleteriousmutationratename, argv[whicharg]);
	whicharg=4;
	//size of a single chromosome, i.e. number of blocks within a single chromosome. so far, we use 50 blocks
	*chromosomesize = atoi(argv[whicharg]);
	strcpy(chromosomesizename, argv[whicharg]);
	whicharg=5;
	//number of chromosomes; remember that this is total number of chromosomes, not ploidy -- all individuals will be diploid; to mimic a human genome, we use 23 chromosomes
	*numberofchromosomes = atoi(argv[whicharg]);
	strcpy(numberofchromosomesname, argv[whicharg]);
	whicharg=6;
	//given that we don't have an exact idea what the beneficial mutation rate is in humans, we use instead a beneficial/deleterious mutation ratio.
	*bentodelmutrate = atof(argv[whicharg]);
	whicharg=7;
	//Mean effect of a beneficial mutation in relative runs. For absolute runs Sb is set to 1.0 in the main function
	*Sbtemp = atof(argv[whicharg]);
	whicharg=8;
	//0 for point; 1 for exponential; 2 for unifrom
	*beneficialdistribution = atoi(argv[whicharg]);
	whicharg=9;
	//0 is root; 1 is single
	*typeofrun = atoi(argv[whicharg]);
	whicharg=10;
	//the slope for the contour line; so far it is used in the relative scheme
	*slopeforcontourline = atof(argv[whicharg]);
	strcpy(slopeforcontourlinename, argv[whicharg]);
	whicharg=11;
	// the seed for the random number generators used in this program *don't forget to be consistent with the seed*
	*randomnumberseed = atoi(argv[whicharg]);
	strcpy(randomnumberseedname, argv[whicharg]);
	whicharg=12;
	//the carrying capacity in a deterministic population with minimal load. That is, the equilibrium population size for a homogeneous population with minimal load
	*K = atoi(argv[whicharg]);
	strcpy(Kname, argv[whicharg]);
	whicharg=13;
	//0 for relative simulation; 1 for absolute
	*relorabs = atoi(argv[whicharg]);
	whicharg=14;
	//the epistasis variable R; R varies between 0 and 1, 1 representing no epistasis and 0 is full epistasis
	*r = atof(argv[whicharg]);
	strcpy(rname, argv[whicharg]);
	whicharg=15;
	//the number of mutations away from worst-case scenario at the beginning of the simulation
	*i_init = atof(argv[whicharg]);
	strcpy(i_initname, argv[whicharg]);
	whicharg=16;
	//the selection coefficient of a beneficial mutation at minimum death rate (i.e. max birth rate). so far we use the arbitrary value 0.01
	*s = atof(argv[whicharg]);
	strcpy(sname, argv[whicharg]);
	whicharg=17;
	//status of tskit API in program, value of 0 tree sequencing is OFF, value of 1 tree sequencing is ON fully, and value of 2 tree sequencing is on after burn-in is completed.
	*tskitstatus = atoi(argv[whicharg]);
	whicharg=18;
	//Indicate whether epistasis is modular or not. 0 for non modular epistasis and 1 for modular epistasis.
	*nonmodormod = atoi(argv[whicharg]);
	whicharg=19;
	//Number of elements per linkage block
	*elementsperlb = atoi(argv[whicharg]);
	strcpy(elementsperlbname, argv[whicharg]);
	whicharg=20;
	//input for snapshot implementation
	*snapshot = atoi(argv[whicharg]);
	whicharg=21;
	//File name of the previous snapshot
	strcpy(prevsnapshotfilename, argv[whicharg]);
	whicharg=22;
	//given that we don't have an exact idea what the mean effect of beneficial mutations are in humans and absolute runs, we use instead a beneficial/deleterious mutation effect ratio. sdmin then scales these parameters
	*SdtoSbratio = atof(argv[whicharg]);
	strcpy(SdtoSbrationame, argv[whicharg]);
	whicharg=23;
	//DFE of deleterious mutations. 0 for point; 1 for exponential
	*deleteriousdistribution = atoi(argv[whicharg]);
	whicharg=24;
	//input the requested size of the raw data file which determines the sampling rate of datapoints
	*rawdatafilesize = atoi(argv[whicharg]);
    whicharg=25;
    //input for reduction in MaxPopSize
    *redinmaxpopsize = atof(argv[whicharg]);
    strcpy(redinmaxpopsizename, argv[whicharg]);
    whicharg=26;
	//status of fixation calcualtion which requires fixation calculation algorithm to run before simplification, value of 0 fixation calculation is OFF, value of 1 fixation calculation ON.
	*calcfixation = atoi(argv[whicharg]);

	if (*beneficialdistribution > 2) {
		printf("[Error] 8th argument (beneficial distribution) must be 0 (point), 1 (exponential) or 2 (uniform). However, by the moment only point mutations are tested for mutaway code \n");
		return -1;
	}
	if (*typeofrun != 0 && *typeofrun != 1 && *typeofrun != 2) {
		printf("[Error] 9th argument (type of run) must be 0 or 1 or 2\n");
		return -1;
	}
	if (*relorabs != 0 && *relorabs != 1) {
		printf("[Error] 13th argument (relative or absolute scheme) must be 0 or 1 \n");
		return -1;
	}
	if (*r < 0 || *r > 1) {
		printf("[Error] 14th argument (epistasis degree R) must be between 0 and 1 \n");
		return -1;
	}
	if (*tskitstatus != 0 && *tskitstatus != 1 && *tskitstatus != 2) {
		printf("[Error] 17th argument (tree sequencing API status) must be either 0 for OFF, 1 for ON, 2 for ON after burn-in \n");
		return -1;
	}
    if (*calcfixation != 0 && *calcfixation != 1) {
		printf("[Error] 26th argument (fixation calculation status) must be either 0 for OFF, 1 for ON \n");
		return -1;
	}
	if (*nonmodormod != 0 && *nonmodormod != 1) {
		printf("[Error] 18th argument (non modular or modular epistasis) must be 0 or 1 \n");
		return -1;
	}
	if (*snapshot != 0 && *snapshot != 1) {
		printf("[Error] 20th argument (is there or not a previous snapshot) must be 0 or 1 \n");
		return -1;
	}
	/* if(*SdtoSbratio != 1.0) {
		printf("[Error] Beaware that this code was originally written for 22nd argument (SdtoSbratio) equal 1.0. Sb = Sd = 1.0. Check that you modify the code before changing this parameter \n");
		return -1;
	} */
	if (*deleteriousdistribution != 0 && *deleteriousdistribution != 1) {
		printf("[Error] 23th argument (deleterious distribution) must be 0 (point) or 1 (exponential). However, by the moment only point mutations are tested for mutaway code \n");
		return -1;
	}
	if (*rawdatafilesize == 0 || *rawdatafilesize > *Nxtimesteps) {
		printf("[Error] 24th argument (size of raw data file) cannot be 0 or larger than # of generations \n");
		return -1;
	}
	if(*r == 1.0 && *bentodelmutrate != 0.0) {
		printf("[Error] Trying to run a simulation with no epistasis and beneficial mutations. This might cause contiuouns pop growth and memory corruption \n");
		return -1;
	} 
	if(*r == 1.0 && *popsize >= *K/((*s) * (*i_init))) {
		printf("[Error] Trying to run a simulation with no epistasis with popsize bigger than kappa. That is a negative per capita birth rate \n");
		return -1;
	}
	if(*r != 1.0) {
        if (*deleteriousdistribution == 0 && *popsize >= (*K)*(1-*r)/(*s)) {
		printf("[Error] Trying to run a simulation with epistasis with popsize bigger than kappa. That is a negative per capita birth rate \n");
		return -1;
	    }else if(*deleteriousdistribution == 1 && *popsize >= -((*K)*log(*r))/(*s)) {
		printf("[Error] Trying to run a simulation with epistasis with popsize bigger than kappa. That is a negative per capita birth rate \n");
		return -1;
    }
	}
	//check that modular elements per linkage block is different than 0 when modular epistasis is used
	if(*nonmodormod == 1 && *elementsperlb == 0) {
		printf("[Error] 19th argument (elements per linkage block) must be other than 0 when 18th argument (modular epistasis) is 1 \n");
		return -1;
	}
    if(*redinmaxpopsize < 0.0){
        printf("[Error] 25th argument (reduction in MaxPopSize per generation) must be a positive integer \n");
        return -1;
    }
    if (*tskitstatus == 0 && *calcfixation != 0) {
		printf("[Error] Fixation calculation is ON while Tree Sequence Recording status is OFF \n");
		return -1;
	}
	return 1;
}
void AssignStringNames(char *beneficialmutationratename, double beneficialmutationrate, char *bendistname, int beneficialdistribution, char *deldistname, int deleteriousdistribution, char *typeofrunname, int typeofrun, char *tskitstatusname, int tskitstatus, char* Sb2name, double Sb2, char *isabsolutename, bool isabsolute, char *iscalcfixationname, bool iscalcfixation) {
	//pointer for the beneficial mutation rate name
	sprintf(beneficialmutationratename, "%1.4f", beneficialmutationrate);

	//pointer for the beneficial distribution name
	if(beneficialdistribution == 0)
		strncpy(bendistname, "point", sizeof("point"));
	else if(beneficialdistribution == 1)
		strncpy(bendistname, "exponential", sizeof("exponential"));
	else if(beneficialdistribution == 2)
		strncpy(bendistname, "uniform", sizeof("uniform"));

	//pointer for the deleterious distribution name
	if(deleteriousdistribution == 0)
		strncpy(deldistname, "point", sizeof("point"));
	else if(deleteriousdistribution == 1)
		strncpy(deldistname, "exponential", sizeof("exponential"));

	//pointer for type of run name
	if(typeofrun == 0)
		strncpy(typeofrunname, "rootSb", sizeof("rootSb"));
	else if(typeofrun == 1)
		strncpy(typeofrunname, "single", sizeof("single"));
    else if(typeofrun == 2)
		strncpy(typeofrunname, "rootNcrit", sizeof("rootNcrit"));

	//pointer for tskitstatus name
	if(tskitstatus == 0) {
		strncpy(tskitstatusname, "OFF", sizeof("OFF"));
	} else if(tskitstatus == 1) {
		strncpy(tskitstatusname, "ON", sizeof("ON"));
	} else if(tskitstatus == 2) {
		strncpy(tskitstatusname, "ON_BURNIN", sizeof("ON_BURNIN"));
	}

	//pointer for Sb name
	sprintf(Sb2name, "%1.4f", Sb2);

	//pointer for isabsolutename
	if(isabsolute)
		strncpy(isabsolutename, "absolute", sizeof("absolute"));
	else
		strncpy(isabsolutename, "relative", sizeof("relative"));

    if(iscalcfixation)
		strncpy(iscalcfixationname, "ON", sizeof("ON"));
	else
		strncpy(iscalcfixationname, "OFF", sizeof("OFF"));
}
