/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<stdio.h>
#include<float.h>
#include<string.h>
#include "pcg_basic.h"
#include<math.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_fit.h>

#define VERBOSE 0
#define VERYVERBOSE 0
#define MISCELLANEOUS 1
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654
#define RANDOMNUMBERSEED 10020
//Random number seed could eventually be an argument input by user.

FILE	*verbosefilepointer;
FILE	*veryverbosefilepointer;
FILE	*miscfilepointer;
FILE	*datafilepointer;

void UpdateLastTwentyGenerations(double *, double);
void DoubleSwap(long double *, long double *);
void DoubleBubbleSort(long double *, int);
double CalculateVarianceInLogFitness(int, long double *, long double);
long double FindFittestWi(long double *, int);
long double Fen_sum(long double *, int);
void Fen_add(long double *, int, long double, int);
long double Fen_range(long double *, int, int);
long double Fen_get(long double *, int);
void Fen_set(long double *, int, long double, int);
void InitializePopulation(long double *, long double *, int , double *, int);
//int ChooseVictimWithoutTree(long double *, int , long double);
int SearchTree(int, int, long double, long double *);
int ChooseParentWithTree(long double *, int, long double);
int ChooseVictim(int);
void RecombineChromosomesIntoGamete(int, int , int , double *, double *, int);
int SampleFromPoisson(float );
void MutateGamete(int , int , double *, float, double);
double CalculateWi(int , int , double *, double *, double , int);
void ReplaceVictim(double *, double *, int, int, int , int , double, long double *, double *, int, long double *, long double *);
//void RecalculateSumOfDeathRates(long double *, int , int , int , struct individual *);

void main(int argc, char *argv[]) {

    printf("\nEntering main function");
    fflush(stdout); //for bug checking

    double Sd;
    Sd = atof(argv[1]);

    int generations;
    generations = atoi(argv[2]);

    int popsize;
    popsize = atoi(argv[3]);

    double deleteriousmutationrate;
    deleteriousmutationrate = atof(argv[4]); //remember that this is the per-locus deleterious mutation rate, not the genome-wide mutation rate.
    
    int chromosomesize;
    chromosomesize = atoi(argv[5]);

    int numberofchromosomes;
    numberofchromosomes = atoi(argv[6]); //remember that this is total number of chromosomes, not ploidy -- all individuals will be diploid.

    double beneficialmutationrate;
    beneficialmutationrate = atof(argv[6]); //remember that this is the per-locus rate, not genome-wide.
    
    double Sb;
    Sb = atof(argv[7]);
    
    int i, j, k;
    double *wholepopulationgenomes;
    int totalpopulationgenomelength;
    int totalindividualgenomelength;
    totalpopulationgenomelength = popsize * numberofchromosomes * 2 * chromosomesize;
    totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
    wholepopulationgenomes = malloc(sizeof(double) * totalpopulationgenomelength);
    long double sumofwis;
    long double *psumofwis = &sumofwis;
    long double *wholepopulationwistree;
    wholepopulationwistree = malloc(sizeof(long double) * popsize);
    int currentparent1;
    int currentparent2;
    int currentvictim;
    long double currentfittestindividualswi;
    double parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
    
    //Following array is to store the variance in log(fitness) for twenty generations at a time for estimating the end of the burn-in phase.
    //The choice of twenty is completely arbitrary and could eventually be an input if it seems worth it.
    double *lasttwentygenerationsvariance;
    double *literallyjustlasttwentygenerations;
    literallyjustlasttwentygenerations = malloc(sizeof(double) * 20);
    lasttwentygenerationsvariance = malloc(sizeof(double) * 20);
    for (k = 0; k < 20; k++) {
        literallyjustlasttwentygenerations[k] = 0.0;
        lasttwentygenerationsvariance[k] = 0.0;
    }
    double slopeofvariance;
    int isburninphaseover = 0;
    int endofburninphase;
    double arbitrarynumber;
    arbitrarynumber = -0.001 * Sd;
    
    double variancesum;
        
    /*Initializes the population as an array of individuals, 
    the current victim as a pointer to an individual,
    and the current parents as pointers to individuals.*/
        
    long double *wholepopulationwisarray;
    wholepopulationwisarray = malloc(sizeof(long double) * popsize);
    //The Fenwick tree does not store each individual's wi, but rather a collection of partial sums.
    //For debugging purposes and data that requires summarizing wis, storing the wis in an array is necessary.

    long double *sortedwisarray;
    sortedwisarray = malloc(sizeof(long double) * popsize);
    //In order to visualize the distribution of fitness in the population,
    //the individual fitnesses need to be sorted, which requires a separate array.
    
    sumofwis = (long double) popsize;
        
    InitializePopulation(wholepopulationwistree, wholepopulationwisarray, popsize, wholepopulationgenomes, totalpopulationgenomelength);
    /*Sets the initial population to have zeroes in all their linkage blocks,
    death rates equal to the baseline wi, and an identifier number.
    It also sums the wis and returns the sum.*/
    
    pcg32_srandom(RANDOMNUMBERSEED, RANDOMNUMBERSEED); // seeds the random number generator.

    verbosefilepointer = fopen("verbose.txt", "w");	//opens the file to which to print verbose data.
    veryverbosefilepointer = fopen("veryverbose.txt", "w"); //opens the file to which to print very verbose data.
    miscfilepointer = fopen("miscellaneous.txt", "w"); //opens the file to which to print miscellaneous data.

    char * datafilename = (char *) malloc(200);
    strcpy(datafilename, "datafor"); //starting the string that will be the name of the data file.
	
    strcat(datafilename, "Sd"); //for adding values of Sd to the data name.
    strcat(datafilename, argv[1]);

    strcat(datafilename, "gens"); //for adding values of generations to the data name.
    strcat(datafilename, argv[2]);

    strcat(datafilename, "popsize"); //for adding values of starting population sizes to the data name.
    strcat(datafilename, argv[3]);

    strcat(datafilename, "mutrate"); //for adding values of mutation rate to the data name (remember that mutation rate is currently the per-locus rate, not per-genome).
    strcat(datafilename, argv[4]);

    strcat(datafilename, "chromsize"); //for adding values of chromosome size to the data name.
    strcat(datafilename, argv[5]);

    strcat(datafilename, "chromnum"); //for adding values of the number of chromosomes to the data name.
    strcat(datafilename, argv[6]);

    datafilepointer = fopen(datafilename, "w"); //opens the file to which to print data to be plotted.
    fprintf(datafilepointer, "Generation,Sum.of.wis,Variance.in.log.fitness\n");
    
    //Following code sets up a separate data file to examine the distribution of fitness in populations.
    /*
     char * fitnessdistributiondatafilename = (char *) malloc(100);
     strcpy(fitnessdistributiondatafilename, "fitnessdistributiondataforSd");
     strcat(fitnessdistributiondatafilename, argv[1]);
     strcat(fitnessdistributiondatafilename, ".txt");
     fitnessdistributiondatafilepointer = fopen(fitnessdistributiondatafilename, "w");
     fprintf(fitnessdistributiondatafilepointer, "Individual.number,Individual.fitness\n");
     */
    
        
    for (i = 0; i < generations; i++) {
        for (j = 0; j < popsize; j++) {
            currentvictim = ChooseVictim(popsize);
            currentparent1 = ChooseParentWithTree(wholepopulationwistree, popsize, sumofwis);
            currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, sumofwis);
            while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
                currentparent2 = ChooseParentWithTree(wholepopulationwistree, popsize, sumofwis);
            }
            
            //Following lines produce a gamete from parent 1 and add deleterious and beneficial mutations to the gamete.
            RecombineChromosomesIntoGamete(currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
            MutateGamete(chromosomesize, numberofchromosomes, parent1gamete, deleteriousmutationrate, Sd);
            MutateGamete(chromosomesize, numberofchromosomes, parent1gamete, beneficialmutationrate, Sb);
            
            //Following lines produce a gamete from parent 2 and add deleterious and beneficial mutations to the gamete.            
            RecombineChromosomesIntoGamete(currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
            MutateGamete(chromosomesize, numberofchromosomes, parent2gamete, deleteriousmutationrate, Sd);
            MutateGamete(chromosomesize, numberofchromosomes, parent2gamete, beneficialmutationrate, Sb);

            ReplaceVictim(parent1gamete, parent2gamete, popsize, currentvictim, chromosomesize, numberofchromosomes, Sd, psumofwis, wholepopulationgenomes, totalindividualgenomelength, wholepopulationwistree, wholepopulationwisarray);
            
            if (VERBOSE) {
                fprintf(verbosefilepointer, "\nSum of wis after generation %d, round %d: %.24Lf", i+1, j+1, sumofwis);
            }
        }
        

        //Following code calculates the variance in log(fitness) of the population.
        //May use an imprecise algorithm -- check before using as data.
        size_t step = 1;
        int m;
        variancesum = CalculateVarianceInLogFitness(popsize, wholepopulationwisarray, sumofwis);
        double c0, cov00, cov01, cov11, sumsq;
        if (isburninphaseover == 0) {
            UpdateLastTwentyGenerations(lasttwentygenerationsvariance, variancesum);
            UpdateLastTwentyGenerations(literallyjustlasttwentygenerations, i+1);
            if (i > 19) {           //to avoid calling the end of the burn-in phase at generation one
                                    //because of setting pre-simulation generations to zeroes
                                    //I just won't start looking for the end of the burn-in phase until 20 generations
                                    //This would be a mild problem if a simulation should end in 20 generations, but that shouldn't ever happen.
                slopeofvariance = 0.0;
                for (m = 0; m < 20; m++) {
                    printf("%f,%f\n", lasttwentygenerationsvariance[m], literallyjustlasttwentygenerations[m]);
                }
                gsl_fit_linear(literallyjustlasttwentygenerations, step, lasttwentygenerationsvariance, step, 20, &c0, &slopeofvariance, &cov00, &cov01, &cov11, &sumsq);
                if (slopeofvariance < arbitrarynumber) {
                    endofburninphase = i;
                    isburninphaseover = 1;
                }
            }
        }        
        
        //This is the main data output, currently the summed fitness and variance in log(fitness) in the population.
        fprintf(datafilepointer, "%d,%Lf,%.18f\n", i+1, sumofwis, variancesum);
        
        //These lines ensure that the magnitude of fitness hasn't declined by too much.
        //At extremely small fitness values, floating-point math becomes imprecise.
        //These lines end the simulation if fitness declines below 10^-10, which should represent a completely degraded population.
        currentfittestindividualswi = FindFittestWi(wholepopulationwisarray, popsize);
        if (currentfittestindividualswi < pow(10.0, -10.0)) {
            fprintf(miscfilepointer, "\nFitness declined to less than 10^-10 during generation %d.", i+1);
            fprintf(datafilepointer, "Fitness declined to catastrophic levels in generation %d.\n", i+1);
            i = generations;
        }
        
        //These lines sort the wis of the population in ascending order and print them as data, to examine the distribution of fitness in the population.
        //The algorithms are currently very inefficient -- if this needs to be performed more than once per simulation, algorithms should be updated.
        /*
         if (i == 1000) {
             for (k = 0; k < popsize; k++) {
                 sortedwisarray[k] = wholepopulationwisarray[k];
             }
             DoubleBubbleSort(sortedwisarray, popsize);
             for (k = 0; k < popsize; k++) {
                 fprintf(fitnessdistributiondatafilepointer, "%d,%.10Lf\n", k+1, sortedwisarray[k]);
             }
         }
         */
        
         
         

    }
        
    //Following code allows for printing full genomes for debugging purposes.
    /*
    fprintf(verbosefilepointer, "\n\nIndividual's genomes:");
    int currentspotingenome, k;
    for (i = 0; i < popsize; i++) {
        fprintf(verbosefilepointer, "\nIndividual %d: ", i+1);
        for (j = 0; j < numberofchromosomes; j++) {
            for (k = 0; k < chromosomesize; k++) {
                currentspotingenome = ((i * totalindividualgenomelength) + k);
                fprintf(verbosefilepointer, "%f ", wholepopulationgenomes[currentspotingenome]);
            }
            fprintf(verbosefilepointer, "\n");
            for (k = 0; k < chromosomesize; k++) {
                currentspotingenome = ((i*totalindividualgenomelength) + (totalindividualgenomelength/2) + k);
                fprintf(verbosefilepointer, "%f ", wholepopulationgenomes[currentspotingenome]);
            }
        }
    }
    */
        
    fclose(verbosefilepointer);
    fclose(veryverbosefilepointer);
    fclose(miscfilepointer);
    fclose(datafilepointer); //closes data files
    free(datafilename);
    //free(fitnessdistributiondatafilename);    
    free(wholepopulationgenomes);
    free(wholepopulationwistree);
    free(wholepopulationwisarray);

    printf("\nOutput is:\n");
    printf("%d", endofburninphase);
    printf("\n");

}

void UpdateLastTwentyGenerations(double lasttwentygenerations[], double newgeneration)
{
    double storage[20];
    int m;
    for (m = 0; m < 19; m++) {
        storage[m] = lasttwentygenerations[m+1];
    }
    storage[19] = newgeneration;
    for (m = 0; m < 20; m++) {
        lasttwentygenerations[m] = storage[m];
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
    variancesum = (variancesum/popsize);
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

void InitializePopulation(long double *wholepopulationwistree, long double *wholepopulationwisarray, int populationsize, double *populationgenomes, int totalpopulationgenomelength) {
	int i, j;
        for (i = 0; i < populationsize; i++) {
            wholepopulationwistree[i] = 1.0; //for relative fitness, all individuals start with probability of being chosen as a parent of 1/N
            wholepopulationwisarray[i] = 1.0;
        }
        //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
	for (i = 0; i < populationsize; i++) {
		j = i + LSB(i+1);
                if (j < populationsize) {
                    wholepopulationwistree[j] += wholepopulationwistree[i];
                }
	}
        for (i = 0; i < totalpopulationgenomelength; i++) {
            populationgenomes[i] = 0.0;
        }
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

int ChooseParentWithTree(long double *wholepopulationwistree, int popsize, long double sumofwis)
{
    long double randomnumberofbirth;
    int newparent = 0;
    
    randomnumberofbirth = (ldexp(pcg32_random(), -32)) * sumofwis;
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize;
    if (leftbound > rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    
    newparent = (SearchTree(leftbound, rightbound, randomnumberofbirth, wholepopulationwistree));
    return newparent;
}

int ChooseVictim(int populationsize)
{
	int randomindividual = pcg32_boundedrand(populationsize);
	return randomindividual;
}

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(int persontorecombine, int chromosomesize, int numberofchromosomes, double *gamete, double *populationgenomes, int totalindividualgenomelength)
{
	int recombinationsite;
	int startchromosome;
	int h, i;
        int startofindividual;
        startofindividual = persontorecombine * totalindividualgenomelength;
	for (h = 0; h < numberofchromosomes; h++) {
		startchromosome = pcg32_boundedrand(2); //generates either a zero or a one to decide to start with chromosome 1 or 2.
		recombinationsite = pcg32_boundedrand(chromosomesize);
		for (i = 0; i < recombinationsite; i++) {
			if (startchromosome == 0) {
				gamete[h*chromosomesize + i] = populationgenomes[startofindividual + (h*chromosomesize) + i];
			}
			else {
				gamete[h*chromosomesize + i] = populationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
			}
		}
                for (i = recombinationsite; i < chromosomesize; i++) {
                        if (startchromosome == 0) {
                                gamete[h*chromosomesize + i] = populationgenomes[startofindividual + totalindividualgenomelength / 2 + (h*chromosomesize) + i];
                        }
                        else {
                                gamete[h*chromosomesize + i] = populationgenomes[startofindividual + (h*chromosomesize) + i];
                        }
		}
	}
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
	} else { 				//for larger means, use rejection method.
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

void MutateGamete(int chromosomesize, int numberofchromosomes, double *gamete, float mutationrate, double mutationeffectsize)
{
	int i;
	float meannumberofmutations = 2.0 * mutationrate * (float) chromosomesize; //Assumes U = 2, could easily be modified to take input values of U.
	int numberofmutations = SampleFromPoisson(meannumberofmutations);
	for(i = 0; i < numberofmutations; i++) {	
		int randomchromosometomutate = pcg32_boundedrand(numberofchromosomes); //if we decide to include heterogenous rates of recombination/mutation, both of these will need to be replaced by a function that weights each linkage block's probability of mutating.
        	int randomblocktomutate = pcg32_boundedrand(chromosomesize);
        	gamete[randomchromosometomutate*chromosomesize + randomblocktomutate] += log(1 + mutationeffectsize); //this is going to be a distribution of effects at some point.
	}
}

//this function seems inefficient, but with recombination and mutation, I'm not sure there's a significantly easier way.
double CalculateWi(int numberofchromosomes, int chromosomesize, double *parent1gamete, double *parent2gamete, double mutationeffectsize, int totalindividualgenomelength)
{
	double newwi = 0.0;
	long double currentlinkageblockssum = 0.0;
	int i;

	for (i = 0; i < (totalindividualgenomelength/2); i++) {
            currentlinkageblockssum += parent1gamete[i];
            currentlinkageblockssum += parent2gamete[i];
	}
	newwi = exp(currentlinkageblockssum);
	return newwi;
}

void ReplaceVictim(double *parent1gamete, double *parent2gamete, int currentpopsize, int currentvictim, int chromosomesize, int numberofchromosomes, double mutationeffectsize, long double *sumofwis, double *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationwistree, long double *wholepopulationwisarray)
{
	int i;
	double newwi;
	
        newwi = CalculateWi(numberofchromosomes, chromosomesize, parent1gamete, parent2gamete, mutationeffectsize, totalindividualgenomelength);
	
        for (i = 0; i < (totalindividualgenomelength/2); i++) {
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + i] = parent1gamete[i];
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
//It would be more efficient to build directly into victim slot, but right now it is possible for the victim to also be a parent, later add an if statement to more efficiently deal with more common case.
	}

	*sumofwis -= wholepopulationwisarray[currentvictim];
	Fen_set(wholepopulationwistree, currentpopsize, newwi, currentvictim);
        wholepopulationwisarray[currentvictim] = (long double) newwi;
	*sumofwis += (long double) newwi;

}
