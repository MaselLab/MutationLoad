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

#define VERBOSE 0
#define VERYVERBOSE 0
#define MISCELLANEOUS 1
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654
#define RANDOMNUMBERSEED 1020
//All of these should probably eventually be arguments input by user.

FILE	*verbosefilepointer;
FILE	*veryverbosefilepointer;
FILE	*miscfilepointer;
FILE	*datafilepointer;

long double Fen_sum(long double *, int);
void Fen_add(long double *, int, long double, int);
long double Fen_range(long double *, int, int);
long double Fen_get(long double *, int);
void Fen_set(long double *, int, long double, int);
void InitializePopulation(long double *, int , long double , int *, int);
//int ChooseVictimWithoutTree(long double *, int , long double);
int SearchTree(int, int, long double, long double *);
int ChooseVictimWithTree(long double *, int, long double);
int ChooseParent(int);
void RecombineChromosomesIntoGamete(int, int , int , int *, int *, int);
int SampleFromPoisson(float );
void MutateGamete(int , int , int *, float );
double CalculateDeathRate(int , int , int *, int *, double , double, int);
void ReplaceVictim(int , int , int, int, int, int, int , int , float , double , long double *, double, int *, int, long double *);
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

	double mutationrate;
	mutationrate = atof(argv[4]); //remember that this is the per-locus mutation rate, not the genome-wide mutation rate.

	long double nulldeathrate;
	nulldeathrate = atof(argv[5]);

	int chromosomesize;
	chromosomesize = atoi(argv[6]);

	int numberofchromosomes;
	numberofchromosomes = atoi(argv[7]); //remember that this is total number of chromosomes, not ploidy -- all individuals will be diploid.

	int i, j;
        int *wholepopulationgenomes;
        int totalpopulationgenomelength;
        int totalindividualgenomelength;
        totalpopulationgenomelength = popsize * numberofchromosomes * 2 * chromosomesize;
        totalindividualgenomelength = numberofchromosomes * 2 * chromosomesize;
        wholepopulationgenomes = malloc(sizeof(int) * totalpopulationgenomelength);
	long double sumofdeathrates;
	long double *psumofdeathrates = &sumofdeathrates;
	long double *wholepopulationdeathratestree;
        wholepopulationdeathratestree = malloc(sizeof(long double) * popsize);
	int currentparent1;
	int currentparent2;
	int currentvictim;
	/*Initializes the population as an array of individuals, 
	the current victim as a pointer to an individual,
	and the current parents as pointers to individuals.*/

	sumofdeathrates = (long double) nulldeathrate * (long double) popsize;
       
        
	InitializePopulation(wholepopulationdeathratestree, popsize, nulldeathrate, wholepopulationgenomes, totalpopulationgenomelength);
	/*Sets the initial population to have zeroes in all their linkage blocks,
	death rates equal to the baseline death rate, and an identifier number.
	It also sums the death rates and returns the sum.*/
    
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

	strcat(datafilename, "nulldeathrate"); //for adding values of nulldeathrate to the data name.
	strcat(datafilename, argv[5]);

	strcat(datafilename, "chromsize"); //for adding values of chromosome size to the data name.
	strcat(datafilename, argv[6]);

	strcat(datafilename, "chromnum"); //for addiving values of the number of chromosomes to the data name.
	strcat(datafilename, argv[7]);

	datafilepointer = fopen(datafilename, "w"); //opens the file to which to print data to be plotted.
	fprintf(datafilepointer, "Generation,Sum.of.death.rates,\n");
    
        
	for (i = 0; i < generations; i++) {
		for (j = 0; j < popsize; j++) {
			currentvictim = ChooseVictimWithTree(wholepopulationdeathratestree, popsize, sumofdeathrates);
			currentparent1 = ChooseParent(popsize);
                        currentparent2 = ChooseParent(popsize);
                        while (currentparent1 == currentparent2) { //probably not ideal, since it'll never break with population sizes of zero or one.
                            currentparent2 = ChooseParent(popsize);
                        }
			ReplaceVictim(i, j, popsize, currentvictim, currentparent1, currentparent2, chromosomesize, numberofchromosomes, mutationrate, Sd, psumofdeathrates, nulldeathrate, wholepopulationgenomes, totalindividualgenomelength, wholepopulationdeathratestree);
			if (VERBOSE) {
                                fprintf(verbosefilepointer, "\nSum of death rates after generation %d, round %d: %.24Lf", i+1, j+1, sumofdeathrates);
			}
		}

		fprintf(datafilepointer, "%d,%Lf,\n", i+1, sumofdeathrates);

		//RecalculateSumOfDeathRates(psumofdeathrates, popsize, i, popsize-1, population);
	}

        long double testsum10 = Fen_sum(wholepopulationdeathratestree, 10);
        
	printf("\nOutput is:\n");

        printf("%Lf", testsum10);
        
	fclose(verbosefilepointer);
	fclose(veryverbosefilepointer);
	fclose(miscfilepointer);
	fclose(datafilepointer); //closes data files
        free(datafilename);
        
        free(wholepopulationgenomes);
        free(wholepopulationdeathratestree);
	
	printf("\n");
}

/*All Fenwick tree functions from Wikipedia page "Fenwick tree" URL:https://en.wikipedia.org/wiki/Fenwick_tree
 This project is licensed under the GNU General Public License version 3.0, 
 * which is compatible with the CC-BY-SA license of Wikipedia text.*/


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

void InitializePopulation(long double *wholepopulationdeathratestree, int populationsize, long double baselinedeathrate, int *populationgenomes, int totalpopulationgenomelength) {
	int i, j;
        for (i = 0; i < populationsize; i++) {
            wholepopulationdeathratestree[i] = baselinedeathrate;
        }
        //this for loop taken from the Fen_init function in sample implementation from 'Fenwick tree' Wikipedia page.
	for (i = 0; i < populationsize; i++) {
		j = i + LSB(i+1);
                if (j < populationsize) {
                    wholepopulationdeathratestree[j] += wholepopulationdeathratestree[i];
                }
	}
        for (i = 0; i < totalpopulationgenomelength; i++) {
            populationgenomes[i] = 0;            
        }
}
/*
int ChooseVictimWithoutTree(long double *wholepopulationdeathratestree, int popsize, long double sumofdeathrates)
{
	long double randomnumberofdeath;
	int potentialvictim = 0;
	double partialsum = wholepopulationdeathrates[0];

	randomnumberofdeath = (ldexp(pcg32_random(), -32)) * sumofdeathrates; //Generates random integer between 0 and 2^32, then 
										//multiplies by 2^-32 to get a float between 0 and 1, 
										//then multiplies by the sum of the death rates. 
	//Replace the next part with a binary search.
        while (partialsum < randomnumberofdeath) {
		potentialvictim++;
		partialsum += wholepopulationdeathrates[potentialvictim];
	}
	if (potentialvictim < popsize) {
		return potentialvictim;
	} else { //this section is to catch bugs where sumofdeathrates has become slightly larger than the actual sum of death rates in the population.
		fprintf(miscfilepointer, "\nError: victim %d chosen on the basis of randomnumberofdeath %Lf indicating 0-1 random number of %Lf in ChooseVictim was outside of the initialized population.", potentialvictim, randomnumberofdeath, randomnumberofdeath/sumofdeathrates);
                potentialvictim = popsize-1;
                return potentialvictim;
	}
}
*/
int SearchTree(int leftbound, int rightbound, long double targetvalue, long double *Fenwicktree)
{
    int middle;
    middle = floor((leftbound+rightbound)/2);
    long double partialsumatmiddle;
    long double partialsumatmiddleminusone;
    partialsumatmiddle = Fen_sum(Fenwicktree, middle);
    partialsumatmiddleminusone = Fen_sum(Fenwicktree, middle-1);
    if(partialsumatmiddle < targetvalue) {
        SearchTree(middle, rightbound, targetvalue, Fenwicktree);
    }
    if(partialsumatmiddle > targetvalue) {
        if(partialsumatmiddleminusone > targetvalue) {
            SearchTree(leftbound, middle, targetvalue, Fenwicktree);
        }
        else {
            return middle;
        }
    }
    return middle;
}

int ChooseVictimWithTree(long double *wholepopulationdeathratestree, int popsize, long double sumofdeathrates)
{
    long double randomnumberofdeath;
    int victim = 0;
    
    randomnumberofdeath = (ldexp(pcg32_random(), -32)) * sumofdeathrates;
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of death rates
    //to get a number between 0 and the sum of death rates.
    
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize-1;
    if (leftbound > rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    
    victim = SearchTree(leftbound, rightbound, randomnumberofdeath, wholepopulationdeathratestree);
    return victim;
}

int ChooseParent(int populationsize)
{
	int randomindividual = pcg32_boundedrand(populationsize);
	return randomindividual;
}

//1 recombination site per chromosome
void RecombineChromosomesIntoGamete(int persontorecombine, int chromosomesize, int numberofchromosomes, int *gamete, int *populationgenomes, int totalindividualgenomelength)
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

void MutateGamete(int chromosomesize, int numberofchromosomes, int *gamete, float mutationrate)
{
	int i;
	float meannumberofmutations = 2.0 * mutationrate * (float) chromosomesize;
	int numberofmutations = SampleFromPoisson(meannumberofmutations);
	for(i = 0; i < numberofmutations; i++) {	
		int randomchromosometomutate = pcg32_boundedrand(numberofchromosomes);
        	int randomblocktomutate = pcg32_boundedrand(chromosomesize);
        	gamete[randomchromosometomutate*chromosomesize + randomblocktomutate] += 1; //this is going to be a distribution of effects at some point.
	}
}

//this function seems inefficient, but with recombination and mutation, I'm not sure there's a significantly easier way.
double CalculateDeathRate(int numberofchromosomes, int chromosomesize, int *parent1gamete, int *parent2gamete, double mutationeffectsize, double nulldeathrate, int totalindividualgenomelength)
{
	double newdeathrate = 0.0;
	int numberofmutations = 0;
	int i;

	for (i = 0; i < (totalindividualgenomelength/2); i++) {
            numberofmutations += parent1gamete[i];
            numberofmutations += parent2gamete[i];
	}
	newdeathrate = (double) nulldeathrate + (mutationeffectsize * (double) numberofmutations);
	return newdeathrate;
}

void ReplaceVictim(int currentgeneration, int currentround, int currentpopsize, int currentvictim, int currentparent1, int currentparent2, int chromosomesize, int numberofchromosomes, float mutationrate, double mutationeffectsize, long double *sumofdeathrates, double nulldeathrate, int *wholepopulationgenomes, int totalindividualgenomelength, long double *wholepopulationdeathratestree)
{
	int i;
	double newdeathrate;
	int parent1gamete[numberofchromosomes*chromosomesize], parent2gamete[numberofchromosomes*chromosomesize];
	RecombineChromosomesIntoGamete(currentparent1, chromosomesize, numberofchromosomes, parent1gamete, wholepopulationgenomes, totalindividualgenomelength);
	RecombineChromosomesIntoGamete(currentparent2, chromosomesize, numberofchromosomes, parent2gamete, wholepopulationgenomes, totalindividualgenomelength);
	MutateGamete(chromosomesize, numberofchromosomes, parent1gamete, mutationrate);
	MutateGamete(chromosomesize, numberofchromosomes, parent2gamete, mutationrate);
	newdeathrate = CalculateDeathRate(numberofchromosomes, chromosomesize, parent1gamete, parent2gamete, mutationeffectsize, nulldeathrate, totalindividualgenomelength);
	for (i = 0; i < (totalindividualgenomelength/2); i++) {
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + i] = parent1gamete[i];
            wholepopulationgenomes[currentvictim*totalindividualgenomelength + totalindividualgenomelength/2 + i] = parent2gamete[i];
//It would be more efficient to build directly into victim slot, but right now it is possible for the victim to also be a parent, later add an if statement to more efficiently deal with more common case.
	}

	*sumofdeathrates -= Fen_get(wholepopulationdeathratestree, currentvictim);
	Fen_set(wholepopulationdeathratestree, currentpopsize, newdeathrate, currentvictim);
	*sumofdeathrates += (long double) newdeathrate;

}

////This function is to check that the procedure for adding and subtracting death rates does not accumulate error in the sum of death rates.
////It should be taken out of the final version, for computational efficiency. 
//void RecalculateSumOfDeathRates(long double *sumofdeathrates, int popsize, int currentgeneration, int currentround, struct individual population[popsize])
//{
//	long double newsumofdeathrates = 0.0;
//	int i;
//	for (i = 0; i < popsize; i++) {
//		newsumofdeathrates += (long double) population[i].deathrate;
//	}
//
//	if (newsumofdeathrates != *sumofdeathrates) {
//		fprintf(miscfilepointer, "\nThe new sum in generation %d, round %d IS NOT equal to the old sum.", currentgeneration+1, currentround+1);
//	        fprintf(miscfilepointer, "\nNew sum of death rates in generation %d, round %d is: %.24Lf", currentgeneration+1, currentround+1, newsumofdeathrates);
//                fprintf(miscfilepointer, "\nOld sum of death rates in generation %d, round %d is: %.24Lf", currentgeneration+1, currentround+1, *sumofdeathrates);
//		fprintf(miscfilepointer, "\nDifference between new and old sum of death rates is: %.24Lf", (newsumofdeathrates - *sumofdeathrates));
//	}
//
//	*sumofdeathrates = newsumofdeathrates;	
//}
