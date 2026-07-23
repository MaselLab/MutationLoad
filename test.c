#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
    double *fitnessArray;
    double *mutatorArray;
    float fitness;
    long double mutationRate;
} Individual;

Individual createIndividual(double *fitnessArray, double *mutatorArray, int totalindividualgenomelength, float minimumMutationRate, float mutationRateGrowthFactor){
    Individual ind;
    ind.fitnessArray = fitnessArray;
    ind.mutatorArray = mutatorArray;

    int i;
    long double currentlinkageblockssum = 0.0;
    long double mutationRateSum = 0.0;
    for (i = 0; i < totalindividualgenomelength; i++){
        currentlinkageblockssum += fitnessArray[i];
        mutationRateSum += mutatorArray[i];
    }
    ind.fitness = exp(currentlinkageblockssum);
    ind.mutationRate = minimumMutationRate*exp(mutationRateGrowthFactor*mutationRateSum);

    return ind;
}

void UpdateIndividual(Individual *ind, int totalindividualgenomelength){
    int i;
    long double currentlinkageblockssum = 0.0;
    long double mutationRateSum = 0.0;

    for (i = 0; i < totalindividualgenomelength; i++){
        currentlinkageblockssum += ind -> fitnessArray[i];
        mutationRateSum += ind -> mutatorArray[i];
    }
    ind -> fitness = exp(currentlinkageblockssum);
    ind -> mutationRate = exp(mutationRateSum);
}

void printIndividual(Individual ind){
    printf("Fitness: %f\n", ind.fitness);
    printf("Mutation Rate: %Lf\n", ind.mutationRate);
}

int main(){
    Individual ind = createIndividual((double[]){0.1, 0.2, 0.3}, (double[]){0.01, 0.02, 0.03}, 3, 0.001, 3);
    printIndividual(ind);
    UpdateIndividual(&ind, 3);
    printIndividual(ind);
    return 0;
}