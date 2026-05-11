#include <stdio.h>

typedef struct{
    double *fitnessArray;
    double *mutatorArray;
    long double fitness;
    long double mutationRate;
} Individual;

Individual createIndividual(double *fitnessArray, double *mutatorArray, int totalindividualgenomelength, int minimumMutationRate, int mutationRateGrowthFactor){
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
    printf("Fitness: %Lf\n", ind.fitness);
    printf("Mutation Rate: %Lf\n", ind.mutationRate);
}

int main(){
    return 0;
}