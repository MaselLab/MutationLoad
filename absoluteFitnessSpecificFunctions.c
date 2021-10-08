int monteCarloStep(int popSize, double sumWi, double* pTimeElapsed, double sumOfDeathRates, struct individual* popArray, int MAX_POP_SIZE, double currentTimeStep) {

    int boolVar;
    int randSeed = 1;

    double deathRate;
    double birthRate;
    double time;
    double mean;
    double deathRateAtPerfection;
    double randomNumber;

    double* pRandomNumber = &randomNumber;

    deathRateAtPerfection = popSize;

    const double b = 10;

    deathRate = sumOfDeathRates;

    //rate of births is calculated using equation used in lab write up
    birthRate = rateOfBirthsCalc(popSize, b, MAX_POP_SIZE, currentTimeStep);


    mean = ((1.0) / (deathRate + birthRate));
    time = ExponentialDerivateOfUnitMeanOne(randSeed);//draws a random number from a distribution with unit mean 1. This occurs because a event is most likely to occur right after a previous event.
    time = (time) * (mean);

    *pTimeElapsed = time + *pTimeElapsed;

    //This is the actual monte carlo step
    boolVar = discoverEvent(deathRate, birthRate, currentTimeStep);

    return boolVar;

}

double rateOfBirthsCalc(int populationSize, double b, int MAX_POP_SIZE, double currentTimeStep) {

    double birthRate;
    double popSizeConvertedToDouble = populationSize * 1.0;

    birthRate = (b) * (popSizeConvertedToDouble) * (1 - (popSizeConvertedToDouble / MAX_POP_SIZE));

    return  birthRate;

}

double ExponentialDerivateOfUnitMeanOne(float idum) {

    float ran1(long* idum);

    float dum;



    do
        dum = ldexp(pcg32_random(), -32);//check if this is how pcg32 is called
    while (dum == 0.0);

    return -log(dum);

}

int discoverEvent(double deathRate, double birthRate, double currentTimeStep) {

    int boolBirth;

    double combinedBirthDeathRate;
    double cutOffPoint;
    double randomNumber;

    randomNumber = ldexp(pcg32_random(), -32);

    //finds what percent of events at this time would be a death event to determine the cuttoff point.
    //if lands on the cutoff point death occurs, techincally skewing towards death but is very minor.
    combinedBirthDeathRate = deathRate + birthRate;
    cutOffPoint = deathRate / combinedBirthDeathRate;

    if (randomNumber > cutOffPoint) {
        boolBirth = 1;
    }
    else if (randomNumber <= cutOffPoint) {
        boolBirth = 0;
    }
    else {

    }

    return boolBirth;

}

int ChooseVictimWithTree(long double* wholepopulationwistree, int popsize, long double sumofwis, long double inverseSumOfWis, int eventNumber)//using pinversesum intesting and is currently unchanged 11/18/2019
{
    long double randomnumberofdeath;
    int newVictim = 0;
    randomnumberofdeath = (long double)ldexp(pcg32_random(), -32) * (inverseSumOfWis);
    //Above line generates a random integer between 0 and 2^32, then multiplies by 2^-32
    //to generate a float between 0 and 1 and then multiplies by the sum of wis
    //to get a number between 0 and the sum of wis.
    int leftbound, rightbound;
    leftbound = 0;
    rightbound = popsize;//change variable name
    if (leftbound >= rightbound) {
        return -1;
        fprintf(miscfilepointer, "\nError: population size is %d.", popsize);
    }
    //Above lines initialize the variables necessary for the SearchTree function and check for an extinct population.
    //the random death is causing a strange number
    newVictim = (SearchTree(leftbound, rightbound, randomnumberofdeath, wholepopulationwistree, eventNumber));//fixed possible error
    return newVictim;
}

int ChooseParent(int populationsize)
{
    int randomindividual = pcg32_boundedrand(populationsize);

    return randomindividual;
}