#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../stat.h"
#include "../vec3.h"
#include "../stats_func.h"
#include "../h2.h"

// Function declerations
void generateMarkovChain(double (*chain)[6], double alpha, double d, int N, gsl_rng * q);
void metropolisStep(double prev[6], double next[6], double alpha, double d, gsl_rng * q);
void setZero(double (*arr)[6], int N);
void randomize(double * v, int sz, double min, double max, gsl_rng * q);

// Global variables
int metropolisCount = 0; // Used for counting the acceptance rate of Metropolis
int metropolisTotal = 0;

// gsl stuff
gsl_rng * init_rng();

int main()
{
    // Output file
    char output_file[255];
    sprintf(output_file, "alphas.dat");     
    FILE * file = fopen(output_file, "w");   
        
    gsl_rng * q = init_rng();
    //Parameters
    int chainLength = 200;
    int equibTime = 100;
    int measuredLength = chainLength - equibTime;
    int numTrials = 100000;
    double d = 0.9; // Stepping parameter

    int k, i, j, trial;
    double alpha_start = 0.05;
    double alpha_end = 0.25;
    double alpha_step = 0.005;
    int nbr_of_alphas = round((alpha_end - alpha_start)/alpha_step) +1;
    
    // Array containing mean local energy and std for different alphas
    double chain[chainLength][6]; // Initialize Markov chain
    double localEnergy[measuredLength];


    double alpha, mean;
    for(k = 0; k < nbr_of_alphas; k++) {
        alpha = alpha_start + k*alpha_step;
        
        for(j = 0; j < measuredLength; ++j)
            localEnergy[j] = 0; // Reset for new alpha measurements

        for (trial = 0; trial < numTrials; ++trial)
        {
            setZero(chain, chainLength);
            randomize(chain[0],6,0,1,q);
            // Sample Markov chain
            generateMarkovChain(chain, alpha, d, chainLength, q);
            
            // Mean energy of the last chainLength - equibTime = measuredLength stepts
            j = 0;
            for (i = equibTime+1; i < chainLength; i++)
                localEnergy[j++] += energy(chain[i], alpha);
        }

        mean = get_mean(localEnergy, measuredLength) / numTrials;
        fprintf(file, "%e\t%e\n", alpha, mean);
        printf("alpha = %.3f: mean = %.3f\n", alpha, mean);
    }
    fclose(file);
    printf("Written to \"%s\".\n", output_file);
    return 0;
}

void randomize(double * v, int sz, double min, double max, gsl_rng * q) {
    int i;
    for(i = 0; i < sz; i++)
        v[i] = min + (max-min) * gsl_rng_uniform(q);
}

void setZero(double (*arr)[6], int N) {
    int i,j;
    for (i = 0; i < N; i++)
        for(j = 0; j < 6; j ++)
            arr[i][j] = 0;
}

void generateMarkovChain(double (*chain)[6], double alpha, double d, int N, gsl_rng * q) {
    int i;

    for(i = 1; i< N; i++) {
        metropolisStep(chain[i-1], chain[i], alpha, d, q);
    }
}

void metropolisStep(double prev[6], double next[6], double alpha, double d, gsl_rng * q) {
    int i;

    double tmp[6];
    double r;
    //int ri = gsl_rng_uniform(q)*6;
    for (i = 0; i < 6; ++i) {
    //for (i = ri; i < ri+1; ++i) {
        r = gsl_rng_uniform(q);
        tmp[i] = prev[i] + d*(r - 0.5);
    }

    double pr = probFunction(tmp, prev, alpha);
    r = gsl_rng_uniform(q);
    
    if (pr > r) {// Accept new step
        for (i = 0; i < 6; ++i)
            next[i] = tmp[i];
    } else {
        for (i = 0; i < 6; ++i)
            next[i] = prev[i];
    }
    
    if (pr >= r)        
        metropolisCount++;
    metropolisTotal++;
}

gsl_rng * init_rng()
{
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
	return q;
}
