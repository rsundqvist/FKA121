#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stat.h"
#include "vec3.h"

// Function declerations
double trialWaveFunction(double * R, double alpha);
double probFunction(double * Rnew, double * Rold, double alpha);
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
    sprintf(output_file, "trajectory.dat", 0);    
        
    gsl_rng * q = init_rng();
    //Parameters
    int chainLength = 100000;
    double alpha = 0.1;
    double d = 1.85;

    // Initialize Markov chain
    double chain[chainLength][6];
    setZero(chain, chainLength);
    randomize(chain[0],6,0,1,q);
    
    // Sample Markov chain
    generateMarkovChain(chain, alpha, d, chainLength, q);
    
    //====================================================================//
    // Simulation complete - print data to file(s)
    //====================================================================//
    int i;
    FILE * file1 = fopen(output_file, "w");
    double x1,y1,z1;
    double x2,y2,z2;
    double theta;
    if (file1 != NULL){
        for (i = 0; i < chainLength; i++) {         
            // Print file1
            x1 = chain[i][0];
            y1 = chain[i][1];
            z1 = chain[i][2];

            x2 = chain[i][3];
            y2 = chain[i][4];
            z2 = chain[i][5];

            double r1[3] = {x1,y1,z1}; // Must redeclare to use this syntax
            double r2[3] = {x2,y2,z2};

            fprintf (file1,"%e %e %e %e %e %e %e\n",
                x1,y1,z1,
                x2,y2,z2,
                angle(r1,r2)
                );
        }
        
        
        // Close file(s)
        fclose(file1);
        printf("Data printed to file: \"%s\".\n", output_file);
    } else {
        printf("File is NULL!\n");
    }

    printf("Acceptance rate: %.3f\n", (double)metropolisCount/metropolisTotal);
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

double trialWaveFunction(double * R, double alpha) {
    // Compute r1 and r2
    double r1[3] = {R[0],R[1],R[2]};
    double r2[3] = {R[3],R[4],R[5]};
    double r1_norm = norm(r1);
    double r2_norm = norm(r2);
    // Compute r12
    double r12_norm = distance(r1, r2);

    return exp(-2*(r1_norm + r2_norm))*(r12_norm/(2*(1 + alpha*r12_norm)));
}

void generateMarkovChain(double (*chain)[6], double alpha, double d, int N, gsl_rng * q) {
    int i;
    double r;

    for(i = 1; i< N; i++) {
        metropolisStep(chain[i-1], chain[i], alpha, d, q);
    }
}

void metropolisStep(double prev[6], double next[6], double alpha, double d, gsl_rng * q) {
    int i;

    double tmp[6];
    double r;
    for (i = 0; i < 6; ++i) {
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

// Computing the probability of accepting a new proposal in the Markov chain
double probFunction(double * Rnew, double * Rold, double alpha){
    return trialWaveFunction(Rnew,alpha)/trialWaveFunction(Rold,alpha);
}
