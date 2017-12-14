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
    char output_file[64];
    sprintf(output_file, "trajectory.dat");    
        
    gsl_rng * q = init_rng();
    //Parameters
    int chainLength = 2200; // # markov steps
    int numTrials = 300;
    double alpha = 0.1; // Trial function parameter
    double d = 0.9; // Stepping parameter
    d = 0.01;

    // Initialize Markov chain
    double chain[chainLength][6];
    double localEnergy[chainLength][numTrials];
    int i;


  	FILE * energy_traj = fopen("energy_traj10.dat", "w");


    int printTimes = numTrials/10;
    for (int trial = 0; trial < numTrials; ++trial)
    {

	    setZero(chain, chainLength);
	    randomize(chain[0],6,0,1,q);
	    // Sample Markov chain
	    generateMarkovChain(chain, alpha, d, chainLength, q);

        for (i = 0; i < chainLength; i++) {     
        	localEnergy[i][trial] = energy(chain[i], alpha);
        }

    	if (trial%printTimes==0) {
    		printf("%.3f\n", (double) trial/numTrials);

        	for (i = 0; i < chainLength; i++) {  
        		fprintf(energy_traj, "%e ", localEnergy[i][trial]);
        	}
        		fprintf(energy_traj, "\n");
    	}
    }


    FILE * file0 = fopen("energy_stats10.dat", "w");
    double mean_E[chainLength], var_E[chainLength];
    for (int i = 0; i < chainLength; ++i)
    {
    	mean_E[i] = get_mean(localEnergy[i], chainLength);
    	var_E[i] = get_variance(localEnergy[i], mean_E[i], chainLength);
    	fprintf(file0, "%e\t%e\n", mean_E[i], var_E[i]);
    }
    fclose(file0);

    //====================================================================//
    // Simulation complete - print data to file(s)
    //====================================================================//
    FILE * file1 = fopen(output_file, "w");
    double x1,y1,z1;
    double x2,y2,z2;
    if (file1 != NULL){
        for (i = 0; i < chainLength; i++) {         
            // Print file1
            x1 = chain[i][0];
            y1 = chain[i][1];
            z1 = chain[i][2];

            x2 = chain[i][3];
            y2 = chain[i][4];
            z2 = chain[i][5];

            double R1[3] = {x1,y1,z1}; // Must redeclare to use this syntax
            double R2[3] = {x2,y2,z2};

            fprintf (file1,"%e %e %e %e %e %e %e\n",
                x1,y1,z1,
                x2,y2,z2,
                angle(R1,R2)
            );
        }
        
        
        // Close file(s)
        fclose(file1);
        printf("Data printed to file: \"%s\".\n", output_file);
    } else {
        printf("File is NULL!\n");
    }

    printf("Acceptance rate: %.3f\n", (double)metropolisCount/metropolisTotal);
    
    int N = chainLength;
    double *values = malloc(sizeof (double[N]));
    for (i = 0; i < N; i++) {
        values[i] = energy(chain[i], alpha);
    }
    statistical_ineff(values, N, 250); // Generate statistical data
    
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
    
    if (pr > r) { // Accept new step
        for (i = 0; i < 6; ++i)
            next[i] = tmp[i];
    } else {      // Reject new step
        for (i = 0; i < 6; ++i)
            next[i] = prev[i];
    }
    
    if (pr >= r) // Track rejection rate
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