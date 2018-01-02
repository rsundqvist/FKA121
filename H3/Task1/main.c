#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../stat.h"

#define FOO -1
#define INITIAL_WALKERS 500
#define MAX_WALKERS INITIAL_WALKERS*5
#define NUMBER_OF_STEPS 10000

// Function declerations
void randomize(double * v, int sz, double min, double max, gsl_rng * q);
int birthAndDeath(double walkers[MAX_WALKERS][3], gsl_rng *q,double a, double b);
void birthWalkers(double walkers[MAX_WALKERS][3], int i, int num, int deadIndicies[10]);
int isDead(double walker[3]);
int get_m(double x, gsl_rng *q, double dTau, double E_t);
int findAvailableIndex(double walkers[MAX_WALKERS][3], int start);
double W(double x, double dTau, double E_t);
double getEnergy(int numWalkers, double E_t, double dTau, double alpha);
int simulate(double energy[NUMBER_OF_STEPS], double walkers[MAX_WALKERS][3], gsl_rng *q,
    double dTau, double dTauSq, double alpha);

// Global variables
int metropolisCount = 0; // Used for counting the acceptance rate of Metropolis
int metropolisTotal = 0;

// gsl stuff
gsl_rng * init_rng();

int main()
{
    // Output file
    char output_file[64];
    sprintf(output_file, "simulation.dat");    
        
    gsl_rng *q = init_rng();
    // Inefficient implementation. Obvious improvement would be to use a linked list, but
    // we didn't feel like implementing one. Missing the C++ STL more than ever...
    double walkers[MAX_WALKERS][3];
    double dTau = 0.1; // Should be in [0.01, 0.1]
    double alpha = 0.01; // Should be in (0,1]
    double dTauSq = sqrt(dTau);
    double energy[NUMBER_OF_STEPS];
    
    // Run simulations
    int trial;
    int numTrials = 100;
    
    double finalEnergy[numTrials];
    double finalNwalkers[numTrials];
    double mean_energy, var_energy;
    double mean_nw, var_nw;
    
    printf("Begin trials: dTau = %.3f, alpha = %.3f\n", dTau, alpha);
    // Simulation numTrials simulations for given values
    for (trial = 0; trial < numTrials; trial++) {
        printf("trial = %d\n", trial);
        
        finalNwalkers[trial] = simulate(energy, walkers, q, dTau, dTauSq, alpha);
        finalEnergy[trial] = energy[NUMBER_OF_STEPS-1];
    }
    mean_energy = get_mean(finalEnergy, numTrials);
    var_energy = get_variance(finalEnergy, mean_energy, numTrials);
    mean_nw = get_mean(finalNwalkers, numTrials);
    var_nw = get_variance(finalNwalkers, mean_nw, numTrials);
    
    printf("mean_energy = %.5f, var_energy = %.5f\n", mean_energy, var_energy);
    printf("mean_nw = %.5f, var_nw = %.5f\n", mean_nw, var_nw);
    printf("Done.\n");
}

int simulate(double energy[NUMBER_OF_STEPS], double walkers[MAX_WALKERS][3], gsl_rng *q,
    double dTau, double dTauSq, double alpha) {
    
    double E_0 = 0.5;
    energy[0] = E_0;
    
    int i, t;
    int numWalkers = INITIAL_WALKERS;
    // Initialize
    for (i = 0; i < INITIAL_WALKERS ; ++i)
    {
        randomize(walkers[i], 3, -1, 1, q);
        walkers[i][1] = 0; // 1D case - set yz-coordinates to zero.
        walkers[i][2] = 0;
    }
    for (i = INITIAL_WALKERS ; i < MAX_WALKERS; ++i)
    {
        walkers[i][0] = NAN; // NaN values signify a "dead" walker.
        walkers[i][1] = 0;
        walkers[i][2] = 0;
    }

    // Begin simulation
    for(t = 1; t < NUMBER_OF_STEPS; t++) { 
        // Update energy
        energy[t] = getEnergy(numWalkers, energy[t-1], dTau, alpha);
        // Update positions
        for (i = 0; i < MAX_WALKERS; ++i)
        {
            if(!isDead(walkers[i]))
            {
                walkers[i][0] += dTauSq*gsl_ran_ugaussian(q);
            }
        }
        // Birth/death process
        numWalkers += birthAndDeath(walkers, q, dTau, energy[t]);
        //printf("t = %d\n", t);
        //printf("t = %d / %d, \t E_k = %.9f, \t N_k = %d \n", t, numberOfSteps, energy[t], numWalkers);
    }
    return numWalkers;
}

int birthAndDeath(double walkers[MAX_WALKERS][3], gsl_rng *q, double dTau, double E_t)
{
    // Keep track of available indices
    int deadIndicies[10] = {
        FOO, FOO, FOO, FOO, FOO,
        FOO, FOO, FOO, FOO, FOO
    };

    int i, j, dbd = 0; // delta births/deaths
    for (i = 0; i < MAX_WALKERS; ++i)
    {
        if (isDead(walkers[i]))
        {
            // Add to list of available walker spots
            for (j = 0; j < 10; ++j)
                if (deadIndicies[j] == FOO) deadIndicies[j] = i;
        }
        else
        {
            int m = get_m(walkers[i][0], q, dTau, E_t);
            if (m == 0)
            {
                // Kill the walker
                walkers[i][0] = NAN;
                for (j = 0; j < 10; ++j)
                    if (deadIndicies[j] == i) deadIndicies[j] = FOO; // Mark as available
                dbd--;
            }
            else if (m > 1)
            {
                m = m-1;
                birthWalkers(walkers, i, m, deadIndicies);
                dbd += m;
            }
        }
    }
    return dbd; // Return the change in number of walkers
}


void birthWalkers(double walkers[MAX_WALKERS][3], int i, int num, int deadIndicies[10])
{
    int cdip = 0; // Current dead indices position
    int index = i;

    double pos = walkers[i][0]; // copy x position

    while (num > 0)
    {
        num--;

        // Suitable position in deadIndices, if there is one
        while(cdip < 10 && deadIndicies[cdip] != FOO)
            cdip++;

        if (cdip < 10)
            index = deadIndicies[cdip];
        else
            index = findAvailableIndex(walkers, 0*(index+1)); // TODO: Don't always search from 0
            

        walkers[index][0] = pos;
    }
}

int findAvailableIndex(double walkers[MAX_WALKERS][3], int start)
{
    int i;
    for (i = start; i < MAX_WALKERS; ++i)
    {
        if (isDead(walkers[i])) return i;
    }

    printf("FAILED TO FIND AVAILABLE INDEX - MAX_WALKERS TOO SMALL!!!\n");
    printf("FAILED TO FIND AVAILABLE INDEX - MAX_WALKERS TOO SMALL!!!\n");
    printf("FAILED TO FIND AVAILABLE INDEX - MAX_WALKERS TOO SMALL!!!\n");
    printf("FAILED TO FIND AVAILABLE INDEX - MAX_WALKERS TOO SMALL!!!\n");
    printf("FAILED TO FIND AVAILABLE INDEX - MAX_WALKERS TOO SMALL!!!\n");
    return -1; // Should never happen!
}

int isDead(double walker[3])
{
    return isnan(walker[0]) || isnan(walker[1]) || isnan(walker[2]);
}

int get_m(double x, gsl_rng *q, double dTau, double E_t)
{
    return (int)(W(x, dTau, E_t) + gsl_rng_uniform(q)); // Eq. 62
}

double W(double x, double dTau, double E_t)
{
    return exp(-dTau*(0.5*x*x - E_t)); // Eq. 60
}


void randomize(double * v, int sz, double min, double max, gsl_rng *q) {
    int i;
    for(i = 0; i < sz; i++)
        v[i] = min + (max-min) * gsl_rng_uniform(q);
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

double getEnergy(int numWalkers, double E_t, double dTau, double alpha) {
    return E_t - alpha/dTau * log((double)numWalkers/INITIAL_WALKERS);
}
