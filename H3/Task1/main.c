#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../stat.h"

#define INITIAL_WALKERS 300
#define MAX_WALKERS INITIAL_WALKERS*5
#define NUMBER_OF_STEPS 100

// Function declerations
void randomize(double * v, int sz, double min, double max, gsl_rng * q);
int birthAndDeath(double walkers[MAX_WALKERS][3], gsl_rng *q,double a, double b);
void birthWalkers(double walkers[MAX_WALKERS][3], int i, int num);
int isDead(double walker[3]);
int get_m(double x, gsl_rng *q, double dTau, double E_t);
int findAvailableIndex(double walkers[MAX_WALKERS][3], int start);
double W(double x, double dTau, double E_t);
double getEnergy(int numWalkers, double E_t, double dTau, double alpha);
int simulate(double energy[NUMBER_OF_STEPS], int walkerCount[NUMBER_OF_STEPS],
    double walkers[MAX_WALKERS][3], gsl_rng *q,
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
    char output_file2[64];
    sprintf(output_file2, "position.dat");    
        
    gsl_rng *q = init_rng();
    // Inefficient implementation. Obvious improvement would be to use a linked list, but
    // we didn't feel like implementing one. Missing the C++ STL more than ever...
    double walkers[MAX_WALKERS][3];
    double dTau = 0.05; // Should be in [0.01, 0.1]
    double alpha = 0.3; // Should be in (0,1]
    double dTauSq = sqrt(dTau);
    double energy[NUMBER_OF_STEPS];
    int walkerCount[NUMBER_OF_STEPS];
    
    // Run simulations
    int trial, i;
    int numTrials = 1000;
    
    double finalEnergy[numTrials];
    double finalNwalkers[numTrials];
    double mean_energy, var_energy;
    double mean_nw, var_nw;
    
    printf("Begin trials: dTau = %.3f, alpha = %.3f\n", dTau, alpha);
    FILE * file = fopen(output_file, "w");
    FILE * file2 = fopen(output_file2, "w");
    
    for (trial = 0; trial < numTrials; trial++) {

        finalNwalkers[trial] = simulate(energy, walkerCount, walkers, q, dTau, dTauSq, alpha);
        finalEnergy[trial] = energy[NUMBER_OF_STEPS-1];

        for (i = 0; i < NUMBER_OF_STEPS; ++i)
        {    
            fprintf (file,"%e %d ", energy[i], walkerCount[i]);
        }
        fprintf (file,"\n");

        for (i = 0; i < MAX_WALKERS; ++i)
        {    
            if (!isDead(walkers[i]))
                fprintf (file2,"%e ", walkers[i][0]);
        }

        if ((trial+1)%(numTrials/10) == 0)
            printf("trial = %d/%d, percent = %.3f\n",
                (trial+1), numTrials, (double)(trial+1)/numTrials);

        fprintf (file2,"\n");
    }
    mean_energy = get_mean(finalEnergy, numTrials);
    var_energy = get_variance(finalEnergy, mean_energy, numTrials);
    mean_nw = get_mean(finalNwalkers, numTrials);
    var_nw = get_variance(finalNwalkers, mean_nw, numTrials);
    
    printf("mean_energy = %.5f, var_energy = %.5f\n", mean_energy, var_energy);
    printf("mean_nw = %.5f, var_nw = %.5f\n", mean_nw, var_nw);
    printf("Done.\n");
}

int simulate(double energy[NUMBER_OF_STEPS], int walkerCount[NUMBER_OF_STEPS],
    double walkers[MAX_WALKERS][3], gsl_rng *q,
    double dTau, double dTauSq, double alpha) {
    
    double E_0 = 0.5;
    energy[0] = E_0;
    
    int i, t;
    int numWalkers = INITIAL_WALKERS;
    walkerCount[0] = numWalkers;
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
        // Update positions 
        for (i = 0; i < MAX_WALKERS; ++i)
        {
            if(!isDead(walkers[i]))
            {
                walkers[i][0] += dTauSq*gsl_ran_ugaussian(q);
            }
        }

        // Birth/death process
        numWalkers += birthAndDeath(walkers, q, dTau, energy[t-1]);
        walkerCount[t] = numWalkers;

        // Update energy
        energy[t] = getEnergy(numWalkers, energy[t-1], dTau, alpha);
    }

    return numWalkers;
}

int birthAndDeath(double walkers[MAX_WALKERS][3], gsl_rng *q, double dTau, double E_t)
{
    int i, dbd = 0; // delta births/deaths
    for (i = 0; i < MAX_WALKERS; ++i)
    {
        if (isDead(walkers[i]))
        {
            // Do nothing.
        }
        else
        {
            int m = get_m(walkers[i][0], q, dTau, E_t);
            if (m == 0)
            {
                walkers[i][0] = NAN; // Kill the walker
                dbd--;
            }
            else if (m > 1)
            {
                m--;
                birthWalkers(walkers, i, m);
                dbd += m;
            }
        }
    }
    return dbd; // Return the change in number of walkers
}


void birthWalkers(double walkers[MAX_WALKERS][3], int i, int num)
{
    int index = i;

    double x = walkers[i][0]; // copy x position
    double y = walkers[i][1]; // z
    double z = walkers[i][2]; // y

    while (num > 0)
    {
        index = findAvailableIndex(walkers, 0);
        walkers[index][0] = x;
        walkers[index][1] = y;
        walkers[index][2] = z;

        num--;
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
    exit(-1);
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
    return exp( -dTau*(0.5*x*x - E_t) ); // Eq. 60
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
