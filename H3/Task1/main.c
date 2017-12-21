#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FOO -1
#define INITIAL_WALKERS 300
#define MAX_WALKERS INITIAL_WALKERS*10

// Function declerations
void randomize(double * v, int sz, double min, double max, gsl_rng * q);
int birthAndDeath(double walkers[MAX_WALKERS][3], gsl_rng *q,double a, double b);
void birthWalkers(double walkers[MAX_WALKERS][3], int i, int num, int deadIndicies[10]);
int isDead(double walker[3]);
int get_m(double x, gsl_rng *q, double dTau, double E_t);
int findAvailableIndex(double walkers[MAX_WALKERS][3], int start);
double W(double x, double dTau, double E_t);
double getEnergy(int numWalkers, double E_t, double dTau, double alpha);

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

    int numberOfSteps = 1000;
    int numWalkers = INITIAL_WALKERS;
    // Inefficient implementation. Obvious improvement would be to use a linked list, but
    // we didn't feel like implementing one. Missing the C++ STL more than ever...
    double walkers[MAX_WALKERS][3];
    double dTau = 0.05; // Should be in [0.01, 0.1]
    double alpha = 1; // Should be in (0,1]
    double dTauSq = sqrt(dTau);
    double energy[numberOfSteps];
    double E_0 = 0.5;
    energy[0] = E_0;

    int i, t;
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

    for(t = 1; t < numberOfSteps; t++) { 
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
        // Update energy
        energy[t] = getEnergy(numWalkers, energy[t-1], dTau, alpha);
        printf("t = %d / %d, \t E_k = %.9f, \t N_k = %d \n", t, numberOfSteps, energy[t], numWalkers);
    }
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
            index = findAvailableIndex(walkers, index+1);
            

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
    return E_t - alpha/dTau * log(numWalkers/INITIAL_WALKERS);
}