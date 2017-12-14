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

double nextAlpha (double R[6], double alpha, int i, double Beta);
double energyGradient(double R[6], double alpha);

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
      
        
    gsl_rng * q = init_rng();
    //Parameters
    int chainLength = 150000;
    int nbr_of_alphas = 100;
    double d = 0.9;    
    double Beta = 0.1; // New alpha scaling parameter

    int k,l, i;
    double alpha = 0.1;
    int nbr_of_chains = 5; // Number for chains for each alpha
    double chain[chainLength][6]; 
    
    // Array containing mean local energy and std for different alphas
    double alphaArray[nbr_of_alphas][3];

    for(k = 0; k < nbr_of_alphas; k++) {
        for(l = 0; l < nbr_of_chains; l++) {
            // Initialize Markov chain
            setZero(chain, chainLength);
            randomize(chain[0],6,0,1,q);
    
            // Sample Markov chain
            generateMarkovChain(chain, alpha, d, chainLength, q);

            // Compute local energy
            int N = chainLength;
            double *values = malloc(sizeof (double[N]));
            for (i = 0; i < N; i++) {
                values[i] = energy(chain[i], alpha);
            }
    
            // Compute energy mean and standard deviation
            double meanE = get_mean(values, N);
            double stdE = get_variance(values, meanE, N);
            free(values); values = NULL;

            // Store relevant values
      
            alphaArray[k][1] += meanE;
            alphaArray[k][2] += stdE;
        }
        alphaArray[k][0] = alpha;
        alphaArray[k][1] /= nbr_of_chains;
        alphaArray[k][2] /= nbr_of_chains;
        printf("k = %d/%d: alpha = %.5f \n", k, nbr_of_alphas, alpha);

        int p = k+1; // γ_p is a scaling factor. p is denotes the iteration number.
        alpha = nextAlpha(chain[chainLength-1], alpha, p, Beta);
    }

    // Print to file
    FILE * file1 = fopen(output_file, "w");
    if (file1 != NULL){
        for (i = 0; i < nbr_of_alphas; i++) {         
            // Print file1
            fprintf (file1,"%e %e %e\n", alphaArray[i][0], alphaArray[i][1], alphaArray[i][2]);
        }
    }
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

double nextAlpha (double R[6], double alpha, int i, double Beta) {
    double A = 1.0;    
    double g = A*pow(i, -Beta); //gamma
    
    printf("    eg = %.5f\n", energyGradient(R, alpha));
    printf("    g = %.5f\n", g);
    return alpha - g*energyGradient(R, alpha);
}

double energyGradient(double R[6], double alpha) {
    double R1[3] = {R[0], R[1], R[2]};
    double R2[3] = {R[3], R[4], R[5]};
    double R1u[3];
    double R2u[3];
    unit(R1, R1u);
    unit(R2, R2u);
    double R12u[3];
    diff(R1u, R2u, R12u);
    
    double R12[3];
    diff(R1, R2, R12);
    
    double r12 = distance(R1, R2);
    double d = 1+alpha*r12;
    
    double num = -2*dot(R12u, R12)*pow(d, 2) + 3*alpha*r12+r12+3;
    double den = pow(d, 5);
    return num/den;
}

