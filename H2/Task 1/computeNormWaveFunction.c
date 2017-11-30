//This program computes the norm for our trial wave function using Metropolis algorithm
#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Function declerations
double trialWavefunction(double * r, double alpha);
void metropolisIntegrate6D(double *ans, double (*gFun)(double[6]), double(*samplePoints)[6], int N);
void metropolisStep(double prev[6], double next[6], double Delta, double (*probFunc)(double[6],double[6]), gsl_rng * q);
void generateMarkovChain(double (*chain)[6], double (*probFunc)(double[6],double[6]), double Delta, int N, gsl_rng * q);

// GSL
gsl_rng * init_rng(); // gsl rng create


int main() {


gsl_rng * q = init_rng();
return 1;
}

double weight(double * r1, double * r2, double norm) {

}

double trialWavefunction(double * r, double alpha) {
    double r1[3];
    r1[0] = r[0];
    r1[1] = r[1];
    r1[2] = r[2];
    double r2[3];
    r2[0] = r[3];
    r2[1] = r[4];
    r2[2] = r[5];    
    double r1_norm = norm(r1);
    double r2_norm = norm(r2); 
    double r12[3];
    diff(r1,r2,r12);
    double r12_norm = norm(r12);
    
    return exp(-2*r1_norm - 2*r2_norm)*exp(r12_norm/(2*(1 + alpha*r12_norm)));
}

void metropolisIntegrate6D(double *ans, double (*gFun)(double[6]), double(*samplePoints)[6], int N) {
    double f;
    double fsum = 0, fsqSum = 0;

    int i;    
    for(i = 0; i < N; i++) {
        f = gFun(samplePoints[i]);
        fsum += f;
        fsqSum += f*f;
    }
    double mu = fsum/N;
    double sigmasq = 1/N * fsqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/N);
}

void generateMarkovChain(double (*chain)[6], double (*probFunc)(double[6],double[6]), double Delta, int N, gsl_rng * q) {
    if (chain[0][0] == 0 && chain[0][0] == 0 && chain[0][0] == 0) {
        printf("The initial position is set to (0, 0, 0). Is this correct?\n");
    }
    
    int i;

    for(i = 1; i< N; i++) {
        metropolisStep(chain[i-1], chain[i], Delta, probFunc, q);
    }
}

void metropolisStep(double prev[6], double next[6], double Delta, double (*probFunc)(double[6],double[6]), gsl_rng * q) {
    int i;
    double tmp[3];
    for (i = 0; i < 3; ++i) {
        //printf("prev[i] = %.5f\n", prev[i]);
        tmp[i] = nextCoordinate(prev[i], Delta, q);
    }
    double pr = probFunc(prev, tmp);
    double r = gsl_rng_uniform(q);
    
    if (pr > r) {// Accept new step
        for (i = 0; i < 3; ++i)
            next[i] = tmp[i];
    } else {
        for (i = 0; i < 3; ++i)
            next[i] = prev[i];
    }
    
    if (pr >= r)        
        count++;
    total++;
    //printf("%.3f (%.6f,%.6f,%.6f) \n\n", pr,next[0],next[1],next[2]);
}
