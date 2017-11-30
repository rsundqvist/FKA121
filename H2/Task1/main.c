#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stat.h"
#include "vec3.h"

// Function declerations
double trialWaveFunction(double * r1, double * r2, double alpha);
double propFunction(double * r1new, double * r1old, double * r2new, double * r2old, double alpha);

int main()
{
    double v1[3] = {1, 0, 1};
    double v2[3] = {1, 2, 3};
    double x[3] = {1, 0, 0};
    double y[3] = {0, 1, 0};
    
    double u[3];
    double c[3];
    double d[3];
    double s[3];
    
    copy(v1, c);
    unit(c, u);
    diff(v1, v2, d);
    sum(v1, v2, s);
    return 0;
}

double trialWaveFunction(double * r1, double * r2, double alpha) {
    double r1_norm = norm(r1);
    double r2_norm = norm(r2);
    double r12[3];
    diff(r1,r2,r12);
    double r12_norm = norm(r12);

    return exp(-2*(r1_norm + r2_norm))*(r12_norm/(2*(1 + alpha*r12_norm)));
}

// Computing the probability of accepting a new proposal in the Markov chain
double propFunction(double * r1new, double * r1old, double * r2new, double * r2old, double alpha){
    return trialWaveFunction(r1new,r2new,alpha)/trialWaveFunction(r1old,r2old,alpha);
}
