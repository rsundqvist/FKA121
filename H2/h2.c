#include "h2.h"
#include "vec3.h"
#include "stats_func.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Computing the probability of accepting a new proposal in the Markov chain
double probFunction(double * Rnew, double * Rold, double alpha) {
    return absWaveFunction(Rnew, alpha)/absWaveFunction(Rold, alpha);
}

double energy(double R[6], double alpha) {
    double R1[3] = {R[0], R[1], R[2]}; // Convert 6-vectors to
    double R2[3] = {R[3], R[4], R[5]}; // two 3-vectors.
    return energy2(R1, R2, alpha);
}

double energy2(double R1[3], double R2[3], double alpha) {
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
    double dsq = d*d;      
    
    // Energy according to eq. (7) in H2b homework description.
    return -4
    + dot(R12u, R12)/( r12 * dsq ) // d^2
    - 1/( r12 * d*dsq ) // d^3
    - 1/( 4 * dsq*dsq ) // d^4
    + 1/r12;
}

double absWaveFunction(double * R, double alpha) {
    double w = d_abs(trialWaveFunction(R, alpha));
    return w*w;
}

double trialWaveFunction(double * R, double alpha) {
    double R1[3] = {R[0],R[1],R[2]}; // Extract R1
    double R2[3] = {R[3],R[4],R[5]}; // Extract R2
    double r1 = norm(R1);
    double r2 = norm(R2);
    double r12= distance(R1, R2); // Compute r12
    return exp( -2*(r1 + r2) + r12/( 2*(1+alpha*r12) ) );
}
