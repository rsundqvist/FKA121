/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589
#define nbr_of_particles 32


void foo(double*, double*);
void calc_E_k(double*, double*, double*, double*, int);
void calc_acc(double*, double*, int, double);

/* Main program */
int main() {
    // Parameters
    double alpha = 0;
    double dt = 0.1;
    double t_max = 25000;
    int nbr_of_timesteps = t_max/dt;

    // displacement, velocity and acceleration
    double q[nbr_of_particles];
    double v[nbr_of_particles];
    double a[nbr_of_particles]; 
    double *p = v; // Since m = 1

    double Q[nbr_of_particles]; // discrete
    double P[nbr_of_particles]; // discrete

    P[0] = sqrt(2*nbr_of_particles);
    double E_k[5][nbr_of_timesteps];

    double omega[nbr_of_particles];
    for (int j = 0; j < nbr_of_particles; j++) {
        omega[j] = 2*sin(j*PI/(nbr_of_particles+1));
    }

    foo(Q, q);
    foo(P, p);

    /* timesteps according to velocity Verlet algorithm */
    for (int i = 1; i < nbr_of_timesteps + 1; i++) {
        calc_acc(a, q, nbr_of_particles, alpha);

        /* v(t+dt/2) */
        for (int j = 0; j < nbr_of_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        } 

        /* q(t+dt) */
        for (int j = 0; j < nbr_of_particles; j++) {
            q[j] += dt * v[j];
        }

        /* a(t+dt) */
        calc_acc(a, q, nbr_of_particles, alpha);

        /* v(t+dt) */
        for (int j = 0; j < nbr_of_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        foo(v, P);
        foo(q, Q);
        //calc_E_k(omega, P, Q, E_k, nbr_of_particles);
    }
}

// The formalae for Q_k and P_k are identical with m = 1
void foo(double *a, double *A)
{
    /* It is useful to construct the transformation matrix outside the main loop */
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    double factor = 1 / ((double) nbr_of_particles + 1);
    for (int i=0; i < nbr_of_particles; i++) {
        for (int j=0; j < nbr_of_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    
    /* Transformation to normal modes A from displacements a.  */
    double sum;
    for (int i = 0; i < nbr_of_particles; i++){
        sum = 0;
        for (int j = 0; j < nbr_of_particles; j++){
            sum += a[j] * trans_matrix[i][j];
        }
        A[i] = sum;
    }  
}
void calc_E_k(double* omega, double* P, double* Q, double* E_k, int size) {
    for (int i = 0; i < size; ++i)
        E_k[i] = 0.5 * (P[i]*P[i] + omega[i]*omega[i]*Q[i]*Q[i]);
}

void calc_acc(double *a, double *q, int size_q, double alpha){
    
     // Boundary conditions
    int i = 0;
    a[i] = q[i+1]-2*q[i]+q[i-1]+alpha*(-q[i+1]+2*q[i]*(q[i+1])); ;
    i = size_q-1;
    a[i] = 2*q[i]+q[i-1]+alpha*(q[i-1]*q[i-1]+2*q[i]*(-q[i-1]));


    for(i = 1; i < size_q-1; i++)
        a[i] = q[i+1]-2*q[i]+q[i-1]+alpha*(q[i-1]*q[i-1]-q[i+1]+2*q[i]*(q[i+1]-q[i-1])); 
}