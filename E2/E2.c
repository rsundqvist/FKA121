/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589
#define nbr_of_particles 32

void Q_mat(double*, double*);

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

    for (int i = 0; i < nbr_of_timesteps; i++) {
        double t = i*dt;

    }


    //TODO
    double Q[nbr_of_particles];
    Q_mat(q, Q);
}

void Q_mat(double *q, double *Q)
{
    //TODO
    /* It is useful to construct the transformation matrix outside the main loop */
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    double factor = 1 / ((double) nbr_of_particles + 1);
    for (int i=0; i < nbr_of_particles; i++) {
        for (int j=0; j < nbr_of_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    
    /* Transformation to normal modes Q from displacements q.  */
    double sum;
    for (int i = 0; i < nbr_of_particles; i++){
        sum = 0;
        for (int j = 0; j < nbr_of_particles; j++){
            sum += q[j] * trans_matrix[i][j];
        }
        Q[i] = sum;
    }  
}

void calc_acc(double *a, double *q, int size_q, double alpha){
    a[0] = 0;
    a[size_q-1] = 0; // Boundary conditions

    for(int i = 1; i < size_q-1; i++)
        a[i] = q[i+1]-2*q[i]+q[i-1]+alpha*(q[i-1]*q[i-1]-q[i+1]+2*q[i]*(q[i+1]-q[i-1])); 
}