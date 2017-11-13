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
void calc_E_k(double omega, double P, double Q, double* E_k, int i);
void calc_acc(double*, double*, int, double);

/* Main program */
int main() {
    FILE *file;
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

    double E_k0[nbr_of_timesteps];
    double E_k1[nbr_of_timesteps];
    double E_k2[nbr_of_timesteps];
    double E_k3[nbr_of_timesteps];
    double E_k4[nbr_of_timesteps];
    int i,j;
    double omega[nbr_of_particles];

    // Initialize values
    E_k0[0] = nbr_of_particles;
    P[0] = sqrt(2*E_k0[0]);
    

    for (j = 0; j < nbr_of_particles; j++) {
        omega[j] = 2*sin(j*PI/(2*nbr_of_particles+2));
    }

    foo(Q, q);
    foo(P, p);
    // timesteps according to velocity Verlet algorithm
    for (i = 1; i < nbr_of_timesteps + 1; i++) {
        if (i%5000 == 0) {
            printf("\tt = %.2f\t\t\%: %.2f\n", i*dt, ((double)i/nbr_of_timesteps));
        }
        
        calc_acc(a, q, nbr_of_particles, alpha);

        // v(t+dt/2)
        for (j = 0; j < nbr_of_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        } 

        // q(t+dt) 
        for (j = 0; j < nbr_of_particles; j++) {
            q[j] += dt * v[j];
        }

        // a(t+dt)
        calc_acc(a, q, nbr_of_particles, alpha);

        // v(t+dt)
        for (j = 0; j < nbr_of_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }

        foo(p, P);
        foo(q, Q);
        
        calc_E_k(omega[0], P[0], Q[0], E_k0, i);
        calc_E_k(omega[1], P[1], Q[1], E_k1, i);
        calc_E_k(omega[2], P[2], Q[2], E_k2, i);
        calc_E_k(omega[3], P[3], Q[3], E_k3, i);
        calc_E_k(omega[4], P[4], Q[4], E_k4, i);
    }

    
    file = fopen("energy.dat","w");
        if (file != NULL){
        printf("%s", "Print to file: ");
        for (i = 0; i < nbr_of_timesteps; i++) {
            //fprintf (file,"%e \t %e \t %e \n", i*dt, E_k0[i], E_k1[i]);
            fprintf (file,"%e \t %e \t %e \t %e \t %e \t %e \n",
                i*dt, E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i]);
        }
        fclose(file);
        printf("energy.dat created!\n");
    } else {
        printf("file is NULL");
    }
}

// The formalae for Q_k and P_k are identical with m = 1
void foo(double *a, double *A)
{
    int i,j;
    /* It is useful to construct the transformation matrix outside the main loop */
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    double factor = 1 / ((double) nbr_of_particles + 1);
    for (i=0; i < nbr_of_particles; i++) {
        for (j=0; j < nbr_of_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    
    /* Transformation to normal modes A from displacements a.  */
    double sum;
    for (i = 0; i < nbr_of_particles; i++){
        sum = 0;
        for (j = 0; j < nbr_of_particles; j++){
            sum += a[j] * trans_matrix[i][j];
        }
        A[i] = sum;
    }  
}

void calc_E_k(double omega, double P, double Q, double* E_k, int i) {
        E_k[i] = 0.5 * (P*P + omega*omega*Q*Q);
}

void calc_acc(double *a, double *q, int size_q, double alpha){
    
     // Boundary conditions
    int i = 0;
    a[i] = q[i+1]-2*q[i]+alpha*(-q[i+1]+2*q[i]*(q[i+1])); ;
    i = size_q-1;
    a[i] = 2*q[i]+q[i-1]+alpha*(q[i-1]*q[i-1]+2*q[i]*(-q[i-1]));

    for(i = 1; i < size_q-1; i++)
        a[i] = q[i+1]-2*q[i]+q[i-1]+alpha*(q[i-1]*q[i-1]-q[i+1]+2*q[i]*(q[i+1]-q[i-1])); 
}
