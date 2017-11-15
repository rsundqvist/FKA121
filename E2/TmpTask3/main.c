/*
E1_main.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2014-10-20.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#define PI 3.141592653589
#define nbr_of_particles 32 /* The number of particles is 3 */

void foo(double * , double *, double (* A) [nbr_of_particles]);

/* Main program */
int main()
{
	/* Declartion of variables */
	double timestep;
	double tmax;
	double nbr_of_timesteps;
	int i,j;
	double timestep_sq,current_time;
	double alpha;
	/* declare file variable */
	FILE *file;

	/* Set variables */
	tmax = 25000;
	timestep = 0.1;
	alpha = 0.1;
	nbr_of_timesteps = tmax/timestep;
	timestep_sq = timestep * timestep;

	/* displacement, velocity and acceleration and energies */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles];
	double P[nbr_of_particles];
	double Q[nbr_of_particles];
	double omega[nbr_of_particles];
	double E_k0[(int)nbr_of_timesteps];
	double E_k1[(int)nbr_of_timesteps];
	double E_k2[(int)nbr_of_timesteps];
	double E_k3[(int)nbr_of_timesteps];
	double E_k4[(int)nbr_of_timesteps];

	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_4 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_5 = malloc((nbr_of_timesteps+1) * sizeof (double));

	  // Transformation matrix
  	double factor;
  	double trans_matrix[nbr_of_particles][nbr_of_particles];
  	factor = 1 / ((double) nbr_of_particles + 1);
  	for (i = 0; i < nbr_of_particles; i++) {
    		for (j = 0; j < nbr_of_particles; j++) {
      			trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
    		}
 	 }
	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	for (i = 1; i < nbr_of_particles; i++) {
		q[i] = 0;
		v[i] = 0;
	}
 	E_k0[0] = nbr_of_particles;
  	P[0] = sqrt(2 * E_k0[0]);
  	foo(P, v, trans_matrix);
  	foo(Q, q, trans_matrix);
	q_1[0] = q[0];
	q_2[0] = q[1];
	q_3[0] = q[2];
	q_4[0] = q[3];
	q_5[0] = q[4];
	
	// Compute omega (eigen frequencies)
	for (j = 0; j < nbr_of_particles; j++) {
  		omega[j] = 2 * sin(j * PI / (2 * nbr_of_particles + 2));
  	}
	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, alpha, nbr_of_particles);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, alpha, nbr_of_particles);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = q[0];
		q_2[i] = q[1];
		q_3[i] = q[2];
		q_4[i] = q[3];
		q_5[i] = q[4];
		/* Update natural coordinates */
		foo(v, P, trans_matrix);
  		foo(q, Q, trans_matrix);
		/* Compute and store energy */
		E_k0[i] = calc_E_k(P[0],omega[0],Q[0]);
		E_k1[i] = calc_E_k(P[1],omega[1],Q[1]);
		E_k2[i] = calc_E_k(P[2],omega[2],Q[2]);
		E_k3[i] = calc_E_k(P[3],omega[3],Q[3]);
		E_k4[i] = calc_E_k(P[4],omega[4],Q[4]);
	}

	/* Print displacement data to output file */
	file = fopen("energy.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e", current_time, E_k0[i], E_k1[i]);	
		fprintf(file, "\n");
	}
	fclose(file);

	file = fopen("displacements.dat","w");
	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i], q_4[i], q_5[i] );	
		fprintf(file, "\n");
	}
	fclose(file);

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	return 0;    
}

/* The formalae for Q_k and P_k are identical with m = 1 */
void foo(double * q, double * Q, double (* trans_matrix) [nbr_of_particles]) {
  int i, j;
  /* Transformation to normal modes Q from displacements q.  */
  double sum;
  for (i = 0; i < nbr_of_particles; i++) {
    sum = 0;
    for (j = 0; j < nbr_of_particles; j++) {
      sum += q[j] * trans_matrix[i][j];
    }
    Q[i] = sum;
  }
}
