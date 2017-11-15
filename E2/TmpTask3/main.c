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

void transformQ(double * q, double * Q,double (*trans_matrix)[nbr_of_particles]);
void inverseTransformQ(double * q, double * Q,double (*trans_matrix)[nbr_of_particles]);

/* Main program */
int main()
{
	/* Declartion of variables */
	double alpha;
	double timestep;
	int i,j,k;
	double timestep_sq,current_time;
	double size_q;
	double tmax;
	int nbr_of_timesteps;

	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles+2];
	double v[nbr_of_particles+2];
	double a[nbr_of_particles+2];
	double Q[nbr_of_particles+2];
	double P[nbr_of_particles+2];
	double *p = v;

	/* Set variables */
	timestep = 0.1;
	tmax = 25000;
	nbr_of_timesteps = (int)tmax/timestep;
	alpha = 0.1;
	size_q = nbr_of_particles + 2;

	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));

	// Determine transformation matrix
    	double factor;
    	double trans_matrix[nbr_of_particles][nbr_of_particles];
	factor = 1 / ((double) nbr_of_particles + 1);
 	for (i=0; i < nbr_of_particles; i++) {
        	for (j=0; j < nbr_of_particles; j++) {
            		trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        	}
    	}
    

	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	

	for (i = 0; i < nbr_of_particles; i++) {
		Q[i] = 0;
		P[i] = 0;
	}
	P[0] = sqrt(2*nbr_of_particles);
	inverseTransformQ(q, Q,trans_matrix);
	inverseTransformQ(p, P,trans_matrix);
	q_1[0] = v[0];
	q_2[0] = v[1];
	q_3[0] = v[2];

	/* -------------------------- Start of Loop -------------------*/


	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, alpha, size_q);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < size_q; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < size_q; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, alpha, size_q);

		/* v(t+dt) */
		for (j = 0; j < size_q; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = v[0];
		q_2[i] = v[1];
		q_3[i] = v[2];
	}

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	fclose(file);

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	return 0;    
}

void transformQ(double * q, double * Q,double (*trans_matrix)[nbr_of_particles]) {
	/* Transformation to normal modes Q from displacements q.  */
	int i,j;
	double sum;
 	for (i = 0; i < nbr_of_particles; i++){
        	sum = 0;
        	for (j = 0; j < nbr_of_particles; j++){
            		sum += q[j+1] * trans_matrix[i][j];
        	}
        	Q[i] = sum;
    }
}

void inverseTransformQ(double * q, double * Q,double (*trans_matrix)[nbr_of_particles]) {
	/* Transformation from normal modes Q to displacements q.  */
	int i,j;
	double sum;
 	for (i = 0; i < nbr_of_particles; i++){
        	sum = 0;
        	for (j = 0; j < nbr_of_particles; j++){
            		sum += Q[j] * trans_matrix[i][j];
        	}
        	q[i+1] = sum;
    }
}
