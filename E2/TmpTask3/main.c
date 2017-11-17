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
double computeModeEnergy(double Qk,double Pk,double omega);

/* Main program */
int main()
{
	/* Declartion of variables */
	double alpha;
	double timestep;
	int i,j,k;
	double timestep_sq,current_time;
	double tmax;
	int nbr_of_timesteps;

	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles];
	double Q[nbr_of_particles];
	double P[nbr_of_particles];
	double omega[nbr_of_particles];
	double *p = v;

	/* Set variables */
	timestep = 0.01;
	tmax = 25000;
	nbr_of_timesteps = (int)tmax/timestep;
	alpha = 0.1;

	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek3 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek4 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek5 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ektot = malloc((nbr_of_timesteps+1) * sizeof (double));

	// Determine transformation matrix
    	double factor;
    	double trans_matrix[nbr_of_particles][nbr_of_particles];
	factor = 1 / ((double) nbr_of_particles + 1);
 	for (i=0; i < nbr_of_particles; i++) {
        	for (j=0; j < nbr_of_particles; j++) {
            		trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        	}
    	}
    
	// Compute omega
	for(i = 0; i < nbr_of_particles; i++) { 
		k = i+1;
		omega[i] = 2.0 * sin(k * PI /(2*nbr_of_particles+2));
	}
	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	

	for (i = 0; i < nbr_of_particles; i++) {
		Q[i] = 0;
		P[i] = 0;
	}
	P[0] = sqrt(2*nbr_of_particles);
	transformQ(Q, q,trans_matrix);
	transformQ(P, p,trans_matrix);
	q_1[0] = a[0];
	q_2[0] = a[1];
	q_3[0] = a[2];

	/* -------------------------- Start of Loop -------------------*/


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
		q_1[i] = a[0];
		q_2[i] = a[1];
		q_3[i] = a[2];


		/* Update normal coordiantes and compute energy */
		transformQ(q, Q,trans_matrix);
		transformQ(p, P,trans_matrix);
		Ek1[i] = computeModeEnergy(Q[0],P[0],omega[0]);
		Ek2[i] = computeModeEnergy(Q[1],P[1],omega[1]);
		Ek3[i] = computeModeEnergy(Q[2],P[2],omega[2]);
		Ek4[i] = computeModeEnergy(Q[3],P[3],omega[3]);
		Ek5[i] = computeModeEnergy(Q[4],P[4],omega[4]);
		/* Total energy */
        for (j = 0; j < nbr_of_particles; j++)
		    Ektot[i] += computeModeEnergy(Q[j], P[j], omega[j]);
	}

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	fclose(file);

	/* Print energy data to output file */
	file = fopen("energy.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e \t %e \t %e \t %e \t %e", current_time, Ek1[i], Ek2[i], Ek3[i], Ek4[i], Ek5[i], Ektot[i]);	
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
            		sum += q[j] * trans_matrix[i][j];
        	}
        	Q[i] = sum;
    }
}


double computeModeEnergy(double Qk,double Pk,double omega) {
	return 0.5*(Pk*Pk + omega*omega*Qk*Qk);
}

