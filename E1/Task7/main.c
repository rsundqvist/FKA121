/*
E1_main.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2014-10-20.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#include "fft_func.h"
#define PI 3.141592653589
//#define nbr_of_timesteps 32767 /* nbr_of_timesteps-1 = power of 2, for best speed */
#define nbr_of_timesteps 262143
//#define nbr_of_timesteps 7 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_particles 3 /* The number of particles is 3 */

/* Main program */
int main()
{
	/* Declartion of variables */
	double timestep;
	int i,j;
	double timestep_sq,current_time;
	double m[nbr_of_particles];
	double kappa;

	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles]; 


	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	long size = (nbr_of_timesteps+1) * sizeof (double);
	double *q_1 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps+1) * sizeof (double));
	double *Ek = malloc((nbr_of_timesteps+1) * sizeof (double)); // Kinetic
	double *Ep = malloc((nbr_of_timesteps+1) * sizeof (double)); // Potential
	double *powspec_data = malloc((nbr_of_timesteps+1) * sizeof (double)); // Power spec
	double *freq = malloc((nbr_of_timesteps+1) * sizeof (double)); // Freq

	/* Set variables */
	timestep = 0.00001;
	kappa = 160/1.6021765;
	timestep_sq = timestep * timestep;
	m[0] = m[2] = 1.0364 * 0.0001 * 15.9994; //mass of oxygen
	m[1] = 1.0364 * 0.0001 * 12.0107; //mass of carbon
	printf("m[0] = %e, m[1] = %e, m[2] = %e\n", m[0], m[1], m[2]);
	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	q[0] = 0.1;
	v[0] = 0;

	for (i = 1; i < nbr_of_particles; i++) {
		q[i] = 0;
		v[i] = 0;
	}
	q_1[0] = q[0];
	q_2[0] = q[1];
	q_3[0] = q[2];
	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, m, kappa, nbr_of_particles);
	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) { //for (j = 0; j < nbr_of_particles; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, nbr_of_particles);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = q[0];
		q_2[i] = q[1];
		q_3[i] = q[2];

		Ep[i] = calc_pe(q, kappa, nbr_of_particles); // Store potential
		Ek[i] = calc_ke(v, nbr_of_particles, m); // Store kinetic
	}

	/* Print displacement data to output file */
	file = fopen("energy.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;

		double Ept = Ep[i];
		double Ekt = Ek[i];
		double Et = Ept + Ekt;

		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, Et, Ept, Ekt);	
		fprintf(file, "\n");
	}
	printf("Energy printed.\n");
	fclose(file);

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	printf("Displacement printed.\n");
	fclose(file);

 	/* make FFT (powerspectrum) */
	powerspectrum(q_1, powspec_data, nbr_of_timesteps);
	powerspectrum_shift(powspec_data, nbr_of_timesteps);
	fft_freq_shift(freq, timestep, nbr_of_timesteps);
	//fft_freq(freq, dt, n);

	/*Save powerspectrum data in file */
    file = fopen("powerspectrum.dat","w");
	for (i = 0; i < nbr_of_timesteps; i++) {
		fprintf (file,"%e \t %e\n", freq[i], powspec_data[i]);
	}
	printf("Power spectrum printed.\n");

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	free(Ek); Ek = NULL;
	free(Ep); Ep = NULL;
	free(powspec_data); powspec_data = NULL;
	free(freq); freq = NULL;
	return 0;    
}
