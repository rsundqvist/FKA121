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

<<<<<<< HEAD
void foo(double * , double *, double (* A) [nbr_of_particles]);
=======
void transformQ(double * q, double * Q,double (*trans_matrix)[nbr_of_particles]);
double computeModeEnergy(double Qk,double Pk,double omega);
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32

/* Main program */
int main()
{
	/* Declartion of variables */
<<<<<<< HEAD
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
	timestep = 0.01;
	alpha = 0.0;
	nbr_of_timesteps = tmax/timestep;
	timestep_sq = timestep * timestep;
	int ir = 100;

	/* displacement, velocity and acceleration and energies */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles];
	double P[nbr_of_particles];
	double Q[nbr_of_particles];
	double omega[nbr_of_particles];
	double E_k0[(int)nbr_of_timesteps/ir];
	double E_k1[(int)nbr_of_timesteps/ir];
	double E_k2[(int)nbr_of_timesteps/ir];
	double E_k3[(int)nbr_of_timesteps/ir];
	double E_k4[(int)nbr_of_timesteps/ir];
	double totalEnergy[(int)nbr_of_timesteps/ir];
	double E0 = nbr_of_particles;
	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((nbr_of_timesteps)/ir * sizeof (double));
	double *q_2 = malloc((nbr_of_timesteps)/ir * sizeof (double));
	double *q_3 = malloc((nbr_of_timesteps)/ir * sizeof (double));
	double *q_4 = malloc((nbr_of_timesteps)/ir * sizeof (double));
	double *q_5 = malloc((nbr_of_timesteps)/ir * sizeof (double));

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
	for (i = 0; i < nbr_of_particles; i++) {
		q[i] = 0;
		v[i] = 0;
	}
 	E_k0[0] = E0;
	for(i = 0; i < nbr_of_particles; i++) {
  		P[i] = sqrt(2 * E0/(nbr_of_particles+1))*sin((i+1)*PI/(i+2));
	}
	foo(P,v,trans_matrix);
	foo(q,Q,trans_matrix);
	q_1[0] = v[0];
	q_2[0] = v[1];
	q_3[0] = v[2];
	q_4[0] = v[3];
	q_5[0] = v[4];
	
	// Compute omega (eigen frequencies)
	for (j = 0; j < nbr_of_particles; j++) {
  		omega[j] = 2 * sin((j+1) * PI / (2 * nbr_of_particles + 2));

	// Compute initial energy
	for(j = 1; j < nbr_of_particles-1; j++)
		totalEnergy[0] += calc_E_k(P[j],omega[j],Q[j]);
  	}
=======
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


>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
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
<<<<<<< HEAD
		if(i%ir==0) {
			int i_log = i/ir;
			q_1[i_log] = v[0];
			q_2[i_log] = v[1];
			q_3[i_log] = v[2];
			q_4[i_log] = v[3];
			q_5[i_log] = v[4];
			/* Update natural coordinates */
			foo(v, P, trans_matrix);
  			foo(q, Q, trans_matrix);
			/* Compute and store energy */
			E_k0[i_log] = calc_E_k(P[0],omega[0],Q[0]);
			E_k1[i_log] = calc_E_k(P[1],omega[1],Q[1]);
			E_k2[i_log] = calc_E_k(P[2],omega[2],Q[2]);
			E_k3[i_log] = calc_E_k(P[3],omega[3],Q[3]);
			E_k4[i_log] = calc_E_k(P[4],omega[4],Q[4]);
			for(j = 1; j < nbr_of_particles-1; j++)
				totalEnergy[i_log] += calc_E_k(P[j],omega[j],Q[j]);
		}
	}

	/* Print displacement data to output file */
	file = fopen("energy.dat","w");

	for (i = 0; i < nbr_of_timesteps/ir; i++) {
		current_time = ir * i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e \t %e \t %e \t %e", current_time, E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i], totalEnergy[i]);	
=======
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
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
		fprintf(file, "\n");
	}
	fclose(file);

<<<<<<< HEAD
	file = fopen("displacements.dat","w");
	for (i = 0; i < nbr_of_timesteps/ir; i++) {
		current_time = ir * i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i], q_4[i], q_5[i] );	
=======
	/* Print energy data to output file */
	file = fopen("energy.dat","w");

	for (i = 0; i < nbr_of_timesteps + 1; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e \t %e \t %e \t %e \t %e", current_time, Ek1[i], Ek2[i], Ek3[i], Ek4[i], Ek5[i], Ektot[i]);	
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
		fprintf(file, "\n");
	}
	fclose(file);

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	return 0;    
}

<<<<<<< HEAD
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
=======
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

>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
