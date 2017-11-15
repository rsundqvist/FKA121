/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589
#define nbr_of_particles 32


void transform(double*, double*, double (*trans_matrix)[]);
void calc_E_k(double* omega, double* P, double* Q, double* E_k, double alpha, int t, int k);
void calc_acc(double*, double*, int, double);
void init_transformation_matrix(double (*trans_matrix)[]);

/* Main program */
int main() {
	FILE *file;
    // Parameters
	double alpha = 0.1;
	double dt = 0.01;
	double t_max = 25000;
	int nbr_of_timesteps = t_max/dt;
	int ir = 100; // Resolution for i. Record every ir:th timestep

    // displacement, velocity and acceleration
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles]; 
    double *p = v; // Since m = 1

    double Q[nbr_of_particles]; // discrete
    double P[nbr_of_particles]; // discrete

    double E_k0 [nbr_of_timesteps/ir];
    double E_k1 [nbr_of_timesteps/ir];
    double E_k2 [nbr_of_timesteps/ir];
    double E_k3 [nbr_of_timesteps/ir];
    double E_k4 [nbr_of_timesteps/ir];
    double x_max[nbr_of_timesteps/ir];
    int i_max[nbr_of_timesteps/ir];
    int i,j;
    double omega[nbr_of_particles];
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    init_transformation_matrix(trans_matrix);

    // Initialize values
    E_k0[0] = nbr_of_particles;
    P[0] = sqrt(2*E_k0[0]);
    transform(P, p, trans_matrix);
    transform(Q, q, trans_matrix);    

    // Eigenfrequencies
    for (j = 0; j < nbr_of_particles; j++) {
    	omega[j] = 2*sin((j+1)*PI/(2*nbr_of_particles+2));
    }

    double sum = 0;
    // timesteps according to velocity Verlet algorithm
    calc_acc(a, q, nbr_of_particles, alpha);
    printf("\n\talpha = %.5f, t_max = %.2f \n", alpha, t_max);
    for (i = 1; i < nbr_of_timesteps + 1; i++) {
    	if (i%(nbr_of_timesteps/10) == 0) {
    		printf("\tt = %.2f\t\t\%: %.2f\n", i*dt, ((double)i/nbr_of_timesteps));
    	}
    	
    	sum += v[2];
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
    	
    	int i_log = i/ir;
    	if (i%ir == 0) {
            transform(q, Q, trans_matrix); 
            transform(p, P, trans_matrix);
            // Update normal coordinates
        	calc_E_k(omega, P, Q, E_k0, alpha, i_log, 0);
        	calc_E_k(omega, P, Q, E_k1, alpha, i_log, 1);
        	calc_E_k(omega, P, Q, E_k2, alpha, i_log, 2);
        	calc_E_k(omega, P, Q, E_k3, alpha, i_log, 3);
        	calc_E_k(omega, P, Q, E_k4, alpha, i_log, 4);

        	for (j = 0; j < nbr_of_particles; j++) {
        		calc_E_k(omega, p, q, x_max, alpha, i_log, j);
        	}
        	double max_abs = 0, max_real = 0;
        	double im = -1;
        	for (j = 0; j < nbr_of_particles; j++) {
        	    double x = a[j];
        	    if (abs(x) > max_abs) {
        	        im = j;
        	        max_abs = abs(x);
        	        max_real = x;
        	    }
        	}
        	x_max[i_log] = max_real;
        	i_max[i_log] = im;
    	}
    }

    
    file = fopen("energy.dat","w");
    if (file != NULL){
    	printf("%s", "Print to file: ");
    	for (i = 0; i < nbr_of_timesteps/ir; i++) {
    		fprintf (file,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %d \n",
    			ir*i*dt, // Time
    			E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i], x_max[i], // Energies
    			i_max[i]); // max value index (whatever is measured)
    	}
    	fclose(file);
    	printf("energy.dat created!\n");
    	
    	// Print to console
    	for (i = 0; i < nbr_of_timesteps/(ir*1000); i++) {
    		printf ("%e \t %e \t %e \t %e \t %e \t %e \t %e \t %d \n",
    			1000*ir*i*dt, // Time
    			E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i], x_max[i], // Energies
    			i_max[i], 0, 0); // max value index (whatever is measured)
    	}
    } else {
    	printf("file is NULL\n");
    }

    printf("sum = %.4f\n",sum);


}

// The formalae for Q_k and P_k are identical with m = 1
void transform(double *q, double *Q, double (*trans_matrix)[nbr_of_particles])
{	
    /* Transformation to normal modes Q from displacements q.  */
	double sum;
	int i, j;
	for (i = 0; i < nbr_of_particles; i++){
		sum = 0;
		for (j = 0; j < nbr_of_particles; j++){
			sum += q[j] * trans_matrix[i][j];
		}
		Q[i] = sum;
	}
}

void init_transformation_matrix(double (*trans_matrix)[nbr_of_particles]) {
	int i,j;
	double factor;	
	factor = 1 / ((double) nbr_of_particles + 1);
	for (i=0; i < nbr_of_particles; i++) {
		for (j=0; j < nbr_of_particles; j++) {
			trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
		}
	}
}

void calc_E_k(double* omega, double* P, double* Q, double* E_k, double alpha, int t, int k) {
	E_k[t] += 0.5 * (P[k]*P[k] + omega[k]*omega[k]*Q[k]*Q[k]);
}

void calc_acc(double *a, double *q, int size_q, double alpha){
   int i;
   q[0] = 0;
   q[size_q-1] = 0;
   for(i = 1; i < size_q-1; i++)
        a[i] = q[i+1] - 2*q[i] + q[i-1] + alpha * (q[i+1]-q[i])*(q[i+1]-q[i]) - alpha * (q[i]-q[i-1])*(q[i]-q[i-1]);
	
   a[0] = q[1] - 2 * q[0] + alpha*(q[1] - q[0])*(q[1] - q[0]) - alpha*q[0];
   a[size_q-1] = - 2 * q[size_q-1] + q[size_q-2] + alpha * q[size_q-1]*q[size_q-1] - alpha*(q[size_q-1] - q[size_q-2])*(q[size_q-1] - q[size_q-2]);

}


