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
	double alpha = 0.01;
	double dt = 0.1;
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
    double E_tot[nbr_of_timesteps/ir];
    int i,j;
    double omega[nbr_of_particles];
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    init_transformation_matrix(trans_matrix);

    // Initialize values
    E_k0[0] = nbr_of_particles;
    P[0] = sqrt(2*E_k0[0]);
    transform(P, p, trans_matrix);
    transform(Q, q, trans_matrix);    

    for (j = 0; j < nbr_of_particles; j++) {
    	omega[j] = 2*sin(j*PI/(2*nbr_of_particles+2));
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
        		calc_E_k(omega, p, q, E_tot, alpha, i_log, j);
        	}
        	double max = 0;
        	for (j = 0; j < nbr_of_particles; j++) {
        	    double v = abs(q[j]);
        	    if (v > max) max = v;
        	}
        	E_tot[i_log] = max;
    	}
    }

    
    file = fopen("energy.dat","w");
    if (file != NULL){
    	printf("%s", "Print to file: ");
    	for (i = 0; i < nbr_of_timesteps/ir; i++) {
            //fprintf (file,"%e \t %e \t %e \n", i*dt, E_k0[i], E_k1[i]);
    		fprintf (file,"%e \t %e \t %e \t %e \t %e \t %e \t %e\n",
    			ir*i*dt, E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i], E_tot[i]);
    	}
    	fclose(file);
    	printf("energy.dat created!\n");
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
	
     // Boundary conditions
	int i;
	q[0] = 0;
	q[size_q-1] = 0;	

	for(i = 0; i < size_q; i++) {
		double qi = q[i], qp = 0, qm = 0;
		if (i > 0) 	  qm = q[i-1];
		if (i < size_q-1) qp = q[i+1];
		
		a[i] = qp - 2*qi + qm + alpha*(qp-qi)*(qp-qi) - alpha*(qi-qm)*(qi-qm);
	}
    //a[i] = q[i+1] - 2*q[i] + q[i-1] + alpha * (q[i+1]-q[i]) * (q[i+1]-q[i]) -  alpha * (q[i]-q[i-1])*(q[i]-q[i-1]);
}


