/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h> 
#include <stdlib.h> 

#define PI 3.141592653589
#define nbr_of_particles 32

<<<<<<< HEAD
void foo(double * , double *, double (* A) [nbr_of_particles]);
void calc_E_k(double * omega, double * P, double * Q, double * E_k, double alpha, int t, int k);
void calc_acc(double * , double * , int, double);

/* Main program */
int main() {
  FILE * file;
  // Parameters
  double alpha = 0.01;
  double dt = 0.01;
  double t_max = 25000;
  int nbr_of_timesteps = t_max / dt;

  // displacement, velocity and acceleration
  double q[nbr_of_particles];
  double v[nbr_of_particles];
  double a[nbr_of_particles];
  double * p = v; // Since m = 1

  double Q[nbr_of_particles]; // discrete
  double P[nbr_of_particles]; // discrete

  int dm = 100;
  double E_k0[nbr_of_timesteps/dm];
  double E_k1[nbr_of_timesteps/dm];
  double E_k2[nbr_of_timesteps/dm];
  double E_k3[nbr_of_timesteps/dm];
  double E_k4[nbr_of_timesteps/dm];
  double E_tot[nbr_of_timesteps/dm];
  int i, j;
  double omega[nbr_of_particles];

  // Transformation matrix
  double factor;
  double trans_matrix[nbr_of_particles][nbr_of_particles];
  factor = 1 / ((double) nbr_of_particles + 1);
  for (i = 0; i < nbr_of_particles; i++) {
    for (j = 0; j < nbr_of_particles; j++) {
      trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
    }
  }

  // Initialize values
  E_k0[0] = nbr_of_particles;
  P[0] = sqrt(2 * E_k0[0]);
  foo(P, p, trans_matrix);
  foo(Q, q, trans_matrix);

  for (j = 0; j < nbr_of_particles; j++) {
    omega[j] = 2 * sin(j * PI / (2 * nbr_of_particles + 2));
  }

  double sum = 0;
  // timesteps according to velocity Verlet algorithm
  calc_acc(a, q, nbr_of_particles, alpha);
  printf("alpha = %.5f, t_max = %.2f \n", alpha, t_max);
  for (i = 1; i < nbr_of_timesteps + 1; i++) {
    if (i % 5000 == 0) {
      printf("\tt = %.2f\t\t\%: %.2f\n", i * dt, ((double) i / nbr_of_timesteps));
    }
    /* v(t+dt/2) */
    for (j = 0; j < nbr_of_particles; j++) {
      v[j] += dt * 0.5 * a[j];
    }
=======

void transform(double*, double*, double (*trans_matrix)[]);
void inverse_transform(double*, double*, double (*trans_matrix)[]);
void calc_E_k(double* omega, double* P, double* Q, double* E_k, double alpha, int t, int k);
void calc_acc(double*, double*, int, double);
void init_transformation_matrix(double (*trans_matrix)[]);

/* Main program */
int main() {
	FILE *file;
    // Parameters
	double alpha = 0;
	double dt = 0.1;
	double t_max = 2500;
	int nbr_of_timesteps = t_max/dt;
	int ir = 100; // Resolution for i. Record every ir:th timestep
	int size_q = nbr_of_particles+2;

    // displacement, velocity and acceleration
	double q[size_q];
	double v[size_q];
	double a[size_q]; 
    double *p = v; // Since m = 1

    double Q[nbr_of_particles]; // discrete
    double P[nbr_of_particles]; // discrete

    double E_k1 [nbr_of_timesteps/ir];
    double E_k2 [nbr_of_timesteps/ir];
    double E_k3 [nbr_of_timesteps/ir];
    double E_k4 [nbr_of_timesteps/ir];
    double E_k5 [nbr_of_timesteps/ir];
    double x_max[nbr_of_timesteps/ir];
    int i_max[nbr_of_timesteps/ir];
    int i,j;
    double omega[nbr_of_particles];
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    init_transformation_matrix(trans_matrix);


    //========================================================================//
    // Initialize
    //========================================================================//

    // Initialize values
    E_k1[0] = nbr_of_particles;
    P[1] = sqrt(2*E_k1[0]);
    transform(P, p, trans_matrix);
    transform(Q, q, trans_matrix);    
>>>>>>> dd2df75adf39020d68e0f3d95083c1432f9e53c4

    /* q(t+dt) */
    for (j = 0; j < nbr_of_particles; j++) {
<<<<<<< HEAD
      q[j] += dt * v[j];
    }

    /* a(t+dt) */
    calc_acc(a, q, nbr_of_particles, alpha);


<<<<<<< HEAD
    /* v(t+dt) */
    for (j = 0; j < nbr_of_particles; j++) {
      v[j] += dt * 0.5 * a[j];
    }
    // Update normal coordinates
    //printf("q_0 =  %.4f \n",q[0]);
    //printf("q_1 =  %.4f \n",q[1]);
    //printf("q_15 =  %.4f \n",q[15]);
    //printf("q_31 =  %.4f \n",q[31]);

    if (i%dm == 0) {
	int m = i/dm;
        calc_E_k(omega, P, Q, E_k0, alpha, m, 0);
        calc_E_k(omega, P, Q, E_k1, alpha, m, 1);
        calc_E_k(omega, P, Q, E_k2, alpha, m, 2);
        calc_E_k(omega, P, Q, E_k3, alpha, m, 3);
        calc_E_k(omega, P, Q, E_k4, alpha, m, 4);

        for (j = 0; j < nbr_of_particles; j++) {
            calc_E_k(omega, P, Q, E_tot, alpha, m, j);
        }
	
    }

    foo(q, Q, trans_matrix);
    foo(p, P, trans_matrix);
  }

  file = fopen("energy.dat", "w");
  if (file != NULL) {
    printf("%s", "Print to file: ");
    for (i = 0; i < nbr_of_timesteps; i++) {
      fprintf(file, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", dm * i * dt, E_k0[i], E_k1[i], E_k2[i], E_k3[i], E_k4[i], E_tot[i]);
    }
    fclose(file);
    printf("energy.dat created!\n");
  } else {
    printf("file is NULL\n");
  }
=======
    	omega[j] = 2*sin(j*PI/(2*nbr_of_particles+2));
=======
    // Eigenfrequencies
    for (j = 1; j < nbr_of_particles-1; j++) {
    	omega[j] = 2*sin(j*PI/(2*(nbr_of_particles+1)));
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
    }

    // timesteps according to velocity Verlet algorithm
    calc_acc(a, q, nbr_of_particles, alpha);
    printf("\n\talpha = %.5f, t_max = %.2f \n", alpha, t_max);
    //========================================================================//
    // Verlet
    //========================================================================//
    for (i = 1; i < nbr_of_timesteps+1; i++) {
    	if (i%(nbr_of_timesteps/10) == 0) {
    		printf("\tt = %.2f\t\t\%: %.2f\n", i*dt, ((double)i/nbr_of_timesteps));
    	}
    	
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
    	
        //====================================================================//
        // Inverse transform
        //====================================================================//
        inverse_transform(q, Q, trans_matrix); 
        inverse_transform(p, P, trans_matrix);
    
    	int i_log = i/ir;
    	if (i%ir == 0) {
            //================================================================//
            // Calculate energy and store
            //================================================================//
        	calc_E_k(omega, P, Q, E_k1, alpha, i_log, 1);
        	calc_E_k(omega, P, Q, E_k2, alpha, i_log, 2);
        	calc_E_k(omega, P, Q, E_k3, alpha, i_log, 3);
        	calc_E_k(omega, P, Q, E_k4, alpha, i_log, 4);
        	calc_E_k(omega, P, Q, E_k5, alpha, i_log, 5);

            // Total energy
        	for (j = 1; j < nbr_of_particles-1; j++) {
        		calc_E_k(omega, p, q, x_max, alpha, i_log, j);
        	}
        	/*
        	double max_abs = 0, max_real = 0;
        	double im = -1;
        	for (j = 0; j < nbr_of_particles; j++) {
        	    double x = v[j];
        	    if (abs(x) > max_abs) {
        	        im = j;
        	        max_abs = abs(x);
        	        max_real = x;
        	    }
        	}
        	x_max[i_log] = max_real;
        	i_max[i_log] = im;*/
    	} // End recording of data/calculation of energy
    } // End for

    
    file = fopen("energy.dat","w");
    if (file != NULL){
    	printf("%s", "Print to file: ");
    	for (i = 0; i < nbr_of_timesteps/ir; i++) {
    		fprintf (file,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %d \n",
    			ir*i*dt, // Time
    			E_k1[i], E_k2[i], E_k3[i], E_k4[i], E_k5[i], x_max[i], // Energies
    			i_max[i]); // max value index (whatever is measured)
    	}
    	fclose(file);
    	printf("energy.dat created!\n");
    	
    	// Print to console
    	for (i = 0; i > nbr_of_timesteps/(ir*1000); i++) {
    		printf ("%e \t %e \t %e \t %e \t %e \t %e \t %e \t %d \n",
    			1000*ir*i*dt, // Time
    			E_k1[i], E_k2[i], E_k3[i], E_k4[i], E_k5[i], x_max[i], // Energies
    			i_max[i], 0, 0); // max value index (whatever is measured)
    	}
    } else {
    	printf("file is NULL\n");
    }
}

<<<<<<< HEAD
    printf("sum = %.4f\n",sum);
>>>>>>> dd2df75adf39020d68e0f3d95083c1432f9e53c4

  printf("sum = %.4f\n", sum);
=======
// ===========================================================================//
// Transformations
// ===========================================================================//
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32

// The formalae for Q_k and P_k are identical with m = 1
void transform(double *q, double *Q, double (*trans_matrix)[nbr_of_particles])
{	
    /* Transformation to normal modes Q from displacements q.  */
	double sum;
	int i, j, bxb = 0;
	for (i = bxb; i < nbr_of_particles-bxb; i++){
		sum = 0;
		for (j = bxb; j < nbr_of_particles-bxb; j++){
			sum += q[j] * trans_matrix[i][j];
		}
		Q[i] = sum;
	}
}

// The formalae for Q_k and P_k are identical with m = 1
<<<<<<< HEAD
<<<<<<< HEAD

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

void calc_E_k(double * omega, double * P, double * Q, double * E_k, double alpha, int t, int k) {
  E_k[t] += 0.5 * (P[k] * P[k] + omega[k] * omega[k] * Q[k] * Q[k]);
}

void calc_acc(double * a, double * q, int size_q, double alpha) {

  // Boundary conditions
  int i;
  q[0] = 0;
  q[size_q - 1] = 0;
=======
void transform(double *q, double *Q, double (*trans_matrix)[nbr_of_particles])
=======
void inverse_transform(double *q, double *Q, double (*trans_matrix)[nbr_of_particles])
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
{	
    /* Transformation to normal modes Q from displacements q.  */
	double sum;
	int i, j, bxb = 0;
	for (i = bxb; i < nbr_of_particles-bxb; i++){
		sum = 0;
		for (j = bxb; j < nbr_of_particles-bxb; j++){
			sum += q[j] * trans_matrix[i][j];
		}
		Q[i] = sum;
	}
}


// ===========================================================================//
// Create matrix
// ===========================================================================//

void init_transformation_matrix(double (*trans_matrix)[nbr_of_particles]) {
	int i,j, bxb = 0;
	double factor;	
	factor = 1 / ((double) nbr_of_particles + 1);
	for (i=bxb; i < nbr_of_particles-bxb; i++) {
		for (j=bxb; j < nbr_of_particles-bxb; j++) {
			trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
		}
	}
}


// ===========================================================================//
// Ok
// ===========================================================================//
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
		
		if (i > 0)        qm = q[i-1];
		if (i < size_q-1) qp = q[i+1];
		
		a[i] = qp - 2*qi + qm + alpha*(qp-qi)*(qp-qi) - alpha*(qi-qm)*(qi-qm);
    }
}
>>>>>>> dd2df75adf39020d68e0f3d95083c1432f9e53c4

  for (i = 0; i < size_q; i++) {
    double qi = q[i], qp = 0, qm = 0;
    if (i > 0) qm = q[i - 1];
    if (i < size_q - 1) qp = q[i + 1];

    a[i] = qp - 2 * qi + qm + alpha * (qp - qi) * (qp - qi) - alpha * (qi - qm) * (qi - qm);
  }
  //a[i] = q[i+1] - 2*q[i] + q[i-1] + alpha * (q[i+1]-q[i]) * (q[i+1]-q[i]) -  alpha * (q[i]-q[i-1])*(q[i]-q[i-1]);
}
