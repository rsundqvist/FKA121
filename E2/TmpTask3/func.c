<<<<<<< HEAD
=======
#include <stdio.h>

>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
/*
E1_func.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2015-10-23.
*/


/*
Function that calculates the acceleration based on the Hamiltonian.
The acceleration is calculated based on the displacements u and then stored in a.
u and a should be vectors of the same size, size_of_u
*/
<<<<<<< HEAD

void calc_acc(double *a, double *q, double alpha,int size_q)
{
   int i;
   q[0] = 0;
   q[size_q-1] = 0;
   /*for(i = 1; i < size_q-1; i++)
        a[i] = q[i+1] - 2*q[i] + q[i-1] + alpha * (q[i+1]-q[i])*(q[i+1]-q[i]) - alpha * (q[i]-q[i-1])*(q[i]-q[i-1]);

	
   a[0] = q[1] - 2 * q[0] + alpha*(q[1] - q[0])*(q[1] - q[0]) - alpha*q[0];
   a[size_q-1] = - 2 * q[size_q-1] + q[size_q-2] + alpha * q[size_q-1]*q[size_q-1] - alpha*(q[size_q-1] - q[size_q-2])*(q[size_q-1] - q[size_q-2]);
*/
   for (i = 1; i < size_q-1; i++) {
     double qi = q[i], qp = 0, qm = 0;
     if (i > 0)
	 qm = q[i - 1];
     if (i < size_q - 1) 
	qp = q[i + 1];
   a[i] = qp - 2 * qi + qm + alpha * (qp - qi) * (qp - qi) - alpha * (qi - qm) * (qi - qm);
   }
}

// Compute the energy of particle with index k at timestep t
double calc_E_k(double Pk, double omegaK, double Qk) {
	return 0.5* (Pk*Pk + omegaK*omegaK *Qk*Qk);
}
=======
void calc_acc(double *a, double *q, double alpha, int sz)
{
    /* Declaration of variables */
    int i;
    sz--;
    /* Calculating the acceleration on the boundaries */
    a[0] = q[1] - 2 * q[0] + alpha * (q[1] - q[0])*(q[1] - q[0]) - alpha*q[0]*q[0];
    a[sz] = q[sz - 1] - 2*q[sz] - alpha * (q[sz] - q[sz-1])* (q[sz] - q[sz-1]) + alpha*q[sz]*q[sz];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < sz ; i++){
        a[i] = q[i+1] - 2*q[i] + q[i-1] + alpha * (q[i+1]-q[i])* (q[i+1]-q[i]) - alpha*(q[i] - q[i-1])*(q[i] - q[i-1]);
        //printf("%.5f \n",a[i]);
    }
}


>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
