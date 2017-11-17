#include <stdio.h>

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


