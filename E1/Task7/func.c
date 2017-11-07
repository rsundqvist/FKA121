/*
E1_func.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2015-10-23.
*/

#include <stdio.h>

/*
Function that calculates the acceleration based on the Hamiltonian.
The acceleration is calculated based on the displacements u and then stored in a.
u and a should be vectors of the same size, size_of_u
*/
void calc_acc(double *a, double *u, double* m, double kappa, int size_of_u)
{
    a[0] = kappa*(u[1] - u[0])/m[0]; // First particle - O
    a[1] = kappa*(u[2] - 2*u[1] + u[0])/m[1]; // Second (middle) particle - C
    a[2] = kappa*(u[1] - u[2])/m[2]; // Third particle - O
    //printf("a = (%e, %e, %e)\n", a[0], a[1], a[2]);
}

/* Function that calculates the potential energy based on the displacements */
double calc_pe(double *u, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    double e = 0;
    /* Calculating the energy on the boundaries */
    //e += kappa*((u[0] - u[1])*(u[0] - u[1])/2+u[0]*u[0]/2);
    //e += kappa*(u[size_of_u - 1])*(u[size_of_u - 1])/2;
    
    /* Calculating the energy of the inner points */
    for (i = 0; i < size_of_u - 1; i++){
        e += kappa*(u[i+1] - u[i])*(u[i+1] - u[i])/2;
    }
    return e;	
}


/* Function that calculates and returns the kinetic energy based on the velocities and masses */
double calc_ke(double *v, int size_of_v, double *m)
{
    /* Declaration of variables */
    int i;
    double e = 0; 
    /* Calculating the energy of the inner points */
    for (i = 0; i < size_of_v; i++){
        e += m[i]*(v[i])*(v[i])/2;
    }
    return e;	
}
