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

void calc_acc(double *a, double *q, double alpha,int size_q)
{
   int i;
   q[0] = 0;
   q[size_q - 1] = 0;

   for (i = 0; i < size_q; i++) {
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
