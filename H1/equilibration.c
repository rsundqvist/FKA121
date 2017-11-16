#include <math.h>
#include <stdlib.h>
#include "equilibration.h"

// http://fy.chalmers.se/~tfsgw/CompPhys/lectures/MD_LectureNotes_171029.pdf
#define K_b 8.617330*0.00001

/*
 * vel: velocities, dt: timestep, Tau_eq: target temperature, Tau_T: time decay constant
 * Tau: current temperature, N: #particles
 */
void equib_temp(double (*vel)[3], double dt, double Tau_eq, Tau_T,
	double Tau, int N) {	
	double alpha_T = 1 + dt/Tau_T * (T_eq - Tau)/Tau; // Eq. 111
	double factor = sqrt(alpha_T);
	
	int j;
	for (int i = 0; i < N; ++i) // v_new = sqrt( alpha_T ) * v_old;
	{
		vel[i][0] = factor * vel[i][0];
		vel[i][1] = factor * vel[i][1];
		vel[i][2] = factor * vel[i][2];
	}
}

/*
 * pos: 3D positions, dt: timestep,
 * Tau: temperature, V: volume, , P: preassure, N: #particles
 */
void equib_preassure(double (*pos)[3], double dt,
	double Tau, double V, double P, int N) {
	double kappa_T = -1/V * (dV/dT)_T; // Eq. 113
	double alpha_P = 1 + kappa_T*dt/Tau_P*(P_eq - P); // Eq. 112
	double factor = cbrt(3, alpha_P);
	
	int j;
	for (int i = 0; i < N; ++i) // v_new = sqrt( alpha_T ) * v_old;
	{
		pos[i][0] = factor * pos[i][0];
		pos[i][1] = factor * pos[i][1];
		pos[i][2] = factor * pos[i][2];
	}
}

/*
 * m: mass per particle, E_k: kinetic energy, N: #particles
 */
double instantaneus_temp (double m, double * E_k, int N) {
	// Eq. 36, page 19
	double sum = 0;
	for (int i = 0; i < N; ++i)
	{
		sum += E_k[i];
	}
	return 2/(3*N*k_B) * sum;
}

/*
 * Tau: temperature, V: volume, W: virial, N: #particles
 */
void preassure (double Tau, double V, double W, int N) {
	// Eq. 109, page 49
	return (N*k_B*T + W)/V;
}
