#include <math.h>
#include <stdlib.h>
#include "equilibration.h"

// http://fy.chalmers.se/~tfsgw/CompPhys/lectures/MD_LectureNotes_171029.pdf
#define k_B 0.0000080000617330

/*
 * vel: velocities, dt: timestep, Tau_eq: target temperature, Tau_T: time decay constant
 * Tau: current temperature, N: #particles
 */
void equib_temp(double (*vel)[3], double dt, double Tau_eq, double Tau_T,
	double Tau, int N) {	
	
	double alpha_T = 1 + dt/Tau_T * (Tau_eq - Tau)/Tau; // Eq. 111
	double factor = sqrt(alpha_T);
	
	int i;
	for (i = 0; i < N; ++i) // v_new = sqrt( alpha_T ) * v_old;
	{
		vel[i][0] = factor * vel[i][0];
		vel[i][1] = factor * vel[i][1];
		vel[i][2] = factor * vel[i][2];
	}
}

/*
 * pos: 3D positions, dt: timestep,
 * Tau: temperature, V: volume, , P: pressure, P_eq: target pressure, N: #particles
 */
void equib_pressure(double (*pos)[3], double dt,
	double Tau_P, double V, double P, double P_eq, int N) {
	
	//double kappa_T = -1/V * (dV/dT)_T; // Eq. 113
	double kappa_T = 0.002219033;
	double alpha_P = 1 + kappa_T*dt/Tau_P*(P_eq - P); // Eq. 112
	double factor = cbrt(alpha_P);
	
	int i;
	for (i = 0; i < N; ++i) // v_new = sqrt( alpha_T ) * v_old;
	{
		pos[i][0] = factor * pos[i][0];
		pos[i][1] = factor * pos[i][1];
		pos[i][2] = factor * pos[i][2];
	}
}

/*
 * m: mass per particle, E_k: kinetic energy, N: #particles
 */
double instantaneus_temp (double E_k, int N) {
	// Eq. 36, page 19
	return 2/(3*N*k_B) * E_k;
}

/*
 * Tau: temperature, V: volume, W: virial, N: #particles
 */
double pressure (double Tau, double V, double W, int N) {
	// Eq. 109, page 49
	return (N*k_B*Tau + W)/V;
}
