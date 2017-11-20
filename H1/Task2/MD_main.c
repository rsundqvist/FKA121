#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#include "utils.h"
#include "equilibration.h"

#define N 256

#define Tau_T 0.5
#define Tau_P 0.5
//#define Tau_eq 773.15
#define Tau_eq 5000
//#define P_eq 0.000000633
#define P_eq 1


void calc_acc(double mass, double (*f)[3], double (*acc)[3]);
double get_kinetic_energy(double mass, double (*vel)[3]);
void setArray3DToZero(double (*arr)[3]);
void randomizeLattice(double (*pos)[3],double a0);

/* Main program */
int main()
{
    //========================================================================//
    // Task specific - H1
    //========================================================================//   
    srand(time(NULL));                      
    double pos[N][3]; // x, y, z
    double vel[N][3];
    double acc[N][3];
    double f[N][3]; // Forces
    double fix = 10000;
    
    // Set arrays to zero
    setArray3DToZero(pos);  
    setArray3DToZero(vel);
    setArray3DToZero(acc);
    setArray3DToZero(f);  

    double a0 = 4.046; // Lattice parameter
    double mass = 0.002796304; // Mass of Al atom
    
    int Nc = 4; // #primitive calls in each direction 
    double L = Nc * a0; // Length of supercell
    double kappa_T = -2;
    double alpha_p;
    double * alpha_pP = &alpha_p;
    double V;
    init_fcc(pos, Nc, a0);
    randomizeLattice(pos, a0); // Introduce some random deviations
    
    //========================================================================//
    // Task specific - H1
    //========================================================================//               
    double T, P, P_prev; // Temperature, Pressure
    P = 0;
    //========================================================================//
    // Setup
    //========================================================================//
    int i, j, i_log;                                                               // i - actual timestep, i_log - logging of timestep data
    double dt = 0.001;
    double t_max = 10;
    int nbr_of_timesteps = t_max/dt;
    int ir = 100; // Resolution for i. Record every ir:th timestep.             // Segfault sensitive.
        
    // Data recording
    double log_data1 [nbr_of_timesteps/ir]; // E_p
    double log_data2 [nbr_of_timesteps/ir]; // E_k
    double log_data3 [nbr_of_timesteps/ir]; // E_tot = E_p + E_k
    double log_data4 [nbr_of_timesteps/ir]; // T, temperature
    double log_data5 [nbr_of_timesteps/ir]; // P, pressure
    
    double log_data6 [nbr_of_timesteps/ir]; // x
    double log_data7 [nbr_of_timesteps/ir]; // y
    double log_data8 [nbr_of_timesteps/ir]; // z
    
    double log_data9 [nbr_of_timesteps/ir]; // x
    double log_data10[nbr_of_timesteps/ir]; // y
    double log_data11[nbr_of_timesteps/ir]; // z
    
    int equilibrate = 500;
    //========================================================================//
    // Verlet
    //========================================================================//
    printf("\nLog resolution: 1 per %d steps, t_max = %.3f \n", ir, t_max);
    for (i = 1; i < nbr_of_timesteps; i++) {
        if (i%(nbr_of_timesteps/20) == 0) { // Print progress - 10%
            printf("\tt = %.2f \t\t %.3f  \n", i*dt, ((double)i/nbr_of_timesteps));
        }
      
        //======================================//
        // Verlet
        //======================================//
        for (j = 0; j < N; j++) { // v(t+dt/2)
            vel[j][0] += dt * 0.5 * acc[j][0];
            vel[j][1] += dt * 0.5 * acc[j][1];
            vel[j][2] += dt * 0.5 * acc[j][2];
        }
        //printf("f = (%2.2f, %2.2f, %2.2f) \n", vel[0][0], vel[0][1], vel[0][2]);
        for (j = 0; j < N; j++) { // q(t+dt)
            pos[j][0] = periodic_boundT( pos[j][0] + dt * vel[j][0], L );
            pos[j][1] = periodic_boundT( pos[j][1] + dt * vel[j][1], L );
            pos[j][2] = periodic_boundT( pos[j][2] + dt * vel[j][2], L );
        }

        //=====================//
        // Accelerations
        //=====================//
        get_forces_AL(f, pos, L, N);  
        calc_acc(mass, f, acc);
        
        for (j = 0; j < N; j++) { // v(t+dt)
            vel[j][0] += dt * 0.5 * acc[j][0];
            vel[j][1] += dt * 0.5 * acc[j][1];
            vel[j][2] += dt * 0.5 * acc[j][2];
        }
        
        double Ek = get_kinetic_energy(mass, vel);
        T = instantaneus_temp (Ek, N);
        P_prev = P;
        V = L*L*L; 
        double W = get_virial_AL(pos, a0, N);
        P = pressure (T, V, W, N);
        if (equilibrate) {
            //printf("fix = %.15f\n", fix);
            equib_pressure(pos, dt, Tau_P, P, P_eq, N, kappa_T, alpha_pP, fix); // Update pressure
            
            a0 *= *alpha_pP; // Rescale cell
            L = Nc * a0;
            /*printf("alpha_p = %.5f\n", alpha_p);
            printf("alpha_pP = %.5f\n", *alpha_pP);
            printf("a0 = %.5f\n", a0);
            printf("L = %.5f\n", L);*/
            
            equib_temp(vel, dt, Tau_eq, Tau_T, T, N); // Update temperature
        }
        
        //====================================================================//
        // Record data (frequency depends on resolution ir);
        //====================================================================//
        if (i%ir == 0) { // Time to log data?
            i_log = i/ir;
            double Ep = get_energy_AL(pos, L, N); // Potential energy
            log_data1[i_log] = Ep;
            log_data2[i_log] = Ek;
            log_data3[i_log] = Ek+Ep;
            
            log_data4[i_log] = T;
            //printf("\tTemp = %.5f \t P = %.5f \n", T,P);
            log_data5[i_log] = P;
            
            log_data6[i_log] = pos[0][0];
            log_data7[i_log] = pos[0][1];
            log_data8[i_log] = pos[0][2];
            
            log_data9[i_log] = pos[69][0];
            log_data10[i_log] = pos[69][1];
            log_data11[i_log] = pos[69][2];
            
            if (equilibrate) {
                equilibrate--;
                printf("Equilibration halt: %d \n", equilibrate);
                printf("\tTemp = %.5f \t P = %.5f \n", T,P);
            }
        }
    }
    printf("\tt = %.2f \t\t %.3f  \n", i*dt, ((double)i/nbr_of_timesteps));
    
    //====================================================================//
    // Simulation complete - print data to file(s)
    // fopen("filename", acc), acc \in {"r", "w", "o"} (read, write, append)
    //====================================================================//
    FILE * file1 = fopen("energy.dat", "w");
    if (file1 != NULL){
        printf("Print to file... ");

        for (i = 0; i < nbr_of_timesteps/ir; i++) { // i are log steps now!
            double t = ir*i*dt;
            
            // Print file1
            fprintf (file1,"%e \t %e \t %e \t %e \t %e \t %e, \t, %e \t %e \t %e \t %e \t %e \t %e \n",
                t, // Time
                log_data1[i], log_data2[i], log_data3[i],
                log_data4[i], log_data5[i],
                log_data6[i], log_data7[i], log_data8[i],
                log_data9[i], log_data10[i], log_data11[i]
                ); // data

            // Print file1
            // ...
        }
        
        
        // Close file(s)
        fclose(file1);
        printf("Data printed!\n");
    } else {
        printf("A file is NULL!\n");
    }

    return 0;
}

void calc_acc(double mass, double (*f)[3], double (*acc)[3]) {
    int i;
    for(i = 0; i < N; i++) {
        acc[i][0] = f[i][0]/mass;  
        acc[i][1] = f[i][1]/mass;
        acc[i][2] = f[i][2]/mass;  
    }
}

// Returns total kinetic energy for all particles
double get_kinetic_energy(double mass, double (*vel)[3]) {
    double energy = 0;
    int i;    
    for(i = 0; i < N; i++) {
        energy += vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2];
    }
    return 0.5 * energy *mass;
}

void setArray3DToZero(double (*arr)[3]) {
    int k;
    for(k = 0; k < N; k++) {
        arr[k][0] = 0;
        arr[k][1] = 0;
        arr[k][2] = 0;
    }
}

void randomizeLattice(double (*pos)[3],double a0) {
    int i;
    double r1,r2,r3;
    for(i = 0; i < N; i++) {
        r1 = (double)rand() / (double)RAND_MAX;
        r2 = (double)rand() / (double)RAND_MAX;
        r3 = (double)rand() / (double)RAND_MAX;
        pos[i][0] += - 0.05 *a0 + 0.1 * a0 * r1;
        pos[i][1] += - 0.05 *a0 + 0.1 * a0 * r2;
        pos[i][2] += - 0.05 *a0 + 0.1 * a0 * r3;
    }
}
