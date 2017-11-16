#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

#define N 256

/* Main program */
int main()
{
    //========================================================================//
    // Task specific - H1.1
    //========================================================================//                         
    double pos[N][3]; // x, y, z
    double vel[N][3];
    double acc[N][3];
    double f[N][3]; // Forces
    
    double a0 = 0.1; // Lattice parameter
    double L = 1; // Length of supercell
    
    int Nc = 50; // #primitive calls in each direction 
    double lattice[4*Nc*Nc*Nc][3];
    init_fcc(lattice, Nc, a0);

    //========================================================================//
    // Setup
    //========================================================================//
    int i, j, i_log;                                                               // i - actual timestep, i_log - logging of timestep data
    double dt = 0.1;
    double t_max = 250;
    int nbr_of_timesteps = t_max/dt;
    int ir = 100; // Resolution for i. Record every ir:th timestep.             // Segfault sensitive.
        
    // Data recording
    double log_data1 [nbr_of_timesteps/ir];
    double log_data2 [nbr_of_timesteps/ir];
    double log_data3 [nbr_of_timesteps/ir];
    
    //========================================================================//
    // Verlet
    //========================================================================//
    printf("\nLog resolution: 1 per %d steps, t_max = %.3f \n", ir, t_max);
    for (i = 1; i < nbr_of_timesteps; i++) {
        if (i%(nbr_of_timesteps/10) == 0) { // Print progress - 10%
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
        for (j = 0; j < N; j++) { // q(t+dt)
            pos[j][0] += dt * vel[j][0];
            pos[j][1] += dt * vel[j][1];
            pos[j][2] += dt * vel[j][2];
        }

        //=====================//
        // Accelerations
        //=====================//
        get_forces_AL(f,pos, L, N);  
        //calc_acc(a, q, m, kappa, nbr_of_particles)
        
        for (j = 0; j < N; j++) { // v(t+dt)
            vel[j][0] += dt * 0.5 * acc[j][0];
            vel[j][1] += dt * 0.5 * acc[j][1];
            vel[j][2] += dt * 0.5 * acc[j][2];
        }
        
        //====================================================================//
        // Record data (frequency depends on resolution ir);
        //====================================================================//
        if (i%ir == 0) { // Time to log data?
            i_log = i/ir;
            double energy = get_energy_AL(pos, L, N); // Potential energy
            
            log_data1[i_log] = energy;
            log_data2[i_log] = 2;
            log_data3[i_log] = 3;
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
            fprintf (file1,"%e \t %e \t %e \t %e     \n",
                t, // Time
                log_data1[i], log_data2[i], log_data3[i]); // data

            // Print file1
            // ...
        }
        
        
        // Close file(s)
        fclose(file1);
        printf("Data printed!\n");
    } else {
        printf("A file is NULL!\n");
    }
    
   
    /*
    int bound = 4;
    for (i = -8; i <= 7; i++) {
        int bounded = periodic_bound(i, bound);
        printf("%d\n", bounded);
    }
    printf("Looks aight.\n");
    */
    
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    //srand(time(NULL));
    double random_value = 0;
    //random_value = (double) rand() / (double) RAND_MAX;
    printf("random_value = %.5f \n", random_value);

    /*
     Descriptions of the different functions in the files initfcc.c and
     alpotential.c are listed below.
    */

    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    return 1;
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
    double energy = get_energy_AL(pos, L, N);
    printf("energy = %.5f", energy);
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    //double virial;
    //virial = get_virial_AL(pos, L, N);
    //printf("virial = %.5f", virial);
    
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    get_forces_AL(f,pos, L, N);  
    return 0; // Prevent err printout
}
