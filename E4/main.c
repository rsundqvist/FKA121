#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 0
#define k_B 0
#define m 0
#define T 0
#define eta 0

//==============//
// Main program //
//==============//
int main()
{
    //========================================================================//
    // Task specific - E4
    //========================================================================//
    double c_0 = exp(-eta*dt);
    double v_th = sqrt(k_B*T/m);        
    double q[N], v[N], a[N]; // pos, vel, acc
    
    //========================================================================//
    // Setup
    //========================================================================//
    gsl_rng *q = init_rng();
    int i, j, i_log;                                                               // i - actual timestep, i_log - logging of timestep data
    double dt = 0.01;
    double t_max = TIME_MAX;
    int nbr_of_timesteps = t_max/dt;
    int ir = 10; // Resolution for i. Record every ir:th timestep.             // Segfault sensitive.
    
    // Data recording
    double log_data1 [nbr_of_timesteps/ir];
    double log_data2 [nbr_of_timesteps/ir];
    double log_data3 [nbr_of_timesteps/ir];
    
    //========================================================================//
    // Verlet
    //========================================================================//
    printf("\nLog resolution: 1 per %d steps, t_max = %.3f \n", ir, t_max);
    
    for (i = 1; i < nbr_of_timesteps; i++) {
        if (100*(i-1)%(nbr_of_timesteps) == 0) { // Print progress
            printf("\tt = %.2f \t\t %.3f  \n", i*dt, ((double)i/nbr_of_timesteps));
            
            if (et > 0) {
                printf("et = %.2f ", et);
                printf("\tT = %.3f \t P = %.10f \n", T,P);
            }
        }
          
        G_1 = gsl_ran_gaussian(q, 1);
        G_2 = gsl_ran_gaussian(q, 1);
    
        //======================================//
        // Verlet
        //======================================//
        for (j = 0; j < N; j++) { // v(t+dt/2)
            v[j] += dt*0.5*a[j] + v_th*sqrt(1-c_0)*G_1;
        }
        //printf("f = (%2.2f, %2.2f, %2.2f) \n", vel[0][0], vel[0][1], vel[0][2]);
        for (j = 0; j < N; j++) { // q(t+dt)
            q[j] += dt * v[j];
        }

        //=====================//
        // Accelerations
        //=====================//
        calc_acc(mass, f, acc);
        
        for (j = 0; j < N; j++) { // v(t+dt)
            v[j] += 0.5*sqrt(c_0)*a[j]*dt + sqrt(c0)*v[k] + v_th*sqrt(1-c0)*G_2;
        }
        
        //====================================================================//
        // Record data
        //====================================================================//
        if (i%ir == 0) { // Time to log data?
            i_log = i/ir;
            log_data1[i_log] = 0;
            log_data2[i_log] = 0;
            log_data3[i_log] = 0;
        }
    }
    printf("\tt = %.2f \t\t %.3f  \n", i*dt, ((double)i/nbr_of_timesteps));
    
    //====================================================================//
    // Simulation complete - print data to file(s)
    //====================================================================//
    FILE * file1 = fopen("e4.dat", "w");
    if (file1 != NULL){
        for (i = 0; i < nbr_of_timesteps/ir; i++) { // i are log steps now!
            double t = ir*i*dt;
            
            // Print file1
            fprintf (file1,"%e \t %e \t %e \t %e \n",
                t, // Time
                log_data1[i], log_data2[i], log_data3[i],
                );

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

gsl_rng * init_rng()
{
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
	return q;
}